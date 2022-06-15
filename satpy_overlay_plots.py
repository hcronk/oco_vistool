"""
satpy_overlay_plots.py

module implementing satpy based plots with OCO-2 data overlays.

intended to integrate with oco_vistool.py
"""

import os, glob, datetime, collections, itertools
import tempfile, shutil
import bz2

import numpy as np

from satpy import Scene 
from satpy.writers import get_enhanced_image
import pyresample
from pyorbital.orbital import get_observer_look

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.io import shapereader

import shapely.geometry as sgeom

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatches

import netCDF4
import h5py

import boto3
import botocore
from botocore.handlers import disable_signing

import fnmatch

from oco_vistool import read_shp
from geo_scan_time import compute_time_offset, approx_scanline_time
from geo_scan_time import get_ABI_timerange_from_filename, get_AHI_timerange_from_filename

import vistool_lib as vl

# to get the location of the city data (a subdir from where the
# code is located).
_code_dir = os.path.dirname(os.path.realpath(__file__))

def _get_view_zenith(scn, obs_lon, obs_lat):
    """
    helper function to compute the view zenith angle to a observed lat/lon
    point (at zero altitude) for the given satellite.
    The satellite position is stored within the satpy Scene object,
    so we get the required info from the scn metadata. Note this only works
    after the true_color composite is loaded.
    """

    obs_datetime = scn.attrs['start_time']
    obs_alt = 0.0

    # orbital parameters are stored in different keys for
    # GOES and Himawari; first branch is the GOES version.
    if 'orbital_parameters' in scn['true_color'].attrs:
        orb_pars = scn['true_color'].attrs['orbital_parameters']
        sat_lon = orb_pars['satellite_nominal_longitude']
        sat_lat = orb_pars['satellite_nominal_latitude']
        sat_alt = orb_pars['satellite_nominal_altitude']
    else:
        sat_lon = scn['true_color'].attrs['satellite_longitude']
        sat_lat = 0.0
        sat_alt = scn['true_color'].attrs['satellite_altitude']

    # in scn object, this is units [m], need [km] for pyorbital function.
    sat_alt *= 1e-3
    # get_observer_look expects array like inputs; if we put sat_lon inside
    # a 1-element list, everything gets broadcast to that shape.
    _, view_el = get_observer_look(
        [sat_lon], sat_lat, sat_alt,
        obs_datetime, obs_lon, obs_lat, obs_alt)
    # convert back to scalar while changing from elevation to zenith.
    view_zenith = 90 - view_el[0]

    return view_zenith


def _get_AHI_times(file_list):
    """
    helper function to get the start and end time from a set of
    Himawari HSD segment files.
    Currently, the only straightforward way to do this appears to be
    a satpy Scene object. I think, without calling the load() function,
    this should not be an expensive operation. However, this does require
    the full datafiles are downloaded.
    """
    
    scn = Scene(reader = 'ahi_hsd', filenames = file_list)
    time_range = scn.start_time, scn.end_time

    return time_range


def get_loc_ABI_files(datetime_utc, data_home, domain, platform, hour_offsets, bands_list, glob_fstr, verbose):
    """
    Helper function for getting the local filesystem GOES ABI files.
    
    Returns a list with relevant filenames.
    
    inputs:
    datetime_utc: datetime object for the needed timeslot
    data_home: origin for the background files (by format: data_home/YYYY/...)
    domain: (C)ONUS or (F)ull disk
    platform: GOES-East (G16) or GOES-West (G17)
    hour_offsets: different offsets to iterate
    bands_list: list of needed bands
    glob_fstr: glob file format
    verbose: T/F - if details are needed
    """
    
    flists = collections.defaultdict(list)
    # iterating over all time options and adding the located files to the list
    for hour_offset, band in itertools.product(hour_offsets, bands_list):
        dt = datetime_utc - datetime.timedelta(hours=hour_offset)

        if verbose:
            print('Searching for data at hour: '+str(dt))
        
        # constructing path to the needed files
        y_dir = dt.strftime('%Y')
        ymdd_dir = dt.strftime('%Y_%m_%d_%j')
        
        # four "?" for the MMSS, and then one more for a character at
        # the end, I am not sure about (frac second?)
        stimestamp = dt.strftime('%Y%j%H') + '????' + '?'
        #print(data_home)
        ddir = os.path.join(data_home, y_dir, ymdd_dir,
                            'abi', 'L1b', 'Rad'+domain)
        if verbose:
            print('composed data dir: '+ddir)

        if not os.access(ddir, os.F_OK):
            print('Data dir not found: '+ddir)
            continue

        glob_str = os.path.join(ddir,glob_fstr.format(
            domain, band, platform, stimestamp))
        if verbose:
            print('glob: '+glob_str)

        located_files = glob.glob(glob_str)
        flists[band] += located_files
    return flists

def get_aws_ABI_files(datetime_utc, domain, platform, hour_offsets, bands_list, glob_fstr, verbose):
    """
    Helper function for getting the AWS S3 GOES ABI filenames.
    
    Returns a list with needed filenames and the AWS S3 bucket object (where the files will be downloaded).
    
    inputs:
    datetime_utc: datetime object for the needed point of time
    domain: (C)ONUS or (F)ull disk
    platform: GOES-East (G16) or GOES-West (G17)
    hour_offsets: different offsets to iterate over
    bands_list: list of needed bands
    glob_fstr: glob file format
    verbose: T/F - if details are needed
    """
    
    # accessing the relevant S3 bucket
    s3 = boto3.resource('s3')
    
    # See https://stackoverflow.com/questions/34865927/can-i-use-boto3-anonymously for reference
    s3.meta.client.meta.events.register('choose-signer.s3.*', disable_signing)
    
    if (platform[-2:] == '16'):
        g_bucket = s3.Bucket('noaa-goes16')
    else:
        g_bucket = s3.Bucket('noaa-goes17')
        
    flists = collections.defaultdict(list)
    
    # iterating over the different hours to locate the needed files
    for hour_offset, band in itertools.product(hour_offsets, bands_list):
        dt = datetime_utc - datetime.timedelta(hours=hour_offset)
        if verbose:
            print('Searching for data at hour: '+str(dt))
        
        # constructing the file path
        y_dir = dt.strftime('%Y')
        dayY_dir = dt.strftime('%j')
        h_dir = dt.strftime('%H')
        
        # four "?" for the MMSS, and then one more for a character at
        # the end, I am not sure about (frac second?)
        stimestamp = dt.strftime('%Y%j%H') + '????' + '?'
        ddir = os.path.join('ABI-L1b-Rad' + domain, y_dir, dayY_dir, h_dir)
        g_path = ddir + '/'
        g_files = g_bucket.objects.filter(Prefix=g_path)

        glob_str = os.path.join(ddir,glob_fstr.format(domain, band, platform, stimestamp))
        if verbose:
            print('glob: '+glob_str)
        
        # access the AWS files by the constructed format
        keys = list()
        for g_file in g_files:
            keys.append(g_file.key)
        located_files = fnmatch.filter(keys, glob_str)
        flists[band] += located_files
    return flists, g_bucket

def download_aws_ABI_files(flist, g_bucket, data_home):
    """
    Helper function for downloading the AWS GOES files given their names.

    Returns a list with downloaded local GOES ABI filenames.
    
    inputs:
    flist: list of needed files paths
    g_bucket: AWS S3 bucket to download the files from
    data_home: data destination for the background files (by format: data_home/ABI-Rad<domain>/...)
    """
    files = list()
    # iterating to and downloading the needed files from the S3 bucket (by the known path)
    for needed_path in flist:
        g_files = g_bucket.objects.filter(Prefix=needed_path)
        for g_file in g_files:
            if not (os.path.isfile(data_home + '/' + g_file.key)):
                if not os.path.isdir(os.path.dirname(data_home + '/' + g_file.key)):
                    os.makedirs(os.path.dirname(data_home + '/' + g_file.key))
                g_bucket.download_file(g_file.key, data_home + '/' + g_file.key)
                files.extend(glob.glob(data_home + '/' + g_file.key))
            else:
                files.extend(glob.glob(data_home + '/' + g_file.key))

    flist = list(set(files))
    return flist

def get_ABI_files(datetime_utc, center_lat,
                  sensor, files_loc, data_home, verbose=False):
    """
    Using the helpers above, accesses the needed GOES ABI files and downloads them (if needed).
    
    Returns a list with accessed GOES filenames and time offsets.
    
    inputs:
    datetime_utc: datetime object for the needed point of time
    center_lat: center latitude of the background image
    sensor: the geostation, device, and domain (GOES16_ABI_C, for example)
    files_loc: if the background images are local or on AWS
    data_home: data origin for the background files (by format: data_home/YYYY/...)
    verbose: if details are needed
    """
    
    # checking the input validity
    domain = sensor[-1]
    platform = sensor.split('_')[0]

    valid_domains = 'C', 'F'
    if domain not in valid_domains:
        raise ValueError('domain must be in '+str(valid_domains))
    valid_platforms = 'GOES16', 'GOES17'
    if platform not in valid_platforms:
        raise ValueError('platform must be in '+str(valid_platforms))
    # code used in filename is G16 or G17 for GOES16, GOES17.
    platform = platform.replace('GOES','G')

    # search within hour of input datetime, and +/-1 hour
    band_list = (1,2,3)
    hour_offsets = (-1,0,1)
    # makes a glob string to match by s-time (I think, start time.)
    glob_fstr = ('OR_ABI-L1b-Rad{0:1s}'+
                 '-M?C{1:02d}_{2:s}_s{3:s}_e*_c*.nc')
    
    # get the needed filenames (and the S3 bucket for AWS)
    if (files_loc == 'aws'):
        flists, g_bucket = get_aws_ABI_files(datetime_utc, domain, platform, hour_offsets, band_list, glob_fstr, verbose)
    else:
        flists = get_loc_ABI_files(datetime_utc, data_home, domain, platform, hour_offsets, band_list, glob_fstr, verbose)
    
    # have a set of 3-hour long lists of files, now just find the one with the
    # smallest time offset to the input datetime_utc.
    flist = []
    time_offset = []
    for band in band_list:
        nfiles = len(flists[band])
        if nfiles == 0:
            print('no files found for ABI band '+str(band))
            continue
        time_offsets_all = []
        for n in range(nfiles):
            # find the approx scan line time for the requested latitude.
            etime, stime = get_ABI_timerange_from_filename(flists[band][n])
            scan_time = approx_scanline_time(stime, etime, center_lat, domain)
            time_offsets_all.append((datetime_utc-scan_time).total_seconds())
        # find minimum time difference
        k = np.argmin(np.abs(time_offsets_all))
        flist.append(flists[band][k])
        time_offset.append(time_offsets_all[k])
    
    # download the files if needed
    if (files_loc == 'aws'):
        flist = download_aws_ABI_files(flist, g_bucket, data_home)
    return flist, time_offset

def get_loc_AHI_files(datetime_utc, data_home, offsets, bands_list):
    """
    Helper function for getting the local filesystem Himawari AHI files.
    
    Returns a list with relevant filenames.
    
    inputs:
    datetime_utc: datetime object for the needed timeslot
    data_home: origin for the background files (by format: data_home/YYYY/...)
    offsets: different time offsets to iterate
    bands_list: list of needed bands
    resolutions_list: list of desired image resolutions, same length as bands_list
        (himawari data has multiple resolution products for each band, in separate files)
    """   
    files = collections.defaultdict(list)

    # round down to nearest 10 minutes.
    rounded_minutes = 10 * (datetime_utc.minute//10)
    rounded_datetime_utc = datetime.datetime(
        datetime_utc.year, datetime_utc.month, datetime_utc.day,
        datetime_utc.hour, rounded_minutes, 0)

    # access and record all the background files by the known path format
    strftime_template = '%Y/%Y_%m_%d_%j/%H%M/HS_H08_%Y%m%d_%H%M'
    for offset, (band, res) in itertools.product(offsets, zip(bands_list, resolutions_list)):
        offset_dt = rounded_datetime_utc + datetime.timedelta(minutes=offset)
        # searches for files with the selected band Bnn, at specified
        # resolution Rnn and use S??10 (ignore the not-segmented file)
        glob_str = os.path.join(
            data_home, (offset_dt.strftime(strftime_template) +
                        '_B{0:02d}_FLDK_R{1:02d}_S??10.DAT'.format(band, res)))
        files_at_offset = glob.glob(glob_str)
        files_at_offset.sort()

        # if data is missing at this time offset, the files_at_offset list
        # will be empty. In that case, don't add to the list.
        if len(files_at_offset) > 0:
            files[band].append(files_at_offset)

    return files

def get_aws_AHI_files(datetime_utc, offsets, bands_list, resolutions_list):
    """
    Helper function for accessing the paths of relevant AWS Himawari background files.

    Returns a dicttionary with AWS S3 Himawari file paths (key: band number, value: list of lists of filepaths by time offset).
    
    inputs:
    datetime_utc: datetime object for the needed timeslot
    offsets: different time offsets to iterate
    bands_list: list of needed bands
    resolutions_list: list of desired image resolutions, same length as bands_list
        (himawari data has multiple resolution products for each band, in separate files)
    """
    files = collections.defaultdict(list)

    # round down to nearest 10 minutes.
    rounded_minutes = 10 * (datetime_utc.minute//10)
    rounded_datetime_utc = datetime.datetime(
        datetime_utc.year, datetime_utc.month, datetime_utc.day,
        datetime_utc.hour, rounded_minutes, 0)

    # connecting to the S3 Himawari bucket
    s3 = boto3.resource('s3')
    
    # See https://stackoverflow.com/questions/34865927/can-i-use-boto3-anonymously for reference
    s3.meta.client.meta.events.register('choose-signer.s3.*', disable_signing)
    hima_bucket = s3.Bucket('noaa-himawari8')
    
    # known AWS Hima filepath format
    strftime_template = 'AHI-L1b-FLDK/%Y/%m/%d/%H%M/'

    for offset, (band, res) in itertools.product(offsets, zip(bands_list, resolutions_list)):
        offset_dt = rounded_datetime_utc + datetime.timedelta(minutes=offset)

        # access the background files by the known path format
        hima_files = hima_bucket.objects.filter(Prefix=offset_dt.strftime(strftime_template))

        # searches for files with the selected band Bnn, at specified
        # resolution Rnn and use S??10 (ignore the not-segmented file)
        glob_str = offset_dt.strftime(strftime_template + 'HS_H08_%Y%m%d_%H%M') + '_B{0:02d}_FLDK_R{1:02d}_S??10.DAT.bz2'.format(band,res)

        # record the AWS files by the needed format
        keys = list()
        for hima_file in hima_files:
            keys.append(hima_file.key)
            
        files_at_offset = fnmatch.filter(keys, glob_str)
        files_at_offset.sort()
           
        # if data is missing at this time offset, the files_at_offset list
        # will be empty. In that case, don't add to the list.
        if len(files_at_offset) > 0:
            files[band].append(files_at_offset)

    return files, hima_bucket    

def download_aws_AHI_files(flist, hima_bucket, data_home):
    """
    Helper function for downloading and decompressing the AWS Himawari files given their AWS S3 paths.

    Returns a list with downloaded local Himawari AHI filenames.
    
    inputs:
    flist: list of needed AWS file paths
    hima_bucket: AWS S3 bucket to download the files from
    data_home: data destination for the background files (by format: data_home/AHI-FLDK/...)
    """
    files = list()

    # iterating to and downloading the needed files from the S3 bucket (by the known path)
    for needed_path in flist:
        hima_files = hima_bucket.objects.filter(Prefix=needed_path)
        # downloading and recording the needed Himawari files by known AWS paths
        for hima_file in hima_files:
            if not (os.path.isfile(data_home + '/' + hima_file.key[:-4]) or os.path.isfile(data_home + '/' + hima_file.key)):
                if not os.path.isdir(os.path.dirname(data_home + '/' + hima_file.key)):
                    os.makedirs(os.path.dirname(data_home + '/' + hima_file.key))
                hima_bucket.download_file(hima_file.key, data_home + '/' + hima_file.key)
            else:
                files.extend(glob.glob(data_home + '/' + hima_file.key[:-4]))
                
        # decompressing the downloaded files
        for filename in glob.glob(data_home + '/' + needed_path):
            zipfile = bz2.BZ2File(filename)
            data = zipfile.read()
            new_filename = filename[:-4] 
            open(new_filename, 'wb').write(data)
            os.remove(filename)
            files.extend(glob.glob(new_filename))
    flist = list(set(files))
    return flist

def get_AHI_files(datetime_utc, center_lat, files_loc, data_home):
    """
    Using the helpers above, accesses the needed Himawari files by date & time and downloads them (if needed).
    
    Returns a list with needed Himawari background filenames and a time offset of the image relative to the request.
    
    inputs:
    datetime_utc: datetime object for the needed point of time
    center_lat: center latitude of the background image
    files_loc: if the background images are local or on AWS
    data_home: origin directory for the background files 
    """

    bands_list = [1, 2, 3, 4]
    offsets = (-10, 0, 10)
    # himawari data contains multiple resolutions for each band;
    # we only want to download the highest resolution for each, as this
    # is what satpy uses to make the composite RGB.
    resolutions_list = [10, 10, 5, 10]

    # accessing the Himawari files paths and bucket from AWS
    if (files_loc == 'aws'):
        files, hima_bucket = get_aws_AHI_files(
            datetime_utc, offsets, bands_list, resolutions_list)

    # accessing the local files
    else:
        files = get_loc_AHI_files(
            datetime_utc, data_home, offsets, bands_list, resolutions_list)

    flist = []
    time_offset = []
    for band in bands_list:
        nfiles = len(files[band])
        if nfiles == 0:
            print('no files found for AHI band '+str(band))
            continue
        time_offsets_all = []
        for n in range(nfiles):
            # use rough approximation: assume the himawari start/end times are exactly
            # the time from the file path, and +10 later.
            # this means we don't have to download the file data.
            # the more accurate variation would use the _get_AHI_times() helper
            # function to get the actual start and end time.
            stime, etime = get_AHI_timerange_from_filename(files[band][n][0])
            scan_time = approx_scanline_time(stime, etime, center_lat, 'F')
            time_offsets_all.append((datetime_utc-scan_time).total_seconds())
        # find minimum time difference
        k = np.argmin(np.abs(time_offsets_all))
        flist.extend(files[band][k])
        time_offset.append(time_offsets_all[k])
    

    # downloading the recorded/accessed files from AWS
    if (files_loc == 'aws'):
        flist = download_aws_AHI_files(flist, hima_bucket, data_home)

    return flist, time_offset

def get_scene_obj(file_list, latlon_extent, sensor, width=750, height=750,
                  tmp_cache=False, resample_method='native_bilinear'):
    """Get Scene object, apply the resample area, to a small box 
    centered on the lat lon point.

    inputs:
    file_list: list of netCDF L1b ABI files, must contain bands 1,2,3
    latlon_extent: extent of displayed region in latlon:
        [min_lon, min_lat, max_lon, max_lat]
    width: number of resampled image pixels in width (x-dimension)
    height: number of resampled image pixels in width (y-dimension)

    tmp_cache: optional keyword, set to True to copy the located files
    to a temporary dir (from tempfile) before loading.
    Note the temporary files are NOT cleaned up.
    Use this option if the data_home access is slow or limited in some
    way, which can be mitigated by copying to a temp file on the local
    filesystem

    resample_method: string keyword to specify the resampling method.
    valid options are:
    nearest: perform nearest neighbor resampling in one step
    bilinear: perform bilinear interpolation in one step
    native_nearest: perform a native interpolation first, to upsample
       lower resolution bands to the highest native resolution; 
       then perform nearest neighbor interpolation to the output grid.
    native_nearest: perform a native interpolation first, to upsample
       lower resolution bands to the highest native resolution; 
       then perform bilinear interpolation to the output grid.

    outputs: the satpy Scene object.

    """

    valid_resample_methods = [
        'nearest', 'bilinear', 'native_nearest', 'native_bilinear']

    if resample_method not in valid_resample_methods:
        raise ValueError('resample_method ' + resample_method +
                         ' is not valid')

    if tmp_cache:
        tdir = tempfile.mkdtemp()
        cached_file_list = []
        for f in file_list:
            src_f = f
            dst_f = os.path.join(tdir, os.path.split(f)[-1])
            shutil.copyfile(src_f, dst_f)
            cached_file_list.append(dst_f)
        scn = Scene(reader=sensor, filenames=cached_file_list)
    else:
        scn = Scene(reader=sensor, filenames=file_list)
    
    scn.load(['true_color'])

    my_area = pyresample.create_area_def(
        'testC', "epsg:3857", width=width, height=height,
        area_extent=latlon_extent, units='degrees')

    if resample_method.startswith('native'):
        tmp_scn = scn.resample(resampler='native')
    else:
        tmp_scn = scn

    # this would split the second string str1_str2,
    # or just return the str if there is no underscore.
    # thus, it should be the resample method after the
    # optional native resampling.
    method = resample_method.split('_')[-1]
    new_scn = tmp_scn.resample(my_area, resampler=method)

    return new_scn


def _setup_axes(fignum, crs, figsize=(10,8), create_colorbar_axis=True):
    """
    Note: this is the OLD method, which is superseded by what is in
    vistool_lib. Leaving it in here for now, until we decided it isn't needed...

    setup gridspec axes

    fignum: MPL figure number
    crs: cartopy transform object

    figsize: optional figure size to send to figure creator.
    create_colorbar_axis: controls whether a colorbar axis is created
        within the MPL figure.

    returns:
    gs: gridspec object
    ax: MPL axis object containing image
    cb_ax: MPL colorbar axis object (If created), otherwise None
    inset_ax: inset axis
    """
    fig = plt.figure(fignum, figsize = figsize)
    fig.clf()

    gs = mpl.gridspec.GridSpec(16, 16)

    ax = plt.subplot(gs[0:-1,2:-2], projection=crs)
    if create_colorbar_axis:
        cb_ax = plt.subplot(gs[0:-1, -1])
    else:
        cb_ax = None

    inset_ax = plt.subplot(gs[6:9, 0:2], projection=ccrs.PlateCarree())
    
    return gs, ax, cb_ax, inset_ax


def annotate_locations(ax, cfg_d):
    """
    add a marker for the ground_site, and city_labels, if set in
    the config dictionary.

    inputs:
    ax: the image display axis
    cfg_d: the configuration dictionary
    """

    if cfg_d['ground_site']:
        print('plotting ground site at: ',
              cfg_d['ground_site'][1], cfg_d['ground_site'][0])
        ax.plot(cfg_d['ground_site'][1], cfg_d['ground_site'][0],
                'w*', markersize=10, transform=ccrs.Geodetic())

    if cfg_d['city_labels']:
        
        populated_places_filename = os.path.join(
            _code_dir , 'natural_earth/ne_10m_populated_places')

        print('loading city label file from:')
        print(_code_dir)
        print(populated_places_filename)

        df = read_shp(populated_places_filename)
        
        relevant_places = df[(df['LATITUDE'] <= cfg_d['lat_ul']) &
                             (df['LATITUDE'] >= cfg_d['lat_lr']) &
                             (df['LONGITUDE'] <= cfg_d['lon_lr']) &
                             (df['LONGITUDE'] >= cfg_d['lon_ul'])]

        for idx, p in relevant_places.iterrows():
            ax.text(p['LONGITUDE'], p['LATITUDE'], p['NAME'],
                    fontsize=7, color=cfg_d['city_labels'],
                    va='bottom', ha='center',
                    transform=ccrs.Geodetic())


def annotate_ll_axes(gs, img_ul, img_lr):
    """
    annotate min/max lat and lon into new, independent axis objects
    in the gridspec.

    inputs:
    gs: the MPL gridspec object
    img_ul, img_lr: the upper left, lower right lat,lon coordinates
        for the displayed image
    """

    ax_minlat = plt.subplot(gs[-1, 1])
    ax_maxlat = plt.subplot(gs[0, 1])
    ax_minlon = plt.subplot(gs[-1, 3])
    ax_maxlon = plt.subplot(gs[-1, -3])

    ax_minlat.axis('off')
    ax_maxlat.axis('off')
    ax_maxlon.axis('off')
    ax_minlon.axis('off')

    fontkw = dict(fontsize=10, fontweight='bold',
                  ha='center', va='center')
    ax_minlat.text(0.5, 1.0, str("%.1f" % img_lr[0]), **fontkw)
    ax_maxlat.text(0.5, 1.0, str("%.1f" % img_ul[0]), **fontkw)
    ax_minlon.text(0.5, 0.0, str("%.1f" % img_ul[1]), **fontkw)
    ax_maxlon.text(0.5, 0.0, str("%.1f" % img_lr[1]), **fontkw)
            
    return (ax_minlat, ax_maxlat, ax_minlon, ax_maxlon)


def make_inset_map(inset_ax, img_ul, img_lr):
    """
    make an inset plot (a small thumbnail image with the coastlines,lakes
    drawn, and a red square showing the LL box for the displayed image)

    inputs:
    inset_ax: MPL axis to use for drawing the map
    img_ul: Upper Left Lat/Lon corner of image
    img_lr: Lower Right Lat/Lon corner of image
    """
    
    inset_extent_x = [img_ul[1], img_lr[1]]
    inset_extent_y = [img_lr[0], img_ul[0]]

    inset_extent_x = [x + 360 if x < 0 else x for x in inset_extent_x]
    inset_extent_y = [y + 180 if y < 0 else y for y in inset_extent_y]

    inset_extent_x[0] -= 20
    inset_extent_y[0] -= 20
    inset_extent_x[1] += 20
    inset_extent_y[1] += 20

    inset_extent_x = [x - 360 if x > 180 else x for x in inset_extent_x]
    inset_extent_y = [y - 180 if y > 90 else y for y in inset_extent_y]

    img_extent = [img_ul[1], img_lr[1], img_lr[0], img_ul[0]]
    extent = [inset_extent_x[0], inset_extent_x[1],
              inset_extent_y[0], inset_extent_y[1]]

    inset_ax.set_extent(extent)

    inset_ax.coastlines()
    inset_ax.add_feature(cfeature.LAKES, edgecolor='black', facecolor='none')
    extent_box = sgeom.box(img_extent[0], img_extent[2], img_extent[1],
                           img_extent[3])
    inset_ax.add_geometries([extent_box], ccrs.PlateCarree(),
                            color='none', edgecolor='red')
    inset_ax.set_aspect('auto')


def plot_scene_obj(ax, scn):
    """
    plot a scene object into a new axis, created in a cleared new figure.

    inputs:
    ax: axis to use for displaying image (via imshow)
    scn: satpy Scene object

    returns:
    MPL image object
    """

    crs = scn['true_color'].attrs['area'].to_cartopy_crs()

    # note, trying to imshow directly on the data attribute of scn,
    # causes exception in MPL:
    # TypeError: Invalid dimensions for image data
    # because shape is (3,600,600), needs to be (600,600,3)?
    #
    # get_enhanced_image returns a trollimage XRImage, the
    # data attribute is an xarray type. Use xarray transpose method to 
    # fix the axes; note we assume certain axis names.

    data_arr = get_enhanced_image(scn['true_color']).data
    data_arr_T = data_arr.transpose('y','x','bands')

    im = ax.imshow(data_arr_T, transform=crs, extent=crs.bounds,
                   origin='upper')

    return im


def nonworldview_overlay_plot(
        cfg_d, ovr_d, odat, out_plot_name=None,
        var_label=None, fignum=10, figsize=(20,20),
        img_xsize=2000):
    """
    Make an overlay plot - GOES/Himawari. 
    This is function to integrate with the vistool plot data flow.

    inputs (all 3 are defined by the corresponding vistool functions)
    cfg_d: configuration dictionary
    ovr_d: overlay information dictionary
    odat: oco2 data dictionary

    optional inputs:
    out_plot_name: string file path to create an output file. The default
        value of None means no output image file is creted.
    var_label: a string to set to the overlay data method (used for the
        colorbar label).
    fignum: matplotlib figure number.
    figsize: matplotlib figure size (width, height) in inches.
    img_xsize: this is used in combination with a vistool_lib function to
        get the number of pixels for the satpy regridded image. This input
        sets the x-size, and the y-size is computed based on the aspect
        ratio of the requested Lat/Lon box.

    returns a dictionary, containing all the axis objects created;
        this is generally most useful for testing or debugging.
        the dictionary contents are:
    image_ax: the MPL axis object containing the image
    cb_ax: the MPL axis object for the colorbar. This will be set to
        invisible if the colorbar was not needed.
    """

    dt = cfg_d['datetime']
    
    latlon_extent = [cfg_d['lon_ul'], cfg_d['lat_lr'],
                     cfg_d['lon_lr'], cfg_d['lat_ul']]
    center_lat = (latlon_extent[1] + latlon_extent[3]) / 2.0
    center_lon = (latlon_extent[0] + latlon_extent[2]) / 2.0

    # accessing the needed files depending on the geostation and files location;
    # downloading if needed
    # also getting the time offset of the background image
    if (cfg_d['sensor'].startswith('GOES')):
        file_list, time_offsets = get_ABI_files(
            dt, center_lat, cfg_d['sensor'], cfg_d['files_loc'], cfg_d['data_home'])
        time_offset = np.mean(time_offsets)/60.0
    else:
        file_list, time_offsets = get_AHI_files(
            dt, center_lat, cfg_d['files_loc'], cfg_d['data_home'])
        time_offset = np.mean(time_offsets)/60.0
        
    if len(file_list) == 0:
        raise ValueError('No files were found for requested date')

    # need to setup xpixels/ypixels here - the pixel size is input
    # to the Scene's area def.
    # note that Rob's method changes ypixels depending on the positions
    # of W,S,E,N., as computed by a Proj object.
    image_size = vl.get_image_size(latlon_extent, "epsg:3857",
                                   img_xsize=img_xsize)

    # getting the scene by the background files
    if (cfg_d['sensor'].startswith('GOES')):
        scn = get_scene_obj(file_list, latlon_extent, 'abi_l1b',
                            width=image_size[0], height=image_size[1],
                            resample_method=cfg_d['resample_method'])
    else:
        scn = get_scene_obj(file_list, latlon_extent, 'ahi_hsd',
                            width=image_size[0], height=image_size[1],
                            resample_method=cfg_d['resample_method'])
    crs = scn['true_color'].attrs['area'].to_cartopy_crs()

    # recompute the time offset, now that the files are loaded.
    # This probably will not change the number for ABI obs, since a
    # precise start/end time is in the filename, but will change
    # things slightly for AHI.
    if (cfg_d['sensor'].startswith('GOES')):
        domain = cfg_d['sensor'][-1]
    else:
        domain = 'F'
    time_offset = compute_time_offset(scn, center_lat, domain, dt)/60.0

    overlay_present = ovr_d is not None
    cbar_needed = overlay_present and ovr_d['cmap'] in plt.colormaps()

    fig, ax, _, cb_ax, layer_cb_ax, fig_scalefactor = vl.setup_axes(
        latlon_extent, crs, fignum=fignum,
        figsize=figsize)
    
    # plotting the retrieved scene
    im = plot_scene_obj(ax, scn)
    
    if overlay_present:
        if ovr_d['var_lims']:
            vmin, vmax = ovr_d['var_lims']
            vl.overlay_data(ax, cb_ax, odat, cmap=ovr_d['cmap'],
                            vmin=vmin, vmax=vmax,
                            var_label=var_label, alpha=ovr_d['alpha'],
                            fig_scalefactor=fig_scalefactor)
        else:
            vl.overlay_data(ax, cb_ax, odat, cmap=ovr_d['cmap'],
                            var_label=var_label, alpha=ovr_d['alpha'],
                            fig_scalefactor=fig_scalefactor)
    else:
        cb_ax.set_visible(False)

    # satpy version never uses the layer colorbar, so turn it off
    layer_cb_ax.set_visible(False)

    todays_date = datetime.datetime.now().strftime('%Y-%m-%d')
    
    # formatting the plot
    if cfg_d['out_plot_title'] == 'auto':
        title_string = vl.create_plot_title_string(cfg_d, ovr_d, odat)
        # add some specific information for Geo backgrounds:
        # view zenith of Geo, and the time offset.
        view_zenith = _get_view_zenith(scn, center_lon, center_lat)
        extra_title_info = (
            ' (Zenith = {:2.0f}$^\\circ$, '.format(view_zenith) +
            '$\\Delta$t = {:4.1f} min.)'.format(time_offset) )
        title_string_lines = title_string.split('\n')
        title_string_lines[0] = title_string_lines[0] + extra_title_info
        title_string = '\n'.join(title_string_lines)
    else:
        title_string = cfg_d['out_plot_title']

    title_size = int(np.round(30 * np.mean(figsize)/20.0))
    ax.set_title(title_string, size=title_size, y=1.01)

    if out_plot_name:
        fig.savefig(out_plot_name)
        print("\nFigure saved at "+ out_plot_name + "\n")

    ax_dict = dict(
        image_ax = ax, cb_ax = cb_ax,
    )

    return ax_dict
