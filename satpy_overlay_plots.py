"""
satpy_overlay_plots.py

module implementing satpy based plots with OCO-2 data overlays.

intended to integrate with oco_vistool.py
"""

import os, glob, datetime, collections, itertools
import tempfile, shutil
import bz2

import numpy as np
import pandas as pd

from satpy import Scene 
from satpy.writers import get_enhanced_image
import pyresample

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import shapely.geometry as sgeom

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatches

import netCDF4
import h5py

import boto3
import botocore

import fnmatch

from oco_vistool import read_shp

# to get the location of the city data (a subdir from where the
# code is located).
_code_dir = os.path.dirname(os.path.realpath(__file__))

def _approx_scanline_time(stime, etime, lat, domain):
    """
    Returns the scanline time, for a particular lat, and domain (C or F).
    see docstring in get_ABI_files for implementation notes.

    returns a single datetime object with the approx scanline time.

    if the lat is out of range, it returns the closest time; in other
    words, out of range will be limited to pick the stime or etime.

    inputs:
    stime: datetime object with start time
    etime: datetime object with end time
    lat: latitude (scalar)
    domain: (C)ONUS or (F)ull disk
    """

    if domain == 'C':
        slat = 50.0
        elat = 15.0
    elif domain == 'F':
        slat = 81.0
        elat = -81.0
    else:
        raise ValueError('invalid domain: '+domain)

    tfrac = (slat - lat) / (slat - elat)

    if tfrac < 0:
        scan_time = stime
    elif tfrac > 1:
        scan_time = etime
    else:
        total_time = (etime-stime).total_seconds()
        scan_time = stime + datetime.timedelta(seconds=total_time*tfrac)

    return scan_time


def get_loc_ABI_files(datetime_utc, data_home, domain, platform, hour_offsets, band_list, glob_fstr, verbose):
    """
    Helper function for getting the local filesystem GOES ABI files.
    
    Returns a list with relevant filenames.
    
    inputs:
    datetime_utc: datetime object for the needed 
    data_home: origin background files directory
    domain: (C)ONUS or (F)ull disk
    platform: GOES-East (G16) or GOES-West (G17)
    hour_offsets: different offsets to iterate
    band_list: list of needed bands
    glob_fstr: glob file format
    verbose: T/F - if details are needed
    """
    
    flists = collections.defaultdict(list)
    # iterating over all time options and adding the located files to the list
    for hour_offset, band in itertools.product(hour_offsets, band_list):
        dt = datetime_utc - datetime.timedelta(hours=hour_offset)

        if verbose:
            print('Searching for data at hour: '+str(dt))
        
        # constructing path to the needed files
        y_dir = dt.strftime('%Y')
        ymdd_dir = dt.strftime('%Y_%m_%d_%j')
        
        # four "?" for the MMSS, and then one more for a character at
        # the end, I am not sure about (frac second?)
        stimestamp = dt.strftime('%Y%j%H') + '????' + '?'
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

def get_aws_ABI_files(datetime_utc, data_home, domain, platform, hour_offsets, band_list, glob_fstr, verbose):
    """
    Helper function for getting the AWS S3 GOES ABI filenames.
    
    Returns a list with needed filenames and the AWS S3 bucket object (where the files will be downloaded).
    
    inputs:
    datetime_utc: datetime object for the needed point of time
    data_home: origin background files directory
    domain: (C)ONUS or (F)ull disk
    platform: GOES-East (G16) or GOES-West (G17)
    hour_offsets: different offsets to iterate over
    band_list: list of needed bands
    glob_fstr: glob file format
    verbose: T/F - if details are needed
    """
    
    # accessing the relevant S3 bucket
    aws_keys = pd.read_csv('Keys.csv', header = 0)
    s3 = boto3.resource('s3', aws_access_key_id = aws_keys['aws_access_key_id'].values[0], aws_secret_access_key = aws_keys['aws_secret_access_key'].values[0])
    
    if (platform[-2:] == '16'):
        g_bucket = s3.Bucket('noaa-goes16')
    else:
        g_bucket = s3.Bucket('noaa-goes17')
        
    flists = collections.defaultdict(list)
    
    # iterating over the different hours to locate the needed files
    for hour_offset, band in itertools.product(hour_offsets, band_list):
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

def download_aws_ABI_files(flist, g_bucket, platform):
    """
    Helper function for downloading the AWS GOES files given their names.

    Returns a list with downloaded local GOES ABI filenames.
    
    inputs:
    flist: list of needed files
    g_bucket: AWS S3 bucket to download the files from
    platform: GOES-East (G16) or GOES-West (G17)
    """
    
    files = list()
    
    # iterating to and downloading the needed files from the S3 bucket (by the known path)
    for needed_path in flist:
        g_files = g_bucket.objects.filter(Prefix=needed_path)
        for g_file in g_files:
            if not (os.path.isfile(platform + '/' + g_file.key)):
                if not os.path.isdir(os.path.dirname(platform + '/' + g_file.key)):
                    os.makedirs(os.path.dirname(platform + '/' + g_file.key))
                g_bucket.download_file(g_file.key, platform + '/' + g_file.key)
                files.extend(glob.glob(platform + '/' + g_file.key))
            else:
                files.extend(glob.glob(platform + '/' + g_file.key))

    flist = list(set(files))
    return flist

def get_ABI_files(datetime_utc, center_lat,
                  sensor, files_loc, data_home = None, verbose=False):
    """
    Using the helpers above, accesses the needed GOES ABI files and downloads them (if needed)
    
    Returns a list with accessed GOES filenames and time offsets.
    
    inputs:
    datetime_utc: datetime object for the needed point of time
    center_lat: center latitude of the background image
    sensor: the geostation, device, and domain (GOES16_ABI_C, for example)
    files_loc: if the background images are local or AWS
    data_home: data origin for the background files if files are stored locally (by format: data_home/YYYY/...)
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
        flists, g_bucket = get_aws_ABI_files(datetime_utc, data_home, domain, platform, hour_offsets, band_list, glob_fstr, verbose)
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
            stimestamp, etimestamp = flists[band][n].split('_')[-3:-1]
            # [1:-1] slice removes 's' or 'e' and the extra number (frac seconds?)
            stime = datetime.datetime.strptime(stimestamp[1:-1], '%Y%j%H%M%S')
            etime = datetime.datetime.strptime(etimestamp[1:-1], '%Y%j%H%M%S')
            scan_time = _approx_scanline_time(stime, etime, center_lat, domain)
            time_offsets_all.append((datetime_utc-scan_time).total_seconds())
        # find minimum time difference
        k = np.argmin(np.abs(time_offsets_all))
        flist.append(flists[band][k])
        time_offset.append(time_offsets_all[k])
    
    # download the files if needed
    if (files_loc == 'aws'):
        flist = download_aws_ABI_files(flist, g_bucket, platform)

    return flist, time_offset

def get_AHI_files(datetime_utc, data_home = None):
    """
    Access the needed Himawari files by date & time and downloads them (if needed).
    
    Returns a list with needed glob-format AHI Himawari files.
    
    inputs:
    datetime_utc: datetime object for the needed point of time
    data_home: data origin for the background files if files are stored locally (by format: data_home/YYYY/...)
    """
    
    # Constructing the files path by the format
    year = str(datetime_utc.year)
    month = str(datetime_utc.month).zfill(2)
    day_m = str(datetime_utc.day).zfill(2)
    day_y = str(datetime_utc.timetuple().tm_yday).zfill(3)
    hour = str(datetime_utc.hour).zfill(2)
    orig_minutes = datetime_utc.minute
    seconds = datetime_utc.second
    minutes = str(int(round(orig_minutes + seconds/60, -1))%60).zfill(2)
    
    bands_list = [1, 2, 3, 4]
    files = list()
    
    # accessing the files locally
    if (data_home is not None):
        for band in bands_list:
            files.extend(glob.glob(data_home + "/" + year + "/" + year + "_" + month + "_" + day_m + "_" + day_y + "/" + hour + 
                                   minutes + "/HS_H08_" + year + month + day_m + "_" + hour + minutes +"_B" + str(band).zfill(2) +
                                  "_FLDK_*.DAT"))  
    # accessing the files and downloading them from AWS 
    else:
        # connecting to the S3 Himawari bucket
        aws_keys = pd.read_csv('Keys.csv', header = 0)
        s3 = boto3.resource('s3', aws_access_key_id = aws_keys['aws_access_key_id'].values[0], aws_secret_access_key = aws_keys['aws_secret_access_key'].values[0])
        hima_bucket = s3.Bucket('noaa-himawari8')
        hima_path = 'AHI-L1b-FLDK/' + year + '/' + month + '/' + day_m +'/' + hour + minutes + '/'
        hima_files = hima_bucket.objects.filter(Prefix=hima_path)
        
        # downloading the needed Himawari files by known AWS paths
        for hima_file in hima_files:
            for band in bands_list:
                if '_B0' + str(band) + '_' in hima_file.key:
                    if not (os.path.isfile('Himawari-08/' + hima_file.key[:-4]) or os.path.isfile('Himawari-08/' + hima_file.key)):
                        if not os.path.isdir(os.path.dirname('Himawari-08/' + hima_file.key)):
                            os.makedirs(os.path.dirname('Himawari-08/' + hima_file.key))
                        hima_bucket.download_file(hima_file.key, 'Himawari-08/' + hima_file.key)
                    else:
                        files.extend(glob.glob('Himawari-08/' + hima_file.key[:-4]))
         
        # decompressing the downloaded files
        for filename in glob.glob('Himawari-08/' + hima_path + '*.bz2'):
            zipfile = bz2.BZ2File(filename)
            data = zipfile.read()
            newfilename = filename[:-4] 
            open(newfilename, 'wb').write(data)
            os.remove(filename)
            files.extend(glob.glob(newfilename))
        files = list(set(files))
    return files

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
        'testC', {'proj':'eqc'}, width=width, height=height,
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


def setup_axes(fignum, crs, figsize=(10,8), create_colorbar_axis=True):
    """
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


def overlay_data(ax, cb_ax, odata, var_label=None, **kw):
    """
    overlay data onto image with vertex lat/lon points.

    inputs:
    ax: axis object where polygons will be drawn
    cb_ax: the colorbar axis object. Note this is only used if
       the cmap is a colormap.
    odata: oco2 overlay data dictionary

    extra kw: cmap (string colormap or color name), vmin, vmax, alpha

    returns:
    None

    method based on
    https://scitools.org.uk/cartopy/docs/latest/gallery/hurricane_katrina.html
    #sphx-glr-gallery-hurricane-katrina-py
    
    the add_geometries() method would accept a list of polygons, but then
    there is no obvious way to color them according to the colorbar - 
    so, instead do them one at a time.
    """

    n_vals = odata['lat'].shape[0]

    if 'cmap' in kw:
        if kw['cmap'] in plt.colormaps():
            cmap_name = kw['cmap']
            use_cmap = True
        elif kw['cmap'] in mpl.colors.cnames.keys():
            color_name = kw['cmap']
            use_cmap = False
        else:
            print(kw['cmap'] + " is not a recognized color or colormap. "+
                  "Data will be displayed in red")
            color_name = 'red'
            use_cmap = False
    else:
        cmap_name = 'jet'
        use_cmap = True

    if 'alpha' in kw:
        alpha = kw['alpha']
    else:
        alpha = 1.0

    if use_cmap:
        if 'vmin' in kw:
            vmin = kw['vmin']
        else:
            vmin = odata['var_data'].min()
        if 'vmax' in kw:
            vmax = kw['vmax']
        else:
            vmax = odata['var_data'].max()

        C_func = mpl.cm.get_cmap(cmap_name)
        N_func = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

    # vertices: loop over footprints and use add_geometries()
    if odata['lat'].ndim == 2:

        for n in range(n_vals):
            poly_pts = np.ma.stack([odata['lon'][n,:], odata['lat'][n,:]]).T
            polygon = sgeom.Polygon(poly_pts)

            if use_cmap:
                facecolor = C_func(N_func(odata['var_data'][n]))
            else:
                facecolor = color_name
            ax.add_geometries([polygon], ccrs.PlateCarree(),
                              facecolor=facecolor, alpha=alpha,
                              edgecolor='none')

    # footprint centers, use scatter()
    else:
        if use_cmap:
            ax.scatter(odata['lon'], odata['lat'], c=odata['var_data'],
                       cmap=cmap_name, vmin=vmin, vmax=vmax,
                       s=9, edgecolor='none', transform=ccrs.PlateCarree())
        else:
            ax.scatter(odata['lon'], odata['lat'], c=color_name,
                       s=9, edgecolor='none', transform=ccrs.PlateCarree())

    if use_cmap:
        cb = mpl.colorbar.ColorbarBase(
            cb_ax, cmap=C_func, orientation = 'vertical', norm = N_func)
        if var_label:
            cb_lab = cb.ax.set_xlabel(var_label, labelpad=8,
                                      fontweight='bold')
            cb_lab.set_fontsize(10)
            cb.ax.xaxis.set_label_position("top")
        for t in cb.ax.yaxis.get_ticklabels():
            t.set_weight("bold")
            t.set_fontsize(12)


def nonworldview_overlay_plot(cfg_d, ovr_d, odat, out_plot_name=None,
                          var_label=None, fignum=10):
    """
    Make an overlay plot - GOES/Himawaei. 
    This is function to integrate with the vistool plot data flow.

    inputs (all 3 are defined by the corresponding vistool functions)
    cfg_d: configuration dictionary
    ovr_d: overlay information dictionary
    odat: oco2 data dictionary

    returns a dictionary, containing all the axis objects created;
        this is generally most useful for testing or debugging.
        the dictionary contents are:
    gridspec: the gridspec object
    image_ax: the MPL axis object containing the image
    cb_ax: the MPL axis object for the colorbar; None if the colorbar
         was not created
    anno_axes: a list containing the 4 dummy axes where the lat/lon corner
         values are displayed
    inset_ax: the MPL axis object with the inset image
    """

    # if there is overlay data present, use the OCO2 obs time
    # to get the data.
    # otherwise, use the datetime object
    if ovr_d:
        dt = datetime.datetime.utcfromtimestamp(np.mean(odat['time']))
    else:
        dt = cfg_d['datetime']

    center_lat = (cfg_d['lat_lr'] + cfg_d['lat_ul']) / 2.0
    
    # getting the needed files depending on the station and files location; downloading if needed
    if (cfg_d['sensor'].startswith('GOES')):
        if (cfg_d['files_loc'] == 'local'):
            file_list, time_offsets = get_ABI_files(
                dt, center_lat, cfg_d['sensor'], cfg_d['files_loc'], cfg_d['data_home'])
        else:
            file_list, time_offsets = get_ABI_files(
                dt, center_lat, cfg_d['sensor'], cfg_d['files_loc'])
        mean_time_offset = np.mean(time_offsets)/60.0
    else:
        if (cfg_d['files_loc'] == 'local'):
            file_list = get_AHI_files(
                dt, cfg_d['data_home'])
        else:
            file_list = get_AHI_files(dt)
    if len(file_list) == 0:
        raise ValueError('No files were found for requested date')

    # convert the LL box corners (in degrees LL) to an extent box
    # [min_lon, min_lat, max_lon, max_lat]
    latlon_extent = [cfg_d['lon_ul'], cfg_d['lat_lr'], 
                     cfg_d['lon_lr'], cfg_d['lat_ul']]
    
    # getting the scene by the background files
    if (cfg_d['sensor'].startswith('GOES')):
        scn = get_scene_obj(file_list, latlon_extent, 'abi_l1b',
                            resample_method=cfg_d['resample_method'])
    else:
        scn = get_scene_obj(file_list, latlon_extent, 'ahi_hsd',
                            resample_method=cfg_d['resample_method'])
    crs = scn['true_color'].attrs['area'].to_cartopy_crs()

    overlay_present = ovr_d is not None
    cbar_needed = overlay_present and ovr_d['cmap'] in plt.colormaps()

    gs, ax, cb_ax, inset_ax = setup_axes(
        fignum, crs, create_colorbar_axis=cbar_needed)
    
    # plotting the retrieved scene
    im = plot_scene_obj(ax, scn)
    
    if overlay_present:
        if ovr_d['var_lims']:
            vmin, vmax = ovr_d['var_lims']
            overlay_data(ax, cb_ax, odat, cmap=ovr_d['cmap'],
                         vmin=vmin, vmax=vmax, var_label=var_label,
                         alpha=ovr_d['alpha'])
        else:
            overlay_data(ax, cb_ax, odat, cmap=ovr_d['cmap'],
                         var_label=var_label, alpha=ovr_d['alpha'])

    make_inset_map(inset_ax, cfg_d['geo_upper_left'], cfg_d['geo_lower_right'])
    todays_date = datetime.datetime.now().strftime('%Y-%m-%d')
    
    # formatting the plot
    if cfg_d['out_plot_title'] == 'auto':
        if cfg_d['sensor'].startswith('GOES'):
            if ovr_d:
                title_string = (
                    'Overlay data from {0:s}' +
                    '\nBackground from {1:s}, '+
                    '\nOverlay time = {2:s},   '+
                    'mean time offset = {3:4.1f} min.,  '+
                    'plot created on {4:s}' )
                title_string = title_string.format(
                    os.path.split(ovr_d['var_file'])[1],
                    os.path.split(file_list[1])[1], 
                    dt.strftime('%Y-%m-%d %H:%M:%S'), mean_time_offset,
                    todays_date)
            else:
                title_string = (
                    'Background from  {0:s}' + 
                    '\n Request time = {1:s},   '+
                    'mean time offset = {2:4.1f} min.,  '+
                    'plot created on {3:s}' )
                title_string = title_string.format(
                    os.path.split(file_list[1])[1],
                    dt.strftime('%Y-%m-%d %H:%M:%S'),
                    mean_time_offset, todays_date)
        else:
            if ovr_d:
                title_string = (
                    'Overlay data from {0:s}' +
                    '\nBackground from {1:s}, '+
                    '\nOverlay time = {2:s},   '+
                    'plot created on {3:s}' )
                title_string = title_string.format(
                    os.path.split(ovr_d['var_file'])[1],
                    os.path.split(file_list[1])[1], 
                    dt.strftime('%Y-%m-%d %H:%M:%S'),
                    todays_date)
            else:
                title_string = (
                    'Background from  {0:s}' + 
                    '\n Request time = {1:s},   '+
                    'plot created on {2:s}' )
                title_string = title_string.format(
                    os.path.split(file_list[1])[1],
                    dt.strftime('%Y-%m-%d %H:%M:%S'), todays_date)
            
    else:
        title_string = cfg_d['out_plot_title']


    ax.set_title(title_string, size='x-small')
    ax.coastlines(resolution='10m', color='lightgray') 
    anno_axes = annotate_ll_axes(
        gs, cfg_d['geo_upper_left'], cfg_d['geo_lower_right'])

    annotate_locations(ax, cfg_d)

    if out_plot_name:
        fig = plt.figure(fignum)
        fig.savefig(out_plot_name, dpi=150)

    ax_dict = dict(
        gridspec = gs, image_ax = ax, cb_ax = cb_ax,
        anno_axes = anno_axes, inset_ax = inset_ax)

    return ax_dict
