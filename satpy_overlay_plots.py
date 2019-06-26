"""
satpy_overlay_plots.py

module implementing satpy based plots with OCO-2 data overlays.

intended to integrate with oco2_modis_vistool.py
"""

import os.path, glob, datetime, collections, itertools
import tempfile, shutil

import numpy as np

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

from oco2_modis_vistool import read_shp

# to get the location of the city data (a subdir from where the
# code is located).
_code_dir = os.path.dirname(os.path.realpath(__file__))

def get_ABI_files(datetime_utc, data_home, domain='C', verbose=False):
    """
    find the file list for an ABI domain, given an input date.
    this will find the images (for bands 1,2,3) closest to the UTC
    date time. If there are no files within 1 hour, throws exception.

    inputs:
    datetime_utc: a python datetime object (in UTC) with the requested time.
    data_home: string containing path to local ABI file archive
    domain: "C" or "F", to use either the CONUS or FD (Full Disk) ABI
        imagery. Defaults to CONUS.
    verbose: set to True to print more information to console
    
    returns:
    file_list: list of strings containing full paths to the located files.
        should contain 3 (for each of bands 1,2,3)
    time_offsets: the time offsets (in seconds) between the input 
        datetime_utc and the start times of the ABI images.

    Implementation notes:
    the data_home is assumed to contain a directory tree with the following
    layout, containing the ABI netCDF4 files:
    data_home
      |--YYYY
           |--YYYY_MM_DD_DOY
              |--abi
                 |--L1b
                    |--RadC
                       |--OR_ABI-L1b-RadC-M3C01_G16_sYYYYDOY*.nc
                       |-- ....
                       |--OR_ABI-L1b-RadC-M3C01_G16_sYYYYDOY*.nc
                          <etc for other bands>
                    |--RadF
                       |--OR_ABI-L1b-RadF-M3C01_G16_sYYYYDOY*.nc
                       |-- ....
                       |--OR_ABI-L1b-RadF-M3C01_G16_sYYYYDOY*.nc
                          <etc for other bands>
                    |--RadM1
                    |--RadM2
           |-- ....
           |--YYYY_MM_DD_DOY

    for example, the first CONUS Band 1 file for 2018:

    data_home/2018/2018_01_01_001/abi/L1b/RadC/
        OR_ABI-L1b-RadC-M3C01_G16_s20180010002196_e20180010004569_c20180010005015.nc

    """

    valid_domains = 'C', 'F'
    if domain not in valid_domains:
        raise ValueError('domain must be in '+str(valid_domains))

    # search within hour of input datetime, and +/-1 hour
    flists = collections.defaultdict(list)
    band_list = (1,2,3)
    hour_offsets = (-1,0,1)

    for hour_offset, band in itertools.product(hour_offsets, band_list):

        dt = datetime_utc - datetime.timedelta(hours=hour_offset)

        if verbose:
            print('Searching for data at hour: '+str(dt))
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

        # makes a glob string to match by s-time (I think, start time.)
        glob_fstr = ('OR_ABI-L1b-Rad' + domain +
                    '-M?C{0:02d}_G16_s{1:s}_e*_c*.nc')
        glob_str = os.path.join(ddir,glob_fstr.format(band, stimestamp))
        if verbose:
            print('glob: '+glob_str)

        located_files = glob.glob(glob_str)

        flists[band] += located_files

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
            stimestamp = flists[band][n].split('_')[-3]
            # removes 's' and the extra number (frac seconds?)
            stimestamp = stimestamp[1:-1]
            file_stime = datetime.datetime.strptime(stimestamp, '%Y%j%H%M%S')
            time_offsets_all.append((datetime_utc-file_stime).total_seconds())
        # find minimum
        k = np.argmin(np.abs(time_offsets_all))
        flist.append(flists[band][k])
        time_offset.append(time_offsets_all[k])
        
    return flist, time_offset


def get_scene_obj(file_list, lat, lon, width_km, height_km,
                  pixels=500, tmp_cache=False, 
                  resample_method='native_bilinear'):
    """Get Scene object, apply the resample area, to a small box 
    centered on the lat lon point.

    inputs:
    file_list: list of netCDF L1b ABI files, must contain bands 1,2,3
    lat: scalar, latitude in deg
    lon: scalar, longitude in deg
    width_km: scalar, size of image (e.g., -N to +N km in longitude/width)
    height_km: scalar, size of image (e.g., -N to +N km in latitude/height)

    tmp_cache: optional keyword, set to True to copy the located files
    to a temporary dir (from tempfile) before loading.
    Note the temporary files are NOT cleaned up.
    Use this option if the data_home access is slow or limited in some
    way, which can be mitigated by copying to a temp file on the local
    filesystem

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
        scn = Scene(reader='abi_l1b', filenames=cached_file_list)
    else:
        scn = Scene(reader='abi_l1b', filenames=file_list)

    scn.load(['true_color'])

    proj_dict = {'proj': 'eqc', 
                 'lat_0': lat, 'lon_0': lon,
                 'a': 6371228.0, 'units': 'm'}

    # the resample's area_extent needs to be in meters (I guess,
    # matching the units in the proj_dict.)
    width_m = width_km * 1e3
    height_m = height_km * 1e3
    my_area = pyresample.AreaDefinition.from_extent(
        'testC', proj_dict, (pixels, pixels),
        (-width_m, -height_m, width_m, height_m))

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


def setup_axes(fignum, crs, figsize=(6,6), create_colorbar_axis=True):
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

    ax = plt.subplot(gs[0:-1,3:-2], projection=crs)
    if create_colorbar_axis:
        cb_ax = plt.subplot(gs[0:-1, -1])
    else:
        cb_ax = None

    inset_ax = plt.subplot(gs[7:9, 0:3], projection=ccrs.PlateCarree())
    
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

    #ax_minlat.set_title(str("%.1f" % img_lr[0]),
    #                    horizontalalignment='left',
    #                    verticalalignment='bottom', **fontkw)
    #ax_maxlat.set_title(str("%.1f" % img_ul[0]),
    #                    horizontalalignment='left',
    #                    verticalalignment='top', **fontkw)
    #ax_minlon.set_title(str("%.1f" % img_ul[1]),
    #                    horizontalalignment='center',
    #                    verticalalignment='top', **fontkw)
    #ax_maxlon.set_title(str("%.1f" % img_lr[1]),
    #                    horizontalalignment='right',
    #                    verticalalignment='top', **fontkw)
                    
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
                       s=2, edgecolor='none', transform=ccrs.PlateCarree())
        else:
            ax.scatter(odata['lon'], odata['lat'], c=color_name,
                       s=2, edgecolor='none', transform=ccrs.PlateCarree())

    if use_cmap:
        cb = mpl.colorbar.ColorbarBase(
            cb_ax, cmap=C_func, orientation = 'vertical', norm = N_func)
        if var_label:
            cb_lab = cb.ax.set_xlabel(var_label, labelpad=8,
                                      fontweight='bold')
            cb_lab.set_fontsize(14)
            cb.ax.xaxis.set_label_position("top")
        for t in cb.ax.yaxis.get_ticklabels():
            t.set_weight("bold")
            t.set_fontsize(12)


def GOES_ABI_overlay_plot(cfg_d, ovr_d, odat, out_plot_name=None,
                          var_label=None, fignum=10):
    """
    make a GOES ABI overlay plot. This is function to integrate with the
    vistool plot data flow.

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
        dt = datetime.datetime.fromtimestamp(np.mean(odat['time']))
    else:
        dt = cfg_d['datetime']
    
    file_list, time_offsets = get_ABI_files(dt, cfg_d['data_home'])

    # convert the LL box corners (in degrees LL) to a center LL (degrees)
    # and the half-height/width (in km) needed for the satpy scene obj.
    center_lat = (cfg_d['lat_ul'] + cfg_d['lat_lr'])/2.0
    center_lon = (cfg_d['lon_ul'] + cfg_d['lon_lr'])/2.0
    d_lat = np.abs(cfg_d['lat_ul'] - cfg_d['lat_lr'])/2.0
    d_lon = np.abs(cfg_d['lon_ul'] - cfg_d['lon_lr'])/2.0
    km_per_deg = 2 * np.pi * 6378.0 / 360.0
    height_km = d_lat * km_per_deg
    width_km = d_lon * km_per_deg * np.cos(np.deg2rad(center_lat))

    scn = get_scene_obj(file_list, center_lat, center_lon,
                        width_km, height_km,
                        resample_method=cfg_d['resample_method'])
    crs = scn['true_color'].attrs['area'].to_cartopy_crs()

    overlay_present = ovr_d is not None
    cbar_needed = overlay_present and ovr_d['cmap'] in plt.colormaps()

    gs, ax, cb_ax, inset_ax = setup_axes(
        fignum, crs, create_colorbar_axis=cbar_needed)

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
    mean_time_offset = np.mean(time_offsets)/60.0
    todays_date = datetime.datetime.now().strftime('%Y-%m-%d')
    ax.set_title('Overlay data from '+
                 os.path.split(ovr_d['var_file'])[1] +
                 '\n background from GOES16-ABI, mean time offset = ' +
                 '{0:4.1f}  min'.format(mean_time_offset) +
                 '\nplot created on ' + todays_date,
                 size='x-small')

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
