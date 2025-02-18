"""
example driver functions to call vistool from other python code.
"""

import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["DASK_NUM_WORKERS"] = "1"
os.environ["DASK_ARRAY__CHUNK_SIZE"] = "16MiB"

from datetime import datetime

import numpy as np

import calendar

from geo_imager_visibility import determine_optimal_geo_satellite
from oco_vistool import load_OCO2_Lite_overlay_data
from satpy_overlay_plots import nonworldview_overlay_plot

def make_geo_image(obs_datetime, latlon_ul, latlon_lr,
                   orbit, target_id, data_rev, out_dir_temp,
                   download_dir = '.', L2_Lite_file=None):
    """
    make an image using true color geostationary data as the background.

    Parameters
    ----------
    obs_datetime : datetime
        Time used to locate Geo imagery. This should be specified to the minute.
        This should be the time of the OCO-3 data collection for this SAM/Target.
    latlon_ul : array-like
        two elements, containing degrees latitude and longitude (in that order) of
        the upper left (e.g., NorthWest) corner of the spatial box to display.
    latlon_lr : array-like
        two elements, containing degrees latitude and longitude (in that order) of
        the lower right (e.g., SouthEast) corner of the spatial box to display.
    orbit : int
        orbit number of associated OCO-3 data. This is used to locate data within
        the Lite file (which can be important for OCO-3 for the uncommon cases
        where there are self-crossings) and for the output filename.
    target_id : str
        the target_id string, e.g. fossilNNNN. This is used only for the output
        plot filename
    data_rev : str
        description string for the data revision. This is used only for the output
        plot filename
    download_dir : str
        directory to use for locally downloaded geo files. By default this is
        the current directory ('.')
    L2_Lite_file : None or str
        If none, then a geostationary image is created from just the datetime and
        latlon/box, and no overlay is drawn. If set to a string, this should be
        a L2 Lite CO2 file containing data to display in an overlay.

    Returns
    -------
    output_plot_file : str
        the string name of the created PNG file.
        If the requested time and latlon box is not contained within any geo
        imager field of regard, then the empty string is returned.

    Example filenames:
    GOESWestC_Lite_B10400Br_r02_none_20200413_5354_cal001.png
    OCO3-GOESWestC_Lite_B10400Br_r02_xco2_bc_qf_20200413_5354_cal001.png

    """
    zenith_limit = 85.0

    lons = [latlon_ul[1], latlon_ul[1], latlon_lr[1], latlon_lr[1]]
    lats = [latlon_ul[0], latlon_lr[0], latlon_lr[0], latlon_ul[0]]

    selected_sensor = determine_optimal_geo_satellite(
        lons, lats, obs_datetime, zenith_limit = zenith_limit )

    # if None is returned, this position is not viewable from the geo imagers.
    if selected_sensor is None:
        return ''

    output_plot_file_fstr = (
        '{sensor:s}{geo_sensor:s}_{data_rev:s}_{var_name:s}_'+
        '{ymd:s}_{orbit:s}_{target_id:s}.png' )

    # a mapping between the geostationary sensor name (needed for satpy_overlay_plots
    # function) and a string token for use in the filename prefix.
    geo_name_mapping = {
        'GOES16_ABI_F':'GOESEastFD', 'GOES16_ABI_C':'GOESEastC',
        'GOES17_ABI_F':'GOESWestFD', 'GOES17_ABI_C':'GOESWestC',
        'GOES18_ABI_F':'GOESWestFD', 'GOES18_ABI_C':'GOESWestC',
        'Himawari-08' : 'Himawari',
        'Himawari-09' : 'Himawari'}
    geo_sensor_prefix = geo_name_mapping[selected_sensor]


    ########
    # Input configuration dictionaries for oco_vistool

    # main configuration: no changes are needed here
    cfg_d = dict(
        out_plot_title = 'auto',
        datetime = obs_datetime,
        lat_ul = latlon_ul[0],
        lat_lr = latlon_lr[0],
        lon_ul = latlon_ul[1],
        lon_lr = latlon_lr[1],
        sensor = selected_sensor,
        files_loc = 'aws',
        data_home = download_dir,
        resample_method = 'native_bilinear')

    # overlay configuration: some changes are possible, to control the way the
    # overlay is drawn on the figure.
    ovr_d = dict(

        # important settings, should not be changed
        sensor = 'OCO-3',
        var_file = L2_Lite_file,
        var_name = 'xco2',
        var_title_string = 'Bias Corrected and Quality Flagged '+r'$X_{CO_2}$',
        lat_name = 'vertex_latitude',
        lon_name = 'vertex_longitude',
        lite_quality = 'good',
        sif_or_co2 = 'CO2',
        lat_ul = latlon_ul[0],
        lat_lr = latlon_lr[0],
        lon_ul = latlon_ul[1],
        lon_lr = latlon_lr[1],

        # important, but could be changed
        cmap = 'viridis',
        alpha = 1.0,

        # possibly unneeded settings
        var_name_only = 'xco2',
        orbit = orbit,

        # unimportant settings, that do not impact the output, but
        # but need to be present with these placeholder values.
        footprint_lims = [1, 8],
        make_background_image = False,
        lat_shift = 0.0,
        lon_shift = 0.0,
        warn = None,
    )

    # end configuration dictionaries
    ########

    # set input overlay data - if no L2_Lite is input, these
    # pass through with None placeholder values. Note some variables
    # change here, that impact the output plot filename.
    if L2_Lite_file:
        odat = load_OCO2_Lite_overlay_data(ovr_d)
        # note that the display  variable limits are set here - this is
        # the 90pct method used in old geo imagery.
        if odat['var_data'].shape[0] == 0:
            print('No overlay data found, not producing image')
            return ''
        ovr_d['var_lims'] = np.percentile(
            odat['var_data'], [10.0, 90.0], interpolation='nearest').tolist()
        sensor = 'OCO3-'
        var_name = 'xco2_bc_qf'
        cbar_name = r'$X_{CO_2}$'+' [ppm]'
    else:
        odat = None
        ovr_d = None
        sensor = ''
        var_name = 'none'
        cbar_name = 'none'

    # build output filename based on parameter values extracted above.
    output_plot_file = output_plot_file_fstr.format(
        sensor = sensor, var_name = var_name,
        geo_sensor = geo_sensor_prefix, data_rev = data_rev,
        ymd = obs_datetime.strftime('%y%m%d'),
        orbit = str(orbit), target_id = target_id)

    try:
        objs = nonworldview_overlay_plot(
            cfg_d, ovr_d, odat, var_label = cbar_name, out_plot_name = None)
    except ValueError:
        print('Image file failed, possibly due to geo imager data gap')
        return ''

    # Time/creation stamp
    objs['image_ax'].text(0.99,0.01,"Created "+str(datetime.now().day)+' '+calendar.month_abbr[datetime.now().month]+' '+str(datetime.now().year), ha='right', va='bottom', transform=objs['image_ax'].transAxes, color='1.0',size=18)

    # could be altered here
    objs['fig'].savefig(out_dir_temp+"/"+output_plot_file) #Give the full path

    return output_plot_file


def sample_goes_run():

    #Specify the full path to the OCO-3 Lite file
    L2_Lite_file = '/data/oco3/scf/product/Lite_B10400Br_r02/2023/10/01/LtCO2/oco3_LtCO2_231001_B10401Br_231128044453s.nc4'

    # info here was manually extracted by looking at:
    # https://ocov3.jpl.nasa.gov/sams/plots.php?sID=35569
    obs_datetime = datetime(2023, 10, 1, 12, 00)
    lat0, lon0 = -34.758, -58.523
    latlon_ul = (lat0 + 1.5, lon0 - 1.5/np.cos(np.deg2rad(lat0)))
    latlon_lr = (lat0 - 1.5, lon0 + 1.5/np.cos(np.deg2rad(lat0)))
    orbit = 24959
    target_id = 'fossil0035'
    data_rev = 'B10401Br_r02'
    download_dir = './tmp'
    out_dir = './' #Now specifying the output directory

    t0 = datetime.now()
    output_plot_file = make_geo_image(
        obs_datetime, latlon_ul, latlon_lr, orbit, target_id, data_rev, out_dir,
        download_dir = download_dir, L2_Lite_file=None)
    print('Made background image : ' + output_plot_file)
    print('Elapsed time: ' + str(datetime.now()-t0))

    t0 = datetime.now()
    output_plot_file = make_geo_image(
        obs_datetime, latlon_ul, latlon_lr, orbit, target_id, data_rev, out_dir,
        download_dir = download_dir, L2_Lite_file=L2_Lite_file)
    print('Made overlay image    : ' + output_plot_file)
    print('Elapsed time: ' + str(datetime.now()-t0))

def sample_himawari_run():

    #Specify the full path to the OCO-3 Lite file
    L2_Lite_file = '/data/oco3/scf/product/Lite_B10400Br_r02/2023/10/01/LtCO2/oco3_LtCO2_231001_B10401Br_231128044453s.nc4'

    # info here was manually extracted by looking at:
    # https://ocov3.jpl.nasa.gov/sams/plots.php?sID=35557
    obs_datetime = datetime(2023, 10, 1, 4, 41)
    lat0, lon0 = 35.8, 119.8
    latlon_ul = (lat0 + 1.5, lon0 - 1.5/np.cos(np.deg2rad(lat0)))
    latlon_lr = (lat0 - 1.5, lon0 + 1.5/np.cos(np.deg2rad(lat0)))
    orbit = 24954
    target_id = 'fossil0078'
    data_rev = 'B10401Br_r02'
    download_dir = './tmp'
    out_dir = './' #Now specifying the output directory

    t0 = datetime.now()
    output_plot_file = make_geo_image(
        obs_datetime, latlon_ul, latlon_lr, orbit, target_id, data_rev, out_dir,
        download_dir = download_dir, L2_Lite_file=None)
    print('Made background image : ' + output_plot_file)
    print('Elapsed time: ' + str(datetime.now()-t0))

    t0 = datetime.now()
    output_plot_file = make_geo_image(
        obs_datetime, latlon_ul, latlon_lr, orbit, target_id, data_rev, out_dir,
        download_dir = download_dir, L2_Lite_file=L2_Lite_file)
    print('Made overlay image    : ' + output_plot_file)
    print('Elapsed time: ' + str(datetime.now()-t0))


if __name__ == "__main__":
    sample_goes_run()
    sample_himawari_run()
