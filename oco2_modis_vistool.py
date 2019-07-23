"""
oco2_modis_vistool.py

The purpose of this script is to pull Aqua-MODIS RGB images from Worldview
using the NASA GIBS API and overlay various OCO-2 data fields for case study
analysis in support of OCO-2 cloud and aerosol screening validation.

GIBS developer documentation:https://wiki.earthdata.nasa.gov/display/GIBS/GIBS+API+for+Developers

Minimum command line call:
python oco2_modis_vistool.py

Output:
Image (.png) named variable_region_date_quality_warn.png by default in specified output directory (default: code directory)

Requirements:

1) json configuration file (default: oco2_modis_vistool_config.json in code directory)

Corresponding author:
Heather Cronk <heather.cronk@colostate.edu>

"""

import warnings
warnings.filterwarnings("ignore")

import os
import sys
from glob import glob
from six import string_types
import collections
import datetime

from OCO2FileOps import *
import h5py

import numpy as np
import math
import pandas as pd

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import shapefile
from shapely.geometry import LineString, Point, Polygon
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec

from osgeo import gdal, osr
import shapely.geometry as sgeom

import xml.etree.ElementTree as ET
import json
import argparse

import re

if sys.version_info < (3,):
    import future
    from builtins import range

class ConfigFile:

    def __init__(self, json_file):
        self.cf = json_file

    def exists(self):
        return glob(self.cf)

    def get_contents(self):
        return(json.load(open(self.cf)))

def update_GIBS_xml(date, xml_file):

    """
    Puts the date of interest into the GIBS XML file
    """

    tree = ET.parse(xml_file)
    root = tree.getroot()

    url = root[0][0].text

    old_date = re.split('/', url)[6]

    new_url = re.sub(old_date, date, url)

    root[0][0].text = new_url
    tree.write(xml_file)

def pull_Aqua_RGB_GIBS(lat_ul, lon_ul, lat_lr, lon_lr, xml_file, tif_file, xsize=1200, ysize=1000):

    """
    Pulls the Aqua RGB imagery from WorldView using GIBS and puts it in specified tif file with associated metadata
    """
    gdal_path = os.popen("which gdal_translate").read().strip()
    cmd = gdal_path + " -of GTiff -outsize "+str(xsize)+" "+str(ysize)+" -projwin "+str(lon_ul)+" "+str(lat_ul)+" "+str(lon_lr)+" "+str(lat_lr)+" "+xml_file+" "+tif_file
    os.system(cmd)

def read_shp(filename):
    """Read shapefile to dataframe w/ geometry.

    Args:
        filename: ESRI shapefile name to be read  (without .shp extension)

    Returns:
        pandas DataFrame with column geometry, containing individual shapely
        Geometry objects (i.e. Point, LineString, Polygon) depending on
        the shapefiles original shape type

    Credit: https://github.com/ojdo/python-tools/blob/master/pandashp.py

    """
    sr = shapefile.Reader(filename)

    cols = sr.fields[:] # [:] = duplicate field list
    if cols[0][0] == 'DeletionFlag':
        cols.pop(0)
    cols = [col[0] for col in cols] # extract field name only
    cols.append('geometry')

    records = [row for row in sr.iterRecords()]

    if sr.shapeType == shapefile.POLYGON:
        geometries = [Polygon(shape.points)
                      if len(shape.points) > 2 else np.NaN  # invalid geometry
                      for shape in sr.iterShapes()]
    elif sr.shapeType == shapefile.POLYLINE:
        geometries = [LineString(shape.points) for shape in sr.iterShapes()]
    elif sr.shapeType == shapefile.POINT:
        geometries = [Point(*shape.points[0]) for shape in sr.iterShapes()]
    else:
        raise NotImplementedError

    data = [r+[g] for r,g in zip(records, geometries)]

    df = pd.DataFrame(data, columns=cols)
    df = df.convert_objects(convert_numeric=True)

    if np.NaN in geometries:
        # drop invalid geometries
        df = df.dropna(subset=['geometry'])
        num_skipped = len(geometries) - len(df)
        warnings.warn('Skipped {} invalid geometrie(s).'.format(num_skipped))
    return df


def _process_overlay_dict(input_dict):
    """
    process OCO-2 overlay data input.
    this is a helper normally called by process_config_dict(), see docstring there.

    In general, this is doing similar work as process_config_dict(), just this function
    acts on the data inside the overlay information from the JSON file (the 1-level
    nested dictionary)
    """
    ovr_d = collections.OrderedDict()

    required_keys = ('file', 'variable', 'lat_name', 'lon_name', 'orbit')
    for k in required_keys:
        if k not in input_dict:
            raise ValueError('Config file, overlay info is missing required key: '+k)

    ovr_d['var_file'] = input_dict['file']
    ovr_d['var_name'] = input_dict['variable']
    # trims off group names (for certain vars) - this is needed at the end
    # for the auto generated filename
    ovr_d['var_name_only'] = re.split('/', ovr_d['var_name'])[-1]
    ovr_d['lat_name'] = input_dict['lat_name']
    ovr_d['lon_name'] = input_dict['lon_name']
    ovr_d['orbit'] = int(input_dict['orbit'])

    ovr_d['band_number'] = input_dict.get('band_number', None)
    # convert empty string to None
    if ovr_d['band_number'] == "":
        ovr_d['band_number'] = None

    # for the moment, this is only needed for one variable (I think),
    # Retrieval/reduced_chi_squared_per_band, and only in V7, because this was 
    # stored as a [N,3] variable.
    if ovr_d['band_number'] is not None:
        if ovr_d['var_name'] == "Retrieval/reduced_chi_squared_per_band":
            err_msg = ("Unexpected band number specification for "+ovr_d['var_name']+
                       ". Options are 1 (0.76 micron), 2 (1.6 micron), or 3 (2.04 micron)")
            try:
                ovr_d['band_number'] = int(ovr_d['band_number'])
            except ValueError:
                raise ValueError(err_msg)
            if ovr_d['band_number'] not in (1,2,3):
                raise ValueError(err_msg)

    if 'variable_plot_lims' in input_dict:
        ovr_d['var_lims'] = input_dict['variable_plot_lims']
        if ovr_d['var_lims']:
            # if this was a string option, then make sure it is one
            # of the two options we have implemented. if not, revert to default.
            if isinstance(ovr_d['var_lims'], string_types):
                if ovr_d['var_lims'] not in ('autoscale_by_orbit', 'autoscale_by_overlay'):
                    print('var_lims "' + ovr_d['var_lims'] +
                          '" is not valid, reverting to "autoscale_by_orbit"')
                    ovr_d['var_lims'] = 'autoscale_by_orbit'
            else:
                try:
                    ovr_d['var_lims'] = float(ovr_d['var_lims'][0]), float(ovr_d['var_lims'][1])
                except:
                    raise ValueError('input variable_plot_lims is invalid')
        else:
            ovr_d['var_lims'] = 'autoscale_by_orbit'
    else:
        ovr_d['var_lims'] = 'autoscale_by_orbit'

    ovr_d['cmap'] = input_dict.get('cmap', 'jet')
    if ovr_d['cmap'] == "":
        ovr_d['cmap'] = 'jet'
    # make sure this is a valid colormap OR color string.
    if ( (ovr_d['cmap'] not in plt.colormaps()) and
         (ovr_d['cmap'] not in mpl.colors.cnames.keys()) ):
        print(ovr_d['cmap']  + " is not a recognized color or colormap. "+
              "Data will be displayed in jet colormap")
        ovr_d['cmap'] = 'jet'

    ovr_d['alpha'] = input_dict.get('transparency', 1)
    if ovr_d['alpha'] == "":
        ovr_d['alpha'] = 1.0

    try:
        ovr_d['alpha'] = float(ovr_d['alpha'])
    except:
        raise ValueError('Input transparency value is invalid. Must be floating point'+
                         'number in range [0,1]')

    if ovr_d['alpha'] < 0 or ovr_d['alpha'] > 1:
        raise ValueError('Input transparency value is invalid. Must be floating point'+
                         'number in range [0,1]')

    try:
        ovr_d['lite_quality'] = input_dict['lite_qf']
    except KeyError :
        print("No quality specifications detected. Output plot will contain all quality soundings")
        ovr_d['lite_quality'] = 'all'
    if not ovr_d['lite_quality']:
        print("No quality specifications detected. Output plot will contain all quality soundings")
        ovr_d['lite_quality'] = 'all'
    if ovr_d['lite_quality'] not in ['all', 'good', 'bad']:
        print("Unexpected quality flag specification. Options are: '', 'all', 'good', or 'bad'")
        print("Exiting")
        sys.exit()

    ovr_d['version_number'] =  re.sub(
        "[A-Za-z_]", "", re.search("_B[0-9a-z]{,5}_",os.path.basename(ovr_d['var_file'])).group())

    if int(ovr_d['version_number']) < 9000:
        try:
            ovr_d['lite_warn_lims'] = input_dict['lite_warn_lims']
        except KeyError:
            print("No warn specifications detected. Output plot will contain all warn levels")
            ovr_d['lite_warn_lims'] = [0, 20]
        if not ovr_d['lite_warn_lims']:
            print("No warn specifications detected. Output plot will contain all warn levels")
            ovr_d['lite_warn_lims'] = [0, 20]
        if ovr_d['lite_warn_lims'][0] > ovr_d['lite_warn_lims'][1]:
            print("Lower warn limit is greater than upper warn limit.")
            print("Exiting")
            sys.exit()
        for lim in ovr_d['lite_warn_lims']:
            if lim not in np.arange(21):
                print("Unexpected warn level specification. Limits must be within [0, 20].")
                print("Exiting")
                sys.exit()
        ovr_d['warn'] = True
    else:
        try:
            ovr_d['lite_warn_lims'] = input_dict['lite_warn_lims']
            if ovr_d['lite_warn_lims']:
                print("Warn levels are not available in B9 and above")
        except:
            pass
        ovr_d['lite_warn_lims'] = None
        ovr_d['warn'] = False

    try:
        ovr_d['footprint_lims'] = input_dict['footprint']
    except KeyError:
        print("No footprint specifications detected. Output plot will contain all footprints")
        ovr_d['footprint_lims'] = [1, 8]
    if not ovr_d['footprint_lims']:
        print("No footprint specifications detected. Output plot will contain all footprints")
        ovr_d['footprint_lims'] = [1, 8]
    if ovr_d['footprint_lims'] == "all":
        ovr_d['footprint_lims'] = [1, 8]
    try:
        if len(ovr_d['footprint_lims']) == 2:
            if ovr_d['footprint_lims'][0] > ovr_d['footprint_lims'][1]:
                print("Lower footprint limit is greater than upper footprint limit.")
                print("Exiting")
                sys.exit()
        elif len(ovr_d['footprint_lims']) > 2:
            print("footprint limits should be a scalar integer or 2-element integer list")
            print("Exiting")
            sys.exit()
    except TypeError:
        # catch for a input scalar integer, put it back in a 1-element list.
        # (the scalar will generate a TypeError on len().)
        ovr_d['footprint_lims'] = [ovr_d['footprint_lims']]
        
    for ft in ovr_d['footprint_lims']:
        if ft not in np.arange(1, 9):
            print("Unexpected footprint specification. Limits must be within [1, 8].")
            print("Exiting")
            sys.exit()

    if len(ovr_d['footprint_lims']) == 2:
        ovr_d['fp_file_tag'] = ("_FP_" + str(ovr_d['footprint_lims'][0]) +
                                "to" + str(ovr_d['footprint_lims'][1]))
    else:
        ovr_d['fp_file_tag'] = "_FP_"+str(ovr_d['footprint_lims'][0])

    if re.search("CO2", os.path.basename(ovr_d['var_file'])):
        ovr_d['sif_or_co2'] = "CO2"
        ovr_d['qf_file_tag'] = "_"+ovr_d['lite_quality']+"_quality"
        if ovr_d['warn']:
            ovr_d['wl_file_tag'] = ("_WL_" + str(ovr_d['lite_warn_lims'][0]) +
                                    "to"+str(ovr_d['lite_warn_lims'][1]))
        else:
            ovr_d['wl_file_tag'] = ""

    elif re.search("SIF", os.path.basename(ovr_d['var_file'])):
        ovr_d['sif_or_co2'] = "SIF"
        ovr_d['qf_file_tag'] = ""
        ovr_d['wl_file_tag'] = ""
    else:
        print("The command line usage of this tool accommodates official CO2 or SIF lite files only.")
        print("Expected filename convention: oco2_LtNNN_YYYYMMDD_B8xxxxx_rxx*.nc, where NNN is CO2 or SIF")
        print("Other data will need to be plotted by importing the tool as a Python module.")
        print("Exiting")
        sys.exit()

    ovr_d['version_file_tag'] = re.search(
        "_B[0-9a-z]{,5}_", os.path.basename(ovr_d['var_file'])).group()[:-1]

    if ovr_d['sif_or_co2'] == "CO2":
        if ovr_d['warn']:
            print("Output plot will include "+ovr_d['lite_quality']+
                  " quality soundings with warn levels within "+
                  str(ovr_d['lite_warn_lims'])+"\n")
        else:
            print("Output plot will include "+ovr_d['lite_quality']+" quality soundings\n")


    return ovr_d


def process_config_dict(input_dict):
    """
    process config file contents.

    this cleans inputs, checking for valid contents, and then creates
    namedtuple objects containing all configuration settings.

    inputs: config_file_dict, the dictionary as returned from the read
    of the JSON dictionary.

    outputs:

    config_dict: a dictionary containing the configuration settings
        related to the base image (MODIS or otherwise)

    contents:
    date: "YYYY-MM-DD"
    straight_up_date: "YYYYMMDD"
    lat_ul
    lon_ul
    lat_lr
    lon_lr
    region: string label
    ground_site: lat, lon of ground site to label
    city_labels: color to use for city labels, or None, ""
    datetime: python datetime object for date.

    overlay_dict: a dictionary containing information related to the
        OCO-2 data overlay for the plot

    contents:
    var_file: filename to use for overlay data (must be OCO2 Lite file)
    var_name: variable name to use
    lat_name: lat variable to use (to control footprint lat, or vertices)
    lon_name: lon variable to use
    orbit: integer orbit number
    var_lims: plot limits, or "autoscale_by_orbit", "autoscale_by_overlay"
    cmap: string colormap name, defaults to 'jet'
    cmap_obj: the MPL linearSegmentedColormap object
    alpha: transparency variable (1=opaque, 0=transparent), float in range [0,1]
    """

    cfg_d = collections.OrderedDict()

    # check for required keys
    required_keys = ('date', 'geo_upper_left', 'geo_lower_right')
    for k in required_keys:
        if k not in input_dict:
            raise ValueError('Config file is missing required key: '+k)

    # first just copy from input dictionary into OrderedDict for config data.
    # note the use of dict.get(key, default) will effectively set default blank or None
    # values.
    cfg_d['date'] = input_dict['date']
    cfg_d['straight_up_date'] = cfg_d['date'].replace("-", "")

    cfg_d['geo_upper_left'] = input_dict['geo_upper_left']
    cfg_d['geo_lower_right'] = input_dict['geo_lower_right']

    try:
        cfg_d['lat_ul'] = float(input_dict['geo_upper_left'][0])
        cfg_d['lon_ul'] = float(input_dict['geo_upper_left'][1])
    except ValueError:
        raise ValueError('geo_upper_left contents are invalid')
    try:
        cfg_d['lat_lr'] = float(input_dict['geo_lower_right'][0])
        cfg_d['lon_lr'] = float(input_dict['geo_lower_right'][1])
    except ValueError:
        raise ValueError('geo_lower_right contents are invalid')

    cfg_d['region'] = input_dict.get('region', '')

    cfg_d['ground_site'] = input_dict.get('ground_site', None)
    if cfg_d['ground_site'] == "":
        cfg_d['ground_site'] = None
    cfg_d['city_labels'] = input_dict.get('city_labels', None)
    if cfg_d['city_labels'] == "":
        cfg_d['city_labels'] = None

    cfg_d['out_plot_dir'] = input_dict.get('out_plot_dir', '')
    cfg_d['out_plot_name'] = input_dict.get('out_plot_name', '')
    cfg_d['out_data_dir'] = input_dict.get('out_data_dir', '')
    cfg_d['out_data_name'] = input_dict.get('out_data_name', '')

    # now do various checks.
    try:
        cfg_d['datetime'] = datetime.datetime.strptime(cfg_d['date'], '%Y-%m-%d')
    except:
        raise ValueError('input field "date" has incorrect format, expecting YYYY-MM-DD')

    if cfg_d['ground_site']:
        cfg_d['ground_site'] = float(cfg_d['ground_site'][0]), float(cfg_d['ground_site'][1])
        if (cfg_d['ground_site'][0] > cfg_d['lat_ul'] or
            cfg_d['ground_site'][0] < cfg_d['lat_lr'] or
            cfg_d['ground_site'][1] > cfg_d['lon_lr'] or
            cfg_d['ground_site'][1] < cfg_d['lon_ul']):
            cfg_d['ground_site'] = None
            print("The ground site is outside the given lat/lon range and "+
                  "will not be included in the output plot.\n")

    if cfg_d['city_labels'] == '':
        cfg_d['city_labels'] = None

    if cfg_d['city_labels'] and cfg_d['city_labels'] not in mpl.colors.cnames.keys():
        print(cities + " is not an available matplotlib color. "+
              "City labels will not be included on the output plot. \n")
        cfg_d['city_labels'] = None

    if not cfg_d['out_plot_dir'] or not glob(cfg_d['out_plot_dir']):
        print("Either there was no output plot location specified or the one specified "+
              "does not exist. Plot output will go in the code directory. \n")
    if not cfg_d['out_data_dir'] or not glob(cfg_d['out_data_dir']):
        print("Either there was no output data location specified or the one specified "+
              "does not exist. Data output will go in the code directory. \n")
    
    if 'oco2_overlay_info' in input_dict:
        ovr_d = _process_overlay_dict(input_dict['oco2_overlay_info'])
        # copy the geo corners
        for k in ('lat_ul', 'lat_lr', 'lon_ul', 'lon_lr'):
            ovr_d[k] = cfg_d[k]
    else:
        ovr_d = collections.OrderedDict()

    return cfg_d, ovr_d


def write_h5_data_subset(out_data_filename, lite_sid_subset,
                         var_lat_subset, var_lon_subset, var_vals_subset,
                         lat_name, lon_name, var_name):
    """
    writes out the h5 subset file, containing the OCO-2 data plotted in the
    overlay.

    If the arrays are masked array types, runs compressed (or compress_rows
    if the arrays are 2-D) on the arrays before writing to h5.

    inputs:
    out_data_filename, the name for the created h5 file
    lite_sid_subset: sounding id array. If this is an empty array, the
        sounding_id dataset will not be created in the output file
    var_lat_subset: subset of latitude, 1 or 2D (if vertices)
    var_lon_subset: subset of longitude, 1 or 2D (if vertices)
    var_vals_subset: data subset
    lat_name: string containing the variable name for latitude
    lon_name: string containing the variable name for longitude
    var_name: string containing the variable name

    output file will contain up to 4 datasets (if sounding_id is included),
    for example, with vertex lat/lon, and using xco2:
    /                        Group
    /sounding_id             Dataset {494}
    /vertex_latitude         Dataset {494, 4}
    /vertex_longitude        Dataset {494, 4}
    /xco2                    Dataset {494}
    """

    if re.search("Masked", str(type(var_lat_subset))):
        if var_lat_subset.ndim > 1:
            lat_data_to_write = np.ma.compress_rows(var_lat_subset)
            lon_data_to_write = np.ma.compress_rows(var_lon_subset)
        else:
            lat_data_to_write = var_lat_subset.compressed()
            lon_data_to_write = var_lon_subset.compressed()
    else:
        lat_data_to_write = var_lat_subset
        lon_data_to_write = var_lon_subset

    if re.search("Masked", str(type(var_vals_subset))):
        if var_vals_subset.ndim > 1:
            var_data_to_write = np.ma.compress_rows(var_vals_subset)
        else:
            var_data_to_write = var_vals_subset.compressed()
    else:
        var_data_to_write = var_vals_subset

    if lite_sid_subset.shape:
        if re.search("Masked", str(type(lite_sid_subset))):
            sid_data_to_write = lite_sid_subset.compressed()
        else:
            sid_data_to_write = lite_sid_subset

    if not lat_name:
        lat_name = "Latitude"
    if not lon_name:
        lon_name = "Longitude"
    if not var_name:
        var_name = "Data"

    outfile = h5py.File(out_data_filename, "w")

    write_ds = outfile.create_dataset(lat_name, data = lat_data_to_write)
    write_ds = outfile.create_dataset(lon_name, data = lon_data_to_write)
    write_ds = outfile.create_dataset(var_name, data = var_data_to_write)
    if lite_sid_subset.shape:
        write_ds = outfile.create_dataset("sounding_id", data = sid_data_to_write)

    outfile.close()

    print("\nData saved at "+out_data_filename)


def _compute_ll_msk(lat_data, lon_data, lat_ul, lat_lr, lon_ul, lon_lr):
    """
    helper to subset the data according to lat/lon upper left, lower right corners.

    prints some console messages if the lat/lon mask removes all the data;
    this would happen if the user input LL box from the JSON does not contain
    the OCO-2 track.
    """

    # get subset of var values etc.
    lat_msk = np.logical_and(lat_data <= lat_ul, lat_data >= lat_lr)
    lon_msk = np.logical_and(lon_data <= lon_lr, lon_data >= lon_ul)
    ll_msk = np.logical_and(lat_msk, lon_msk)

    # Ensure that for vertex input, if any of the vertices are outside the lat/lon
    # limits, that the footprint will be masked out.
    if lat_data.ndim == 2:
        ll_msk = np.all(ll_msk, axis=1)

    # if there are no common points, print some info back to console
    # to suggest where the LL box needs to be to capture some OCO2 data.
    if not np.any(ll_msk):
        lon_subset_idx = lon_msk.nonzero()[0]
        print("No OCO-2 data found in the input lat/lon box, with coords:")
        print("lat/lon UL: ", lat_ul, lon_ul)
        print("lat/lon LR: ", lat_lr, lon_lr)
        if lon_subset_idx.shape[0] == 0:
            print("No OCO-2 data found within the lon range")
        else:
            min_idx = lon_subset_idx.min()
            max_idx = lon_subset_idx.max()
            print("Lite indices where the longitude is between " +
                  str(lon_ul) + " and " + str(lon_lr) +": " +
                  str(min_idx) + "-" + str(max_idx))
            print("Latitude range for those indices: " +
                  str(lat_data[min_idx]) + "-" + str(lat_data[max_idx]))
            print("Latitude range given: " + str(lat_lr) + "-" + str(lat_ul))
        
    return ll_msk
    

def load_OCO2_Lite_overlay_data(ovr_d):
    """
    loads OCO2 Lite data, given the processed overlay dictionary as input

    returns a dictionary containing the OCO-2 data subsetted to the overlay
    lat/lon corners. contains the following fields:

    sif_or_co2: "CO2" or "SIF", denoting the type of the data file
    data_long_name: string long name for the variable (from the netCDF attributes)
    data_units: string units for the variable (from the netCDF attributes)
    data_fill: the missing data fill value
    data: the subsetted data values, with shape (N,)
    lat, lon: the subsetted lat/lon, with shape (N,) (for footprint lat/lon), 
        or shape (N,4) (for footprint vertex lat/lon)
    """

    # data dictionary that will be constructed.
    dd = collections.OrderedDict()

    ### Prep OCO-2 Variable ###

    h5 = h5py.File(ovr_d['var_file'], "r")

    try:
        data_obj = h5[ovr_d['var_name']]
    except:
        print(ovr_d['var_name']+" DNE in "+ovr_d['var_file'])
        print("Check that the variable name includes any necessary group paths. Ex: /Preprocessors/dp_abp")
        print("Exiting")
        sys.exit()

    data = data_obj[:]

    if ovr_d['sif_or_co2'] == "CO2":
        try:
            dd['data_long_name'] = data_obj.attrs.get('long_name')[0]
        except:
            print("Problem reading long name for " + ovr_d['var_name'])
            dd['data_long_name'] = ""
        try:
            dd['data_units'] = data_obj.attrs.get('units')[0]
        except:
            print("Problem reading units for " + ovr_d['var_name'])
            dd['data_units'] = ""
        try:
            dd['data_fill'] = data_obj.attrs.get('missing_value')[0]
        except:
            print("Problem reading missing value for " + ovr_d['var_name'])
            dd['data_fill'] = ""
    if ovr_d['sif_or_co2'] == "SIF":
        dd['data_long_name'] = re.split("/", ovr_d['var_name'])[-1]
        try:
            dd['data_units'] = data_obj.attrs.get('unit')
        except:
            print("Problem reading units for " + ovr_d['var_name'])
            dd['data_units'] = ""
        try:
            dd['data_fill'] = data_obj.attrs.get('missing_value')[0]
        except:
            print("Problem reading missing value for " + ovr_d['var_name'])
            dd['data_fill'] = ""

    # the string variables changed in B7 vs B8 and later; B7 stored them
    # in a way that it loads as a bytes type in python, later builds loads
    # as a string.
    for v in ('data_long_name', 'data_units'):
        if type(dd[v]) is not str:
            dd[v] = dd[v].decode('utf-8')

    try:
        lat_data = h5[ovr_d['lat_name']][:]
    except:
        print(ovr_d['lat_name']+" DNE in "+ovr_d['var_file'])
        print("Check that the variable name includes any necessary group paths. Ex: SoundingGeometry/sounding_latitude")
        print("Exiting")
        sys.exit()

    try:
        lon_data = h5[ovr_d['lon_name']][:]
    except:
        print(ovr_d['lon_name']+" DNE in "+ovr_d['var_file'])
        print("Check that the variable name includes any necessary group paths. Ex: SoundingGeometry/sounding_longitude")
        print("Exiting")
        sys.exit()
    h5.close()

    if lat_data.ndim != lon_data.ndim:
        print(ovr_d['lat_name']+" and "+ovr_d['lon_name']+" have different dimensions. Exiting")
        sys.exit()

    if ovr_d['var_name'] == "Retrieval/reduced_chi_squared_per_band":
        if not ovr_d['band_number']:
            print(ovr_d['var_name'] + " is stored per band. "+
                  "Please select a band number (1=0.76 micron, 2=1.6 micron, 3=2.04 micron). Exiting")
            sys.exit()
        else:
            data = data[:,ovr_d['band_number']-1]

    # use helper object to get the remaining variables (sid, warnlevel, QF, orbit, footprint)
    if ovr_d['sif_or_co2'] == "CO2":
        lite_file = LiteCO2File(ovr_d['var_file'])
        lite_file.open_file()
        if ovr_d['warn']:
            lite_warn = lite_file.get_warn()
        lite_qf = lite_file.get_qf()
    else:
        lite_file = LiteSIFFile(ovr_d['var_file'])
        lite_file.open_file()
    lite_time = lite_file.get_time()
    lite_sid = lite_file.get_sid()
    lite_orbit = lite_file.get_orbit()
    lite_footprint = lite_file.get_footprint()
    lite_file.close_file()

    # create a list of masks to logical_and together
    # orbit, WL, QF are easy, LatLon is more complicated, so there is a helper function.
    # LatLon is done last, we need to check for the autoscale limits there
    # before applying the subset: the autoscale_by_orbit will use all data from
    # the orbit (and WL, QF, etc)
    # Note: consider checking these masks as they are computed, to inform the user
    # if one or more removes all available observations.
    msk_list = []

    if ovr_d['orbit']:
        orbit_msk = lite_orbit == ovr_d['orbit']
        msk_list.append(orbit_msk)
        
    if ovr_d['sif_or_co2'] == "CO2":

        if ovr_d['lite_quality'] == 'good':
            qf_msk = lite_qf == 0
        elif ovr_d['lite_quality'] == 'bad':
            qf_msk = lite_qf == 1
        else:
            qf_msk = np.ones(lite_qf.shape, np.bool)
        msk_list.append(qf_msk)

        if ovr_d['warn']:
            warn_msk = np.logical_and(
                lite_warn >= ovr_d['lite_warn_lims'][0],
                lite_warn <= ovr_d['lite_warn_lims'][1])
            msk_list.append(warn_msk)

    if len(ovr_d['footprint_lims']) == 2:
        fp_msk = np.logical_and(
            lite_footprint >= ovr_d['footprint_lims'][0],
            lite_footprint <= ovr_d['footprint_lims'][1])
    else:
        fp_msk = lite_footprint == ovr_d['footprint_lims'][0]
    msk_list.append(fp_msk)

    combined_msk = msk_list[0].copy()
    for msk in msk_list[1:]:
        combined_msk = np.logical_and(combined_msk, msk)

    # apply the combined mask to this point (includes all
    # but latlon), to get the by orbit data range.
    tmp_data = data[combined_msk]
    dd['orbit_var_lims'] = np.ma.min(tmp_data), np.ma.max(tmp_data)

    # now apply the latlon mask, and proceed to do the
    # subsetting.
    ll_msk = _compute_ll_msk(lat_data, lon_data,
                             ovr_d['lat_ul'], ovr_d['lat_lr'],
                             ovr_d['lon_ul'], ovr_d['lon_lr'])
    combined_msk = np.logical_and(combined_msk, ll_msk)

    dd['var_data'] = data[combined_msk]
    # note the ellipsis will work equally well for 1D or 2D (vertex) data
    dd['lat'] = lat_data[combined_msk, ...]
    dd['lon'] = lon_data[combined_msk, ...]
    dd['sounding_id'] = lite_sid[combined_msk]
    dd['time'] = lite_time[combined_msk]

    return dd


def do_modis_overlay_plot(
    geo_upper_left, geo_lower_right, date,
    var_lat, var_lon, var_vals, var_vals_missing=None, lite_sid=np.empty([]),
    var_lims=None, interest_pt=None,
    cmap='jet', alpha=1, lat_name=None, lon_name=None, var_name=None,
    out_plot="vistool_output.png", out_data="vistool_output.h5", var_label=None, cities=None):

    lat_ul = geo_upper_left[0]
    lon_ul = geo_upper_left[1]
    lat_lr = geo_lower_right[0]
    lon_lr = geo_lower_right[1]

    ### Pull Aqua-MODIS RGB from GIBS ###

    update_GIBS_xml(date, xml_file)

    print("Pulling RGB")
    try:
        pull_Aqua_RGB_GIBS(lat_ul, lon_ul, lat_lr, lon_lr,  xml_file, code_dir+'/intermediate_RGB.tif')
    except:
        print("Problem pulling RGB. Check that the geolocation bounds specified in the configuration file are for the upper left hand corner and the lower right hand corner")


    ### Pull in and prep RGB tif file ###

    ds = gdal.Open(code_dir+'/intermediate_RGB.tif')

    data = ds.ReadAsArray()
    gt = ds.GetGeoTransform()
    proj = ds.GetProjection()

    inproj = osr.SpatialReference()
    inproj.ImportFromWkt(proj)

    width = ds.RasterXSize
    height = ds.RasterYSize

    #Calculate lat/lon lims of RGB
    minx = gt[0]
    miny = gt[3] + width*gt[4] + height*gt[5]
    maxx = gt[0] + width*gt[1] + height*gt[2]
    maxy = gt[3]


    #Check if color is single or map

    if cmap in plt.colormaps():
        color_or_cmap = "cmap"
    elif cmap in mpl.colors.cnames.keys():
        color_or_cmap = "color"
    else:
        print(cmap + " is not a recognized color or colormap. Data will be displayed in red")
        cmap = 'red'
        color_or_cmap = "color"

    fig_x = abs(gt[1])
    fig_y = abs(gt[5])

    while fig_x < 1 and fig_y < 1:
        fig_x *= 10
        fig_y *= 10

    while fig_x < 5 or fig_y < 5:
        if fig_x < 10 and fig_y < 10:
            fig_x *= 1.1
            fig_y *= 1.1
        else:
            break

    if color_or_cmap == "cmap":
        fig_x += 2

    # note - if var_vals.shape is zero, the var_lims are not needed, since
    # the colorbar and scatter or polygon collection will not be plotted.
    if var_lims is None and var_vals.shape[0] > 0:
            var_lims = [var_vals.min(), var_vals.max()]

    ### Plot prep ###
    states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none')


    ### Plot the image ###
#    if var_vals.shape:
#        fig = plt.figure(figsize=(fig_x + 1,fig_y))
#    else:
#        fig = plt.figure(figsize=(fig_x,fig_y))

    fig = plt.figure(figsize=(fig_x + 1,fig_y))
    gs =  gridspec.GridSpec(16, 16)

    img = plt.imread(code_dir+'/intermediate_RGB.tif')
    img_extent = (minx, maxx, miny, maxy)

    ax = plt.subplot(gs[0:-1, 3:-2], projection=ccrs.PlateCarree())
    ax_pos = ax.get_position()
    ax.imshow(img, origin='upper', transform=ccrs.PlateCarree(), extent=img_extent, aspect='auto')
    ax.coastlines(resolution='10m', color='black', linewidth=1)
    ax.add_feature(states_provinces, edgecolor='black', linewidth=1)
    ax.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=1)

    if cities is not None:

        populated_places_filename = code_dir+'/natural_earth/ne_10m_populated_places'

        df = read_shp(populated_places_filename)

        relevant_places = df[(df['LATITUDE'] <= maxy) &
                          (df['LATITUDE']>= miny) &
                          (df['LONGITUDE']<= maxx) &
                          (df['LONGITUDE']>= minx)]

        for idx, p in relevant_places.iterrows():

            #print p['NAME'], p['LATITUDE'], p['LONGITUDE']
            ax.text(p['LONGITUDE'], p['LATITUDE'], p['NAME'], fontsize=7, color=cities, va='bottom', ha='center', transform=ccrs.Geodetic())

    if interest_pt is not None:
        ax.plot(interest_pt[1], interest_pt[0], 'w*', markersize=10, transform=ccrs.Geodetic())

    ylocs, ylabels = plt.yticks()
    xlocs, xlabels = plt.xticks()

    ax_minlat = plt.subplot(gs[-1, 1])
    ax_maxlat = plt.subplot(gs[0, 1])
    ax_minlon = plt.subplot(gs[-1, 3])
    ax_maxlon = plt.subplot(gs[-1, -2])

    ax_minlat.axis('off')
    ax_maxlat.axis('off')
    ax_maxlon.axis('off')
    ax_minlon.axis('off')

    ax_minlat.set_title(str("%.1f" % ylocs[0]), horizontalalignment='left', verticalalignment='bottom', fontsize=10, fontweight='bold')
    ax_maxlat.set_title(str("%.1f" % ylocs[-1]), horizontalalignment='left', verticalalignment='top', fontsize=10, fontweight='bold')
    ax_minlon.set_title(str("%.1f" % xlocs[0]), horizontalalignment='center', verticalalignment='top', fontsize=10, fontweight='bold')
    ax_maxlon.set_title(str("%.1f" % xlocs[-1]), horizontalalignment='right', verticalalignment='top', fontsize=10, fontweight='bold')
    
    patches = []

    if var_lat.ndim == 2:
        zip_it = np.ma.dstack([var_lon, var_lat])

    if color_or_cmap == "cmap" and var_vals.shape[0] > 0:
    
        # vertex points for each footprint
        if var_lat.ndim == 2 and var_vals.ndim == 1:
            for row in range(zip_it.shape[0]):
                polygon = mpatches.Polygon(zip_it[row,:,:])
                patches.append(polygon)

            p = mpl.collections.PatchCollection(patches, cmap=cmap, alpha=alpha, edgecolor='none')
            p.set_array(var_vals)
            p.set_clim(var_lims[0], var_lims[1])
            ax.add_collection(p)

        # or just plot the central (sounding) lat lon.
        else:
            ax.scatter(var_lon, var_lat, c=var_vals,
                       cmap=cmap, edgecolor='none', s=2, vmax=var_lims[1], vmin=var_lims[0])
            
        cb_ax1 = plt.subplot(gs[0:-1, -1])
        norm = mpl.colors.Normalize(vmin = var_lims[0], vmax = var_lims[1])
        cmap_obj = mpl.cm.get_cmap(cmap)
        cb1 = mpl.colorbar.ColorbarBase(cb_ax1, cmap=cmap_obj, orientation = 'vertical', norm = norm)
        if var_label:
            cb1_lab = cb1.ax.set_xlabel(var_label, labelpad=8, fontweight='bold')
            cb1_lab.set_fontsize(14)
            cb1.ax.xaxis.set_label_position("top")
        for t in cb1.ax.yaxis.get_ticklabels():
            t.set_weight("bold")
            t.set_fontsize(12)

    if color_or_cmap == "color" and var_vals.shape[0] > 0:

        # vertex points for each footprint
        if var_lat.ndim == 2 and var_vals.ndim == 1:
            for row in range(zip_it.shape[0]):
                polygon = mpatches.Polygon(zip_it[row,:,:], color=cmap)
                patches.append(polygon)
            p = mpl.collections.PatchCollection(patches, alpha=alpha, edgecolor='none',
                                                match_original=True)
            ax.add_collection(p)

        # or just plot the central (sounding) lat lon.
        else:
            ax.scatter(var_lon, var_lat, c=cmap, edgecolor='none', s=2)


    # during testing, it appears that sometimes the scatter
    # could cause MPL to shift the axis range - I think because one
    # scatter point goes slightly out of the display range.
    # so, here force it back to the original domain.
    ax.axis(img_extent)

    inset_extent_x = [minx, maxx]
    inset_extent_y = [miny, maxy]

    inset_extent_x = [x + 360 if x < 0 else x for x in inset_extent_x]
    inset_extent_y = [y + 180 if y < 0 else y for y in inset_extent_y]

    inset_extent_x[0] -= 20
    inset_extent_y[0] -= 20
    inset_extent_x[1] += 20
    inset_extent_y[1] += 20

    inset_extent_x = [x - 360 if x > 180 else x for x in inset_extent_x]
    inset_extent_y = [y - 180 if y > 90 else y for y in inset_extent_y]

    inset_ax = plt.subplot(gs[7:9, 0:3], projection=ccrs.PlateCarree())
    inset_ax.set_extent([inset_extent_x[0], inset_extent_x[1], inset_extent_y[0], inset_extent_y[1]])

    inset_ax.coastlines()
    inset_ax.add_feature(cfeature.LAKES, edgecolor='black', facecolor='none')
    extent_box = sgeom.box(minx, miny, maxx, maxy)
    inset_ax.add_geometries([extent_box], ccrs.PlateCarree(), color='none', edgecolor='red')
    inset_ax.set_aspect('auto')

    fig.savefig(out_plot, dpi=150, bbox_inches='tight')
    print("\nFigure saved at "+out_plot)
    #print("code directory in subroutine:", code_dir)
    os.remove(code_dir+'/intermediate_RGB.tif')


### Static Defnitions

code_dir = os.path.dirname(os.path.realpath(__file__))
xml_file = code_dir+'/GIBS_Aqua_MODIS_truecolor.xml'

if __name__ == "__main__":

    ### Dynamic Definitions: get information from config file ###

    parser = argparse.ArgumentParser(description="Get configuration file")
    parser.add_argument(
        'config_file_loc', type=str,
        default=code_dir+'/oco2_modis_vistool_config.json', nargs='?',
        help="Name of config file (default: oco2_modis_vistool_config.json in code directory)")
    args = parser.parse_args()
    config_file = ConfigFile(args.config_file_loc)

    if config_file.exists():
        input_dict = config_file.get_contents()
        for k in input_dict.keys():
            input_dict[k.lower()] = input_dict.pop(k)
    else:
        print('The expected configuration file '+ args.config_file_loc + ' DNE in ' + code_dir)
        print('Exiting')
        sys.exit()

    cfg_d, ovr_d = process_config_dict(input_dict)

    # with no overlay data requested, simply generate the MODIS image by itself and exit.
    if not ovr_d:

        if not cfg_d['out_plot_name']:
            if cfg_d['region']:
                out_plot_name = ("MODISimagery_" + cfg_d['region'] +
                                 "_" + cfg_d['straight_up_date'] + ".png")
            else:
                out_plot_name = ("MODISimagery_" + cfg_d['straight_up_date'] + ".png")
            out_plot_name = os.path.join(cfg_d['out_plot_dir'], out_plot_name)
        else:
            out_plot_name = os.path.join(cfg_d['out_plot_dir'], 
                                         cfg_d['out_plot_name'])

        do_modis_overlay_plot(cfg_d['geo_upper_left'],
                              cfg_d['geo_lower_right'],
                              cfg_d['date'], np.array([]), np.array([]), np.array([]), np.array([]),
                              interest_pt=cfg_d['ground_site'], cmap='black',
                              out_plot=out_plot_name, cities=cfg_d['city_labels'])

        sys.exit()


    odat = load_OCO2_Lite_overlay_data(ovr_d)
    
    # here, handle the var limit options.
    # if specific limits were input, via a 2-element list,
    # that will be passed directly to do_modis_overlay_plot().
    # if autoscaling by orbit, find those min/max now, and create the 2 element list.
    if ovr_d['var_lims'] == 'autoscale_by_orbit':
        ovr_d['var_lims'] = odat['orbit_var_lims']
    # otherwise, convert to None, so then do_modis_overlay_plot() will then derive an autoscale
    # min/max according to the points within the overlay.
    if ovr_d['var_lims'] == 'autoscale_by_overlay':
        ovr_d['var_lims'] = None

    ### Plot prep ###

    # create a compact colorbar label, including var name and units.
    if odat['data_long_name']:
        oco2_data_long_name = re.sub("_", " ", odat['data_long_name'])
        cbar_strings = re.split(' ', oco2_data_long_name)
        cbar_cap_strings = ""

        for i, s in enumerate(cbar_strings):
            cap_s = s
            if not s[0].isupper():
                cap_s = s.capitalize()
            if i % 2 != 0:
                cbar_cap_strings = cbar_cap_strings + cap_s +"\n"
            else:
                cbar_cap_strings = cbar_cap_strings + cap_s +" "
        cbar_cap_strings = cbar_cap_strings[:-1]
        cbar_name = cbar_cap_strings+'\n('+odat['data_units']+')'
    else:
        cbar_name = ""

    # construct filenames for the plot and optional h5 output file.
    if not cfg_d['out_plot_name']:
        out_plot_name = ovr_d['var_name_only']
        if cfg_d['region']:
            out_plot_name = out_plot_name + "_" + cfg_d['region']
        out_plot_name = (out_plot_name + "_" +
                         cfg_d['straight_up_date'] + 
                         ovr_d['version_file_tag']+
                         ovr_d['qf_file_tag']+
                         ovr_d['wl_file_tag']+
                         ovr_d['fp_file_tag']+".png")
        out_plot_name = os.path.join(cfg_d['out_plot_dir'], out_plot_name)
    else:
        out_plot_name = os.path.join(cfg_d['out_plot_dir'], 
                                     cfg_d['out_plot_name'])

    if not cfg_d['out_data_name']:
        out_data_name = ovr_d['var_name_only']
        if cfg_d['region']:
            out_data_name = out_data_name + "_" + cfg_d['region']
        out_data_name = (out_data_name + "_" +
                         cfg_d['straight_up_date'] + 
                         ovr_d['version_file_tag']+
                         ovr_d['qf_file_tag']+
                         ovr_d['wl_file_tag']+
                         ovr_d['fp_file_tag']+".h5")
        out_data_name = os.path.join(cfg_d['out_data_dir'], out_data_name)
    else:
        out_data_name = os.path.join(cfg_d['out_data_dir'], 
                                     cfg_d['out_data_name'])

    do_modis_overlay_plot(cfg_d['geo_upper_left'], cfg_d['geo_lower_right'],
                          cfg_d['date'], odat['lat'], odat['lon'], odat['var_data'],
                          var_vals_missing=odat['data_fill'],
                          lite_sid=odat['sounding_id'],
                          var_lims=ovr_d['var_lims'], interest_pt=cfg_d['ground_site'],
                          cmap=ovr_d['cmap'], alpha=ovr_d['alpha'],
                          lat_name=ovr_d['lat_name'], lon_name=ovr_d['lon_name'],
                          var_name=ovr_d['var_name'],
                          out_plot=out_plot_name, out_data=out_data_name,
                          var_label=cbar_name, cities=cfg_d['city_labels'])

    ### Write data to hdf5 file ###
    if odat['var_data'].shape[0] == 0:
        print("\nThere is no valid data for the given subset criteria.")
    else:
        write_h5_data_subset(
            out_data_name, odat['sounding_id'],
            odat['lat'], odat['lon'], odat['var_data'],
            ovr_d['lat_name'], ovr_d['lon_name'], ovr_d['var_name'])
