"""
oco2_modis_vistool.py

The purpose of this script is to pull Aqua-MODIS RGB images from Worldview
using the NASA GIBS API and overlay various OCO-2 data fields for case study 
analysis in support of OCO-2 cloud and aerosol screening validation.

GIBS developer documentation:https://wiki.earthdata.nasa.gov/display/GIBS/GIBS+API+for+Developers

Minimum command line call:
python oco2_modis_vistool.py

Output:
Image (.png) named _ in specified output directory (default: code directory)

Requirements:

1) json configuration file (default: oco2_modis_vistool_config.json in code directory)
2) This script is meant to run on ocomaster at Colorado State University. 
   To run it elsewhere, the data locations will need to be updated 
   and modules will need to be verified.

Authors:
Natalie Tourville <natalie.tourville@colostate.edu>
Heather Cronk <heather.cronk@colostate.edu>

Revision history:
08/2016 (N.D. Tourville; H.Q. Cronk): Initial version

"""

import warnings
warnings.filterwarnings("ignore")

import os
import sys
from glob import glob

from OCO2FileOps import *
import h5py

import numpy as np

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib as mpl
import matplotlib.pyplot as plt

from osgeo import gdal, osr
from shapely.ops import transform as geom_transform

import xml.etree.ElementTree as ET
import json
import argparse

import re

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
    
    cmd = "/usr/bin/gdal_translate -of GTiff -outsize "+str(xsize)+" "+str(ysize)+" -projwin "+str(lon_ul)+" "+str(lat_ul)+" "+str(lon_lr)+" "+str(lat_lr)+" "+xml_file+" "+tif_file

    os.system(cmd)

### Static Defnitions

code_dir = os.path.dirname(os.path.realpath(__file__))
xml_file = code_dir+'/GIBS_Aqua_MODIS_truecolor.xml'


### Dynamic Definitions: get information from config file ###

parser = argparse.ArgumentParser(description="Get configuration file")
parser.add_argument('config_file_loc', type=str, default=code_dir+'/oco2_modis_vistool_config.json', nargs='?',
                    help="Name of config file (default: oco2_modis_vistool_config.json in code directory)")
args = parser.parse_args()
config_file = ConfigFile(args.config_file_loc)

if config_file.exists():
    orbit_info_dict = config_file.get_contents()
else:
    print 'The expected configuration file '+ args.config_file_loc + ' DNE in ' + code_dir
    print 'Exiting'
    sys.exit()


date = orbit_info_dict['date']
straight_up_date = date.replace("-", "")
lat_ul = orbit_info_dict['geo_upper_left'][0]
lon_ul = orbit_info_dict['geo_upper_left'][1]
lat_lr = orbit_info_dict['geo_lower_right'][0]
lon_lr = orbit_info_dict['geo_lower_right'][1]
region = orbit_info_dict['region']
overlay_info_dict = orbit_info_dict['oco2_overlay_info']
var_file = overlay_info_dict['file']
if not glob(var_file):
    print var_file+" does not exist."
    print "Exiting"
    sys.exit()
var_name = overlay_info_dict['variable']
var_lims = overlay_info_dict['variable_plot_lims']
lat_name = overlay_info_dict['lat_name']
lon_name = overlay_info_dict['lon_name']
orbit_int = overlay_info_dict['orbit']

#delta_lat = lat_ul - lat_lr
#delta_lon = lon_ul - lon_lr
#if delta_lon > 180:
#    delta_lon -= 360
#delta_lon = abs(delta_lon)
#
#print delta_lat
#print delta_lon
##sys.exit()

if re.search('oco2_Lt', var_file):
    lite = True
    print "\nLite overlay file detected. Checking for QF and Warn specs..."
    
    try:
        lite_quality = overlay_info_dict['lite_QF']
    except:
        print "No quality specifications detected. Output plot will contain all quality soundings"
	lite_quality = 'all'
    if not lite_quality:
        print "No quality specifications detected. Output plot will contain all quality soundings"
	lite_quality = 'all'
    if lite_quality not in ['', 'all', 'good', 'bad']:
        print "Unexpected quality flag specification. Options are '', 'all', 'good', or 'bad'"
	print "Exiting"
	sys.exit()
    
    try:
        lite_warn_lims = overlay_info_dict['lite_warn_lims']
    except:
        print "No warn specifications detected. Output plot will contain all warn levels"
	lite_warn_lims = [0, 20]
    if not lite_warn_lims:
        print "No warn specifications detected. Output plot will contain all warn levels"
	lite_warn_lims = [0, 20]
    if lite_warn_lims[0] > lite_warn_lims[1]:
        print "Lower warn limit is greater than upper warn limit."
	print "Exiting"
	sys.exit()
    for lim in lite_warn_lims:
        if lim not in np.arange(21):
            print "Unexpected warn level specification. Limits must be within [0, 20]."
	    print "Exiting"
	    sys.exit()
    
print "Output plot will include "+lite_quality+" quality soundings with warn levels within "+str(lite_warn_lims)+"\n"

try:
    interest_pt = orbit_info_dict['ground_site']
except:
    interest_pt = []
try:
    output_dir = orbit_info_dict['output_dir']
except:
    output_dir = code_dir
if not output_dir or not glob(output_dir):
    print "Either there was no output location specified or the one specified does not exist. Output will go in the code directory"
    output_dir = code_dir

### Pull Aqua-MODIS RGB from GIBS ###

update_GIBS_xml(date, xml_file)

print "Pulling RGB"
try:
    pull_Aqua_RGB_GIBS(lat_ul, lon_ul, lat_lr, lon_lr,  xml_file, code_dir+'/intermediate_RGB.tif')
except:
    print "Problem pulling RGB. Check that the geolocation bounds specified in the configuration file are for the upper left hand corner and the lower right hand corner"


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


### Prep OCO-2 Variable ###

h5 = h5py.File(var_file)
try:
    oco2_data = h5[var_name][:]
except:
    print var_name+" DNE in "+var_file
    print "Check that the variable name includes any necessary group paths. Ex: /Preprocessors/dp_abp"
    print "Exiting"
    sys.exit()
try:
    lat_data = h5[lat_name][:]
except:
    print lat_name+" DNE in "+var_file
    print "Check that the variable name includes any necessary group paths. Ex: SoundingGeometry/sounding_latitude"
    print "Exiting"
    sys.exit()
try:
    lon_data = h5[lon_name][:]
except:
    print lon_name+" DNE in "+var_file
    print "Check that the variable name includes any necessary group paths. Ex: SoundingGeometry/sounding_longitude"
    print "Exiting"
    sys.exit()
h5.close()

if lite:
    
    qf_file_tag = "_all_quality"

    lite_file = LiteFile(var_file)
    lite_file.open_file()
    lite_lat = lite_file.get_lat()
    lite_lon = lite_file.get_lon()
    lite_sid = lite_file.get_sid()
    lite_xco2 = lite_file.get_xco2()
    lite_warn = lite_file.get_warn()
    lite_qf = lite_file.get_qf()
    lite_orbit = lite_file.get_orbit()
    lite_file.close_file()


    orbit_subset = np.where(lite_orbit == orbit_int)
    lite_lat = lite_lat[orbit_subset]
    lite_lon = lite_lon[orbit_subset]
    lite_sid = lite_sid[orbit_subset]
    lite_qf = lite_qf[orbit_subset]
    lite_xco2 = lite_xco2[orbit_subset]
    lite_warn = lite_warn[orbit_subset]
    
    lite_lat_subset_mask = set(np.where(np.logical_and(lite_lat <= maxy, lite_lat >= miny))[0])
    lite_lon_subset_mask = set(np.where(np.logical_and(lite_lon <= maxx, lite_lon >= minx))[0])
	
    lite_latlon_subset_mask = list(lite_lat_subset_mask.intersection(lite_lon_subset_mask))
        
    lite_lat = lite_lat[lite_latlon_subset_mask]
    lite_lon = lite_lon[lite_latlon_subset_mask]
    lite_sid = lite_sid[lite_latlon_subset_mask]
    lite_qf = lite_qf[lite_latlon_subset_mask]
    lite_xco2 = lite_xco2[lite_latlon_subset_mask]
    lite_warn = lite_warn[lite_latlon_subset_mask]

    print "Number of Lite soundings:", len(lite_sid)
    
    
    if lite_quality == 'good':
    
        quality_mask = np.where(lite_qf == 0)
	
	lite_lat = lite_lat[quality_mask]
	lite_lon = lite_lon[quality_mask]
	lite_xco2 = lite_xco2[quality_mask]
	lite_warn = lite_warn[quality_mask]
	
	qf_file_tag = "_good_quality"
	
    if lite_quality == 'bad':
    
        quality_mask = np.where(lite_qf == 1)
	
	lite_lat = lite_lat[quality_mask]
	lite_lon = lite_lon[quality_mask]
	lite_xco2 = lite_xco2[quality_mask]
	lite_warn = lite_warn[quality_mask]
	
	qf_file_tag = "_bad_quality"
	
    warn_mask = np.where(np.logical_and(lite_warn <= lite_warn_lims[1], lite_warn >= lite_warn_lims[0]))[0]
    
    lite_lat = lite_lat[warn_mask]
    lite_lon = lite_lon[warn_mask]
    lite_xco2 = lite_xco2[warn_mask]
    lite_warn = lite_warn[warn_mask]

    wl_file_tag = "_WL_"+str(lite_warn_lims[0])+"to"+str(lite_warn_lims[1])


### Plot prep ###
states_provinces = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='50m',
    facecolor='none')

outfile = output_dir+'/'+region+"_"+straight_up_date+qf_file_tag+wl_file_tag+".png"

### Plot the image ###
fig = plt.figure(figsize=(5,10))

img = plt.imread(code_dir+'/intermediate_RGB.tif')
img_extent = (minx, maxx, miny, maxy)

ax = plt.axes(projection=ccrs.PlateCarree())
ax.imshow(img, origin='upper', transform=ccrs.PlateCarree(), extent=img_extent)
ax.coastlines(resolution='10m', color='black', linewidth=1)
ax.add_feature(states_provinces, edgecolor='black', linewidth=1)
ax.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=1)
if interest_pt:
    ax.plot(interest_pt[1], interest_pt[0], 'w*', markersize=10, transform=ccrs.Geodetic())

ax.scatter(lite_lon, lite_lat, c=lite_xco2, cmap='jet', edgecolor='none', s=2, vmax=var_lims[1], vmin=var_lims[0])

cb_ax1 = fig.add_axes([.85, .3, .04, .4])
norm = mpl.colors.Normalize(vmin = var_lims[0], vmax = var_lims[1])
cb1 = mpl.colorbar.ColorbarBase(cb_ax1, cmap='jet', orientation = 'vertical', norm = norm)
cb1_lab = cb1.ax.set_xlabel('XCO2\n(ppm)')
cb1_lab.set_fontsize(9)
cb1.ax.xaxis.set_label_position("top")
cb1.ax.tick_params(labelsize=8)

#print "\nSaving figure. You may see a warning here but it does not affect anything"
fig.savefig(outfile, dpi=150, bbox_inches='tight')
print "\nFigure saved at "+outfile
os.remove(code_dir+'/intermediate_RGB.tif')
print "All Done!"
