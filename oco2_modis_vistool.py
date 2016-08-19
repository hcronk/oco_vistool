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

import os
import sys
from glob import glob

from OCO2FileOps import *

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
    
    cmd = "gdal_translate -of GTiff -outsize "+str(xsize)+" "+str(ysize)+" -projwin "+str(lon_ul)+" "+str(lat_ul)+" "+str(lon_lr)+" "+str(lat_lr)+" "+xml_file+" "+tif_file

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
lat_ul = orbit_info_dict['geo_upper_left'][0]
lon_ul = orbit_info_dict['geo_upper_left'][1]
lat_lr = orbit_info_dict['geo_lower_right'][0]
lon_lr = orbit_info_dict['geo_lower_right'][1]
region = orbit_info_dict['region']
try:
    interest_pt = orbit_info_dict['ground_site']
except:
    interest_pt = []
try:
    output_dir = orbit_info_dict['output_dir']
except:
    output_dir = code_dir
if not output_dir or not glob(output_dir):
    print "output dir DNE, using code dir instead"
    output_dir = code_dir

### Pull Aqua-MODIS RGB from GIBS ###
update_GIBS_xml(date, xml_file)

try:
    pull_Aqua_RGB_GIBS(lat_ul, lon_ul, lat_lr, lon_lr,  xml_file, code_dir+'/intermediate_RGB.tif')
except:
    print "Problem pulling RGB. Check that the geolocation bounds specified in the configuration file are for the upper left hand corner and the lower right hand corner"


 




