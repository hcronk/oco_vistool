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

Authors:
Heather Cronk <heather.cronk@colostate.edu>
Natalie Tourville <natalie.tourville@colostate.edu>

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
import math
import pandas as pd

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import shapefile
from shapely.geometry import LineString, Point, Polygon
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.path as mpath
import matplotlib.patches as mpatches

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

def do_modis_overlay_plot(
    geo_lower_right, geo_upper_left, date, 
    var_lat, var_lon, var_vals, var_lims=None, interest_pt=None, 
    cmap='jet', alpha=1,
    outfile=None, var_label=None, cities=None):

    if var_lims is None:
        var_lims = [var_vals.min(), var_vals.max()]

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
    
    if cmap in plt.cm.datad.keys():
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
    
    if var_vals.shape:
	# get subset of var values etc.
	latlon_subset_mask = np.logical_and(
            np.logical_and(var_lat <= maxy, var_lat >= miny), 
            np.logical_and(var_lon <= maxx, var_lon >= minx) )
        
	#Ensure that if any of the vertices are outside the lat/lon limits, all are masked
	latlon_subset_mask_1d = np.all(latlon_subset_mask, axis=1)
	latlon_subset_mask_2d = np.dstack([latlon_subset_mask_1d, latlon_subset_mask_1d, latlon_subset_mask_1d, latlon_subset_mask_1d])[0,:,:]
    
	var_lon_subset = np.ma.masked_where(latlon_subset_mask_2d == False, var_lon)
	var_lat_subset = np.ma.masked_where(latlon_subset_mask_2d == False, var_lat)
	var_vals_subset = np.ma.masked_where(latlon_subset_mask_1d == False, var_vals)
	
	zip_it = np.ma.dstack([var_lon_subset, var_lat_subset])

	
#This is for 1D lat/lon. Update this to give the appropriate info for the 2D lat/lon arrays!
#	if var_lon_subset.size == 0 or var_lat_subset.size == 0:
#            lat_subset_idx = set(np.where(np.logical_and(var_lat <= maxy, var_lat >= miny))[0])
#	    lon_subset_idx = set(np.where(np.logical_and(var_lon <= maxx, var_lon >= minx))[0])
#	    latlon_subset_idx = list(lat_subset_idx.intersection(lon_subset_idx))
#            print("\nThe lat/lon ranges given have no common points for the OCO-2 ground track")
#	    #print("Indices where the latitude is between " + str(miny) + " and " + str(maxy) +": " + str(min(lat_subset_idx)) + "-" + str(max(lat_subset_idx)))
#            print("Indices where the longitude is between " + str(minx) + " and " + str(maxx) +": " + str(min(lon_subset_idx)) + "-" + str(max(lon_subset_idx)))
#	    print("Latitude range for those indices: " + str(var_lat[min(lon_subset_idx)]) + "-" + str(var_lat[max(lon_subset_idx)]))
#            print("Latitude range given: " + str(miny) + "-" + str(maxy))
#	    print("Indices of intersection:", latlon_subset_idx)
#	    print("Exiting")
#	    os.remove(code_dir+'/intermediate_RGB.tif')
#	    sys.exit()

    ### Plot prep ###
    states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none')
			  
				 
    ### Plot the image ###
    fig = plt.figure(figsize=(fig_x,fig_y))

    img = plt.imread(code_dir+'/intermediate_RGB.tif')
    img_extent = (minx, maxx, miny, maxy)

    ax = plt.axes(projection=ccrs.PlateCarree())
    ax_pos = ax.get_position()
    ax.imshow(img, origin='upper', transform=ccrs.PlateCarree(), extent=img_extent)
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
    
    if var_vals.shape:
	g1 = ax.gridlines(draw_labels=True, alpha = 0.5)
	g1.xlabels_top = False
	g1.ylabels_right = False
	g1.xlabel_style = {'size': 7}
	g1.ylabel_style = {'size': 7}
	g1.xformatter = LONGITUDE_FORMATTER
	g1.yformatter = LATITUDE_FORMATTER
    
        patches = []

	if color_or_cmap == "cmap":
	    for row in xrange(zip_it.shape[0]):
	        polygon = mpatches.Polygon(zip_it[row,:,:]) 
	        patches.append(polygon)
	
#	    ax.scatter(var_lon_subset, var_lat_subset, c=var_vals_subset, 
#        	       cmap=cmap, edgecolor='none', s=2, vmax=var_lims[1], vmin=var_lims[0])
		       
            p = mpl.collections.PatchCollection(patches, cmap=cmap, alpha=alpha, edgecolor='none')
	    p.set_array(var_vals_subset)
	    p.set_clim(var_lims[0], var_lims[1])
	    ax.add_collection(p)
	    
	    cb_ax1 = fig.add_axes([ax_pos.x1, ax_pos.y0 + .2, .04, .4])
	    norm = mpl.colors.Normalize(vmin = var_lims[0], vmax = var_lims[1])
	    cb1 = mpl.colorbar.ColorbarBase(cb_ax1, cmap=cmap, orientation = 'vertical', norm = norm)
	    cb1_lab = cb1.ax.set_xlabel(var_label, labelpad=8)
	    cb1_lab.set_fontsize(7)
	    cb1.ax.xaxis.set_label_position("top")
	    cb1.ax.tick_params(labelsize=6)
	if color_or_cmap == "color":
#    	    ax.scatter(var_lon_subset, var_lat_subset, c=cmap, edgecolor='none', s=2)
            for row in xrange(zip_it.shape[0]):
	        polygon = mpatches.Polygon(zip_it[row,:,:], color=cmap) 
	        patches.append(polygon)
	    p = mpl.collections.PatchCollection(patches, alpha=alpha, edgecolor='none', match_original=True)
	    ax.add_collection(p)
		

    fig.savefig(outfile, dpi=150)
    print("\nFigure saved at "+outfile)
    #print("code directory in subroutine:", code_dir)
    os.remove(code_dir+'/intermediate_RGB.tif')


### Static Defnitions

code_dir = os.path.dirname(os.path.realpath(__file__))
xml_file = code_dir+'/GIBS_Aqua_MODIS_truecolor.xml'

if __name__ == "__main__":

    ### Dynamic Definitions: get information from config file ###

    parser = argparse.ArgumentParser(description="Get configuration file")
    parser.add_argument('config_file_loc', type=str, default=code_dir+'/oco2_modis_vistool_config.json', nargs='?',
                	help="Name of config file (default: oco2_modis_vistool_config.json in code directory)")
    args = parser.parse_args()
    config_file = ConfigFile(args.config_file_loc)

    if config_file.exists():
	orbit_info_dict = config_file.get_contents()
    else:
	print('The expected configuration file '+ args.config_file_loc + ' DNE in ' + code_dir)
	print('Exiting')
	sys.exit()


    date = orbit_info_dict['date']
    straight_up_date = date.replace("-", "")
    lat_ul = orbit_info_dict['geo_upper_left'][0]
    lon_ul = orbit_info_dict['geo_upper_left'][1]
    lat_lr = orbit_info_dict['geo_lower_right'][0]
    lon_lr = orbit_info_dict['geo_lower_right'][1]
    try:    
	region = orbit_info_dict['region']
    except:
	region = ""
    overlay_info_dict = orbit_info_dict['oco2_overlay_info']
    if overlay_info_dict:  
	var_file = overlay_info_dict['file']
	if not glob(var_file):
            print(var_file+" does not exist.")
            print("Exiting")
            sys.exit()
	var_name = overlay_info_dict['variable']
	var_plot_name = re.split('/', var_name)[-1]
        try:
	    var_lims = overlay_info_dict['variable_plot_lims']
	except:
	    var_lims = []
	lat_name = overlay_info_dict['lat_name']
	lon_name = overlay_info_dict['lon_name']
	try: 
            orbit_int = overlay_info_dict['orbit']
	except: 
            orbit_int = False    

	if re.search('oco2_Lt', var_file):
            lite = True
            #print("\nLite overlay file detected. Checking for QF and Warn specs...")

            try:
        	lite_quality = overlay_info_dict['lite_QF']
            except:
        	print("No quality specifications detected. Output plot will contain all quality soundings")
        	lite_quality = 'all'
            if not lite_quality:
        	print("No quality specifications detected. Output plot will contain all quality soundings")
        	lite_quality = 'all'
            if lite_quality not in ['', 'all', 'good', 'bad']:
        	print("Unexpected quality flag specification. Options are: '', 'all', 'good', or 'bad'")
        	print("Exiting")
        	sys.exit()

            try:
        	lite_warn_lims = overlay_info_dict['lite_warn_lims']
            except:
        	print("No warn specifications detected. Output plot will contain all warn levels")
        	lite_warn_lims = [0, 20]
            if not lite_warn_lims:
        	print("No warn specifications detected. Output plot will contain all warn levels")
        	lite_warn_lims = [0, 20]
            if lite_warn_lims[0] > lite_warn_lims[1]:
        	print("Lower warn limit is greater than upper warn limit.")
        	print("Exiting")
        	sys.exit()
            for lim in lite_warn_lims:
        	if lim not in np.arange(21):
                    print("Unexpected warn level specification. Limits must be within [0, 20].")
                    print("Exiting")
                    sys.exit()

	print("Output plot will include "+lite_quality+" quality soundings with warn levels within "+str(lite_warn_lims)+"\n")

	try:
	    cmap = overlay_info_dict['cmap']
	except:
	    cmap = ""
	if not cmap:
	    cmap = "jet"
	
	try:
	    alpha = overlay_info_dict['transparency']
	except:
	    alpha = ""
	if not alpha:
	    alpha = 1
	if alpha < 0 or alpha > 1:
	    print("Unexpected transparency specification. Value must be within [0, 1]. Output plot will contain full-color.\n")
	    alpha = 1
    try:
        interest_pt = orbit_info_dict['ground_site']
	if not interest_pt:
	    interest_pt = None
	if (interest_pt[0] > lat_ul or interest_pt[0] < lat_lr or interest_pt[1] > lon_lr or interest_pt[1] < lon_ul):
	    interest_pt = None
	    print("The ground site is outside the given lat/lon range and will not be included in the output plot.\n")
    except:
        interest_pt = None
    try:
        cities = orbit_info_dict['city_labels']
	if cities and cities not in mpl.colors.cnames.keys():
	    print(cities + " is not an available matplotlib color. City labels will not be included on the output plot. \n")
	    cities = None
	if not cities:
	    cities = None
    except:
        cities = None
    try:
	output_dir = orbit_info_dict['output_dir']
    except:
	output_dir = code_dir
    if not output_dir or not glob(output_dir):
	print("Either there was no output location specified or the one specified does not exist. Output will go in the code directory. \n")
	output_dir = code_dir

    try:
	outfile_name = orbit_info_dict['outfile']
    except:
	outfile_name = ""
    
    if not overlay_info_dict:
        
        if not outfile_name:
	    if region:
        	outfile_name = "MODISimagery_"+region+"_"+straight_up_date+".png"
	    else:
        	outfile_name = "MODISimagery_"+straight_up_date+".png"
	outfile = output_dir+"/"+outfile_name
	
	do_modis_overlay_plot(orbit_info_dict['geo_lower_right'], 
                              orbit_info_dict['geo_upper_left'],
                              date, np.array([]), np.empty([]), np.empty([]), 
                              interest_pt=interest_pt, cmap='black',
			      outfile=outfile, cities=cities)
        sys.exit()
			  
    

    ### Prep OCO-2 Variable ###

    h5 = h5py.File(var_file)
    try:
	oco2_data_obj = h5[var_name]
	oco2_data = h5[var_name][:]
	oco2_data_long_name = oco2_data_obj.attrs.get('long_name')[0]
	oco2_data_units = oco2_data_obj.attrs.get('units')[0]
    except:
	print(var_name+" DNE in "+var_file)
	print("Check that the variable name includes any necessary group paths. Ex: /Preprocessors/dp_abp")
	print("Exiting")
	sys.exit()
    try:
	lat_data = h5[lat_name][:]
    except:
	print(lat_name+" DNE in "+var_file)
	print("Check that the variable name includes any necessary group paths. Ex: SoundingGeometry/sounding_latitude")
	print("Exiting")
	sys.exit()
    try:
	lon_data = h5[lon_name][:]
    except:
	print(lon_name+" DNE in "+var_file)
	print("Check that the variable name includes any necessary group paths. Ex: SoundingGeometry/sounding_longitude")
	print("Exiting")
	sys.exit()
    h5.close()

    if lite:

	qf_file_tag = "_all_quality"

	lite_file = LiteFile(var_file)
	lite_file.open_file()
	lite_lat = lite_file.get_lat()
	lite_lon = lite_file.get_lon()
	lite_vert_lat = lite_file.get_vertex_lat()
	lite_vert_lon = lite_file.get_vertex_lon()
	lite_sid = lite_file.get_sid()
	lite_xco2 = lite_file.get_xco2()
	lite_warn = lite_file.get_warn()
	lite_qf = lite_file.get_qf()
	lite_orbit = lite_file.get_orbit()
	lite_file.close_file()
	
        if orbit_int:
	    orbit_subset = np.where(lite_orbit == orbit_int)
	    lite_lat = lite_lat[orbit_subset]
	    lite_lon = lite_lon[orbit_subset]
	    lite_vert_lat = lite_vert_lat[orbit_subset, :]
	    lite_vert_lon = lite_vert_lon[orbit_subset, :]
	    lite_sid = lite_sid[orbit_subset]
	    lite_qf = lite_qf[orbit_subset]
	    lite_xco2 = lite_xco2[orbit_subset]
	    lite_warn = lite_warn[orbit_subset]
	    oco2_data = oco2_data[orbit_subset]
	    lat_data = lat_data[orbit_subset]
	    lon_data = lon_data[orbit_subset]
	    
        if lite_quality == 'good':
    
            quality_mask = np.where(lite_qf == 0)
	    qf_file_tag = "_good_quality"

	    lite_lat = lite_lat[quality_mask]
	    lite_lon = lite_lon[quality_mask]
	    lite_vert_lat = lite_vert_lat[0, quality_mask, :]
	    lite_vert_lon = lite_vert_lon[0, quality_mask, :]
	    lite_xco2 = lite_xco2[quality_mask]
	    lite_warn = lite_warn[quality_mask]
	    oco2_data = oco2_data[quality_mask]
	    lat_data = lat_data[quality_mask]
	    lon_data = lon_data[quality_mask]
	
        if lite_quality == 'bad':

            quality_mask = np.where(lite_qf == 1)
	    qf_file_tag = "_bad_quality"

	    lite_lat = lite_lat[quality_mask]
	    lite_lon = lite_lon[quality_mask]
	    lite_vert_lat = lite_vert_lat[0, quality_mask, :]
	    lite_vert_lon = lite_vert_lon[0, quality_mask, :]
	    lite_xco2 = lite_xco2[quality_mask]
	    lite_warn = lite_warn[quality_mask]
	    oco2_data = oco2_data[quality_mask]
	    lat_data = lat_data[quality_mask]
	    lon_data = lon_data[quality_mask]
	
	warn_mask = np.where(np.logical_and(lite_warn <= lite_warn_lims[1], lite_warn >= lite_warn_lims[0]))[0]

	lite_lat = lite_lat[warn_mask]
	lite_lon = lite_lon[warn_mask]
	lite_vert_lat = lite_vert_lat[0, warn_mask, :]
	lite_vert_lon = lite_vert_lon[0, warn_mask, :]
	lite_xco2 = lite_xco2[warn_mask]
	lite_warn = lite_warn[warn_mask]
	oco2_data = oco2_data[warn_mask]
	lat_data = lat_data[warn_mask]
	lon_data = lon_data[warn_mask]
	

	wl_file_tag = "_WL_"+str(lite_warn_lims[0])+"to"+str(lite_warn_lims[1])
	
    if not var_lims:
	var_lims = [np.min(oco2_data), np.max(oco2_data)]
        vmax = int(math.ceil(var_lims[1]))
        vmin = int(math.floor(var_lims[0]))
    
    vmax = var_lims[1]
    vmin = var_lims[0]

    ### Plot prep ###

    oco2_data_long_name = re.sub("_", " ", oco2_data_long_name)
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
    cbar_name = cbar_cap_strings+'\n('+oco2_data_units+')'

    if not outfile_name:
	if region:
            outfile_name = var_plot_name+"_"+region+"_"+straight_up_date+qf_file_tag+wl_file_tag+".png"
	else:
            outfile_name = var_plot_name+"_"+straight_up_date+qf_file_tag+wl_file_tag+".png"
    outfile = output_dir+"/"+outfile_name

#    do_modis_overlay_plot(orbit_info_dict['geo_lower_right'], 
#                          orbit_info_dict['geo_upper_left'],
#                          date, lat_data, lon_data, oco2_data, 
#                          var_lims=[vmin,vmax], interest_pt=interest_pt, cmap=cmap,
#			  outfile=outfile, var_label=cbar_name, cities=cities)

    do_modis_overlay_plot(orbit_info_dict['geo_lower_right'], 
                          orbit_info_dict['geo_upper_left'],
                          date, lite_vert_lat, lite_vert_lon, oco2_data, 
                          var_lims=[vmin,vmax], interest_pt=interest_pt, 
			  cmap=cmap, alpha=alpha,
			  outfile=outfile, var_label=cbar_name, cities=cities)
