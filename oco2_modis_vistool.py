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
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec

from osgeo import gdal, osr
import shapely.geometry as sgeom

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
    geo_upper_left, geo_lower_right, date,
    var_lat, var_lon, var_vals, var_vals_missing=None, lite_sid=np.empty([]),
    orbit_start_idx=0, var_lims=None, interest_pt=None, 
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
        
        if var_lat.ndim == 2:
            
            N = var_lat.shape[1]
                
            #Ensure that if any of the vertices are outside the lat/lon limits, all are masked
            latlon_subset_mask_1d = np.all(latlon_subset_mask, axis=1)
            latlon_subset_mask_2d = np.dstack([latlon_subset_mask_1d] * N)[0,:,:]
            
            var_lon_subset = np.ma.masked_where(latlon_subset_mask_2d == False, var_lon)
            var_lat_subset = np.ma.masked_where(latlon_subset_mask_2d == False, var_lat)
            if var_vals.ndim == 2:
                var_vals_subset = np.ma.masked_where(latlon_subset_mask_2d == False, var_vals)
            else: 
                var_vals_subset = np.ma.masked_where(latlon_subset_mask_1d == False, var_vals)
            if lite_sid.shape:
                lite_sid_subset = np.ma.masked_where(latlon_subset_mask_1d == False, lite_sid)
            
            if var_lon_subset.count() == 0 or var_lat_subset.count() == 0:
                print("\nWARNING: The lat/lon ranges given have no common points for the OCO-2 ground track")
                try:
                    lat_subset_idx = set(np.where(np.logical_and(var_lat <= maxy, var_lat >= miny))[0])
                    lon_subset_idx = set(np.where(np.logical_and(var_lon <= maxx, var_lon >= minx))[0])
                    latlon_subset_idx = list(lat_subset_idx.intersection(lon_subset_idx))
                    #print("Indices where the latitude is between " + str(miny) + " and " + str(maxy) +": " + str(min(lat_subset_idx)) + "-" + str(max(lat_subset_idx)))
                    print("Along-track indices where the longitude is between " + str(minx) + " and " + str(maxx) +": " + str(min(lon_subset_idx)+orbit_start_idx) + "-" + str(max(lon_subset_idx)+orbit_start_idx))
                    print("Latitude range for those indices: " + str(min(var_lat[min(lon_subset_idx)])) + "-" + str(min(var_lat[max(lon_subset_idx)])))
                    print("Latitude range given: " + str(miny) + "-" + str(maxy))
                    print("Along-track indices of intersection:", latlon_subset_idx)
                    #print("Exiting")
                    #os.remove(code_dir+'/intermediate_RGB.tif')
                    #sys.exit()
                except:
                    pass
                out_data = False
                var_vals = np.empty([])
             
        else:
            var_lon_subset = var_lon[latlon_subset_mask]
            var_lat_subset = var_lat[latlon_subset_mask]
            var_vals_subset = var_vals[latlon_subset_mask]
            if lite_sid.shape:
                lite_sid_subset = lite_sid[latlon_subset_mask]
                
        zip_it = np.ma.dstack([var_lon_subset, var_lat_subset])
        
    else:
            var_lon_subset = var_lon[latlon_subset_mask]
            var_lat_subset = var_lat[latlon_subset_mask]
            var_vals_subset = var_vals[latlon_subset_mask]
            if lite_sid.shape:
                lite_sid_subset = lite_sid[latlon_subset_mask]
                
            if var_lon_subset.size == 0 or var_lat_subset.size == 0:
                print("\nWARNING: The lat/lon ranges given have no common points for the OCO-2 ground track")
                try:
                    lat_subset_idx = set(np.where(np.logical_and(var_lat <= maxy, var_lat >= miny))[0])
                    lon_subset_idx = set(np.where(np.logical_and(var_lon <= maxx, var_lon >= minx))[0])
                    latlon_subset_idx = list(lat_subset_idx.intersection(lon_subset_idx))
                    #print("Indices where the latitude is between " + str(miny) + " and " + str(maxy) +": " + str(min(lat_subset_idx)) + "-" + str(max(lat_subset_idx)))
                    print("Indices where the longitude is between " + str(minx) + " and " + str(maxx) +": " + str(min(lon_subset_idx)+orbit_start_idx) + "-" + str(max(lon_subset_idx)+orbit_start_idx))
                    print("Latitude range for those indices: " + str(var_lat[min(lon_subset_idx)]) + "-" + str(var_lat[max(lon_subset_idx)]))
                    print("Latitude range given: " + str(miny) + "-" + str(maxy))
                    print("Indices of intersection:", latlon_subset_idx)
                    #print("Exiting")
                    #os.remove(code_dir+'/intermediate_RGB.tif')
                    #sys.exit()
                except:
                    pass
                out_data = False
                var_vals = np.empty([])
    
    min_segment_lat = np.ma.min(var_lat_subset)
    max_segment_lat = np.ma.max(var_lat_subset)
    min_segment_lon = np.ma.min(var_lon_subset)
    max_segment_lon = np.ma.max(var_lon_subset)
    
    print min_segment_lat, max_segment_lat
    print min_segment_lon, max_segment_lon
    
    if var_vals_missing:
        var_vals_subset = np.ma.masked_where(var_vals_subset == var_vals_missing, var_vals_subset)
        if var_vals_subset.count() == 0:
            print("\nThere is no valid data for the given subset criteria.")
            out_data = False
            var_vals = np.empty([])
    
    if var_lims is None:
        var_lims = [var_vals_subset.min(), var_vals_subset.max()]
    
    ### Write data to hdf5 file ###
    
    if out_data:
        
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
        
        if re.search("Masked", str(type(var_lat_subset))):
            if var_vals_subset.ndim > 1:
                var_data_to_write = np.ma.compress_rows(var_vals_subset)
            else:
                var_data_to_write = var_vals_subset.compressed()
        else:
            var_data_to_write = var_vals_subset
        
        outfile = h5py.File(out_data, "w")
        if not lat_name:
            lat_name = "Latitude"
        if not lon_name:
            lon_name = "Longitude"
        if not var_name:
            var_name = "Data"
        
        write_ds = outfile.create_dataset(lat_name, data = lat_data_to_write)
        write_ds = outfile.create_dataset(lon_name, data = lon_data_to_write)
        write_ds = outfile.create_dataset(var_name, data = var_data_to_write)
        if lite_sid.shape:
            if re.search("Masked", str(type(lite_sid_subset))):
                write_ds = outfile.create_dataset("sounding_id", data = lite_sid_subset.compressed())
            else:
                write_ds = outfile.create_dataset("sounding_id", data = lite_sid_subset)
        outfile.close()
        print("\nData saved at "+out_data)

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

    ax_minlat.set_title(str(ylocs[0]), horizontalalignment='left', verticalalignment='bottom', fontsize=10, fontweight='bold')
    ax_maxlat.set_title(str(ylocs[-1]), horizontalalignment='left', verticalalignment='top', fontsize=10, fontweight='bold')
    ax_minlon.set_title(str(xlocs[0]), horizontalalignment='center', verticalalignment='top', fontsize=10, fontweight='bold')
    ax_maxlon.set_title(str(xlocs[-1]), horizontalalignment='right', verticalalignment='top', fontsize=10, fontweight='bold')
    
    if var_vals.shape:
    
        patches = []

        if color_or_cmap == "cmap":
            if var_lat.ndim == 2 and var_vals.ndim == 1: 
                for row in xrange(zip_it.shape[0]):
                    polygon = mpatches.Polygon(zip_it[row,:,:]) 
                    patches.append(polygon)
                    
                p = mpl.collections.PatchCollection(patches, cmap=cmap, alpha=alpha, edgecolor='none')
                p.set_array(var_vals_subset)
                p.set_clim(var_lims[0], var_lims[1])
                ax.add_collection(p)

            else:
                ax.scatter(var_lon_subset, var_lat_subset, c=var_vals_subset, 
                           cmap=cmap, edgecolor='none', s=2, vmax=var_lims[1], vmin=var_lims[0])

            cb_ax1 = plt.subplot(gs[0:-1, -1])
            norm = mpl.colors.Normalize(vmin = var_lims[0], vmax = var_lims[1])
            cb1 = mpl.colorbar.ColorbarBase(cb_ax1, cmap=cmap, orientation = 'vertical', norm = norm)
            if var_label:
                cb1_lab = cb1.ax.set_xlabel(var_label, labelpad=8, fontweight='bold')
                cb1_lab.set_fontsize(14)
                cb1.ax.xaxis.set_label_position("top")
            for t in cb1.ax.yaxis.get_ticklabels():
                t.set_weight("bold")
                t.set_fontsize(12)
        if color_or_cmap == "color":
            for row in xrange(zip_it.shape[0]):
                polygon = mpatches.Polygon(zip_it[row,:,:], color=cmap) 
                patches.append(polygon)
            p = mpl.collections.PatchCollection(patches, alpha=alpha, edgecolor='none', match_original=True)
            ax.add_collection(p)
    
    #ax.plot([min_segment_lon, max_segment_lat], [max_segment_lon, min_segment_lat], c='white', transform=ccrs.Geodetic()) 
    
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
orbit_start_idx = 0
orbit_end_idx = 0

if __name__ == "__main__":

    ### Dynamic Definitions: get information from config file ###

    parser = argparse.ArgumentParser(description="Get configuration file")
    parser.add_argument('config_file_loc', type=str, default=code_dir+'/oco2_modis_vistool_config.json', nargs='?',
                        help="Name of config file (default: oco2_modis_vistool_config.json in code directory)")
    args = parser.parse_args()
    config_file = ConfigFile(args.config_file_loc)

    if config_file.exists():
        orbit_info_dict = config_file.get_contents()
        for k in orbit_info_dict.keys():
            orbit_info_dict[k.lower()] = orbit_info_dict.pop(k)
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
    try:
        overlay_info_dict = orbit_info_dict['oco2_overlay_info']
    except:
        overlay_info_dict = {}
    if overlay_info_dict:
        for k in overlay_info_dict.keys():
            overlay_info_dict[k.lower()] = overlay_info_dict.pop(k)
        var_file = overlay_info_dict['file']
        if not glob(var_file):
            print(var_file+" does not exist.")
            print("Exiting")
            sys.exit()
        var_name = overlay_info_dict['variable']
        var_plot_name = re.split('/', var_name)[-1]
        try: 
            band_number = overlay_info_dict['band_number']
            band_number = int(band_number)
        except: 
            band_number = False
        if var_name == "Retrieval/reduced_chi_squared_per_band" and band_number not in [1, 2, 3, False]:
                print("Unexpected band number specification. Options are 1 (0.76 micron), 2 (1.6 micron), or 3 (2.04 micron)")
                print(var_name + " is stored per band. Please select an appropriate band number. Exiting")
                sys.exit()

        if 'variable_plot_lims' in overlay_info_dict:
            var_lims = overlay_info_dict['variable_plot_lims']
            if var_lims:
                # if this was a string option, then make sure it is one 
                # of the two options we have implemented. if not, revert to default.
                if type(var_lims) in (str, unicode):
                    if var_lims not in ('autoscale_by_orbit', 'autoscale_by_overlay'):
                        print('var_lims "' + var_lims + '" is not valid, reverting to "autoscale_by_orbit"')
                        var_lims = 'autoscale_by_orbit'
            else:
                # empty list, tuple, None, etc 
                var_lims = 'autoscale_by_orbit'
        else:
            var_lims = 'autoscale_by_orbit'

        lat_name = overlay_info_dict['lat_name']
        lon_name = overlay_info_dict['lon_name']
        try: 
            orbit_int = overlay_info_dict['orbit']
        except: 
            orbit_int = False

        if re.search('oco2_Lt', var_file):
            lite = True
            #print("\nLite overlay file detected. Checking for QF and Warn specs...")
            
            if re.search("CO2", os.path.basename(var_file)):
                sif_or_co2 = "CO2"
                qf_file_tag = "_all_quality"
            elif re.search("SIF", os.path.basename(var_file)):
                sif_or_co2 = "SIF"
                qf_file_tag = ""
                wl_file_tag = ""
            else:
                print("The command line usage of this tool accommodates official CO2 or SIF lite files only.") 
                print("Expected filename convention: oco2_LtNNN_YYYYMMDD_B8xxxxx_rxx*.nc, where NNN is CO2 or SIF")
                print("Other data will need to be plotted by importing the tool as a Python module.")
                print("Exiting")
                sys.exit()

            try:
                lite_quality = overlay_info_dict['lite_qf']
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
            
            try:
                footprint_lims = overlay_info_dict['footprint']
            except:
                print("No footprint specifications detected. Output plot will contain all footprints")
                footprint_lims = [1, 8]
            if not footprint_lims:
                print("No footprint specifications detected. Output plot will contain all footprints")
                footprint_lims = [1, 8]
            if footprint_lims == "all":
                footprint_lims = [1, 8]
            if len(footprint_lims) == 2:
                if footprint_lims[0] > footprint_lims[1]:
                    print("Lower footprint limit is greater than upper footprint limit.")
                    print("Exiting")
                    sys.exit()
            for ft in footprint_lims:
                if ft not in np.arange(1, 9):
                    print("Unexpected footprint specification. Limits must be within [1, 8].")
                    print("Exiting")
                    sys.exit()
            
            if sif_or_co2 == "CO2":
                print("Output plot will include "+lite_quality+" quality soundings with warn levels within "+str(lite_warn_lims)+"\n")

        try:
            cmap = overlay_info_dict['cmap']
        except:
            cmap = ""
        if not cmap:
            cmap = "jet"

        try:
            alpha = float(overlay_info_dict['transparency'])
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
        out_plot_dir = orbit_info_dict['out_plot_dir']
    except:
        out_plot_dir = code_dir
    if not out_plot_dir or not glob(out_plot_dir):
        print("Either there was no output location specified or the one specified does not exist. Output will go in the code directory. \n")
        out_plot_dir = code_dir
    try:
        out_plot_name = orbit_info_dict['out_plot_name']
    except:
        out_plot_name = ""
    
    
    try:
        out_data_dir = orbit_info_dict['out_data_dir']
    except:
        out_data_dir= code_dir
    if not out_data_dir or not glob(out_data_dir):
        print("Either there was no output location specified or the one specified does not exist. Output will go in the code directory. \n")
        out_data_dir = code_dir
    try:
        out_data_name = orbit_info_dict['out_data_name']
    except:
        out_data_name = ""
    
    
    if not overlay_info_dict:
        
        if not out_plot_name:
            if region:
                out_plot_name = "MODISimagery_"+region+"_"+straight_up_date+".png"
            else:
                out_plot_name = "MODISimagery_"+straight_up_date+".png"
        out_plot_name = os.path.join(out_plot_dir, out_plot_name)

        do_modis_overlay_plot(orbit_info_dict['geo_upper_left'], 
                              orbit_info_dict['geo_lower_right'],
                              date, np.array([]), np.empty([]), np.empty([]), np.empty([]),
                              interest_pt=interest_pt, cmap='black',
                              out_plot=out_plot_name, cities=cities)
                                  
        sys.exit()
  
    

    ### Prep OCO-2 Variable ###

    h5 = h5py.File(var_file)
    if var_name:
        try:
            oco2_data_obj = h5[var_name]
        except:
            print(var_name+" DNE in "+var_file)
            print("Check that the variable name includes any necessary group paths. Ex: /Preprocessors/dp_abp")
            print("Exiting")
            sys.exit()
        oco2_data = h5[var_name][:]
    else:
        oco2_data_obj = h5['xco2']
        oco2_data = np.ones_like(oco2_data_obj[:])
    if sif_or_co2 == "CO2":
        oco2_data_long_name = oco2_data_obj.attrs.get('long_name')[0]
        oco2_data_units = oco2_data_obj.attrs.get('units')[0]
        oco2_data_fill = oco2_data_obj.attrs.get('missing_value')[0]
    if sif_or_co2 == "SIF":
        oco2_data_long_name = re.split("/", var_name)[-1]
        try:
            oco2_data_units = oco2_data_obj.attrs.get('unit').decode('utf-8')
        except:
            oco2_data_units = None
        try:
            oco2_data_fill = oco2_data_obj.attrs.get('missing_value')[0]
        except:
            oco2_data_fill = None
      
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
    
    if lat_data.ndim != lon_data.ndim:
        print(lat_name+" and "+lon_name+" have different dimensions. Exiting")
        sys.exit()
    
    if var_name == "Retrieval/reduced_chi_squared_per_band":
        if not band_number:
            print(var_name + " is stored per band. Please select a band number (1=0.76 micron, 2=1.6 micron, 3=2.04 micron). Exiting")
            sys.exit()
        else:
            oco2_data = oco2_data[:,band_number-1]
    
    if lite:
        
        if sif_or_co2 == "CO2":
            lite_file = LiteCO2File(var_file)
            lite_file.open_file()
            lite_warn = lite_file.get_warn()
            lite_qf = lite_file.get_qf()
        else:
            lite_file = LiteSIFFile(var_file)
            lite_file.open_file()
        lite_sid = lite_file.get_sid()
        lite_orbit = lite_file.get_orbit()
        lite_footprint = lite_file.get_footprint()
        lite_file.close_file()
        
        if orbit_int:
            orbit_subset = np.where(lite_orbit == orbit_int)
            orbit_start_idx = orbit_subset[0][0]
            orbit_end_idx = orbit_subset[0][-1]
            if sif_or_co2 == "CO2":
                lite_qf = lite_qf[orbit_subset]
                lite_warn = lite_warn[orbit_subset]
            lite_sid = lite_sid[orbit_subset]
            oco2_data = oco2_data[orbit_subset]
            lite_footprint = lite_footprint[orbit_subset]   
            if lat_data.ndim == 2:
                lat_data = np.squeeze(lat_data[orbit_subset, :])
                lon_data = np.squeeze(lon_data[orbit_subset, :])
            else:
                lat_data = lat_data[orbit_subset]
                lon_data = lon_data[orbit_subset]
            
        if sif_or_co2 == "CO2":
            if lite_quality == 'good':

                quality_mask = np.where(lite_qf == 0)
                qf_file_tag = "_good_quality"
                
                lite_qf = lite_qf[quality_mask]
                lite_warn = lite_warn[quality_mask]
                lite_sid = lite_sid[quality_mask]
                oco2_data = oco2_data[quality_mask]
                lite_footprint = lite_footprint[quality_mask]
                if lat_data.ndim == 2:
                    lat_data = np.squeeze(lat_data[quality_mask, :])
                    lon_data = np.squeeze(lon_data[quality_mask, :])
                else:
                    lat_data = lat_data[quality_mask]
                    lon_data = lon_data[quality_mask]

            if lite_quality == 'bad':

                quality_mask = np.where(lite_qf == 1)
                qf_file_tag = "_bad_quality"
                
                lite_qf = lite_qf[quality_mask]
                lite_warn = lite_warn[quality_mask]
                lite_sid = lite_sid[quality_mask]
                oco2_data = oco2_data[quality_mask]
                lite_footprint = lite_footprint[quality_mask]
                if lat_data.ndim == 2:
                    lat_data = np.squeeze(lat_data[quality_mask, :])
                    lon_data = np.squeeze(lon_data[quality_mask, :])
                else:
                    lat_data = lat_data[quality_mask]
                    lon_data = lon_data[quality_mask]

            warn_mask = np.where(np.logical_and(lite_warn <= lite_warn_lims[1], lite_warn >= lite_warn_lims[0]))[0]
            lite_qf = lite_qf[warn_mask]
            lite_warn = lite_warn[warn_mask]
            lite_sid = lite_sid[warn_mask]
            oco2_data = oco2_data[warn_mask]
            lite_footprint = lite_footprint[warn_mask]
            if lat_data.ndim == 2:
                lat_data = np.squeeze(lat_data[warn_mask, :])
                lon_data = np.squeeze(lon_data[warn_mask, :])
            else:
                lat_data = lat_data[warn_mask]
                lon_data = lon_data[warn_mask]

            wl_file_tag = "_WL_"+str(lite_warn_lims[0])+"to"+str(lite_warn_lims[1])
        
        if len(footprint_lims) == 2:
            footprint_mask = np.where(np.logical_and(lite_footprint <= footprint_lims[1], lite_footprint >= footprint_lims[0]))[0]
            fp_file_tag = "_FP_"+str(footprint_lims[0])+"to"+str(footprint_lims[1]) 
        else:
            footprint_mask = np.where(lite_footprint == footprint_lims)
            fp_file_tag = "_FP_"+str(footprint_lims[0])

        if sif_or_co2 == "CO2":
            lite_qf = lite_warn[footprint_mask]
            lite_warn = lite_warn[footprint_mask]
        lite_sid = lite_sid[footprint_mask]
        oco2_data = oco2_data[footprint_mask]
        lite_footprint = lite_footprint[footprint_mask]
        if lat_data.ndim == 2:
            lat_data = np.squeeze(lat_data[footprint_mask, :])
            lon_data = np.squeeze(lon_data[footprint_mask, :])
        else:
            lat_data = lat_data[footprint_mask]
            lon_data = lon_data[footprint_mask]   
    
    # here, handle the var limit options. 
    # if specific limits were input, via a 2-element list, 
    # that will be passed directly to do_modis_overlay_plot().
    # if autoscaling by orbit, find those min/max now, and create the 2 element list.
    if var_lims == 'autoscale_by_orbit':
        var_lims = [np.min(oco2_data), np.max(oco2_data)]
    # otherwise, convert to None, so then do_modis_overlay_plot() will then derive an autoscale 
    # min/max according to the points within the overlay.
    if var_lims == 'autoscale_by_overlay':
        var_lims = None

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

    if not out_plot_name:
        if region:
            out_plot_name = var_plot_name+"_"+region+"_"+straight_up_date+qf_file_tag+wl_file_tag+fp_file_tag+".png"
        else:
            out_plot_name = var_plot_name+"_"+straight_up_date+qf_file_tag+wl_file_tag+fp_file_tag+".png"
    out_plot_name = os.path.join(out_plot_dir, out_plot_name)
    
    if not out_data_name:
        if region:
            out_data_name = var_plot_name+"_"+region+"_"+straight_up_date+qf_file_tag+wl_file_tag+fp_file_tag+".h5"
        else:
            out_data_name = var_plot_name+"_"+straight_up_date+qf_file_tag+wl_file_tag+fp_file_tag+".h5"
    out_data_name = os.path.join(out_data_dir, out_data_name)
    
    do_modis_overlay_plot(orbit_info_dict['geo_upper_left'],
                          orbit_info_dict['geo_lower_right'], 
                          date, lat_data, lon_data, oco2_data, oco2_data_fill, lite_sid,
                          orbit_start_idx, var_lims=var_lims, interest_pt=interest_pt, 
                          cmap=cmap, alpha=alpha,lat_name=lat_name, lon_name=lon_name, var_name=var_name,
                          out_plot=out_plot_name, out_data=out_data_name, var_label=cbar_name, cities=cities)
