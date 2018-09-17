# oco2_modis_vistool.py

The purpose of this script is to pull Aqua-MODIS RGB images from Worldview
using the NASA GIBS API and overlay various OCO-2 data fields for case study 
analysis in support of OCO-2 cloud and aerosol screening validation.

It can be used as a command line tool to plot any data in the OCO-2 CO2 or SIF Lite file, filtered by warn level, quality flag, and/or footprint number where applicable.  
It can also be called as a function from within another Python program for expanded use with other datasets.

## GETTING STARTED

1) Clone this code \
   `git clone http://github.com/hcronk/oco2_modis_vistool`
   
   We suggest obtaining the package via a `git clone` rather than downloading the tarball because
   the tool is still under active development. A clone will enable you to quickly update your copy 
   to the current stable version using a `git pull` command on your end, rather than downloading 
   and comparing files.  

2) Make sure you have the necessary system and Python requirements.
    #### Option 1: DIY   
    
    ##### System Requirements:
    GDAL version 1.9.1 or greater with cURL support enabled.  
        
    To check if cURL is enabled for GDAL, use the command `gdalinfo --format WMS`  
    If cURL is enabled you will see information about the WMS format, if not, 
    you will get an error message and you will need to reconfigure GDAL to support cURL.  
    **NOTE:** If this comes with your Anaconda distribution of Python, you may need
    to set an environmental variable called GDAL_DATA in your bash profile that points to
    gdal within your Anaconda folder.
    (Ex: `export GDAL_DATA="/local/home/hcronk/anaconda/share/gdal"`)

    ##### Python requirements:
    Python 2.7 or higher
      
    For those new to Python, it is easiest to get the Anaconda distribution: 
    https://www.continuum.io/download and then use the Anaconda utility conda 
    to fetch and install packages.


    The conda installs required on top of what comes with Anaconda 5.2 for both Python 2.7 and Python 3.6 are as follows:
    + conda install -c conda-forge cartopy
    + conda install -c conda-forge gdal=2.2.4
    + conda install -c conda-forge ncurses
    + conda install -c conda-forge hdf5
    + conda install -c conda-forge future

    #### Option 2: Docker Container
    Dockerfiles and docker-compose files to set up an appropriate Python 2.7 or Python 3.6 container are included in the repo    

3) Test run from the command line  
   `python oco2_modis_vistool.py`


4) Look for the output PNG and HDF5 in the code directory

5) Update the JSON configuration file with information from your scene of interest. Options include:
   #### Required:
     **date** (string, format:"YYYY-MM-DD"):  
         The date of your case study. This is used to get the correct RGB. 

     **geo_upper_left** (list, format:[lat, lon]):  
         The latitude and longitude coordinates of the upper left hand corner of your area of 
         interest.

     **geo_lower_right** (list, format:[lat, lon]):  
         The latitude and longitude coordinates of the lower right hand corner of your area of 
         interest.

   #### Optional:
	**region** (string, default: none):  
    	The region of your case study. This is used to build the default output filename in the 
        command line execution.

     **cmap** (string, default: "jet"):  
       The desired colormap or color for the overlay data.  
       Built-in colormap options can be found here:
       http://matplotlib.org/examples/color/colormaps_reference.html  
       Built-in single color options can be found here:
       http://matplotlib.org/examples/color/named_colors.html

     **oco2_overlay_info** (dict, default: none):  
         The information for the OCO-2 variable to overlay on the MODIS RGB.  
		**{**  
		**file** (string):  
        		This is *required* for a variable overlay.  
			The full path to the file that contains the variable to overlay.  
		**variable** (string, format:"Preprocessors/dp_adp"):  
        		This is *required* for a variable overlay.  
			The full path to the variable within the file specified.  
		**band_number** (string or integer, format: "1" or 1):  
			This is optional for variables that are stored per-band in the same array.  
		**variable_plot_lims** (list, format:[min, max], 
			default: [integer floor of variable's minimum value, integer ceiling of variable's 
            		maximum value]):  
            		The maximum and minimum values used for plotting and color assignment. Any data values 
            		beyond the min and max will show up on the plot as the minimum or maximum color, 
            		respectively.  
		**lat_name** (string, format:"latitude"):  
			This is *required* for a variable overlay for command line execution.  
            		The name of the latitude array associated with the variable. Used for subsetting and 
            		plotting.  
		**lon_name** (string, format:"longitude"):  
			This is *required* for a variable overlay for command line execution.  
			The name of the longitude array associated with the variable. Used for subsetting and 
            		plotting.  
		**orbit** (integer, format:11016):  
			This is optional for a variable overlay.  
            		The OCO-2 orbit number. Used for subsetting.  
		**lite_QF** (string, format:"good", default: "all"):  
			The quality flag applied to the CO2 Lite File variable. Used for filtering.  
            		Valid options: "good", "bad", "all"  
		**lite_warn_lims** (list, format:[0,15], default: [0,20]):  
			The range of acceptable warn levels applied to the CO2 Lite File variable. Used for 
            		filtering.  
            		Valid options: integers from 0 to 20  
		**footprint** (1 or 2 element list, format: [2:4] or [1], default: [1:8]):  
           		The footprint number applied to the CO2 or SIF Lite File variable. Used for filtering.  
		**transparency** (integer between 0 and 1, format: 0, default: 1):  
           		The level of transparency applied to the OCO-2 data where 1 is opaque and 0 is 
           		transparent.  
		**}**  
       
	**ground_site** (list, format:[lat lon], default: none):  
    	The latitude and longitude coordinates of a ground site of interest. It will show up as 
        a white star on the output plot.

	**city_labels** (string, default: none):  
    	This is optional to display city names on the output plot in the color specified in this 
        string. They are on a 10m scale, and work best if you are displaying a small area.  
        Python built-in named colors are available here: 
        http://matplotlib.org/examples/color/named_colors.html

	**out_plot_dir** (string, default: the directory that you are running the code in):  
    	The full output directory path for the output plot.

	**out_plot_name** (string, default example: xco2_Amazon_20150701_all_quality_WL_0to20.png):  
    	The name of the output plot.

	**out_data_dir** (string, default: the directory that you are running the code in):  
    	The full output directory path for the output data file.

	**out_data_name** (string, default example: xco2_Amazon_20150701_all_quality_WL_0to20.h5):  
    	The name of the output data file.


## Other Information

GIBS developer documentation:  https://wiki.earthdata.nasa.gov/display/GIBS/GIBS+API+for+Developers


Minimum command line call:  
`python oco2_modis_vistool.py`


Minimum function call:  
```
from oco2_modis_vistool import do_modis_overlay_plot

do_modis_overlay_plot([lat_of_upper_lefthand_corner, lon_of_upper_lefthand_corner],
                      [lat_of_lower_righthand_corner, lon_of_lower_righthand_corner], 
		      date, lat_data, lon_data, variable_data)
```  


## A Note about Questions, Suggestions, and Issues

For the benefit of the entire user community, it is requested that you submit questions, suggestions, and issues via the GitHub issue tracker so that everyone is able to see the work that is in progress. 
https://github.com/hcronk/oco2_modis_vistool/issues

---
Developed and maintained by:  
Heather Cronk  
heather.cronk@colostate.edu
