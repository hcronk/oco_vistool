# oco_vistool.py

The purpose of this script is to pull different RGB images (depending on the user's choice from Encoding.csv) from Worldview using the NASA GIBS API and overlay various OCO-2 data fields for 
case study analysis in support of OCO-2 cloud and aerosol screening validation.

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

    ##### Python requirements:
    Python 2.7 or higher
      
    For those new to Python, it is easiest to get the Anaconda distribution: 
    https://www.continuum.io/download and then use the Anaconda utility conda 
    to fetch and install packages.


    The conda installs required on top of what comes with Anaconda 5.2 for both Python 2.7 and Python 3.6 are as follows:
    + conda install -c conda-forge cartopy
    + conda install -c conda-forge ncurses
    + conda install -c conda-forge hdf5
    + conda install -c conda-forge future

    #### Option 2: Docker Container
    Dockerfiles and docker-compose files to set up an appropriate Python 2.7 or Python 3.6 container are included in the repo    

3) Test run from the command line  
   `python oco_vistool.py`


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
     **sensor** (string):
         The sensor name that is desired to be used for the background image. Supported sensors at the moment are: 
         "Worldview", "GOES16_ABI_C", "GOES16_ABI_F", "GOES17_ABI_C", "GOES17_ABI_F".
   #### Required if the sensor field, required above, is "Worldview" (otherwise, not needed):
     **layer** (integer, format: integer between 0 and maximum code number in Encoding.csv):
         The code of the background layer to be used by the vistool. Each layer is given its code (encoded) 
         in the Encoding.csv file in the code directory. The file contains the relation: layer name (by the
         NASA Worldview standards) - code (generated unique number). So, when adding new
         layers to the CSV, use the above standards.

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
         The information for the OCO-2 variable to overlay on the Background image.  
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

	**out_plot_name** (string, default example: MODIS_Aqua_CorrectedReflectance_TrueColor_xco2_Amazon_20150701_B7302b_QF_good_WL_0to20.png):  
    	The name of the output plot.

	**out_data_dir** (string, default: the directory that you are running the code in):  
    	The full output directory path for the output data file.

	**out_data_name** (string, default example: xco2_Amazon_20150701_B7302b_QF_good_WL_0to20.h5):  
    	The name of the output data file.

	**out_background_name** (string, default example: MODIS_Aqua_CorrectedReflectance_TrueColor_Amazon_20150701.png):  
    	The name of the output background file, created if no overlay data is present or make_background_image is enabled.


## Other Information

GIBS developer documentation:  https://wiki.earthdata.nasa.gov/display/GIBS/GIBS+API+for+Developers

#### WMTS Capabilities File and Its Parser
There is a script in the code directory called "parse_wmts.py" which parses the wmts.xml capabilities file. The
reason of this Python parser is to get the relation, layer name - layer specifications XML url, for all NASA Worldview
layers (as of June 2021) into a separate CSV file. This relation CSV file is used for plotting and styling the output image.
The parser used a copy of the wmts capabilities file which was downloaded into the code directory from: https://gibs.earthdata.nasa.gov/wmts/epsg4326/best/wmts.cgi?SERVICE=WMTS&REQUEST=GetCapabilities

Minimum command line call:  
`python oco_vistool.py`


Minimum function call:  
```
from oco_vistool import do_overlay_plot

do_overlay_plot([lat_of_upper_lefthand_corner, lon_of_upper_lefthand_corner],
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
