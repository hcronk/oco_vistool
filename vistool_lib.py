from datetime import datetime
import time

import numpy as np

import cartopy.crs as ccrs
from cartopy.io import shapereader
from shapefile import ShapefileException
from cartopy.mpl.geoaxes import GeoAxes

from pyproj import Proj

import shapely.geometry as sgeom
from shapely.strtree import STRtree

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1 import make_axes_locatable

import calendar

def get_image_size(latlon_extent, proj_string, img_xsize = 2000):
    """
    compute the image size using R Nelson method.
    This is mainly for testing, I don't think we will try to
    integrate this directly into the code.

    inputs:
    latlon_extent: 4 element array like with West, South, East, North
        values describing the lat/lon box being displayed.
    proj_string: the pyProj projection string; this should match the
        the project actually used in the data.
        ususally this should be "epsg:3857"
    img_xsize: number of x-dimension pixels in the image. the number
        of ypixels is derived based on the lat/lon box.

    returns:
    imgsize: a 2 element tuple with (xpixels, ypixels).
        the xpixels value will be equal to the img_xsize.
    """

    W, S, E, N = latlon_extent

    m1 = Proj(proj_string, preserve_units=True)
    W_xpos, S_ypos = m1(W,S)
    E_xpos, N_ypos = m1(E,N)
    xpixels = img_xsize
    ypixels = int((N_ypos - S_ypos) / (E_xpos - W_xpos) * xpixels)
    imgsize = xpixels, ypixels

    return imgsize


def get_overlay_boundingbox(odat, target_latlon, edge_width=0.25):
    """
    compute the image size using R Nelson method.
    This is mainly for testing, I don't think we will try to
    integrate this directly into the code.

    inputs:
    odat: the overlay data dictionary (see oco_vistool.py)
    target_latlon: the targeted latlon position, as a 2 element
        array-like. This is used to ensure the overlay bounding box
        includes both the overlay data and targeted position.
    edge_width: degrees latitude of "padding" to include in the
        bounding box. note the amount of longitude is automatically
        scaled by dividing the latitude edge width by the cosine of lat.

    returns a bounding box:
    [min_lon, min_lat, max_lon, max_lat]
    that encloses the overlay data

    There must be at least one valid data point within the overlay
    data (odat) for this to work.

    """

    target_lat, target_lon = target_latlon
    min_lat = min(ovr_d['lat'].min(), target_lat)
    min_lon = min(ovr_d['lon'].min(), target_lon)
    max_lat = max(ovr_d['lat'].max(), target_lat)
    max_lon = max(ovr_d['lon'].max(), target_lon)

    cos_lat = np.cos(np.deg2rad(target_lat))
    # copying method from R. Nelson to compute default LL box.
    box_N = np.round(max_lat + buffer_width, 1)
    box_S = np.round(min_lat - buffer_width, 1)
    box_E = np.round(max_lon + buffer_width/cos_lat, 1)
    box_W = np.round(min_lon - buffer_width/cos_lat, 1)

    box = [box_W, box_S, box_E, box_N]

    return box

def add_coastlines(latlon_extent, ax):
    """
    new method for adding coastlines, to replace the cartopy builtin
    method. This is dependent on the GSHHS data, "h"igh resolution
    shapefiles. We rely on the cartopy shapereader to manage the
    GSHHS datafiles.

    Copies the method created by Rob Nelson for the OCO3 quicklook
    imagery.

    inputs:
    latlon_extent: W,S,E,N of the bounding box in a 4-element arraylike.
    ax: the axis to use for drawing coastlines.
    """

    shp1 = shapereader.Reader(shapereader.gshhs('h', 1))
    shp2 = shapereader.Reader(shapereader.gshhs('h', 2))
    shp3 = shapereader.Reader(shapereader.gshhs('h', 3))

    coastlines_all = [[geoms1 for geoms1 in shp1.geometries()],
                      [geoms2 for geoms2 in shp2.geometries()],
                      [geoms3 for geoms3 in shp3.geometries()]]
    coastlines = [item for sublist in coastlines_all for item in sublist]

    # STRtree is used to efficiently discard most of the coastline data,
    # which is a global dataset. In effect, this loop identifies the
    # small subset of coastline data that resides in the latlon box.
    W, S, E, N = latlon_extent
    plotting_box = sgeom.Polygon([[W,S],[E,S],[E,N],[W,N],[W,S]])
    tree = STRtree([plotting_box])
    coasts_in_box = []
    for i in range(len(coastlines)):
        if len(tree.query(coastlines[i])) > 0: coasts_in_box.append(coastlines[i])
    ax.add_geometries(coasts_in_box, ccrs.PlateCarree(),
                       facecolor='None', edgecolor='k')


def _add_background_map(latlon_extent, ax1):
    """
    for testing only: copy Rob Nelson's method for adding the static
    background map imagery from ESRI ArcGIS server.
    This hides two imports which no other functions need (PIL, urllib).
    """
    from PIL import Image
    import urllib

    W, S, E, N = latlon_extent

    #Background map
    m1 = Proj("epsg:3857", preserve_units=True)
    W_xpos, S_ypos = m1(W,S)
    E_xpos, N_ypos = m1(E,N)
    xpixels = 2000
    ypixels = int((N_ypos - S_ypos) / (E_xpos - W_xpos) * xpixels)
    W_3857, S_3857 = m1(W,S)
    E_3857, N_3857 = m1(E,N)
    url = ('http://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/' +
           f'export?bbox={W_3857},{S_3857},{E_3857},{N_3857}&bboxSR=3857&imageSR=3857' +
           f'&size={xpixels},{ypixels},&dpi=96&format=png32&transparent=true&f=image')
    try: ESRI = np.array(Image.open(urllib.request.urlopen(url)))
    except: #Sometimes this fails randomly, so try again
        time.sleep(20)
        print('Bad response from ESRI, sleeping for 20...')
        ESRI = np.array(Image.open(urllib.request.urlopen(url)))
    im1 = ax1.imshow(ESRI,extent=ax1.get_extent(),origin="upper")


def setup_axes(latlon_extent, crs, fignum=1,
               figsize=(20,20), background_map=False):
    """
    setup axes to match the Quicklook imagery.
    copying code from R. Nelson

    inputs:
    latlon_extent: W,S,E,N of the bounding box in a 4-element arraylike.
    crs: the projection to use with the plotting axes. Normally derived from
        the satpy Scene object. Typically this will be "epsg:3857"
    fignum: matplotlib figure number.
    figsize: matplotlib figure size (width, height) in inches.
    background_map: set to True to use the ESRI supplied background map.
        generally this should be false, since this is designed to use
        Geo imager data. For tests only.

    returns:
    fig : the created matplotlib figure object
    ax1 : the axis for the display imagery and overlay
    ax2 : the inset axis
    cax : colorbar axis
    fig_scalefactor : a float scalefactor to control size of certain
        annotations, relative to the (20,20) "standard" fig size.
    """

    W, S, E, N = latlon_extent

    fig = plt.figure(fignum, figsize=figsize)
    fig.clf()
    ax1 = plt.axes(projection=crs)
    ax1.set_extent([W,E,S,N], ccrs.PlateCarree())

    # an attempt to make some of the annotations scale with figure size.
    # figsize 20,20 is the "default", which is what is used for the actual
    # quicklook imagery on the SAM site.
    fig_scalefactor = 0.5*sum(figsize) / 20.0
    inset_marker_size = int(np.round(100 * fig_scalefactor))
    grid_labelsize = int(np.round(24 * fig_scalefactor))
    inset_ax_size = 2.0 * fig_scalefactor
    cbax_width = 0.4 * fig_scalefactor

    #Load the coastlines
    add_coastlines(latlon_extent, ax1)

    #Grid
    gl = ax1.gridlines(draw_labels=True, color="0.75")
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': grid_labelsize}
    gl.ylabel_style = {'size': grid_labelsize}

    if background_map:
        add_background_map(latlon_extent, ax1)

    #Colorbar
    divider = make_axes_locatable(ax1)
    cbax = divider.append_axes(
        "right", size=cbax_width, pad=0.25*fig_scalefactor, axes_class=plt.Axes)
    layer_cbax = divider.append_axes(
        "bottom", size=cbax_width, pad=1*fig_scalefactor, axes_class=plt.Axes)

    #Globe inset
    ax2 = inset_axes(
        ax1,width=inset_ax_size, height=inset_ax_size,
        loc="upper right", axes_class=GeoAxes,
        axes_kwargs=dict(map_projection=ccrs.Orthographic(((E+W)/2.),((N+S)/2.))))
    ax2.set_global()
    ax2.scatter(((W+E)/2.),((N+S)/2.),c='r',s=inset_marker_size,
                marker='*',transform=ccrs.PlateCarree())
    ax2.stock_img(name='ne_shaded')

    #Coastlines in the globe inset
    ax2.coastlines(color="0.25")

    #Time stamp
    plt.text(0.99,0.01,"Created "+str(datetime.now().day)+' '+calendar.month_abbr[datetime.now().month]+' '+str(datetime.now().year)+"\nSource: NASA/JPL-Caltech",ha='right', va='bottom', transform=ax1.transAxes,color='1.0',size=18)

    return fig, ax1, ax2, cbax, layer_cbax, fig_scalefactor


def create_plot_title_string(cfg_d, ovr_d, odat):
    """
    create a default plot title string, based on the various input data.
    Since we are trying to squeeze in a lot of different pieces of information
    into the title, this requires basically all the input data (the config and
    overlay dictionaries, as well as the actual overlay data.)

    """

    title = ''
    if ovr_d is not None:
        title += ovr_d['sensor'] + ' ' + ovr_d['var_title_string'] + '\n'
    title += 'Background from ' + cfg_d['sensor'] + '\n'

    if odat is not None:
        title += odat['operation_mode'] + ' Mode'
        if odat['target_id'] == 'none' and odat['target_name'] == 'none':
            title += '\n'
        else:
            title += ", " + odat['target_id'] + ', "' + odat['target_name'] + '"\n'

    if odat is not None:
        title += odat['data_version'] + '\n'

    if odat is not None:
        dt = datetime.utcfromtimestamp(odat['time'][0])
        time_string = dt.strftime('%H:%M UTC %-d %b %Y')
    else:
        time_string = cfg_d['datetime'].strftime('%H:%M UTC %-d %b %Y')

    title += time_string

    if ovr_d is not None:
        title += ', Orbit ' + str(ovr_d['orbit'])

    return title


def overlay_data(ax, cb_ax, odata, var_label=None,
                 fig_scalefactor=1.0, **kw):
    """
    overlay data onto background image.
    if vertices are available, shapely polygon obects are used;
    otherwise mpl.scatter is used.
    This function also takes care of the colorbar associated to
    the overlaid data.

    inputs:
    ax: axis object where polygons will be drawn
    cb_ax: the colorbar axis object. Note this is only used if
       the cmap is a colormap.
    odata: oco2 overlay data dictionary

    var_label: string containing a label to add to the colorbar,
         typically the name and units of the variable being overlaid.
    fig_scalefactor: floating point number to scale the annotation
         text sizes in the colorbar.

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
        cbar = plt.colorbar(mpl.cm.ScalarMappable(norm=N_func, cmap=C_func),
                            extend='both',cax=cb_ax)
        cbar.set_label(var_label, size=28*fig_scalefactor,
                       rotation=270, labelpad=35*fig_scalefactor)
        cbar.ax.tick_params(labelsize=22*fig_scalefactor)
        cbar.ax.yaxis.get_offset_text().set_fontsize(22*fig_scalefactor)
