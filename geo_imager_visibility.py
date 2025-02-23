from datetime import datetime

import numpy as np

import pyresample
from satpy.readers._geos_area import get_area_extent, get_area_definition
from pyorbital.orbital import get_observer_look

def _GOES_fixed_data():
    # helper to return fixed data used for all GOES AreaDef creation
    # copied attributes directly from "goes_imager_projection" variable
    d = dict(
        perspective_point_height = 35786023.0,
        semi_major_axis = 6378137.0,
        semi_minor_axis = 6356752.31414,
        sweep_angle_axis = 'x',
        ncols  = 10848,
        nlines = 10848,
        x = [-0.151858,  0.151858],
        y = [ 0.151858, -0.151858],
        spatial_resolution = '1km at nadir',
    )
    
    return d

def _GOES_East_data(mode='F'):
    # helper to return mode-specific values that are
    # specific to GOES-East
    d = dict(
        longitude_of_projection_origin = -75.0,
        orbital_slot = 'GOES-East',
    )
    # other values, now specific to the image (CONUS or FD)
    if mode == 'F':
        # the fixed numbers will be OK
        pass
    elif mode == 'C':
        d['ncols']  = 5000
        d['nlines'] = 3000
        d['x'] = [-0.101346, 0.038626]
        d['y'] = [ 0.128226, 0.044254]
    else:
        raise ValueError('Bad mode, must be "C" or "F"')
    return d

def _GOES_West_data(mode='F'):
    # helper to return mode-specific values that are
    # specific to GOES-West
    d = dict(
        longitude_of_projection_origin = -137.0,
        orbital_slot = 'GOES-West',
    )
    # other values, now specific to the image (CONUS or FD)
    if mode == 'F':
        # the fixed numbers will be OK
        pass
    elif mode == 'C':
        d['ncols']  = 5000
        d['nlines'] = 3000
        d['x'] = [-0.069986, 0.069986]
        d['y'] = [ 0.128226, 0.044254]
    else:
        raise ValueError('Bad mode, must be "C" or "F"')
    return d


def GOES_AreaDef(sat, mode='F'):
    """
    obtain the pyresample AreaDefinition for GOES.

    Parameters
    ----------
    sat : str
        either "W" or "E", specifying GOES-West or -East.
    mode : str
        either "F" or "C", specifying Full-Disk or CONUS

    Returns
    -------
    d : pyresample.geometry.AreaDefinition
        the Area Definition object needed by satpy/pyresample/etc.
    
    """
    d = _GOES_fixed_data()
    if sat == 'W':
        _d = _GOES_West_data(mode=mode)
    elif sat == 'E':
        _d = _GOES_East_data(mode=mode)
    else:
        raise ValueError('sat must be "E" or "W"')
    d.update(_d)

    # rest is copied from satpy.readers.abi_base.NC_ABI_BASE._get_areadef_fixedgrid()
    # with almost no modification
    a = d['semi_major_axis']
    b = d['semi_minor_axis']
    h = d['perspective_point_height']
    lon_0 = d['longitude_of_projection_origin']
    sweep_axis = d['sweep_angle_axis']
    
    x_l = d['x'][0]
    x_r = d['x'][-1]
    y_l = d['y'][-1]
    y_u = d['y'][0]
    x_half = (x_r - x_l) / (d['ncols'] - 1) / 2.
    y_half = (y_u - y_l) / (d['nlines'] - 1) / 2.
    area_extent = (x_l - x_half, y_l - y_half, x_r + x_half, y_u + y_half)
    area_extent = tuple(np.round(h * val, 6) for val in area_extent)

    proj_dict = {'proj': 'geos',
                 'lon_0': float(lon_0),
                 'a': float(a),
                 'b': float(b),
                 'h': h,
                 'units': 'm',
                 'sweep': sweep_axis}

    area_def = pyresample.geometry.AreaDefinition(
        d['orbital_slot'],
        d['spatial_resolution'],
        'abi_fixed_grid',
        proj_dict,
        d['ncols'],
        d['nlines'],
        np.asarray(area_extent))

    return area_def


def Himawari_AreaDef():
    """
    obtain the pyresample AreaDefinition for Himawari.

    Parameters
    ----------
    None

    Returns
    -------
    d : pyresample.geometry.AreaDefinition
        the Area Definition object needed by satpy/pyresample/etc.
    
    """
    
    pd = {'cfac': 40932549,
          'lfac': 40932549,
          'coff': 5500.5,
          'loff': 5500.5,
          'a': 6378137.0,
          'h': 35785863.0,
          'b': 6356752.3,
          'ssp_lon': 140.7,
          'nlines': 11000,
          'ncols': 11000,
          'scandir': 'N2S',
          'a_name': 'FLDK',
          'a_desc': 'AHI FLDK area',
          'p_id': 'geosh8'}

    area_ext = get_area_extent(pd)
    
    area_def = get_area_definition(pd, area_ext)

    return area_def


def get_sensor_zenith(sat_lon, sat_lat, sat_alt, obs_lon, obs_lat, obs_alt):
    """
    get_sensor_zenith: return the zenith angle from the observer point
    (a lon/lat/alt on Earth) to a satellite observer (also specified
    by lon/lat/alt)
    This is primarily a helper function for get_sensor_zenith_bysat().

    All inputs are assumed to be scalars.

    Parameters
    ----------
    sat_lon : float
        longitude of subsatellite point in degrees
    sat_lat : float
        latitude of subsatellite point in degrees
    sat_alt : float
        altitude of satelliate in km
    obs_lon : float
        longitude of observed point in degrees
    obs_lat : float
        latitude of observed point in degrees
    obs_alt : float
        altitude of observed point in km

    Returns
    -------
    sat_zen : float
        the zenith angle of the satellite from the observed point.
        
    """
    
    #
    # implementation note:
    # pyorbital.get_observer_look needs the time as an input, even though
    # this calculation does not actually require it (see long answer below).
    # short answer: we can just pick any time we want as a dummy value,
    # and proceed.
    #
    # If you check the calculation with different times, you do get the
    # same result regardless of the time (to within fp error).
    # After looking through pyorbital code, I think the issue is that they
    # refer to a T.S. Kelso algorithm that describes a calculation
    # with the satellite position in ECI, and the observer position in
    # ECEF; in that algorithm, the observer position is converted to
    # ECI. The pyorbital calculation follows this algorithm so it converts
    # the satellite position (which we actually already have in ECEF)
    # and converts it to ECI.
    #

    dummy_utc_time = datetime(2000,1,1)
    # get_obsever_look also assumes array-like inputs (certain calculations
    # have np.where applied). If we put the first value in a 1-element
    # array, the rest of the calc gets internally broadcasted.
    sat_azi, sat_el = get_observer_look(
        [sat_lon], sat_lat, sat_alt, dummy_utc_time,
        obs_lon, obs_lat, obs_alt)
    # output is 1-element ndarray (see above), so convert back to scalar.
    sat_zen = 90 - sat_el[0]

    return sat_zen


def get_sensor_zenith_bysat(obs_lon, obs_lat, sat, mode=None):
    """
    get_sensor_zenith_bysat:
    return the zenith angle from the observer point (a lon/lat on Earth,
    assumed to be at sea level) to the given geostationary satellite.
    Input lat/lon is assumed to be scalar (i.e. we compute this for only
    a single point)

    Parameters
    ----------
    obs_lon : float
        longitude of observed point in degrees
    obs_lat : float
        latitude of observed point in degrees
    sat : str
        string specifying the satellite. Options: 'himawari',
        'GOESW', or 'GOESE'
    mode : str
        'F' or 'C' for (F)ull Disk or (C)ONUS. Only valid
        for GOES, as we do not attempt to use more than 1 mode
        from Himawari.

    Returns
    -------
    sat_zen : float
        the zenith angle of the satellite from the observed point.
        If the requested lat/lon is not visible by the requested satellite,
        np.nan is returned.

    """

    if sat not in ('GOESE', 'GOESW', 'himawari'):
        raise ValueError('sat must be one of "GOESE", "GOESW", or "himawari"')

    if sat == 'himawari':
        area_def = Himawari_AreaDef()
    else:
        area_def = GOES_AreaDef(sat[-1], mode)
    
    try:
        col, row = area_def.get_array_indices_from_lonlat(obs_lon, obs_lat)
        visible = True
    except ValueError:
        visible = False

    if visible:
        sat_alt = area_def.proj_dict['h']
        sat_lat = 0.0
        sat_lon = area_def.proj_dict['lon_0']
        sensor_zenith = get_sensor_zenith(
            sat_lon, sat_lat, sat_alt/1000.0, obs_lon, obs_lat, 0.0)
    else:
        sensor_zenith = np.nan

    return sensor_zenith


def is_visible(obs_lon, obs_lat, sat, mode=None, zenith_limit=None):
    """
    is_visible:
    determine if a certain observed lat/lon is visible by the geostationary
    satellite. Input lat/lon is assumed to be scalar (i.e. a single location)

    Parameters
    ----------
    obs_lon : float
        longitude of observed point in degrees
    obs_lat : float
        latitude of observed point in degrees
    sat : str
        string specifying the satellite. Options: 'himawari',
        'GOESW', or 'GOESE'
    mode : str
        'F' or 'C' for (F)ull Disk or (C)ONUS. Only valid
        for GOES, as we do not attempt to use more than 1 mode
        from Himawari.
    zenith_limit : None or float
        If a floating point number is input, then this is applied as a
        threshold value for determining visibility; e.g. if 80 is input,
        then the point is visible only if the zenith angle is less than
        80 degrees.
        If None, no such check is performed

    Returns
    -------
    visible : bool
        the zenith angle of the satellite from the observed point.
        If the requested lat/lon is not visible by the requested satellite,
        np.nan is returned.
    """

    if sat not in ('GOESE', 'GOESW', 'himawari'):
        raise ValueError('sat must be one of "GOESE", "GOESW", or "himawari"')

    if sat == 'himawari':
        area_def = Himawari_AreaDef()
    else:
        area_def = GOES_AreaDef(sat[-1], mode)

    try:
        col, row = area_def.get_array_indices_from_lonlat(obs_lon, obs_lat)
        visible = True
    except ValueError:
        visible = False

    if visible and zenith_limit:
        sat_alt = area_def.proj_dict['h']
        sat_lat = 0.0
        sat_lon = area_def.proj_dict['lon_0']
        sensor_zenith = get_sensor_zenith(
            sat_lon, sat_lat, sat_alt/1000.0, obs_lon, obs_lat, 0.0)
        if sensor_zenith >= zenith_limit:
            visible = False

    return visible


def determine_optimal_geo_satellite(lons, lats, obs_datetime, zenith_limit=None):
    """
    determine_optimal_geo_satellite:
    Determine which geostationary satellite imager is the most optimal for
    the given observed region. The region is specified by a set of lat/lon
    describing the region's corner points or outline polygon.

    The basic method is check the visibility of the surface location for
    each geo satellite. In the case the location can be viewed by more than
    one satellite (for example, the overlap region between GOES-East and -West),
    the satellite with the smaller viewing zenith angle is selected.

    Since a bounding polygon is input, the geo satellite is selected only
    if all points are visible. In other words, a region that is only partly
    visible will be treated as not visible.

    For GOES, preference is always given to the CONUS view over the Full Disk.

    Parameters
    ----------
    lons : array-like of float
        longitudes, in degrees, of the region - this should contain
        a set of corner points, or other bounding polygon.
    lats : array-like of float
        longitudes, in degrees, of the region - this should contain
        a set of corner points, or other bounding polygon.
    obs_datetime : datetime
        time of observation as python datetime. This is used to select between
        satellites when there are operational swaps, e.g.,
        Himawari-08 to 09. The datetime needs only to specifiy the date, the
        time of day is not required.
    zenith_limit : None or float
        set to a floating point value to enforce a viewing zenith limit.
        For example, setting this to 80.0 will mean the region needs to
        be entirely visible by the geo satellite with view zenith angles
        less than 80.0 degrees.

    Returns
    -------
    selected_sensor : str or None
        a string sensor name, among the following choices (these are
        the string ID names needed for the vistool functions)
        Himawari-08
        Himawari-09
        GOES16_ABI_C, GOES16_ABI_F
        GOES17_ABI_C, GOES17_ABI_F
        GOES18_ABI_C, GOES18_ABI_F

       If none of the geo satellites can view the requested location,
       None is returned.
    """

    llpairs = list(zip(lons, lats))

    visible_sensors = []
    sensor_zeniths = []

    center_lat = np.mean(lats)
    center_lon = np.mean(lons)

    checks = [is_visible(*ll, 'himawari', zenith_limit=zenith_limit) for ll in llpairs]

    if np.all(checks):
        sensor_zeniths.append(
            get_sensor_zenith_bysat(center_lon, center_lat, 'himawari', mode='F'))
        # last complete day in Himawari-08 archive is 2022 12 12.
        # newer data from Himawari-09.
        if obs_datetime < datetime(2022, 12, 13):
            visible_sensors.append('Himawari-08')
        else:
            visible_sensors.append('Himawari-09')

    checksF = [is_visible(*ll, 'GOESW', mode='F', zenith_limit=zenith_limit) for ll in llpairs]

    if np.all(checksF):
        
        checksC = [is_visible(*ll, 'GOESW', mode='C', zenith_limit=zenith_limit) for ll in llpairs]
        if np.all(checksC):
            sensor_zeniths.append(
                get_sensor_zenith_bysat(center_lon, center_lat, 'GOESW', mode='C'))
            # GOES-18 became GOES-West on January 10, 2023.
            # see https://registry.opendata.aws/noaa-goes/
            if obs_datetime < datetime(2023, 1, 10):
                visible_sensors.append('GOES17_ABI_C')
            else:
                visible_sensors.append('GOES18_ABI_C')
        else:
            sensor_zeniths.append(
                get_sensor_zenith_bysat(center_lon, center_lat, 'GOESW', mode='F'))
            # GOES-18 became GOES-West on January 10, 2023.
            # see https://registry.opendata.aws/noaa-goes/
            if obs_datetime < datetime(2023, 1, 10):
                visible_sensors.append('GOES17_ABI_F')
            else:
                visible_sensors.append('GOES18_ABI_F')

    checksF = [is_visible(*ll, 'GOESE', mode='F', zenith_limit=zenith_limit) for ll in llpairs]

    if np.all(checksF):

        checksC = [is_visible(*ll, 'GOESE', mode='C', zenith_limit=zenith_limit) for ll in llpairs]
        if np.all(checksC):
            sensor_zeniths.append(
                get_sensor_zenith_bysat(center_lon, center_lat, 'GOESE', mode='C'))
            visible_sensors.append('GOES16_ABI_C')
        else:
            sensor_zeniths.append(
                get_sensor_zenith_bysat(center_lon, center_lat, 'GOESE', mode='F'))
            visible_sensors.append('GOES16_ABI_F')

    if len(visible_sensors) > 0:
        min_zenith_idx = np.argmin(sensor_zeniths)
        selected_sensor = visible_sensors[min_zenith_idx]
    else:
        selected_sensor = None

    return selected_sensor
