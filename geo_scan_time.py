import datetime

def approx_scanline_time(stime, etime, lat, domain):
    """
    Returns the scanline time, for a particular lat, and domain (C or F,
    for (C)ONUS or (F)ull disk.)

    Note this is a rough approximation, assuming the actual collection
    time for a given scan line is a linear interpolation between the
    start and end time in the file. It also assumes the north-most
    scan line is collected first (Should be true for ABI, AHI).
    This linear interpolation is not exactly true, because the actual
    scan pattern is much more complicated. Doing the correct
    calculation would be very difficult, and this approximation should
    be accurate to a minute or two.

    returns a single datetime object with the approx scanline time.

    if the lat is out of range, it returns the closest time; in other
    words, out of range will be limited to pick the stime or etime.

    inputs:
    stime: datetime object with start time
    etime: datetime object with end time
    lat: latitude (scalar)
    domain: (C)ONUS or (F)ull disk
    """

    if domain == 'C':
        slat = 50.0
        elat = 15.0
    elif domain == 'F':
        slat = 81.0
        elat = -81.0
    else:
        raise ValueError('invalid domain: '+domain)

    tfrac = (slat - lat) / (slat - elat)

    if tfrac < 0:
        scan_time = stime
    elif tfrac > 1:
        scan_time = etime
    else:
        total_time = (etime-stime).total_seconds()
        scan_time = stime + datetime.timedelta(seconds=total_time*tfrac)

    return scan_time


def compute_time_offset(scn, lat, domain, datetime_utc):
    """
    compute the time offset from a given datetime in UTC, relative
    to the observation time of a geostationary image.

    uses _get_scanline_time() to estimate the actual observation time
    of the time in the image.

    Inputs:
    scn : a satpy Scene object containing the geostationary image data
    lat : the center lat of the requested region
    domain : (C)ONUS or (F)ull disk
    datetime_utc : a python datetime object, assumed in UTC, for the
        requested time

    returns:
    time_offset : the offset between the input time and the observation time.
        In other words, datetime_utc - obs_time; reported in seconds.
    """

    stime = scn.start_time
    etime = scn.end_time
    
    scan_time = approx_scanline_time(stime, etime, lat, domain)

    time_offset = (datetime_utc - scan_time).total_seconds()

    return time_offset


def get_AHI_timerange_from_filename(filename):
    """
    get the approximate time range (start and end time of the image
    collection period) from the AHI filename.

    This is roughly speaking the 'planned' time, and assumes the
    start time is the 10-minute quantized time and the end time
    is +10 minutes later.

    Returns:
    stime: the start time as a python datetime object
    etime: the end time as a python datetime object
    """

    filename_only = filename.split('/')[-1]
    stime = datetime.datetime.strptime(filename_only[7:20], '%Y%m%d_%H%M')
    etime = stime + datetime.timedelta(minutes=10)

    return stime, etime


def get_ABI_timerange_from_filename(filename):
    """
    get the start and end times from the ABI image from the filename.

    Returns:
    stime: the start time as a python datetime object
    etime: the end time as a python datetime object
    """

    stimestamp, etimestamp = filename.split('_')[-3:-1]
    # [1:-1] slice removes 's' or 'e' and the extra number (frac seconds?)
    stime = datetime.datetime.strptime(stimestamp[1:-1], '%Y%j%H%M%S')
    etime = datetime.datetime.strptime(etimestamp[1:-1], '%Y%j%H%M%S')

    return stime, etime
