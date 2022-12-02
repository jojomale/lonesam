# Johanna Lehr, 17-Dec-2021

"""
Helper routines for base, sds-db modules.
"""


import os
import numpy as np
from scipy.signal import get_window

from obspy.core import UTCDateTime as UTC
from obspy.signal.filter import bandpass

import logging

from . import dqclogging
# Create the global logger
logger = dqclogging.create_logger()
module_logger = logging.getLogger(logger.name+'.util')


def _create_hdf5_attribs(f, stationcode, starttime, endtime, 
                        freqs, winlen_seconds, proclen_seconds):
    """
    Set meta data in hdf5-file as attributes.

    Parameters
    ----------------
    f : h5py.File
        open, editable h5py.File-object.
    stationcode : str
        string of format 
        "{network}.{station}.{location}.{channel}"
    starttime : UTCDateTime
        starttime of data in file
    endtime : UTCDateTime
        endtime of data in file. to be precise, it is the
        starttime of the last processing window, i.e.
        the actual endtime-proclen_seconds 
    freqs : 2-tuple
        min, max frequency between which amplitude was filtered
    winlen_seconds : int
        length of window over which amplitudes were extracted
        and psds computed, in seconds
    proclen_seconds : int
        length of time series that was processed at once, 
        in seconds. E.g. 1 day of data
    """
    f.attrs.create('stationcode', np.string_(stationcode))
    f.attrs.create('starttime',  np.array(starttime.timetuple()[:6]))
    f.attrs.create('endtime', np.array(endtime.timetuple()[:6]))
    f.attrs.create('amplitude_frequencies', np.array(freqs))
    f.attrs.create('seconds_per_window', np.array(winlen_seconds))
    f.attrs.create('seconds_per_proclen', np.array(proclen_seconds))


def process_stream(st, inv, starttime, endtime, 
                    sampling_rate=None):
    """
    Standard processing steps to perform on seismic data
    before further analysis. Returns single trace.

    Includes:
    - remove sensitivity
    - merge data and fill gaps with ``np.nan``
    - trim, pad with ``np.nan``
    - demean by subtracting ``np.nanmean``

    
    Parameters
    ----------------
    st : obspy.core.Stream
        raw seismic data. Should be only one station!
    inv : obspy.core.Inventory
        meta data of station, should include response
        information 
    starttime : obspy.core.UTCDateTime
        starttime to which trace is trimmed
    endtime : obspy.core.UTCDateTime
        endtime to which trace is trimmed

    Returns
    ------------
    tr : obspy.core.Trace
        Single trace.

    Warning
    ---------
    If more than one station-channel combination is present,
    we will return only the first one. However, there can be 
    multiple traces of the same station-channel.

    Notes
    ------
    If data is obtained from a datacenter, it may contain 
    gaps. If so, the stream contains multiple traces, which
    we ``st.merge()`` into a single trace, filling the gaps 
    with Nan.

    """
    st.remove_sensitivity(inv)
    #st.merge(fill_value=np.nan)
    resample(st, sampling_rate)
    st.merge(fill_value=np.nan)
    st.trim(starttime, endtime, pad=True, fill_value=np.nan, 
            nearest_sample=False)
    if len(st) > 1:
        module_logger.warning("More than 1 trace in stream %s!"+ 
            " Using first trace only" % st)
    tr = st[0]
    
    # Demean ignoring gaps
    tr.data = tr.data - np.nanmean(tr.data)
    return tr


def resample(st, sampling_rate):
    """
    Resample traces in stream to common sampling rate 
    if necessary.

    We use `tr.interpolate()` to increase and
    `tr.resample()` to decrease sampling rate.
    """
    for tr in st:
        if tr.stats.sampling_rate < sampling_rate:
            module_logger.info("Up-sampling of {} from {:g} Hz to {} Hz".format(
                tr.id, tr.stats.sampling_rate, sampling_rate
            ))
            tr.interpolate(sampling_rate)
            
        elif tr.stats.sampling_rate > sampling_rate:
            module_logger.info("Down-sampling {} from {:g} Hz to {} Hz".format(
                tr.id, tr.stats.sampling_rate, sampling_rate
            ))
            tr.resample(sampling_rate, no_filter=False)
           


def merge_different_samplingrates(st):
    """
    Merge stream, fill gaps with nans.
    If traces have different sampling rates, resample to highest.
    """
    try:
        st.merge(fill_value=np.nan)
    except Exception:
        sr = max([tr.stats.sampling_rate for tr in st])
        module_logger.warning("Found different sampling rates. " + 
                    "Resampling to highest ({:g} Hz).".format(sr))
        st.resample(sr)
        st.merge(fill_value=np.nan)


def get_adjacent_frames(tr, starttime, nf, winlen_samples):
    """
    Reshape vector into (``nf``, ``winlen_samples``). These
    frames do not overlap!.

    We slice the trace at ``starttime`` and add as many samples
    as necessary to fill the targeted shape.

    One frame contains the data over which amplitude and psd
    are computed.

    Parameters
    ---------------
    tr : obspy.core.Trace
        seismic data trace
    starttime : obspy.core.UTCDateTime
        where to start in the trace
    nf : int
        number of frames
    winlen_samples : int
        number of samples in one frame
    """
    #print(starttime)
    #print(tr)
    ntot = int(nf*winlen_samples)
    data = tr.slice(starttime, endtime=None).data[:ntot]
    return data.reshape((nf, winlen_samples))



def get_overlapping_tapered_frames(tr, starttime, nf, winlen_samples,
                           taper_samples):
    """
    Splits the vector up into (overlapping) Tukey windows.
    
    A Tukey window combines a boxcar-window with cosine-tapered
    flanks at the egdes.    

    Frames containing any Nans are set entirely to Nan.
    Loosely based on `obspy.signal.util.enframe` but faster.

    We slice the trace at ``starttime`` and add as many samples
    as necessary to fill the targeted shape.

    One frame contains the data over which amplitude and psd
    are computed.

    Parameters
    ---------------
    tr : obspy.core.Trace
        seismic data trace
    starttime : obspy.core.UTCDateTime
        where to start in the trace
    nf : int
        number of frames
    winlen_samples : int
        number of samples in one frame
    taper_samples : int
        number of samples for the tapered region. Should be
        long enough to accommodate potential filter effects.
    """
    sr = tr.stats.sampling_rate
    
    # Samples in window including tapers
    nwin = int(winlen_samples + 2*taper_samples)
    
    # Total number of samples of trace to process
    proclen_samples = int(nf * winlen_samples + 2*taper_samples)
    
    # Cut out the needed data
    x = tr.slice(starttime-taper_samples/sr).data[:proclen_samples]
    
    # Ratio of tapers to total window size
    a =  2*taper_samples / nwin
    win = get_window(('tukey', a), nwin, fftbins=False)
    
    # From obspy.signal.enframe()
    #nx = len(x)
    #nwin = len(win)
    if (len(win) == 1):
        length = win
    else:
        length = nwin
    #nf = int(np.fix((nx - length + winlen_samples) // winlen_samples))
    # f = np.zeros((nf, length))
    indf = winlen_samples * np.arange(nf)
    f = x[np.expand_dims(indf, 1) + 
          np.expand_dims(np.arange(length), 0)]
    f = f * win
    f[np.any(np.isnan(f), axis=1),:] = np.nan
    #no_win, _ = f.shape
    return f, taper_samples



def get_amplitude(tr, starttime, fmin, fmax, overlap, 
        winlen_samples, nf, percentile=75):
    """
    Extract percentile of amplitude in consecutive
    time windows.

    The seismic data in ``tr`` is bandpass-filtered before
    extracting the amplitudes.

    The percentile is computed over ``winlen_samples``.
    Windows are created using either 
    ``get_overlapping_tapered_frames`` or
    ``get_adjacent_frames``


    Parameters
    --------------
    tr : Trace
        time series of seismic data to process
    starttime : UTCDateTime
        starttime in trace
    fmin, fmax : floats
        min and max frequency of bandpass
    overlap : int
        length of overlap in seconds. Only used if nans are
        present in data. See Note.
    winlen_samples : int
        Size of windows (frames) over which percentile is 
        computed, in samples.
    nf : int
        number of frames into which data is divided
    percentile : int, float [75]
        percentile of amplitude, which we get for each frame.


    Note
    ----------
    If trace is free of nans, we can simply filter the whole trace at 
    once and use reshape to get the windows.
    If there are nans, the obspy filter function only filters
    the data up to the first nan. Thus you can loose almost the entire
    trace because of a single nan in the early part. To avoid this, we
    check for nans and if there are some, we filter the data per
    window. This means however, we have to ensure again some overlap
    and the windowing will be slower.
    """
    sr = tr.stats.sampling_rate
    
    if np.any(np.isnan(tr.data)):
        module_logger.info('Found nans in %s' % tr)
        taper_samples = int(overlap*sr)
        
        data, taper_samples = get_overlapping_tapered_frames(tr, 
                            starttime, nf, winlen_samples,
                           taper_samples)
        data = bandpass(data, fmin, fmax, sr)
        data = data[:,taper_samples:-taper_samples]
                                    
    else:
        tr.filter('bandpass', freqmin=fmin, freqmax=fmax)
        data = get_adjacent_frames(tr, starttime+overlap, nf, winlen_samples)
        
    prctl = np.nanpercentile(data, percentile, axis=1)
    return prctl




def iter_years(startdate, enddate):
    """
    Iterator that returns date of first and last day
    of each year between ``startdate`` and ``enddate``.
    For years of ``startdate`` and ``enddate``, respective
    dates are returns.

    Intended to use for processing data in batches of 1 year, 
    i.e. output files will cover 1y of data.


    Parameters
    ------------
    startdate : UTCDateTime
    enddate : UTCDateTime



    Example
    ---------
    .. code-block:: python
        
        from data_quality_control import util
        from obspy.core import UTCDateTime as UTC

        sdate = UTC("2018-05-20")
        edate = UTC("2021-10-30")
        for s, e in util.iter_years(sdate, edate):
            print(s, e)

        >>> 2018-05-20T00:00:00.000000Z 2018-12-31T00:00:00.000000Z
        >>> 2019-01-01T00:00:00.000000Z 2019-12-31T00:00:00.000000Z
        >>> 2020-01-01T00:00:00.000000Z 2020-12-31T00:00:00.000000Z
        >>> 2021-01-01T00:00:00.000000Z 2021-10-30T00:00:00.000000Z

    """
    # Init loop counters
    _startdate = startdate
    _enddate = UTC("{:d}-12-31".format(_startdate.year))
    while _startdate < enddate:
        # Check if end-of-year date is larger than acual enddate
        if _enddate > enddate: 
            _enddate = enddate
            module_logger.debug("Reset enddate to %s" % _enddate)
        # Process 1 year (or less)
        module_logger.info("\nProcessing %s - %s" % (_startdate, _enddate))
        yield _startdate, _enddate

        _startdate = UTC("{:d}-01-01".format(_startdate.year+1))
        _enddate = UTC("{:d}-12-31".format(_startdate.year))
        # module_logger.debug("New _startdate = %s" % _startdate)
        # module_logger.debug("New _enddate = %s" % _enddate)


def get_end_of_month(stime):
    """
    Return last day of month in ``stime``
    as ``UTCDateTime.date``
    """
    t = next_month(stime)
    return UTC((t - 1).date)


def next_month(stime):
    """
    Return first day of month after ``stime`` as
    ``UTCDateTime``. Respects turn of year.
    """
    quot, rem = np.divmod(stime.month + 1, 12)
    if rem==0:
        rem = 1
    stime = UTC("{:d}-{:02d}-01".format(stime.year+quot, rem)) 
    return stime
    

def iter_month(startdate, enddate):
    """
    Iterator that returns date of first and last day
    of each month between ``startdate`` and ``enddate``.
    For months of ``startdate`` and ``enddate``, respective
    dates are returns.

    Intended to use for processing data in batches of 1 month, 
    i.e. output files will cover 1 month of data.


    Parameters
    ------------
    startdate : UTCDateTime
    enddate : UTCDateTime


    Example
    ---------
    .. code-block:: python
        
        from data_quality_control import util
        from obspy.core import UTCDateTime as UTC

        sdate = UTC("2020-05-20")
        edate = UTC("2021-10-30")
        for s, e in util.iter_month(sdate, edate):
            print(s.date, e.date)

        >>> 2020-05-20 2020-05-31
        >>> 2020-06-01 2020-06-30
        >>> 2020-07-01 2020-07-31
        >>> 2020-08-01 2020-08-31
        >>> 2020-09-01 2020-09-30
        >>> 2020-10-01 2020-10-31
        >>> 2020-11-01 2020-12-31
        >>> 2021-01-01 2021-01-31
        >>> 2021-02-01 2021-02-28
        >>> 2021-03-01 2021-03-31
        >>> 2021-04-01 2021-04-30
        >>> 2021-05-01 2021-05-31
        >>> 2021-06-01 2021-06-30
        >>> 2021-07-01 2021-07-31
        >>> 2021-08-01 2021-08-31
        >>> 2021-09-01 2021-09-30
        >>> 2021-10-01 2021-10-25

    """
    # Init loop counters
    _startdate = startdate
    _enddate = get_end_of_month(startdate)
    while _startdate < enddate:
        # Check if end-of-year date is larger than acual enddate
        if _enddate > enddate: 
            _enddate = enddate
            module_logger.debug("Reset enddate to %s" % _enddate)
        # Process 1 year (or less)
        module_logger.info("\nProcessing %s - %s" % (_startdate, _enddate))
        yield _startdate, _enddate

        _startdate = next_month(_startdate) 
        _enddate = get_end_of_month(_startdate)
        



def iter_timeinc(startdate, enddate, inc, timelevel):
    """
    General iterator that yields starttime, endtime of 
    every ``inc`` time interval.

    Use for ``inc`` < 1 month.

    Parameters
    -----------
    startdate : UTCDateTime
    enddate : UTCDateTime
    inc : int
        duration of interval in seconds
    timelevel : int (3, 4, 5, 6)
        Up to which index of the timetuple should be used.
        E.g. for ``inc=3600``, ``timelevel=4``
        (== iter_hours()).
        For ``inc=12*3600`` also use ``timelevel=4``.
        But for ``inc=24*3600``, use ``timelevel=3`` 
        (== iter_days()).

    """
    _startdate = UTC(*np.array(startdate.timetuple())[:timelevel])
    _enddate = _startdate + inc
    while _startdate < enddate:
        # Check if end-of-year date is larger than acual enddate
        if _enddate > enddate: 
            _enddate = enddate
            module_logger.debug("Reset enddate to %s" % _enddate)
        # Process 1 year (or less)
        module_logger.info("\nProcessing %s - %s" % (_startdate, _enddate))
        yield _startdate, _enddate

        _startdate = _startdate + inc
        _enddate = _startdate + inc


def iter_days(starttime, endtime):
    """
    Returns starttime and endtime of every day (24h) 
    between ``starttime`` and ``endtime``.

    Wrapper for 
    ``iter_timeinc(starttime, endtime, inc=24*3600, timelevel=3)``
    """
    return iter_timeinc(starttime, endtime, 24*3600, 3)
        

def iter_hours(starttime, endtime):
    """
    Returns starttime and endtime of every hour (3600s) 
    between ``starttime`` and ``endtime``.

    Wrapper for 
    ``iter_timeinc(starttime, endtime, inc=3600, timelevel=4)``
    """
    return iter_timeinc(starttime, endtime, 3600, 4)


# Format strings for output-filenames, 
# depending on `fileunit` of Processor
FNAME_BASE = "{outdir}"+os.path.sep+"{network}.{station}.{location}.{channel}_"
FNAME_FMTS = {
    "year" :    FNAME_BASE + "{year:04d}.hdf5",
    "month" :   FNAME_BASE + "{year:04d}-{month:02d}.hdf5",
    "day" :     FNAME_BASE + "{year:04d}-{month:02d}-{day:02d}.hdf5",
    "hour" :    FNAME_BASE + "{year:04d}-{month:02d}-{day:02d}-{hour:02d}.hdf5",
    }

TIME_ITERATORS = {
    "year": iter_years,
    "month": iter_month,
    "day" : iter_days,
    "hour": iter_hours
}


# Mapper of characters to define numpy.datetime data types,
# depending on time increment
datetime_flags = {"Y" : 24*3600*365,
       "M" : 24*3600*28,
       "D" : 24*3600,
       "h" : 3600,
       "m" : 60,
       "s" : 1}


def choose_datetime_inc(dt):
    """
    Determine a human-friendly time increment
    for dt (in seconds) and numpy.datetime 
    character.

    Used to create ticklabels for time axis.

    Example
    ---------------
    - ``dt = 3600`` --> ``k = "H", inc = 1``
    - ``dt = 7200`` --> ``k = "H", inc = 2``
    - ``dt = 900`` --> ``k = "m", inc = 15``

    Parameters
    ------------
    dt : int
        time increment in seconds

    Returns
    ------------
    k : str ["M", "D", "h", "m", "s"]
        character specifying numpy.datetime data type
        ``dtype='datetime64[{}]'.format(k)``
    inc : int
        time increment in unit of ``k``

    """
    for k, v in datetime_flags.items():
        if v <= dt :
            break
    inc = int(dt / v)
    return k, inc
