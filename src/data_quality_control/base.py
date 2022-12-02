"""
This module provides low-level functionalities to extract
amplitudes and power-spectral densities (PSDs) from seismic
data.

It is intended to provide the most flexible interface in
terms of data source and selecting time intervals for window 
size, processing length and file units.

The processing workflow is managed by the 
``GenericProcessor``-class, which can used as base for
more customized classes. E.g. the ``sds_db.SDSProcessor``
is taylored to using the sds-client of obspy.
"""

from datetime import timedelta
import time

import os.path

import numpy as np
from scipy.signal import welch

import matplotlib.pyplot as plt
# Use try just in case older versions don't know the style
try:
    plt.style.use('tableau-colorblind10')
except FileNotFoundError:
    pass


from obspy.core import UTCDateTime as UTC

import plotly.graph_objects as go

import h5py

from . import util, dqclogging

import logging
# Create the global logger
logger = dqclogging.create_logger()
module_logger = logging.getLogger(logger.name+'.base')


STATIONCODE = "{network}.{station}.{location}.{channel}"


default_processing_params = dict(
    overlap = 60,
    amplitude_frequencies = (4,14),
    nperseg = 2048,
    winlen_seconds = 3600,
    proclen_seconds = 24*3600,
    sampling_rate = 100
)



class GenericProcessor():
    """
    Organize processing of data sets and save results
    to disk.

    Parameters
    -------------------
    network : str
        network code
    station : str
        station code
    location : str
        location code
    channel : str
        channel code
    dataclient : obspy.clients
        client to obtain data from 
        (``dataclient.get_waveforms()``)
    invclient : obspy.clients
        client to get meta data from
        (``invclient.get_stations()``)
    outdir : str, path-like
        directory to which results are written
    preprocessing : func
        Function, that performs basic preprocessing
        of seismic data. Must take obspy.Stream as argument
        and return obspy.Trace.
    fileunit : {"hour", "day", "month", "year"}, ["year"]
        Interval after which a new output file is started.
        Determines also format string for file-name
        E.g.: ``fileunit="year"`` means all results for data
        from 1 year (01-January - 31-December) will be stored
        in one file. If start/end time for processing are different
        from 01-January/31-December files will be truncated
        accordingly.
    **procparams : 
        processing parameters as keyword arguments

    
    Note
    ----------
    Format of filenames will be determined by ``fileunit`` as
    
    .. code-block:: python

        FNAME_BASE = "{outdir}/{network}.{station}.{location}.{channel}_"  
        FNAME_FMTS = {
            "year" :    FNAME_BASE + "{year:04d}.hdf5",
            "month" :   FNAME_BASE + "{year:04d}-{month:02d}.hdf5",
            "day" :     FNAME_BASE + "{year:04d}-{month:02d}-{day:02d}.hdf5",
            "hour" :    FNAME_BASE + "{year:04d}-{month:02d}-{day:02d}-{hour:02d}.hdf5",
            }

    ``year, month, day, hour`` come from datetime object.
    
    
    Note
    ------
    ``expand_nslc()`` can be customized for a specific 
    data base or data client. If implemented, this may
    allow for wildcards in 
    ``station, network, location, channel``
    
    """

    def __init__(self, network, station, location, channel, 
                dataclient, invclient,
                outdir='.',  preprocessing=None,
                fileunit="year",
                **procparams):
        
        self.network = network
        self.station = station
        self.location = location
        self.channel = channel
        self.outdir = outdir
        self.client = dataclient
        self.invclient = invclient
        self.fileunit = fileunit
        self.logger = logging.getLogger(module_logger.name+
                            '.'+"GenericProcessor")
        self.logger.setLevel(logging.DEBUG)

        if preprocessing is None:
            preprocessing = util.process_stream
        self.preprocessing = preprocessing
        
        if "procparams" in procparams:
            self.proc_params = procparams['procparams']
            self.proc_params.update(**procparams)
        else:
            self.proc_params = ProcessingParameters(**procparams)

        ## Derived attributes:
        self.win_per_proclen = int(self.proc_params.proclen_seconds / 
                            self.proc_params.winlen_seconds)
        
        self.fname_fmt = util.FNAME_FMTS[self.fileunit]


    def __repr__(self) -> str:
        """
        Print attributes and their values
        """
        nsc = "{}.{}.{}.{}".format(self.network, self.station, 
                                    self.location, self.channel)
        d = "Data is sent to {}".format(self.outdir)

        sds = "Data client:\n{}".format(
            "; ".join(["{}={}".format(k, str(v)) for k,v 
                            in self.client.__dict__.items()]))
        inv = "Inventory client {}".format(
            "; ".join(["{}={}".format(k, str(v)) for k,v 
                            in self.invclient.__dict__.items()]))
        proc = "Processing settings:\n{}".format(self.proc_params)
        pprc = "Preprocessing of seismic data: {}".format(self.preprocessing.__name__)
        wpd = "Windows per proclen: {:d}".format(self.win_per_proclen)
        #s = ["{} = {}".format(k, str(v)) for k, v in self.__dict__.items()]
        return "\n".join([nsc, d, sds, inv, proc, pprc, wpd])




    def iter_time(self, starttime, endtime):
        """
        Iterator that yields start/end of intervals
        between ``starttime``, ``endtime`` depending
        on ``fileunit``.

        Example
        ---------
        If ``fileunit='year'``:

        .. code-block :: python

            gp = GenericProcessor(... fileunit="year")
            stime = UTC("2018-05-14")
            etime = UTC("2021-09-20")
            for s, e in gp.iter_time(stime, etime):
                print(s.date, e.date)
            >>>> 2018-05-14 2018-12-31
            >>>> 2019-01-01 2019-12-31
            >>>> 2020-01-01 2020-12-31
            >>>> 2021-01-01 2021-09-20

        """
        return util.TIME_ITERATORS[self.fileunit](starttime, endtime)


    def iter_nslc(self):
        """
        Iterator that yields each combination
        of network, station, location, channel from
        ``self._networks, self._stations, self._locations, self._channels``.
        These lists are set by ``expand_nslc()``
        """
        for n in self._networks:
            for s in self._stations:
                for l in self._locations:
                    for c in self._channels:
                        yield n, s, l, c


    def expand_nslc(self, starttime=None, endtime=None):
        """
        Expand network, station, location, channel 
        settings into lists.

        This may need to be customized to the data base.
        It sets lists of respective codes as attributes 
        ``_networks, __stations, _locations, _channels``
        """
        for attrib in ["network", "station", "channel", "location"]:
            val = self.__getattribute__(attrib)
            if isinstance(val, str):
                val = [val]
            self.__setattr__("_"+attrib+"s", val)



    def process(self, starttime, endtime, force_new_file=False):
        """
        Execute processing.

        This is the method you most likely want to call!

        Iterates through processing interval defined by ``fileunit``
        and calls NSCProcessor().process() for each NSLC combination.
        
        Manages creation of output files.
        """
        
        T0 = time.time()
        self.starttime = UTC(starttime)
        self.endtime = UTC(endtime)


        for _starttime, _endtime in self.iter_time(starttime, endtime):
            if any([self.fileunit == x for x in ["year", "month", "day"]]):
                _endtime = _endtime + 24*3600 - self.proc_params.proclen_seconds

            self.expand_nslc(_starttime, _endtime)
            for n, s, l, c in self.iter_nslc():
                nscproc = NSCProcessor(n, s, c, l, self.client,
                                self.invclient, procparams=self.proc_params)
                #print(nscproc)
                output = nscproc.process(_starttime, _endtime)
                output.trim_nan()

                fout = self.get_ofile(nscproc, _starttime, force_new_file)
                with fout:
                    output.insert_in_file(fout)


        walltime = timedelta(seconds=time.time()-T0)
        self.logger.info("Finished. Took %s h" % walltime)


    def create_ofile(self, nscprocessor, starttime):
        """
        Create output file.
        """
        ofilename = self.fname_fmt.format(
                        outdir=self.outdir,
                year=starttime.year, 
                month=starttime.month, 
                day=starttime.day,
                hour=starttime.hour,
                **nscprocessor.nsc_as_dict())
        
        fstime, fetime = self.allocate_file_start_end(starttime) 
        n_proclen = int((fetime + self.proc_params.proclen_seconds - fstime) / 
                        self.proc_params.proclen_seconds)

        self.logger.info("Creating output file %s" % ofilename)
        self.logger.info("Starttime=%s, endtime=%s, n_proclen=%s" %
                    (fstime, fetime, n_proclen))

        f = h5py.File(ofilename, "w")
        

        f.create_dataset("amplitudes", 
                shape=(n_proclen, self.win_per_proclen), 
                            fillvalue=np.nan)
        nfreqs = self.proc_params.nperseg // 2 + 1
        f.create_dataset("psds", 
                        shape=(n_proclen, self.win_per_proclen, nfreqs), 
                            fillvalue=np.nan)
        f.create_dataset("frequency_axis", shape=(nfreqs,), 
                            fillvalue=np.nan)

        util._create_hdf5_attribs(f, nscprocessor.stationcode,
                                fstime, fetime,
                                self.proc_params.amplitude_frequencies,
                                self.proc_params.winlen_seconds,
                                self.proc_params.proclen_seconds
                                )
    
        return f

    
    def get_ofile(self, nscprocessor, starttime, force_new_file=False):
        """
        Open or create new output file for results from an
        ``NSCProcessor``

        Parameters
        -----------
        nscprocessor : NSCProcessor
        starttime : UTCDateTime
            nominal starttime of file
        force_new_file : bool [False]
            If True, always create new, fresh file, i.e.
            override existing ones.
        
        Returns
        -------------
        f : h5py.File
            open file handler
        """
        
        ofilename = self.fname_fmt.format(
                        outdir=self.outdir,
                year=starttime.year, 
                month=starttime.month, 
                day=starttime.day,
                hour=starttime.hour,
                **nscprocessor.nsc_as_dict())

        if force_new_file:
            return self.create_ofile(nscprocessor, starttime)

        try:
            f = h5py.File(ofilename, "r+")
        except FileNotFoundError:
            f = self.create_ofile(nscprocessor, starttime)
        return f


    def allocate_file_start_end(self, starttime):
        """
        Determine datetime of start/end time of data
        covered in output file depending on ``fileunit``.

        Example
        ---------
        - if ``fileunit="year", starttime=UTC("2021-04-10")``: 
            ``stime=2021-01-01``,
            ``etime = 2021-12-31``
        - if ``fileunit="month, starttime=UTC("2021-04-10")``:
            ``stime = 2021-04-01``,
            ``etime = 2021-04-30``

        """
        if self.fileunit == "year":
            stime = UTC(starttime.year, 1,1,0,0,0)
            etime = (UTC(starttime.year, 12, 31, 0,0,0) + 
                    24*3600 - self.proc_params.proclen_seconds)
        elif self.fileunit == "month":
            stime = UTC(starttime.year, starttime.month, 1, 0,0,0)
            etime = (util.get_end_of_month(stime) + 
                    24*3600 - self.proc_params.proclen_seconds)
        elif self.starttime == "day":
            stime = UTC(starttime.year, starttime.month, starttime.day,
                        0, 0, 0)
            etime = stime + 3600*24 - self.proc_params.proclen_seconds
        elif self.fileunit == "hour":
            stime = UTC(starttime.year, starttime.month, starttime.day,
                        starttime.hour, 0, 0)
            etime = stime + 3600
        return stime, etime




class ProcessingParameters():
    """
    Bundle processing parameters.

    Parameters
    -----------------
    overlap : int [60]
        length of taper (overlap) if overlapping frames
        are necessary. In seconds.
    amplitude_frequencies : 2-tuple [(4,14)]
        min, max frequency of bandpass, applied before
        processing amplitudes
    nperseg : int [2048]
        samples per segment for FFT in Welch-PSD. 
        Passed to ``scipy.signal.welch()``
    winlen_seconds : int [3600]
        length of time window (=frame) over which amplitude
        and psd are computed. In seconds
    proclen_seconds : int [24*3600]
        length of data processed at once. In seconds. 
        Default is 1 day, corresponding to length of data
        in 1 sds-file
    """
    def __init__(self, **kwargs)-> None:
        # overlap, amplitude_frequencies, nperseg,
        #                 winlen_seconds, proclen_seconds
        for k, v in default_processing_params.items():
            if k in kwargs:
                self.__setattr__(k, kwargs[k])
            else:
                self.__setattr__(k, v)

    def __repr__(self) -> str:
        """
        Print attributes and their values
        """
        s = ["{} = {}".format(k, str(v)) for k, v in self.__dict__.items()]
        return "\n".join(s)

    def get_dict(self):
        """
        Get attributes as dict.
        """
        return self.__dict__
    

    def update(self, **kwargs):
        """
        Update an attribute
        """
        for k, v in kwargs.items():
            if k in self.__dict__:
                self.__setattr__(k, v)



class NSCProcessor():
    """
    Process one specific network-station-location-channel
    combination.

    Processing means extracting percentile of amplitude
    and power spectral density over some time window.

    Parameters
    -------------------
    netw : str
        network code
    stat : str
        station code
    chan : str
        channel code
    loc : str
        location code
    dataclient : obspy.clients
        client to obtain data from 
        (``dataclient.get_waveforms()``)
    invclient : obspy.clients
        client to get meta data from
        (``invclient.get_stations()``)
    **procparams : 
        processing parameters as keyword arguments
    """
    def __init__(self, netw, stat, chan, loc,
                    dataclient, invclient, **procparams):
        self.network = netw
        self.station = stat
        self.channel = chan
        self.location = loc
        self.dataclient = dataclient
        self.invclient = invclient

        self.logger = logging.getLogger(module_logger.name+
                            '.'+"NSCProcessor")
        self.logger.setLevel(logging.DEBUG)

        if "procparams" in procparams:
            self.processing_params = procparams['procparams']
            self.processing_params.update(**procparams)
        else:
            self.processing_params = ProcessingParameters(**procparams)
        
        self.stationcode = STATIONCODE.format(**self.nsc_as_dict())


    def __repr__(self) -> str:
        """
        Print attributes and their values
        """
        s = ["{} = {}".format(k, str(v)) for k, v in self.__dict__.items()]
        return "\n".join(s)


    def nsc_as_dict(self):
        """
        Return network, station, location, channel as dict.
        """
        d = dict(
            network = self.network,
            station = self.station,
            location = self.location,
            channel = self.channel            
        )
        return d


    def process(self, starttime, endtime, 
                preprocessing=None):
        """
        Compute noise amplitudes and spectra of raw data 
        within time range.

        Between start/end time, we load waveform data in 
        chunks of ``proclen_seconds``. Amplitude and PSD
        are computed for each frame of length 
        ``winlen_seconds`` that chunk. So the number of
        frames per chunk is 
        ``nf = proclen_seconds / winlen_seconds``

        Empty streams are replaced by ``np.nan``.

        Parameters
        -----------
        startdate : UTCDateTime 
        enddate : UTCDateTime
        preprocessing : func [None]
            function that takes obspy.Stream and returns 
            obspy.Trace. If ``None``, we use 
            ``dataqc.util.process_stream()`` which merges 
            and trims data, fills gaps with ``np.nan`` and
            removes mean.

        Returns
        ----------
        output : BaseProcessedData
            An object, which bundles amplitude and psd
            matrix, frequency axis for psd and processing
            parameters. 
            Shape of amplitude will be ``(n_proclens, nf)``,
            ``n_proclens`` being the number of chunks 
            between start/end time.
            Shape of psd will be ``(n_proclens, nf, nperseg//2+1)``.
            The last axis is the number of frequencies in psd.
        """
        self.starttime = starttime
        self.endtime = endtime
        if preprocessing is None:
            preprocessing = util.process_stream
        fmin, fmax = self.processing_params.amplitude_frequencies
        AMP = []
        PXX = []
        # frequency_axis = []
        start_after = 0  # counter for missing frames at beginning
        output = BaseProcessedData(starttime, endtime, 
                        self.stationcode,
                        self.processing_params.amplitude_frequencies, 
                        self.processing_params.winlen_seconds,
                        self.processing_params.proclen_seconds)
        starttime = starttime-self.processing_params.overlap
        inv = self.invclient.get_stations(
            starttime=starttime, endtime=endtime, level='response',
            **self.nsc_as_dict())
        self.logger.info("Processing %s" % self.stationcode)
        while starttime <= self.endtime - self.processing_params.overlap:
            endtime = (starttime + 
                        self.processing_params.proclen_seconds + 
                        2*self.processing_params.overlap)
            self.logger.debug("%s - %s" % (starttime, endtime))
            
            st = self.dataclient.get_waveforms(starttime=starttime, endtime=endtime, 
                                    **self.nsc_as_dict())
            try:
                tr = preprocessing(st, inv, starttime, endtime, 
                        self.processing_params.sampling_rate)
            # No data in trace:
            except IndexError:
                # self.logger.debug("No data for %s" % UTC((starttime + overlap).date))
                self.logger.info("No data for %s" % starttime)
                starttime = starttime + self.processing_params.proclen_seconds
                
                # Shape of output is determined by seismic data properties like
                # sampling rate. While we don't have data, we cannot know the
                # shape so we only count how many outputs we missed.
                if len(AMP) == 0:
                    start_after += 1
                    self.starttime = starttime + self.processing_params.overlap
                    output.startdate = UTC((starttime + self.processing_params.overlap))
                else:
                    AMP.append(np.ones(AMP[-1].shape)*np.nan)
                    PXX.append(np.ones(PXX[-1].shape)*np.nan)
                continue

            # Get some numbers
            #sr = tr.stats.sampling_rate
            nf = int(self.processing_params.proclen_seconds/
                     self.processing_params.winlen_seconds)
            #proclen_samples = proclen * sr
            winlen_samples = int(self.processing_params.winlen_seconds * self.processing_params.sampling_rate)
            
            # Spectra
            data = util.get_adjacent_frames(tr, 
                    starttime+self.processing_params.overlap, 
                    nf, winlen_samples)
            frequency_axis, P = welch(data, fs=self.processing_params.sampling_rate, 
                    nperseg=self.processing_params.nperseg, axis=1)
                
            # Amplitude
            prctl = util.get_amplitude(tr, 
                        starttime+self.processing_params.overlap, 
                        fmin, fmax,
                        self.processing_params.overlap, 
                        winlen_samples, nf)
            
            AMP.append(prctl) #amp[1:-1])
            PXX.append(P) # pxx[1:-1,:])

            starttime = starttime + self.processing_params.proclen_seconds

        if len(AMP) > 0:
            output.set_data(np.array(AMP),
                            np.array(PXX),
                            np.array(frequency_axis))
        return output



class BaseProcessedData():
    """
    Handle processed amplitudes and spectra together
    with some meta data. Read/write HDF5

    Parameters
    ------------------
    startdate : UTCDateTime [None]
        startdate of data 
    enddate : UTCDateTime [None],
    stationcode : str ["...."]
        station code of format 
        "{network}.{station}.{location}.{channel}" 
    amplitude_frequencies : 2-tuple [(None,None)]
        min, max frequency between which data was filtered
        for amplitude extraction
    seconds_per_window : int [None]
        seconds per frame over which amplitude and psd were
        computed
    proclen_seconds : int [None]
        length of data in seconds that were processed at once
    """
    def __init__(self,startdate=None, enddate=None,
                stationcode="....", 
                amplitude_frequencies=(None,None),
                seconds_per_window=None, proclen_seconds=None):
        self.amplitudes = None
        self.psds = None
        self.frequency_axis = None
        self.stationcode = stationcode
        self.amplitude_frequencies = amplitude_frequencies
        self.seconds_per_window = seconds_per_window
        if startdate is not None:
            startdate = UTC(startdate)
        self.startdate = startdate
        if enddate is not None:
            enddate = UTC(enddate) 
        self.enddate = enddate
        self.proclen_seconds = proclen_seconds
        
        self.logger = logging.getLogger(module_logger.name+
                            '.'+"BaseProcessedData")
        self.logger.setLevel(logging.DEBUG)

    def get_nslc(self):
        """
        Split station code into network, station, 
        location, channel
        """
        return self.stationcode.split('.')


    def from_file(self, fname):
        """
        Read properties and data from HDF5 file.
        """
        with h5py.File(fname, 'r') as fin:
            self.amplitude_frequencies = fin.attrs['amplitude_frequencies']
            self.seconds_per_window = fin.attrs['seconds_per_window']
            self.startdate = UTC(*fin.attrs['starttime'])
            self.enddate = UTC(*fin.attrs['endtime'])
            self.stationcode = fin.attrs['stationcode']
            self.amplitudes = np.array(fin['amplitudes'])
            self.frequency_axis = np.array(fin['frequency_axis'])
            self.psds = np.array(fin['psds'])
            self.proclen_seconds = fin.attrs['seconds_per_proclen']
        return self

    
    def to_file(self, outdir="."):
        """
        Write content to HDF5 file. 

        Filename is created from stationcode, startdate
        and enddate.
        """
        fname = "{}_{}_{}.hdf5".format(
                self.stationcode, self.startdate.date, 
                  self.enddate.date)
        fname = os.path.join(outdir, fname)
        if not self.has_data():
            self.logger.warning("No data for file %s. Skipping" % fname)
            return

        self.logger.info("Writing data to %s" % fname)

        with h5py.File(fname, 'w') as fout:
            fout.create_dataset('amplitudes', data=self.amplitudes)
            fout.create_dataset('psds', data=self.psds)
            fout.create_dataset('frequency_axis', data=self.frequency_axis)
            
            util._create_hdf5_attribs(fout, self.stationcode, 
                                self.startdate, self.enddate,
                                self.amplitude_frequencies,
                                self.seconds_per_window,
                                self.proclen_seconds)
          

    def insert_in_file(self, fout):
        """
        Insert data in existing hdf5-file.

        From our own start/end times and those in file,
        and the frame and processing lengths, we determine
        the indices at which to insert the new data.

        Parameters
        -------------
        fout : h5py.File
            file-object, ready for write in which to insert 
            the data
        
        Note
        -------
        Raises ValueError if
        - we have no data attributes
        - starttime in file > self.startdate

        Warning
        -----------
        Overrides existing data in file.
        """
        if not self.has_data():
            self.logger.warn("output has no data to insert")
            return
            #raise RuntimeWarning("output has no data to insert")

        fstime = UTC(*fout.attrs["starttime"])

        # We only insert if file starts at or before data
        # Or should we cut of the extending part on our own?
        if fstime > self.startdate:
            msg = "Targeted file starts after data. Cannot insert."
            self.logger.error(msg)
            raise ValueError(msg)

        i = ((self.startdate - fstime) / 
                    self.proclen_seconds)
        j = ((self.enddate + self.proclen_seconds - fstime) / 
                        self.proclen_seconds)
       
        i, j = int(i), int(j)
        self.logger.debug("starttime: %s" % self.startdate)
        self.logger.debug("endtime: %s" % self.enddate)
        self.logger.debug("Amplitude matrix shape: %s" % 
            ", ".join([str(s) for s in self.amplitudes.shape]))
        self.logger.debug("Target shape: (%s)" % 
            ", ".join([str(s) for s in fout["amplitudes"][i:j,:].shape]))
        self.logger.debug("Total shape of target %s" %
            ", ".join([str(s) for s in fout["amplitudes"][:].shape]))
        self.logger.debug("Targeted index range %s:%s" % (i,j))

        fout["amplitudes"][i:j,:] = self.amplitudes
        fout["psds"][i:j,:,:] = self.psds
        fout["frequency_axis"][:] = self.frequency_axis

    
    def set_data(self, amplitudes, psds, freqs):
        """
        Add amplitudes, psds and frequency axis data.

        Raises RuntimeWarning if shapes are inconsistent, 
        e.g. ``psds.shape[-1] != freqs.size``

        Parameters
        --------------
        amplitudes : ndarray, 2d 
            amplitude data of shape ``(n_proclen, n_winlen)``
        psds : ndarray, 3d
            power spectral densities of shape
            ``(n_proclen, n_winlen, n_frequencies)``
        freqs : ndarray, 1d
            frequency axis corresponding to psds of shape
            ``(n_frequencies,)``
        """
        # make sure shapes of data are consistent
        if amplitudes.shape != psds.shape[:-1]:
            msg = ("Added data has inconsistent shapes!"+
                "amplitudes: %s; psds %s" % (amplitudes.shape, psds.shape))
            self.logger.warn(msg)
            raise RuntimeWarning(msg)
        if psds.shape[-1] != freqs.size:
            msg = ("Added frequencies and psds have inconsisent size!" +
                    "psds: {}".format(psds.shape[-1] + 
                    ", freqs: {}".format(freqs.size))
                    )
            self.logger.warn(msg)       
            raise RuntimeWarning(msg)

        # check if proclen fits
        nwin = amplitudes.shape[-1]
        inp_proclen = nwin*self.seconds_per_window
        
        if not np.isclose(inp_proclen, self.proclen_seconds):
            msg = ("Inconsistent length of processing slice." + 
                "{:g} s in added data ".format(inp_proclen) +
                "vs {:g} s allocated".format(self.proclen_seconds)
            )
            self.logger.warn(msg)
            raise RuntimeWarning(msg)

        self.amplitudes = amplitudes
        self.psds = psds
        self.frequency_axis = freqs


        

    def has_data(self):
        """
        Check if amplitude and/or spectral data are available
        """
        if self.amplitudes is None and self.psds is None:
            return False
        else:
            return True

    
    def trim_nan(self):
        """
        Remove entire nan-rows (days) at both ends and 
        adjust start/enddate.
        """
        # If no data, we can exit right away
        if not self.has_data(): 
            return

        n = 0
        while len(self.amplitudes) > 0:
            if np.all(np.isnan(self.amplitudes[0,:])):
                self.amplitudes = np.delete(self.amplitudes, 0, axis=0)
                n = n+1
            else:
                break
        self.psds = self.psds[n:,:,:]

        m = 0
        while len(self.amplitudes) > 0:
            if np.all(np.isnan(self.amplitudes[-1,:])):
                self.amplitudes = np.delete(self.amplitudes, -1, axis=0)
                m = m+1
            else:
                break
        if m > 0:
            self.psds = self.psds[:-m,:,:]

        # Check if any data is left
        if len(self.amplitudes) == 0:
            self.amplitudes = None
            self.psds = None
        # if yes, adjust start/enddate
        else:
            nw = self.amplitudes.shape[1]
            self.startdate = self.startdate + n*nw*self.seconds_per_window
            self.enddate = self.enddate - m*nw*self.seconds_per_window
            

    def extend_from_file(self, file):
        """
        Read new result data from HDF5-file and add to
        existing.

        Combines ``.from_file()`` and ``.extend()``

        Parameters
        --------------
        file : str, path-like
            Name of hdf5-file. Must have appropriate structure,
            obviously.
        """
        # If there is no data yet, we can simply read it from file
        if (self.amplitudes is None and
            self.psds is None):
            self.from_file(file)
            return
        
        # If we already have some data, we need to insert the
        # new data at the right place and fill potential gaps.
        # Read the new data
        self += BaseProcessedData().from_file(file)
        # new = BaseProcessedData()
        # new.from_file(file)
        # self.extend(new)
        
        
    def __iadd__(self, new):
        """
        Use += notation to add more data
        """
        if not isinstance(new, BaseProcessedData):
            raise TypeError("You can only add another BaseProcessedData")
        self.extend(new)
        return self
        

    def extend(self, new):
        """
        Add a new output object.

        Adjusts start/endtimes, removes leading/trailing
        Nans. Should only work if processing parameters,
        frequency axis are compatible.
        """
        # Get total number of days to get new array sizes
        
        tmin = min(self.startdate, new.startdate)
        tmax = max(self.enddate, new.enddate)
        #days = timedelta(seconds=tmax-tmin).days+1

        nrows = int((tmax-tmin+self.proclen_seconds) / self.proclen_seconds)
        
        # If shapes, processing parameters, frequency axis or 
        # stationcodesare inconsistent, we get an error here
        self._check_sanity(new)
        
        # Initialize data containers for merged data
        amps_shp, psds_shp = (list(self.amplitudes.shape), 
                                list(self.psds.shape))
        amps_shp[0] = nrows #days
        psds_shp[0] = nrows #days
        
        new_amps = np.ones(amps_shp)*np.nan
        new_psds = np.ones(psds_shp)*np.nan
        
        # Insert the data. We ignore overlaps here. New data
        # overwrites existing data if they overlap
        for d in [self, new]:
            #i = timedelta(seconds=d.startdate-tmin).seconds
            i = int((d.startdate - tmin) / d.proclen_seconds )
            n = len(d.amplitudes)
            
            new_amps[i:i+n,:] = d.amplitudes
            new_psds[i:i+n,:,:] = d.psds
            
        self.amplitudes = new_amps
        self.psds = new_psds
        self.startdate = UTC(tmin.date)
        self.enddate = UTC(tmax.date)
        self.trim_nan()


    def _check_sanity(self, new):
        """
        Check consistency of new and existing data.

        Runs:
        - ``_chech_shapes``
        - ``_check_frequencies``
        - ``_check_codes``

        Throws IOError is any one of them fails.
        """
        # Any of these methods raises IOError
        # if inconsistencies are found
        self._check_shapes(new)
        self._check_frequencies(new)
        self._check_codes(new)
        return True
    
        
    def _check_shapes(self, new):
        """
        Checks if number of processing windows
        and frequency axis of psd and amplitude are
        compatible.
        Throws IOError if not.
        """
        amps_shp = self.amplitudes.shape
        psds_shp = self.psds.shape
        if (amps_shp[-1] == new.amplitudes.shape[-1] and
            psds_shp[1:] == new.psds.shape[1:]):
            return amps_shp, psds_shp
        else:
            self.logger.error("Data in files have inconsistent shapes!")
            raise IOError("Data in files have inconsistent shapes!")
            #return False, False

    
    def _check_frequencies(self, new):
        """
        Check if frequency axis are similar using np.isclose().
        Throws IOError if not.
        """
        if not np.all(np.isclose(self.frequency_axis, new.frequency_axis)):
            self.logger.error("Frequency axis are different!")
            raise IOError("Frequency axis are different!")
        else:
            return True

    
    def _check_codes(self, new):
        """
        Check if station codes are identical. 
        Throws IOError if not.
        """
        if self.stationcode != new.stationcode:
            self.logger.error("Station codes are different! " +
                "Cannot extend existing data.")
            raise IOError("Station codes are different! " +
                "Cannot extend existing data.")
        else:
            return True

    
    def _plotly_amplitudes(self):
        z = self.amplitudes #np.clip(AMP, None,  a_max=5.5)
        sh_0, sh_1 = z.shape
        y, x = np.linspace(0, sh_0-1, sh_0), np.linspace(0, sh_1-1, sh_1)
        fig = go.Figure(data=[go.Surface(z=z, x=x, y=y, name='amplitude', cmin=2, cmax=None)])
        fig.update_layout(title='75%-amplitude', autosize=True,
                          width=800, height=500,
                          scene=dict(aspectmode='manual', aspectratio=dict(x=1, y=2, z=0.5))
                          #margin=dict(l=65, r=50, b=65, t=90)
                         )
        fig.show()


    def plot_amplitudes(self, tax=None, ax=None, 
        labelinc=1):
        """
        Plot amplitude matrix using matplotlib.

        Parameters
        ---------------
        tax : ndarray
            labels of time axis. If not given, labels
            are created from start/endtime
        ax : matplotlib.axes
            Axes to use. If not given, new figure is 
            created.

        Returns
        ----------
        ax : matplotlib.axes
            Axes containing the plot.
        
        """
        if ax is None:
            fig, ax = plt.subplots(1,1)

        dtflag, dtinc = util.choose_datetime_inc(self.proclen_seconds)


        if tax is None:
            tax = np.arange(self.startdate, 
                            self.enddate+self.proclen_seconds, 
                            dtinc,
                            dtype='datetime64[{}]'.format(dtflag))
        yticks = np.arange(0, len(tax), labelinc)

        im = ax.imshow(self.amplitudes)
        ax.set_title("Amplitude data matrix");
        ax.set_xlabel('windows')
        ax.set_ylabel('proclen');
        ax.set_yticks(yticks)
        ax.set_yticklabels(tax);
        plt.colorbar(im, ax=ax)
        return ax


    def plot_psds(self, func=None, tax=None, ax=None, 
            N_freqlabels=8):
        """
        Plot PSD matrix as (time, frequency) using 
        matplotlib.

        First two axes (proclens, winwlens) are flattened.

        Parameters
        ---------------
        func: function
            Apply to PSD data before plotting, e.g. np.log.
        tax : ndarray
            labels of time axis. If not given, labels
            are created from start/endtime
        ax : matplotlib.axes
            Axes to use
        N_freqlabels : int
            Number of tick labels of frequency axis

        Returns
        ----------
        ax : matplotlib.axes
            Axes containing the plot.
        
        """
        
        if ax is None:
            fig, ax = plt.subplots(1,1)
        dtflag, dtinc = util.choose_datetime_inc(
                            self.proclen_seconds)

        i, j, k = self.psds.shape
        M = self.psds.reshape(i*j, k)
        if func is not None:
            M = func(M)
        if tax is None:
            tax = np.arange(self.startdate, 
                            self.enddate+self.proclen_seconds,
                            dtinc,
                        dtype='datetime64[{}]'.format(dtflag))
        yticks = np.arange(0,len(M), j)
        cax = ax.imshow(M)
        ax.set_yticks(yticks);
        ax.set_yticklabels(tax);
        plt.colorbar(cax, ax=ax) 

        xticks_inc = int(k / N_freqlabels)
        xticks = np.arange(0, self.frequency_axis.size, xticks_inc)
        xticklabels = np.round(self.frequency_axis, 2)[::xticks_inc]
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels)    
        ax.set_xlabel('Hz')

        ax.set_title("Power spectral densities")
        return ax

    
    def __repr__(self):
        """
        String representation of properties.
        """
        nsc = "Results for {}\n".format( self.stationcode) + 20*"-"
        s = "Starttime: {}".format(self.startdate)
        e = "Enddate: {}".format(self.enddate)
        d = "Days: {:d}".format(self.amplitudes.shape[0])
        shp1 = "Amplitude shape = {}".format(self.amplitudes.shape)
        shp2 = "PSD shape = {}".format(self.psds.shape)
        pr1 = "Seconds per window = {:g}".format(self.seconds_per_window)
        pr2 = "Amplitude for {:g} - {:g} Hz".format(*self.amplitude_frequencies)

        return "\n".join([nsc, s,e,d,shp1,shp2, pr1, pr2])



