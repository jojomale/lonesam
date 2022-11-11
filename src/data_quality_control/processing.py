"""
**OLD !!!**

Initial classes and routines for extracting amplitudes and
computing PSDs from seismic data.

"""
import configparser
from datetime import timedelta
import time
from glob import glob
import os.path
from types import prepare_class
import numpy as np

from scipy.signal import welch, get_window

from obspy.signal.filter import bandpass
from obspy.clients.filesystem.sds import Client
from obspy.clients.fdsn import RoutingClient
from obspy.core import UTCDateTime as UTC
# from obspy.signal import util

import plotly.graph_objects as go

import h5py

import logging
logger = logging.getLogger('processing')
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)  # set level
cformatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                            datefmt='%y-%m-%d %H:%M:%S')
ch.setFormatter(cformatter)
if not logger.hasHandlers():
    logger.addHandler(ch)

#from .processing import ProcessingParameters as PP

STATIONCODE = "{network}.{station}.{location}.{channel}"

default_settings = dict(network = '*',
    station = '*',
    channel = '*',
    overlap = 60,
    amplitude_frequencies = (4,14),
    nperseg = 2048,
    winlen_in_s = 3600,
    proclen = 24*3600,
    sds_root = os.path.abspath('.'),
    inventory_routing_type = "eida-routing",
    sds_client_dict = {},
    outdir = ".")

default_processing_params = dict(
    overlap = 60,
    amplitude_frequencies = (4,14),
    nperseg = 2048,
    winlen_seconds = 3600,
    proclen_seconds = 24*3600,
)

def get_default_config(outname=None):
    """
    Return default settings as ConfigParser.

    Write to outname if not None.
    """
    config_from_dict = configparser.ConfigParser()

    sec = "PROCESSING PARAMETERS"
    config_from_dict[sec] = {}
    for k, v in default_settings.items():
        print(k, v)
        if isinstance(v, str):
            v = v
        elif isinstance(v, int) or isinstance(v, float):
            v = str(v)
        elif isinstance(v, list) or isinstance(v, tuple):
            v = ','.join([str(vi) for vi in v])
        elif isinstance(v, dict):
            config_from_dict.read_dict({k.upper(): v})
            continue
        config_from_dict.set(sec, k, v)
        
    if outname is not None:
        with open(outname, 'w') as fp:
            config_from_dict.write(fp)
    return config_from_dict



class ProcessingParameters():
    def __init__(self, **kwargs)-> None:
        # overlap, amplitude_frequencies, nperseg,
        #                 winlen_seconds, proclen_seconds
        for k, v in default_processing_params.items():
            if k in kwargs:
                self.__setattr__(k, kwargs[k])
            else:
                self.__setattr__(k, v)
        # self.overlap = overlap
        # self.amplitude_frequencies = amplitude_frequencies
        # self.nperseg = nperseg
        # self.winlen_seconds = winlen_seconds
        # self.proclen_seconds = proclen_seconds

    def __repr__(self) -> str:
        """
        Print attributes and their values
        """
        s = ["{} = {}".format(k, str(v)) for k, v in self.__dict__.items()]
        return "\n".join(s)

    def get_dict(self):
        return self.__dict__
    

    def update(self, **kwargs):
        for k, v in kwargs.items():
            if k in self.__dict__:
                self.__setattr__(k, v)


class NSCProcessor():
    def __init__(self, netw, stat, chan, loc,
                    dataclient, invclient, **procparams):
        self.network = netw
        self.station = stat
        self.channel = chan
        self.location = loc
        self.dataclient = dataclient
        self.invclient = invclient

        if "procparams" in procparams:
            self.processing_params = procparams['procparams']
            self.processing_params.update(**procparams)
        else:
            self.processing_params = ProcessingParameters(**procparams)

        # self.processing_params = ProcessingParameters(**procparams)
        # if isinstance(procparams, PP):
        #     self.processing_params = procparams
        # elif isinstance(procparams, dict):
        #     self.processing_params = ProcessingParameters(**procparams)
        # else:
        #     raise ValueError("procparams must be dict or ProcessingParameters." + 
        #                     "Not {}".format(type(procparams)))
        self.stationcode = STATIONCODE.format(**self.nsc_as_dict())


    def __repr__(self) -> str:
        """
        Print attributes and their values
        """
        s = ["{} = {}".format(k, str(v)) for k, v in self.__dict__.items()]
        return "\n".join(s)


    def nsc_as_dict(self):
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
        for a specific network, station, channel combination 
        within time range.

        Parameters
        -----------
        startdate : UTCDateTime 
        enddate : UTCDateTime
        
        """
        self.starttime = starttime
        self.endtime = endtime
        if preprocessing is None:
            preprocessing = process_stream
        fmin, fmax = self.processing_params.amplitude_frequencies
        AMP = []
        PXX = []
        frequency_axis = []
        start_after = 0  # counter for missing frames at beginning
        output = BaseProcessedData(starttime, endtime, 
                        self.stationcode,
                        self.processing_params.amplitude_frequencies, 
                        self.processing_params.winlen_seconds,)
        starttime = starttime-self.processing_params.overlap
        inv = self.invclient.get_stations(
            starttime=starttime, endtime=endtime, level='response',
            **self.nsc_as_dict())
        logger.info("Processing %s" % self.stationcode)
        while starttime <= self.endtime - self.processing_params.overlap:
            endtime = (starttime + 
                        self.processing_params.proclen_seconds + 
                        2*self.processing_params.overlap)
            logger.debug("%s - %s" % (starttime, endtime))
            
            st = self.dataclient.get_waveforms(starttime=starttime, endtime=endtime, 
                                    **self.nsc_as_dict())
            try:
                tr = preprocessing(st, inv, starttime, endtime)
            # No data in trace:
            except IndexError:
                # logger.debug("No data for %s" % UTC((starttime + overlap).date))
                logger.info("No data for %s" % starttime)
                starttime = starttime + self.processing_params.proclen_seconds
                
                # Shape of output is determined by seismic data properties like
                # sampling rate. While we don't have data, we cannot know the
                # shape so we only count how many outputs we missed.
                if len(AMP) == 0:
                    start_after += 1
                    self.starttime = starttime + self.processing_params.overlap
                    output.startdate = UTC((starttime + self.processing_params.overlap).date)
                else:
                    AMP.append(np.ones(AMP[-1].shape)*np.nan)
                    PXX.append(np.ones(PXX[-1].shape)*np.nan)
                continue

            # Get some numbers
            sr = tr.stats.sampling_rate
            nf = int(self.processing_params.proclen_seconds/
                     self.processing_params.winlen_seconds)
            #proclen_samples = proclen * sr
            winlen_samples = int(self.processing_params.winlen_seconds * sr)
            
            # Spectra
            data = get_adjacent_frames(tr, starttime+self.processing_params.overlap, nf, 
                                    winlen_samples)
            frequency_axis, P = welch(data, fs=sr, 
                    nperseg=self.processing_params.nperseg, axis=1)
                
            # Amplitude
            prctl = get_amplitude(tr, starttime+self.processing_params.overlap, 
                                fmin, fmax,
                                self.processing_params.overlap, winlen_samples, nf)
            
            AMP.append(prctl) #amp[1:-1])
            PXX.append(P) # pxx[1:-1,:])

            starttime = starttime + self.processing_params.proclen_seconds

        if len(AMP) > 0:
            output.set_data(np.array(AMP),
                            np.array(PXX),
                            np.array(frequency_axis))
        # self.amps = np.array(AMP)
        # self.psds = np.array(PXX)
        # self.frequency_axis = np.array(frequency_axis)
        # self._start_after = start_after
        # self.n_proclens = len(self.amps)

        # return self.amps, self.psds, self.frequency_axis, self._start_after
        return output





class SDSDataBaseProcessor():
    """
    Class to manage processing of SDS data base

    NSC can contain wildcards
    """
    OFILENAME_FMT = "{outdir}"+os.path.sep+"{network}.{station}.{channel}_{year}.hdf5"
    PROCLEN_SECONDS = 24*3600
    PROCLENS_PER_FILE = 366

    def __init__(self, network, station, channel, 
                sds_root, inventory_routing_type, 
                outdir='.', sdsdict={}, preprocessing=None,
                **procparams):
        
        self.network = network
        self.station = station
        self.channel = channel
        self.outdir = outdir
        self.client = Client(sds_root, **sdsdict)
        self.invclient = RoutingClient(inventory_routing_type)
        
        if preprocessing is None:
            preprocessing = process_stream
        self.preprocessing = preprocessing
        
        procparams['proclen_seconds'] = self.PROCLEN_SECONDS
        if "procparams" in procparams:
            self.proc_params = procparams['procparams']
            self.proc_params.update(**procparams)
        else:
            self.proc_params = ProcessingParameters(**procparams)

        ## Derived attributes:
        self.win_per_day = int(self.PROCLEN_SECONDS / 
                            self.proc_params.winlen_seconds)
        
        
    def __repr__(self) -> str:
        """
        Print attributes and their values
        """
        nsc = "{}.{}.{}".format(self.network, self.station, self.channel)
        d = "Data is sent to {}".format(self.outdir)

        sds = "SDS-Client for data:\n{}".format(
            "; ".join(["{}={}".format(k, str(v)) for k,v 
                            in self.client.__dict__.items()]))
        inv = "Inventory from Routing client {}".format(self.invclient._url)
        proc = "Processing settings:\n{}".format(self.proc_params)
        pprc = "Preprocessing of seismic data: {}".format(self.preprocessing.__name__)
        wpd = "Windows per day: {:d}".format(self.win_per_day)
        #s = ["{} = {}".format(k, str(v)) for k, v in self.__dict__.items()]
        return "\n".join([nsc, d, sds, inv, proc, pprc, wpd])




    def get_nsc(self, starttime, endtime):
        """
        Expand network, station, channel attributes (which can contain
        wildcards) to determine all unique networks, stations and channels 
        available within in time range.

        Parameters
        --------------
        starttime : UTCDatetime
            begin of time range, we only use the year
        endtime : UTCDateTime
            end of time range, we only use the year

        """

        networks = set()
        stations = set()
        channels = set()
        years = np.arange(starttime.year, endtime.year+1)
        for y in years:
            logger.debug("%s" % str(y))
            for n in glob(os.path.join(self.client.sds_root, 
                                        str(y), self.network)):
                nw = os.path.split(n)[-1]
                networks.add(nw)
                logger.debug("%s" % nw)
                for s in glob(os.path.join(n, self.station)):
                    sn = os.path.split(s)[-1]
                    stations.add(sn)
                    logger.debug("%s" % sn)
                    for c in glob(os.path.join(s, self.channel)+
                                    '.{}'.format(self.client.sds_type)):
                        ch = os.path.split(c)[-1].split('.')[0]
                        channels.add(ch)
                        logger.debug("%s" % ch)
                        
        return networks, stations, channels



    def process(self, startdate, enddate):
        """
        Process data between given time range and write
        results to disk.

        We generally assume, that all other parameters for 
        the database search (network, station, channel, etc)
        and processing parameters have been set before during
        initialization. Nevertheless, they can be overridden
        here.

        Parameters
        -------------
        startdate : str, UTCDateTime
            begin of time range for analysis. Can be UTCDateTime
            object or string readable by UTCDateTime
        enddate : str, UTCDateTime
            end of time range, same as startdate
        outdir : str, path-like
            directory where to write results to
        processing_kwargs :
            see RawDataProcessor
        """
        T0 = time.time()
        self.startdate = UTC(startdate)
        self.enddate = UTC(enddate)

        #self._check_set_attr(processing_kwargs)

        # # else override outdir attrib if given here,
        # # Could include it in processing_kwargs
        # if outdir is not None:
        #     self.outdir = outdir

        # Init loop counters
        _startdate = startdate
        _enddate = UTC("{:d}-12-31".format(_startdate.year))
        while _startdate < enddate:
            # Check if end-of-year date is larger than acual enddate
            if _enddate > enddate: 
                _enddate = enddate
                logger.debug("Reset enddate to %s" % _enddate)

            # Process 1 year (or less)
            logger.info("\nProcessing %s - %s" % (_startdate, _enddate))
            networks, stations, channels = self.get_nsc(_startdate, _enddate)

            for n in networks:
                for s in stations:
                    for c in channels:
                        nscproc = NSCProcessor(n, s, c, "", self.client,
                                self.invclient, procparams=self.proc_params)
                        #print(nscproc)
                        output = nscproc.process(_startdate, _enddate)
                        output.trim_nan()

                        fout = self.get_ofile(nscproc, _startdate.year)
                        with fout:
                            i = output.startdate.julday-1
                            j = output.enddate.julday
                            fout["amplitudes"][i:j,:] = output.amplitudes
                            fout["psds"][i:j,:,:] = output.psds
                            fout["frequency_axis"][:] = output.frequency_axis

            # Proceed loop
            # If we work with dates here, we may avoid problems with leap years
            _startdate = UTC("{:d}-01-01".format(_startdate.year+1))
            _enddate = UTC("{:d}-12-31".format(_startdate.year))
            logger.debug("New _startdate = %s" % _startdate)
            logger.debug("New _enddate = %s" % _enddate)
        walltime = timedelta(seconds=time.time()-T0)
        logger.info("Finished. Took %s h" % walltime)




    def create_ofile(self, nscprocessor, year):
        ofilename = self.OFILENAME_FMT.format(outdir=self.outdir,
                year=year, **nscprocessor.nsc_as_dict())
        
        year = int(year)
        f = h5py.File(ofilename, "w")
        f.create_dataset("amplitudes", 
                shape=(self.PROCLENS_PER_FILE, self.win_per_day), 
                            fillvalue=np.nan)
        nfreqs = self.proc_params.nperseg // 2 + 1
        f.create_dataset("psds", 
                        shape=(self.PROCLENS_PER_FILE, self.win_per_day, nfreqs), 
                            fillvalue=np.nan)
        f.create_dataset("frequency_axis", shape=(nfreqs,), 
                            fillvalue=np.nan)

        _create_hdf5_attribs(f, nscprocessor.stationcode,
                                UTC(year, 1, 1, 0, 0, 0),
                                UTC(year, 12, 31, 0, 0, 0),
                                self.proc_params.amplitude_frequencies,
                                self.proc_params.winlen_seconds)
        # f.attrs.create('stationcode', nscprocessor.stationcode)
        # f.attrs.create('starttime',
        #     np.array([int(year), 1, 1, 0, 0, 0]))
        # f.attrs.create('endtime',
        #     np.array([int(year), 12, 31, 0, 0, 0]))
        # f.attrs.create('amplitude_frequencies', 
        #                 np.array(self.proc_params.amplitude_frequencies))
        # f.attrs.create('seconds_per_window', 
        #                 self.proc_params.winlen_seconds)
        return f

    
    def get_ofile(self, nscprocessor, year):
        ofilename = self.OFILENAME_FMT.format(outdir=self.outdir,
                year=year, **nscprocessor.nsc_as_dict())
        try:
            f = h5py.File(ofilename, "r+")
        except FileNotFoundError:
            f = self.create_ofile(nscprocessor, year)
        return f







class RawDataProcessor():
    """
    Compute amplitude levels and power spectral densities from
    raw seismic data in an sds-filesystem and store in HDF5-files.

    Parameters
    -------------
    network : str ['*']
        network(s), can contain Unix-style wildcards, 
        e.g. "G*" or "G[KR]"
    station : str ['*']
        stations(s) to process, can contain Unix-style
        wildcards
    channel : ['*']
        channels(s), see above
    overlap : int [60]
        seconds by which daily files are extended to
        accomodate filter effects etc. Also used, if 
        overlapping windows are necessary during 
        processing of a single file due to data gaps.
    amplitude_frequencies : tuple of 2 floats [(4,14)]
        lower and upper frequency for bandpass filter,
        which is applied before amplitude analysis
    nperseg : int [2048]
        number of samples in segment for spectral estimation,
        passed to scipy.signal.welch
    winlen_in_s : [3600]
        window length in seconds over which amplitude and 
        spectra are computed
    proclen : [24*3600]
        length of data stream to process at once. Should be 
        the nominal length of the files in data base. 
        Usually 1 day.
    outdir : str, pathlike ["."]
        where to put processed data
    sds_root : str, path-like ['.']
        root directory of sds-filesystem, passed to
        obspy.clients.filesystem.sds.Client.
    inventory_routing_type : ["eida-routing"]
        routing client over which to get inventories,
        one of "eida-routing" or "iris-federator",
        passed to obspy.clients.fdsn.RoutingClient
    sds_client_dict : dict [{}]
        additional keyword arguments for 
        obspy.clients.filesystem.sds.Client
    """
    def __init__(self, **kwargs):
      
        # Set the intended attributes from default_settings
        for k, v in default_settings.items():
            self.__setattr__(k, v)

        # Override attributes from configfile if given
        #print(kwargs)
        try:
            v = kwargs.pop('configfile')
            if v is None:
                raise KeyError
            self.from_config(v)
            self.configfile = v
        except KeyError:
            pass
            
        # Override attributes directly given as argument
        self._check_set_attr(kwargs)

        self.sdsclient = Client(self.sds_root)
        self.invclient = RoutingClient(self.inventory_routing_type)
        self.processed_codes = set()


    def _check_set_attr(self, kwargs):
        # Override attributes directly given as argument
        for k, v in kwargs.items():
            try:
                self.__getattribute__(k)
                logger.debug("Overriding %s to %s" %
                    (k, str(v)))
            except AttributeError:
                raise AttributeError("Unknown keyword {}.".format(k))
            self.__setattr__(k, v)

    
    
    def from_config(self, configfile):

        if not os.path.isfile(configfile):
            raise FileNotFoundError("No config-file %s " % 
                    configfile)

        logger.info("Reading settings from %s" % configfile)
        config = configparser.ConfigParser()
        config.read(configfile)
        
        sec = "PROCESSING PARAMETERS"
        # Get strings
        str_kws = [k for k, v in default_settings.items() 
           if isinstance(v, str)]
        
        for k in str_kws:
            v = config[sec].get(k) 
            self.__setattr__(k, v)
        
        # Get number parameters (we expect only ints)
        num_kws = [k for k,v in default_settings.items()
                   if any([isinstance(v, int), 
                          isinstance(v, float)])]
        for k in num_kws:
            v = config[sec].getint(k)
            self.__setattr__(k, v)
        
        # Get tuples
        k = 'amplitude_frequencies'
        v = [float(i) for i in 
             config[sec][k].split(',')]
        self.__setattr__(k, v)
        
        sec = "SDS_CLIENT_DICT"
        sds_client_dict = {}
        for k in ['sds_type', 'format']:
            v = config[sec].get(k)
            if v is not None:
                sds_client_dict[k] = v
        
        for k in ["fileborder_samples", "fileborder_seconds"]:
            v = config[sec].getint(k)
            if v is not None:
                sds_client_dict[k] = v
        self.__setattr__('sds_client_dict', sds_client_dict)
        return config
        
        
    def print(self):
        """
        Print attributes and their values
        """
        for k, v in self.__dict__.items():
            print(k, ':', v)


    def process_stream(self, st, inv, starttime, endtime):
        """
        Basic processing of seismic stream and returns trace

        Applies

        - remove sensitivity to convert from counts to m/s
        - merge traces in stream, fill gaps with Nans
        - trims to starttime, endtime (should be ``proclen``)
          and extends with Nans if necessary 

        """
        return process_stream(st, inv, starttime, endtime)


    def process(self, startdate, enddate, outdir=None,
                ignore_years=False,
                    **processing_kwargs ):
        """
        Process data between given time range and write
        results to disk.

        We generally assume, that all other parameters for 
        the database search (network, station, channel, etc)
        and processing parameters have been set before during
        initialization. Nevertheless, they can be overridden
        here.

        Parameters
        -------------
        startdate : str, UTCDateTime
            begin of time range for analysis. Can be UTCDateTime
            object or string readable by UTCDateTime
        enddate : str, UTCDateTime
            end of time range, same as startdate
        outdir : str, path-like
            directory where to write results to
        ignore_years : bool
            if True, we don't start new output files at change of
            year. Only recommended for short periods crossing 
            change of year.
        processing_kwargs :
            see RawDataProcessor
        """
        T0 = time.time()
        self.startdate = UTC(startdate)
        self.enddate = UTC(enddate)

        self._check_set_attr(processing_kwargs)

        # else override outdir attrib if given here,
        # Could include it in processing_kwargs
        if outdir is not None:
            self.outdir = outdir

        # Init loop counters
        _startdate = startdate
        if ignore_years:
            _enddate = enddate
        else:
            _enddate = UTC("{:d}-12-31".format(_startdate.year))
        while _startdate < enddate:
            # Check if end-of-year date is larger than acual enddate
            if _enddate > enddate: 
                _enddate = enddate
                logger.debug("Reset enddate to %s" % _enddate)

            # Process 1 year (or less)
            logger.info("\nProcessing %s - %s" % (_startdate, _enddate))
            networks, stations, channels = self.get_nsc(_startdate, _enddate)

            for n in networks:
                for s in stations:
                    for c in channels:
                        data = self.process_nsc(_startdate, _enddate,
                                                    n, s, c)
                        
                        data.trim_nan()
                        data.to_file(self.outdir)

            # Proceed loop
            if ignore_years: break

            # If we work with dates here, we may avoid problems with leap years
            _startdate = UTC("{:d}-01-01".format(_startdate.year+1))
            _enddate = UTC("{:d}-12-31".format(_startdate.year))
            logger.debug("New _startdate = %s" % _startdate)
            logger.debug("New _enddate = %s" % _enddate)
        walltime = timedelta(seconds=time.time()-T0)
        logger.info("Finished. Took %s h" % walltime)


    def get_nsc(self, starttime, endtime):
        """
        Expand network, station, channel attributes (which can contain
        wildcards) to determine all unique networks, stations and channels 
        available within in time range.

        Parameters
        --------------
        starttime : UTCDatetime
            begin of time range, we only use the year
        endtime : UTCDateTime
            end of time range, we only use the year
        """

        networks = set()
        stations = set()
        channels = set()
        years = np.arange(starttime.year, endtime.year+1)
        for y in years:
            logger.debug("%s" % str(y))
            for n in glob(os.path.join(self.sdsclient.sds_root, 
                                        str(y), self.network)):
                nw = os.path.split(n)[-1]
                networks.add(nw)
                logger.debug("%s" % nw)
                for s in glob(os.path.join(n, self.station)):
                    sn = os.path.split(s)[-1]
                    stations.add(sn)
                    logger.debug("%s" % sn)
                    for c in glob(os.path.join(s, self.channel)+
                                    '.{}'.format(self.sdsclient.sds_type)):
                        ch = os.path.split(c)[-1].split('.')[0]
                        channels.add(ch)
                        logger.debug("%s" % ch)
                        
        return networks, stations, channels


    def process_nsc(self, startdate, enddate, 
                    network, station, channel,
                    **processing_kwargs ):
        """
        Compute noise amplitudes and spectra of raw data 
        for a specific network, station, channel combination 
        within time range.

        Parameters
        -----------
        startdate : UTCDateTime 
        enddate : UTCDateTime
        network : str
            network to process, here **no wildcards** should be
            used
        station : str
        channel : str


        Returns
        ---------
        output : BaseProcessedData
            an object, that contains resulting amplitudes, spectra,
            station code and relevant processing parameters


        Warning
        ---------
        The ``network``, ``station`` and ``channel`` are different
        from the ones in the attrbutes if these contain wildcards. 
        They are rather a specific combination of n,s,c available 
        in the database given the wildcards.
        """
        station_dict = {'network': network,
                        'station': station,
                        'location': '',
                        'channel': channel}
        self._check_set_attr(processing_kwargs)
        fmin, fmax = self.amplitude_frequencies
        output = BaseProcessedData(startdate, enddate, 
                        ".".join(list(station_dict.values())),
                        self.amplitude_frequencies, self.winlen_in_s)
        
        AMP = []
        PXX = []
        starttime = startdate-self.overlap
        inv = self.invclient.get_stations(
            starttime=startdate, endtime=enddate, level='response',
            **station_dict)
        logger.info("Processing %s" % output.stationcode)
        while starttime <= enddate - self.overlap:
            endtime = starttime + self.proclen + 2*self.overlap
            logger.debug("%s - %s" % (starttime, endtime))
            
            st = self.sdsclient.get_waveforms(starttime=starttime, endtime=endtime, 
                                    **station_dict)
            try:
                tr = self.process_stream(st, inv, starttime, endtime)
            # No data in trace:
            except IndexError:
                # logger.debug("No data for %s" % UTC((starttime + overlap).date))
                logger.info("No data for %s" % starttime)
                starttime = starttime + self.proclen
                # Reset starttime if there was no data so far.
                if len(AMP) == 0:
                    output.startdate = UTC((starttime + self.overlap).date)
                # Append Nan-data to cover gap
                else:
                    AMP.append(np.ones(AMP[-1].shape)*np.nan)
                    PXX.append(np.ones(PXX[-1].shape)*np.nan)
                continue

            # Get some numbers
            sr = tr.stats.sampling_rate
            nf = int(self.proclen/self.winlen_in_s)
            #proclen_samples = proclen * sr
            winlen_samples = int(self.winlen_in_s * sr)
            
            # Spectra
            data = get_adjacent_frames(tr, starttime+self.overlap, nf, 
                                    winlen_samples)
            freq, P = welch(data, fs=sr, nperseg=self.nperseg, axis=1)
                
            # Amplitude
            prctl = get_amplitude(tr, starttime+self.overlap, fmin, fmax,
                                self.overlap, winlen_samples, nf)
            
            AMP.append(prctl) #amp[1:-1])
            PXX.append(P) # pxx[1:-1,:])

            starttime = starttime + self.proclen
        
        if len(AMP) > 0:
            output.set_data(np.array(AMP),
                            np.array(PXX),
                            freq)
        return output


    def get_station_dict(self):
        station_dict = {'network': self.network,
                        'station': self.station,
                        'location': '',
                        'channel': self.channel}
        return station_dict


    def update(self, startdate=UTC(UTC.now().date)-2*24*3600, 
                    enddate=UTC(UTC.now().date)-24*3600,
                    ignore_years=False):
        """
        Update existing database with new data
        """
        T0 = time.time()
        # Find latest entry in raw data base
        # Find latest data in processed database
        logger.info("Updating...")
        networks, stations, channels = self.get_nsc(startdate, enddate)

        # Cases
        # time range > 1 year
        # enddate.year != startdate.year --> need to create new file

        # Init loop counters
        _startdate = startdate
        if ignore_years:
            _enddate = enddate
        else:
            _enddate = UTC("{:d}-12-31".format(_startdate.year))
        while _startdate < enddate:
            # Check if end-of-year date is larger than acual enddate
            if _enddate > enddate: 
                _enddate = enddate
                logger.debug("Reset enddate to %s" % _enddate)

            # Process 1 year (or less)
            logger.info("\nProcessing %s - %s" % (_startdate, _enddate))
            networks, stations, channels = self.get_nsc(_startdate, _enddate)

            for n in networks:
                for s in stations:
                    for c in channels:
                        # Check if there is data for this time range

                        data = self.process_nsc(_startdate, _enddate,
                                                    n, s, c)
                        
                        data.trim_nan()
                        data.to_file(self.outdir)

            # Proceed loop
            if ignore_years: break

            # If we work with dates here, we may avoid problems with leap years
            _startdate = UTC("{:d}-01-01".format(_startdate.year+1))
            _enddate = UTC("{:d}-12-31".format(_startdate.year))
            logger.debug("New _startdate = %s" % _startdate)
            logger.debug("New _enddate = %s" % _enddate)
        walltime = timedelta(seconds=time.time()-T0)
        logger.info("Finished. Took %s h" % walltime)



class BaseProcessedData():
    """
    Handle processed amplitudes and spectra together
    with some meta data. Read/write HDF5
    """
    def __init__(self,startdate=None, enddate=None,
                stationcode="....", 
                amplitude_frequencies=(None,None),
                seconds_per_window=None):
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
            logger.warning("No data for file %s. Skipping" % fname)
            return

        logger.info("Writing data to %s" % fname)

        with h5py.File(fname, 'w') as fout:
            fout.create_dataset('amplitudes', data=self.amplitudes)
            fout.create_dataset('psds', data=self.psds)
            fout.create_dataset('frequency_axis', data=self.frequency_axis)
            
            _create_hdf5_attribs(fout, self.stationcode, 
                                self.startdate, self.enddate,
                                self.amplitude_frequencies,
                                self.seconds_per_window
                                )
            # fout.attrs.create('stationcode', self.stationcode)
            # fout.attrs.create('starttime',
            #     np.array(self.startdate.timetuple()[:6]))
            # fout.attrs.create('endtime',
            #     np.array(self.enddate.timetuple()[:6]))
            # fout.attrs.create('amplitude_frequencies', 
            #                 np.array(self.amplitude_frequencies))
            # fout.attrs.create('seconds_per_window', 
            #                 self.seconds_per_window)

    
    def set_data(self, amplitudes, psds, freqs):
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
        # Get total number of days to get new array sizes
        tmin = min(self.startdate, new.startdate)
        tmax = max(self.enddate, new.enddate)
        days = timedelta(seconds=tmax-tmin).days+1
        
        # If shapes, processing parameters, frequency axis or 
        # stationcodesare inconsistent, we get an error here
        self._check_sanity(new)
        
        # Initialize data containers for merged data
        amps_shp, psds_shp = (list(self.amplitudes.shape), 
                                list(self.psds.shape))
        amps_shp[0] = days
        psds_shp[0] = days
        
        new_amps = np.ones(amps_shp)*np.nan
        new_psds = np.ones(psds_shp)*np.nan
        
        # Insert the data. We ignore overlaps here. New data
        # overwrites existing data if they overlap
        for d in [self, new]:
            i = timedelta(seconds=d.startdate-tmin).days
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
        """
        # Any of these methods raises IOError
        # if inconsistencies are found
        self._check_shapes(new)
        self._check_frequencies(new)
        self._check_codes(new)
        return True
    
        
    def _check_shapes(self, new):
        amps_shp = self.amplitudes.shape
        psds_shp = self.psds.shape
        if (amps_shp[-1] == new.amplitudes.shape[-1] and
            psds_shp[1:] == new.psds.shape[1:]):
            return amps_shp, psds_shp
        else:
            logger.error("Data in files have inconsistent shapes!")
            raise IOError("Data in files have inconsistent shapes!")
            #return False, False

    
    def _check_frequencies(self, new):
        if not np.all(np.isclose(self.frequency_axis, new.frequency_axis)):
            logger.error("Frequency axis are different!")
            raise IOError("Frequency axis are different!")
        else:
            return True

    
    def _check_codes(self, new):
        if self.stationcode != new.stationcode:
            logger.error("Station codes are different! " +
                "Cannot extend existing data.")
            raise IOError("Station codes are different! " +
                "Cannot extend existing data.")
        else:
            return True

    
    def _plot_amplitudes(self):
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
        
    
    def __repr__(self):
        nsc = "Results for {}\n".format( self.stationcode) + 20*"-"
        s = "Starttime: {}".format(self.startdate)
        e = "Enddate: {}".format(self.enddate)
        d = "Days: {:d}".format(self.amplitudes.shape[0])
        shp1 = "Amplitude shape = {}".format(self.amplitudes.shape)
        shp2 = "PSD shape = {}".format(self.psds.shape)
        pr1 = "Seconds per window = {:g}".format(self.seconds_per_window)
        pr2 = "Amplitude for {:g} - {:g} Hz".format(*self.amplitude_frequencies)

        return "\n".join([nsc, s,e,d,shp1,shp2, pr1, pr2])




def process_stream(st, inv, starttime, endtime):
    st.remove_sensitivity(inv)
    st.merge(fill_value=np.nan)
    st.trim(starttime, endtime, pad=True, fill_value=np.nan, 
            nearest_sample=False)
    if len(st) > 1:
        logger.warning("More than 1 trace in stream %s! Using first trace only" % 
                st)
    tr = st[0]
    
    # Demean ignoring gaps
    tr.data = tr.data - np.nanmean(tr.data)
    return tr




def get_data(dataclient, starttime, overlap, proclen,
             inv, network, station, location, channel,):
    
    starttime = starttime - overlap
    endtime = starttime + proclen + 2*overlap
    st = dataclient.get_waveforms(network, station, 
                             location, channel,starttime, endtime)
    st.remove_sensitivity(inv)
    st.merge(fill_value=np.nan)
    st.trim(starttime, endtime, pad=True, fill_value=np.nan, 
            nearest_sample=False)
    if len(st) > 1:
        raise RuntimeWarning("More than 1 trace in stream!")
    return st[0]    
    

def get_adjacent_frames(tr, starttime, nf, winlen_samples):
    """
    Reshape vector into (``nf``, ``winlen_samples``)
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
    
    Frames containing any Nans are set entirely to Nan.
    Loosely based on `obspy.signal.util.enframe`
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
        logger.info('Found nans in %s' % tr)
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


def _create_hdf5_attribs(f, stationcode, starttime, endtime, 
                        freqs, winlen_seconds):

    f.attrs.create('stationcode', np.string_(stationcode))
    f.attrs.create('starttime',  np.array(starttime.timetuple()[:6]))
    f.attrs.create('endtime', np.array(endtime.timetuple()[:6]))
    f.attrs.create('amplitude_frequencies', np.array(freqs))
    f.attrs.create('seconds_per_window', np.array(winlen_seconds))



