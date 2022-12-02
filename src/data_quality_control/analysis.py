#import configparser
from datetime import timedelta, time
from glob import glob
from pathlib import Path
import numpy as np

#from scipy.signal import welch, get_window

#from obspy.signal.filter import bandpass
#from obspy.clients.filesystem.sds import Client
#from obspy.clients.fdsn import RoutingClient
from obspy.core import UTCDateTime as UTC
# from obspy.signal import util

import plotly.graph_objects as go

import matplotlib.pyplot as plt
plt.style.use('tableau-colorblind10')

import h5py

from . import processing, base, util, dqclogging

import logging
# Create the global logger
logger = dqclogging.create_logger()
module_logger = logging.getLogger(logger.name+'.analysis')



class Analyzer():
    def __init__(self, 
                 datadir, stationcode, fileunit="year",
                #stime="00:00", etime="23:59:59:999999"
                ):
        #self.sdate = UTC(sdate).date
        #self.edate = UTC(edate).date
        #self._update_time(stime, etime)
        #self._update_datetime()
        self.datadir = datadir
        self.stationcode = stationcode
        self.fileunit  = fileunit
        self.iter_time = util.TIME_ITERATORS[self.fileunit]
        
        # Get fmtstr of data files
        fmtstr_base, sep, fmtstr_time = util.FNAME_FMTS[self.fileunit].rpartition("_")
        self.fmtstr = (fmtstr_base.format(
                        outdir=self.datadir, **self.nslc_as_dict()) + 
                        sep + fmtstr_time)
        self.logger = logging.getLogger(module_logger.name+
                            '.'+"Analyzer")
        self.logger.setLevel(logging.DEBUG)
        
        #self.files = self.get_filenames()
        
    def get_available_datafiles(self):
        """
        Return list with all available HDF5-filenames for 
        self.stationcode in self.datadir
        """
        #return glob(self.fmtstr.rpartition("_")[0] + "*.hdf5")
        return [str(f) for f in 
                Path(self.datadir).glob(self.stationcode+"*.hdf5")]


    def get_available_timerange(self):
        """
        Get list of available files and read startdate of first 
        and enddate of last file in sorted list.
        """
        flist = self.get_available_datafiles()
        flist.sort()
        startdate = base.BaseProcessedData().from_file(flist[0]).startdate
        enddate = base.BaseProcessedData().from_file(flist[-1]).enddate
        return startdate, enddate


    def _get_filenames(self):
        """
        Get filenames within set time range.

        `self.starttime` and `self.endtime` must be set!
        """
        
        logger.info("Looking for data file %s" % self.fmtstr)
        files = []
        
        for starttime, endtime in self.iter_time(self.starttime, self.endtime):
            files.append(self.fmtstr.format(year=starttime.year, 
                                        month=starttime.month, 
                                        day=starttime.day,
                                        hour=starttime.hour))
   
        return sorted(files)
        

    def nslc_as_dict(self):
        d = {k: v for k, v in zip(["network", "station", "location", "channel"], 
                                  self.stationcode.split("."))}
        return d
    
    
    def _tstr2time(self, t):
        return time(*[int(s) for s in t.split(':')])
        
            
    def _update_datetime(self):
        self.starttime = UTC("{}T{}".format(self.sdate, self.stime))
        self.endtime = UTC("{}T{}".format(self.edate, self.etime))


    def _update_time(self, stime, etime):
        if stime:
            self.stime = self._tstr2time(stime)
        if etime:
            self.etime = self._tstr2time(etime)
        self._update_datetime()
    

    def set_time(self, stime, etime):
        self.starttime = stime
        self.endtime = etime


    def get_all_available_data(self):
        startdate, enddate = self.get_available_timerange()
        self.get_data(startdate, enddate)


    def _iter_open_hdf5(self, fnamelist=None):
        """
        Generator that returns open h5py.File object for
        each filename in fnamelist. Closes file before
        yielding next file and before error is raised.

        If `fnamelist=None`, we use `self.files`. Causes
        error if not set.
        """
        if fnamelist is None:
            fnamelist = self.files

        for fname in fnamelist:
            logger.debug("Opening file %s" % fname)
            try:
                val = h5py.File(fname, 'r')
                # Return file object
                yield val
                # Close before proceding
                val.close()
            # Always close file before we 
            # present the error
            except Exception as e:
                val.close()
                logger.exception("Error while opening file %s\n" % fname + e)
                raise

            
    def get_data(self, starttimes, endtime=None):
        """
        Load amplitudes, psds and metadata from HDF5-files for
        specified times.

        `starttimes` can be given as array-like object containing
        UTCDateTimes or a single UTCDateTime. In the first case
        Psds are loaded only for the windows corresponding to 
        the times in `starttimes`. In the second case, `endtime`
        must be given as UTCDateTime. Psds are loaded for the 
        entire range between start and end time.

        If list of starttimes is given, starttime and endtime properties
        of the `Analyzer` are set to minimum and maximum times
        in list.

        Amplitudes are always loaded for the entire range between
        `starttime` and `endtime`.

        Parameters
        --------------
        starttimes : array-like or UTCDateTime
            If UTCDateTime, we expect also `endtime` as UTCDateTime.
            PSDs are loaded for the time range between start and endtime.
            if array-like, `endtime` is ignored. Psds are loaded for 
            windows listed in starttimes only.
        endtime : None or UTCDateTime
            only required if `starttimes` is UTCDateTime. End of time
            range for data selection.
        """
        if isinstance(starttimes, (list, np.ndarray, tuple)):
            stime = min(starttimes)
            etime = max(starttimes)
            self.timeax_psd = [np.datetime64(t) for t in starttimes]
        elif isinstance(starttimes, UTC) and isinstance(endtime, UTC):
            etime = endtime
            stime = starttimes
            starttimes = [stime, etime]
            self.timeax_psd = None
        else:
            raise UserWarning("Need to give either list of times or" + 
                            "start and endtime of time range. " + 
                            "Times must be obspy.UTCDateTimes.")

        self.set_time(stime, etime)
        self.files = self._get_filenames()

        DATA = base.BaseProcessedData()
        for fname in self.files:
            logger.debug("Loading %s" % fname)
            DATA.extend_from_file(fname)
        
        inds_amp, inds_psd = self._get_data_indices(DATA, starttimes)

        self.amps = DATA.amplitudes[inds_amp,:]
        self.psds = DATA.psds.reshape((-1, DATA.psds.shape[-1]))[inds_psd]
        self.freqax = DATA.frequency_axis
        self.proclen_seconds = DATA.proclen_seconds
        self.winlen_seconds = DATA.seconds_per_window
        self.nwin = self.amps.shape[1]
        
        if not self.timeax_psd:
            k, inc = util.choose_datetime_inc(self.winlen_seconds)
            self.timeax_psd = np.arange(stime, etime, inc, dtype="datetime64[{}]".format(k))

        return DATA


    def _get_data_indices(self, DATA, times):
        # Amplitude indices: we only select whole processing units (e.g. 1 day)
        # Get indices of data slices
        i = int((self.starttime - DATA.startdate) / 
                     DATA.proclen_seconds)
        j = int((self.endtime + DATA.proclen_seconds - DATA.startdate) / 
                         DATA.proclen_seconds)
        inds_amp = slice(i, j)

        # PSD indices
        inds = [int((t - DATA.startdate) / DATA.seconds_per_window )
                for t in times] 
        if len(inds) == 2:
            inds_psd = slice(*inds)
        else:
            # i, j = np.unravel_index(inds, DATA.amplitudes.shape)
            inds_psd = inds

        return inds_amp, inds_psd


    def infostr(self):
        t = (self.stationcode + "<br>" +
            "{} - {}<br>".format(self.sdate, self.edate) +
            "{} - {}".format(self.stime, self.etime))
        return t


    def plot_spectrogram(self, ax=None, **kwargs):
        if ax is None:
            fig, ax = plt.subplots(1,1)
        else:
            fig = ax.get_figure()
        
        pmesh = ax.pcolormesh(self.timeax_psd, self.freqax, 
                            np.log10(self.psds.T)*10, 
                            shading="auto",
                            **kwargs)
        fig.autofmt_xdate()
        plt.colorbar(pmesh, ax=ax, 
                    label=r'power spectral density, dB($\frac{m^2}{s^2\cdot Hz}$)')        
        ax.set_xlabel("time")
        ax.set_ylabel("frequency, Hz")
        return fig


    def plot3d(self):
        return self.plot3d_amplitudes(), self.plot3d_psds()


    def plot3d_amplitudes(self, func=None):

        title = ("Hourly 75%-amplitude<br>" + 
                    self.infostr()
            )

        if func:
            z = func(self.amps)
        else:
            z = self.amps
        dateax, timeax = self._get_time_axis()
        
        # Numpy-datetime can give you a **really** hard time to convert
        # between the different increments....
        xticks = [str(timedelta(
            **{np_td2datetime_td_keywords[str(timeax.dtype)] : int(np.int64(s))})) 
               for s in timeax]
        
        char = str(timeax.dtype)[-2]
        timeax = np.array(timeax, dtype=np.datetime64(None, char))
        #print(xticks, timeax)
        fig = self._plotly_3dsurface(timeax, dateax, z,
                        name="amplitudes")

        #xticks = [str(s).split("T")[-1] for s in timeax]
        
        fig.update_layout(title=title,
            scene=dict(
                xaxis=dict(title='Time', ticktext=xticks, tickvals=timeax),
                yaxis=dict(title='Date'),
                zaxis=dict(title="m/s")
            )
        )
        return fig

        

    def plot3d_psds(self, func=None):

        title = ("Hourly power spectral density\n" + 
                    self.infostr()
            )
        if func:
            z = func(self.psds)
        else:
            z = self.psds
        try:
            funcname = func.__name__+"(", ")"
        except AttributeError:
            funcname = "", ""
        nwin = z.shape[1]
        z = z.reshape((z.shape[0]*z.shape[1], z.shape[2]))
        dateax, timeax = self._get_time_axis()
        datetimeax = dateax[:,None] + timeax[None,:]
        y = datetimeax.ravel()
        x = self.freqax
        fig = self._plotly_3dsurface(x, y, z, name="psds")
        fig.update_layout(title=title, 
                        scene=dict(
                            xaxis=dict(title='Frequency, Hz'),
                            yaxis=dict(title='Datetime'),
                            zaxis=dict(title="psd, {}m^2/s^2/Hz{}".format(*funcname)
                                        )
                                )
                            ),
                        
        return fig


    def plot_psd(self):
        pass

    


    def _plotly_3dsurface(self,x,y, z, name=None, cmin=None, cmax=None):
        #sh_0, sh_1 = z.shape
        #y, x = np.linspace(0, sh_0-1, sh_0), np.linspace(0, sh_1-1, sh_1)
        fig = go.Figure(data=[go.Surface(z=z, x=x, y=y, name=name, 
                                            cmin=cmin, cmax=cmax)])
        fig.update_layout(autosize=True,
                          width=800, height=500,
                          scene=dict(aspectmode='manual',
                                     aspectratio=dict(x=1, y=2, z=0.5))
                          #margin=dict(l=65, r=50, b=65, t=90)
                         )
        #fig.show()
        return fig


    def _get_time_axis(self):
        
        dtflag, dtinc = util.choose_datetime_inc(self.proclen_seconds)

        dateax = np.arange(self.starttime, 
                            self.endtime+self.proclen_seconds,
                            dtinc,
                        dtype='datetime64[{}]'.format(dtflag))
        
        dtflag, dtinc = util.choose_datetime_inc(self.winlen_seconds)
        dur = self.nwin*self.winlen_seconds
        timeax = np.arange(0, int(dur/util.datetime_flags[dtflag])+dtinc, 
                           dtinc, 
                           dtype='datetime64[{}]'.format(dtflag))
                                                
        timeax = np.arange(0, self.nwin*dtinc, np.timedelta64(dtinc, dtflag))
        #print(timeax)
        return dateax, timeax

       

np_td2datetime_td_keywords = {'timedelta64[{}]'.format(v[0]) : v.lower() for 
                              v in ["minutes", "hours", "Days", "Months", "Years"]}



