#import configparser
from datetime import timedelta #, time
#from glob import glob
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

from . import base, util, dqclogging

import logging
# Create the global logger
logger = dqclogging.create_logger()
module_logger = logging.getLogger(logger.name+'.analysis')



class Analyzer():
    def __init__(self, 
                 datadir, nscl_code, fileunit="year",
                #stime="00:00", etime="23:59:59:999999"
                ):
        #self.sdate = UTC(sdate).date
        #self.edate = UTC(edate).date
        #self._update_time(stime, etime)
        #self._update_datetime()
        self.datadir = datadir
        self.stationcode = nscl_code
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
        self.logger.info("Looking for pattern " + str(Path(self.datadir).joinpath(
                        self.stationcode+"_"+util.FNAME_WILDCARD[self.fileunit]+".hdf5")))
        return [str(f) for f in 
                Path(self.datadir).glob(self.stationcode+"_"+util.FNAME_WILDCARD[self.fileunit]+".hdf5")]


    def get_available_timerange(self):
        """
        Get list of available files and read startdate of first 
        and enddate of last file in sorted list.
        """
        flist = self.get_available_datafiles()
        flist.sort()
        data = base.BaseProcessedData().from_file(flist[0])
        data.trim_nan()
        startdate = data.startdate

        data = base.BaseProcessedData().from_file(flist[-1])
        data.trim_nan()
        enddate = data.enddate
        return startdate, enddate


    def _get_filenames(self, starttime, endtime):
        """
        Get filenames within time range.
        """
        
        logger.info("Looking for data file %s" % self.fmtstr)
        files = []
        
        for _starttime, _endtime in self.iter_time(starttime, endtime):
            files.append(self.fmtstr.format(year=_starttime.year, 
                                        month=_starttime.month, 
                                        day=_starttime.day,
                                        hour=_starttime.hour))
   
        return sorted(files)
        

    def nslc_as_dict(self):
        d = {k: v for k, v in zip(["network", "station", "location", "channel"], 
                                  self.stationcode.split("."))}
        return d
    
    
    # def _tstr2time(self, t):
    #     return time(*[int(s) for s in t.split(':')])
        
            
    # def _update_datetime(self):
    #     self.starttime = UTC("{}T{}".format(self.sdate, self.stime))
    #     self.endtime = UTC("{}T{}".format(self.edate, self.etime))


    # def _update_time(self, stime, etime):
    #     if stime:
    #         self.stime = self._tstr2time(stime)
    #     if etime:
    #         self.etime = self._tstr2time(etime)
    #     self._update_datetime()
    

    def set_time(self, stime, etime):
        self.starttime = stime
        self.endtime = etime
        self.logger.info("New start and end time: {} and {}".format(
            self.starttime, self.endtime))


    def get_all_available_data(self):
        startdate, enddate = self.get_available_timerange()
        self.get_data(startdate, enddate)


    def _iter_open_hdf5(self, fnamelist=None):
        """
        Generator that returns open h5py.File object for
        each filename in fnamelist. Closes file before
        yielding next file and before error is raised.

        Useful for testing or if direct access to the hdf5-files
        is needed.

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

        
        if isinstance(starttimes, UTC):
            if not isinstance(endtime, UTC):
                endtime = UTC()
            starttimes = [starttimes, endtime]
            timerange = True
        elif isinstance(starttimes, (list, np.ndarray, tuple)):
            timerange = False
        else:
            raise UserWarning("Need to give either list of times or" + 
                            "start and endtime of time range. " + 
                            "Times must be obspy.UTCDateTimes.")

        self.logger.debug("timerange is set to {}".format(timerange))
        starttimes = sorted(starttimes)
        etime = starttimes[-1]
        stime = starttimes[0]
            
        self.files = self._get_filenames(stime, etime)

        DATA = base.BaseProcessedData()
        #self.logger.debug("{} - {}".format(str(DATA.startdate), str(DATA.enddate)))
        for fname in self.files:
            logger.debug("Loading %s" % fname)
            DATA.extend_from_file(fname)

        ## I don't know why but Nans are not always trimmed correctly, leading
        ## to false dates in DATA.
        DATA.trim_nan()
        self.logger.info("Available time range in data: {}-{}".format(
            DATA.startdate, DATA.enddate
        ))
        
        ## Reduce start/end time to those in DATA if out of available range
        if stime > DATA.enddate + DATA.proclen_seconds:
            msg = ("Requested time range starts after data is available." + 
                        "Data available from {} - {}").format(
                            *self.get_available_timerange())
            self.logger.exception(msg)
            raise RuntimeError(msg)

        if etime <= DATA.startdate:
            msg = ("Requested time range ends before data is available." + 
                            "Data available from {} - {}").format(
                            *self.get_available_timerange())
            self.logger.exception(msg)
            raise RuntimeError(msg)
        
        if  stime < DATA.startdate:      
            stime = DATA.startdate
            self.logger.info("Adjusting starttime to available: {}".format(stime))

        if DATA.enddate + DATA.proclen_seconds < etime:
            # if DATA.enddate + DATA.proclen_seconds <= stime:
            #     msg = ("Requested time range starts after available data." + 
            #                     "Data ends at {}").format(DATA.enddate)
            #     self.logger.exception(msg)
            #     raise RuntimeError(msg)
                
            # else:
            etime = DATA.enddate + DATA.proclen_seconds
            self.logger.info("Adjusting endtime to available: {}".format(etime))
        self.set_time(stime, etime)

        # print(len(starttimes))
        # print(timerange)

        inds_amp, inds_psd, timeax_psd = self._get_data_indices(DATA, starttimes, timerange)
        #self.logger.debug("Indices amplitude: {}".format(str(inds_amp)))
        #self.logger.debug("Indices PSD: {}".format(str(inds_psd)))
        
        ## It might make more sense to add the DATA directly
        self.amps = DATA.amplitudes[inds_amp,:]
        self.psds = DATA.psds.reshape((-1, DATA.psds.shape[-1]))[inds_psd]
        #self.logger.debug(self.psds.shape)
        self.freqax = DATA.frequency_axis
        self.proclen_seconds = DATA.proclen_seconds
        self.winlen_seconds = DATA.seconds_per_window
        self.nwin = self.amps.shape[1]
        self.timeax_psd = timeax_psd
        self.amplitude_frequency_range = DATA.amplitude_frequencies
        return DATA


    def _get_indices_timeax_timerange(self, DATA):
        times = [self.starttime, self.endtime]
        inds_amp, inds_psd = util._get_data_indices(DATA, times)
        k, inc = util.choose_datetime_inc(DATA.seconds_per_window)
        timeax_psd = np.arange(self.starttime, self.endtime, inc, 
                            dtype="datetime64[{}]".format(k))
        return inds_amp, slice(*inds_psd), timeax_psd

    def _get_indices_timeax_timelist(self, DATA, starttimes):
        self.logger.debug("len(input starttimes): {:d}".format(
            len(starttimes)))
        starttimes = [t for t in starttimes if 
                t >= self.starttime and t < self.endtime]
        #print(starttimes)
        self.logger.debug("len(starttimes) within available time: {:d}".format(
            len(starttimes)))
        timeax_psd = np.array([np.datetime64(t) for t in starttimes])
        self.logger.debug("Returning indices for timelist")
        inds_amp, inds_psd = util._get_data_indices(DATA, starttimes)
        self.logger.debug("Inds_psd: {:g}".format(len(inds_psd)))
        self.logger.debug("max ind {:g}".format(max(inds_psd)))
        return inds_amp, inds_psd, timeax_psd

    def _get_data_indices(self, DATA, starttimes, timerange):
        if timerange:
            return self._get_indices_timeax_timerange(DATA)
        else:
            return self._get_indices_timeax_timelist(DATA, starttimes)


    # def infostr(self):
    #     t = (self.stationcode + "<br>" +
    #         "{} - {}<br>".format(self.sdate, self.edate) +
    #         "{} - {}".format(self.stime, self.etime))
    #     return t


    def plot_spectrogram(self, ax=None, **kwargs):
        if ax is None:
            fig, ax = plt.subplots(1,1)
        else:
            fig = ax.get_figure()
        
        pmesh = ax.pcolormesh(self.timeax_psd, self.freqax, 
                            np.log10(self.psds.T*1e9**2), 
                            shading="auto",
                            **kwargs)
        fig.autofmt_xdate()
        plt.colorbar(pmesh, ax=ax, 
                    #label=r'power spectral density, dB($\frac{m^2}{s^2\cdot Hz}$)'
                    label=r'power spectral density, $\log_{10}(\frac{nm^2}{s^2\cdot Hz}$)')        
        ax.set_xlabel("time")
        ax.set_ylabel("frequency, Hz")
        return fig


    def plot3d(self):
        return self.plot3d_amplitudes(), self.plot3d_psds()


    def plot3d_amplitudes(self, func=None):

        
        if func:
            z = func(self.amps)
        else:
            z = self.amps
        dateax, timeax = self._get_time_axis()

        title = ("Hourly 75%-amplitude<br>" + 
            "{} - {}<br>".format(min(dateax), max(dateax)) +
            "{} - {}".format(min(timeax), max(timeax))
            )

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

        

    def plot3d_psds(self, func=None, zlabel=None):
        """
        
        Notes
        ----------
        Latex rendering for me works in title but not on axis labels.
        """

        if func:
            z = func(self.psds)
            if zlabel is None:
                try:
                    funcname = func.__name__+"(", ")"
                except AttributeError:
                    funcname = "", ""
                zlabel = "psd, {}m^2/s^2/Hz{}".format(*funcname)
        else:
            z = np.log10(self.psds*1e9**2)
            zlabel = r"$\alpha \log_{10}\frac{nm^2}{s^2Hz}$"
        
        #nwin = z.shape[1]
        #z = z.reshape((z.shape[0]*z.shape[1], z.shape[2]))
        #dateax, timeax = self._get_time_axis()
        #datetimeax = dateax[:,None] + timeax[None,:]
        #y = datetimeax.ravel()
        y = self.timeax_psd
        x = self.freqax
        fig = self._plotly_3dsurface(x, y, z, name="psds")

        title = ("Hourly power spectral density<br>" + 
           "{} - {}<br>".format(min(y), max(y))
           )
        #title = r'$\frac{\alpha^2}{\beta}$'

        fig.update_layout(title=title, 
                        scene=dict(
                            xaxis=dict(title='Frequency, Hz'),
                            yaxis=dict(title='Datetime'),
                            zaxis=dict(title=zlabel#"psd, {}m^2/s^2/Hz{}".format(*funcname)
                                        )
                                )
                            )
        
                        
        return fig


    # def plot_psd(self):
    #     pass

    


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


    def __repr__(self) -> str:
        s1 = "Analyzer for station {}".format(self.stationcode)
        datadir = "Datadir: {}".format(str(self.datadir))
        fileunit = "HDF5-file covers 1 {}".format(self.fileunit)
        fmtstr = "Filename pattern: {}".format(self.fmtstr)
        loglevel = "Loglevel: {}".format(self.logger.level)

        try:
            l1 = "I have data for {} - {}".format(self.starttime, self.endtime)
            shp1 = "Amplitude shape = {}".format(self.amps.shape)
            shp2 = "PSD shape = {}".format(self.psds.shape)
            pr1 = "Seconds per window = {:g}".format(self.winlen_seconds)
            pr2 = "Amplitude for {:g} - {:g} Hz".format(
                *self.amplitude_frequency_range)

            s2 = "\n".join([l1, shp1, shp2, pr1, pr2])
        except AttributeError as e:
            s2 = ""
            self.logger.debug(e)
        
        return "\n".join([s1, datadir, fileunit, fmtstr, loglevel, s2])


       

np_td2datetime_td_keywords = {'timedelta64[{}]'.format(v[0]) : v.lower() for 
                              v in ["minutes", "hours", "Days", "Months", "Years"]}



