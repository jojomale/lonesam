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

wildcards = ["?", "*"]


class Analyzer(base.BaseProcessedData):
    def __init__(self, 
                 datadir, nslc_code, fileunit="year",
                ):
        super().__init__(stationcode=nslc_code)
        self.datadir = datadir
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


    def nslc_as_dict(self):
        d = {k: v for k, v in zip(["network", "station", "location", "channel"], 
                                  self.stationcode.split("."))}
        return d

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
        self.logger.debug("Looking for earliest available time")
        data = base.BaseProcessedData().from_file(flist[0])
        data.trim_nan()
        startdate = data.startdate

        self.logger.debug("Looking for latest available time")
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


    def _check_times(self, starttimes, endtime):
        """
        Determines if given times indicate time range
        or time list.

        Raises UserWarning if input is wrong.
        """
        if isinstance(starttimes, UTC):
            if not isinstance(endtime, UTC):
                endtime = UTC()
            starttimes = [starttimes, endtime]
            self.timerange = True
        elif isinstance(starttimes, (list, np.ndarray, tuple)):
            self.timerange = False
        else:
            raise UserWarning("Need to give either list of times or" + 
                            "start and endtime of time range. " + 
                            "Times must be obspy.UTCDateTimes.")

        self.logger.debug("timerange is set to {}".format(self.timerange))
        return sorted(starttimes)


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
        
        if any([char in self.stationcode for char in wildcards]):
            raise RuntimeError("Station code {} ".format(self.stationcode) + 
                "contains wildcard characters. " + 
                "Can only get data for defined netw.stat.loc.chan!")
        
        self.logger.info("Renewing data")
        self._reset()
        starttimes = self._check_times(starttimes, endtime)
        etime = starttimes[-1]
        stime = starttimes[0]
            
        self.files = self._get_filenames(stime, etime)
        #print("Before loading", self.stationcode)
        for fname in self.files:
            logger.debug("Loading %s" % fname)
            self.extend_from_file(fname)
            #print("During load", self.stationcode)

        self.trim_nan()
        #self.fill_days()
        self._check_if_requested_times_are_available(stime, etime)
        
        self._check_shape_vs_time()
        self.timeax_psd = self._get_psd_datetimeax()

        if not self.timerange:
            self.filter_psds_for_times(starttimes)
        

    def trim(self, starttime=None, endtime=None, 
                fill_value=None):
        """
        Remove data outside of start/end time.
        """
        self.logger.debug("Trimming data if necessary")

        if starttime is None:
            starttime = self.startdate
        if endtime is None:
            endtime = self.enddate

        nbeg = int((starttime - self.startdate) / 
                self.seconds_per_window)
        if nbeg >= 0:
            new_amps = self.amplitudes[nbeg:]
            new_psds = self.psds[nbeg:,:]
        elif fill_value is not None:
            fill = np.ones((-nbeg, self.frequency_axis.size))*fill_value
            new_amps = np.hstack((fill[:,0], self.amplitudes))
            new_psds = np.vstack((fill, self.psds))
        else:
            new_amps = self.amplitudes
            new_psds = self.psds
            starttime = self.startdate
        
        nend = int((self.enddate - endtime) / 
                self.seconds_per_window)
        if nend >= 0:
            new_amps = new_amps[:-nend]
            new_psds = new_psds[:-nend,:]
        elif fill_value is not None:
            fill = np.ones((-nend, self.frequency_axis.size))*fill_value
            new_amps = np.hstack((new_amps, fill[:,0]))
            new_psds = np.vstack((new_psds, fill))
        else:
            endtime = self.enddate

        self.set_data(new_amps, new_psds, self.frequency_axis)
        self.set_time(starttime, endtime)
        self.timeax_psd = self._get_psd_datetimeax()
        self._check_shape_vs_time()
        self.logger.debug("New shapes: amps={}, psds={}".format(
                self.amplitudes.shape, self.psds.shape
        ))
        


    def _reset(self):
        for attr in ["amplitude_frequencies",
                    "seconds_per_window",
                    "startdate", "enddate", 
                    "amplitudes", "psds", "frequency_axis", 
                    ]:
            try:
                self.__setattr__(attr, None)
            except AttributeError:
                continue


    def _check_if_requested_times_are_available(self, stime, etime):
        ## Reduce start/end time to those in self if out of available range
        if stime > self.enddate:
            self._reset()
            msg = ("Requested time range starts after data is available." + 
                        "Data available from {} - {}").format(
                            *self.get_available_timerange())
            self.logger.exception(msg)
            raise RuntimeError(msg)

        if etime <= self.startdate:
            self._reset()
            msg = ("Requested time range ends before data is available." + 
                            "Data available from {} - {}").format(
                            *self.get_available_timerange())
            self.logger.exception(msg)
            raise RuntimeError(msg)
        
        if  stime < self.startdate:      
            stime = self.startdate
            self.logger.info("Adjusting starttime to available: {}".format(stime))

        if self.enddate < etime:
            etime = self.enddate
            self.logger.info("Adjusting endtime to available: {}".format(etime))
        


    def filter_psds_for_times(self, timelist):
        self.logger.debug("len(input timelist): {:d}".format(
            len(timelist)))
        timelist = [t for t in timelist if 
                t >= self.startdate and t < self.enddate]
        #print(starttimes)
        self.logger.debug("len(timelist) within available time: {:d}".format(
            len(timelist)))
        timeax = np.array([np.datetime64(t) for t in timelist])
        time_idx = [int((t - self.startdate)/self.seconds_per_window )
                    for t in timelist]
        self.timeax_psd = timeax
        self.psds = self.psds[time_idx,:]



    def __repr__(self) -> str:
        s1 = "Analyzer for station {}".format(self.stationcode)
        datadir = "Datadir: {}".format(str(self.datadir))
        fileunit = "HDF5-file covers 1 {}".format(self.fileunit)
        fmtstr = "Filename pattern: {}".format(self.fmtstr)
        loglevel = "Loglevel: {}".format(self.logger.level)

        try:
            l1 = "I have data for {} - {}".format(self.startdate, self.enddate)
            shp1 = "Amplitude shape = {}".format(self.amplitudes.shape)
            shp2 = "PSD shape = {}".format(self.psds.shape)
            pr1 = "Seconds per window = {:g}".format(self.seconds_per_window)
            pr2 = "Amplitude for {:g} - {:g} Hz".format(
                *self.amplitude_frequencies)

            s2 = "\n".join([l1, shp1, shp2, pr1, pr2])
        except AttributeError as e:
            s2 = "No data attached."
            self.logger.debug(e)
        
        return "\n".join([s1, datadir, fileunit, fmtstr, loglevel, s2])


    def plot_spectrogram(self, ax=None, func=None, 
            colorbarlabel="",
             **kwargs):
        """
        Plot power spectral densities as spectrogram.

        We use matplotlibs `pcolormesh` to plot the PSDs
        in a spectrogram-like way, i.e. time on x-axis, 
        frequency on y-axis and PSD as color.

        Parameters
        -----------
        ax : matplotlib.axes
            axes to plot in
        func : callable
            modifies PSDS as `func(self.psds.T)` before plotting.
            If not given, we plot `np.log10(self.psds.T * 1e9**2)`
            which is `log10(nm^2/s^2/Hz`.
        colorbarlabel : str [""]
            set colorbar label. Useful in combination with `func`
            to set correct units for color scale.
        kwargs :  
            keyword arguments passed to pcolormesh. 
            `vmax` can be callable to set `vmax=vmax(Z)`.
            If not given, we use: 
            `cmap=plt.cm.afmhot`,
            `shading=auto`,
            `vmax=np.nanmax(Z)`


        """

        if not "cmap" in kwargs:
            kwargs["cmap"] = plt.cm.afmhot
        if not "shading" in kwargs:
            kwargs["shading"] = "auto"
        
        
   
        if ax is None:
            fig, ax = plt.subplots(1,1)
        else:
            fig = ax.get_figure()

        if self.timerange:
            tax = self.timeax_psd
        else:
            tax = np.arange(self.timeax_psd.size)
            nticks = 10
            dtick = tax.size // nticks
            xticks = tax[::dtick]
            xticklabels = self.timeax_psd[::dtick]

        if func:
            Z = func(self.psds.T)
            colorbarlabel = colorbarlabel
        else:
            Z = np.log10(self.psds.T*1e9**2)
            colorbarlabel = r'power spectral density, $\log_{10}(\frac{nm^2}{s^2\cdot Hz})$'
        
        if not "vmax" in kwargs:
            kwargs["vmax"] = 0.9*np.nanmax(Z)
        elif callable(kwargs["vmax"]):
            kwargs["vmax"] = kwargs["vmax"](Z)
        else:
            kwargs["vmax"] = kwargs["vmax"]
        self.logger.debug("Kwargs passed to pcolormesh are {}".format(kwargs))
        
        pmesh = ax.pcolormesh(tax, self.frequency_axis, Z, 
                            **kwargs)

        if not self.timerange:
            self.logger.debug("Setting xticks")
            ax.set_xticks(xticks)
            ax.set_xticklabels(xticklabels)
        fig.autofmt_xdate()

        plt.colorbar(pmesh, ax=ax, 
                    label=colorbarlabel
                    )        
        ax.set_xlabel("time")
        ax.set_ylabel("frequency, Hz")
        return fig


    def plot3d(self):
        return self.plot3d_amplitudes(), self.plot3d_psds()


    def plot3d_amplitudes(self, func=None):
        z = self.reshape_amps_to_days()
        if func:
            z = func(z)

        dateax, timeax = self._get_date_and_time_axis_for_amplitude_matrix()

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
      
        y = self.timeax_psd
        x = self.frequency_axis
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




class Interpolator(Analyzer):
    def __init__(self, datadir, nslc_code, fileunit="year"):
        super().__init__(datadir, nslc_code, fileunit)
        
    def interpolate(self):
        TSTA, TEND = self.get_available_timerange()
        
        for tsta, tend in self.iter_time(TSTA, TEND):
            print(tsta, tend)
        # for f in files:
        #     data = base.BaseProcessedData().from_file(f)
        #     print(f)
        #     print(data)
        #     print()


       

np_td2datetime_td_keywords = {'timedelta64[{}]'.format(v[0]) : v.lower() for 
                              v in ["minutes", "hours", "Days", "Months", "Years"]}



