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
    """
    Handler to view and manage processed amplitudes and PSDs.

    Loads data for desired time range and provides plotting
    routines.
    
    Parameters
    ------------------ 
    datadir : :py:class:`pathlib.Path`
        directory of HDF5-files
    nslc_code : str
        code of format {network}.{station}.{location}.{channel}
    fileunit : {"year", "month", "day", "hour"}
        flag indicating the expected time range per file. 
        Determines the search pattern for files.

    Attributes
    ----------------
    iter_time(starttime, endtime) : func
        Generator providing start and endtime of HDF5-files
        within requested timerange
        depending on :py:attr:`.fileunit`
        
        See 
        
        - ``fileunit="year"``: :py:func:`data_quality_control.util.iter_years` ,
        - ``fileunit="month"``::py:func:`data_quality_control.util.iter_month` , 
        - ``fileunit="day"``::py:func:`data_quality_control.util.iter_days` ,
        - ``fileunit="hour"``::py:func:`data_quality_control.util.iter_hours`

    """

    def __init__(self, 
                 datadir, nslc_code, fileunit="year",
                ):
        super().__init__(stationcode=nslc_code)
        self.logger = logging.getLogger(module_logger.name+
                            '.'+"Analyzer")
        self.logger.setLevel(logging.DEBUG)

        self.datadir = datadir
        self.fileunit  = fileunit
        self.iter_time = util.TIME_ITERATORS[self.fileunit]
        

    @property        
    def fmtstr(self):
        """
        Name pattern of hdf5-files.
        """
        fmtstr_base, sep, fmtstr_time = util.FNAME_FMTS[
                            self.fileunit].rpartition("_")
        return (fmtstr_base.format(
                outdir=self.datadir, **self.nslc_as_dict()) + 
                sep + fmtstr_time)


    def nslc_as_dict(self):
        """
        Return station code as dictionary.
        """
        d = {k: v for k, v in zip(["network", "station", "location", "channel"], 
                                  self.stationcode.split("."))}
        return d


    def get_available_datafiles(self):
        """
        Return list with all available HDF5-filenames for 
        :py:attr:`stationcode` in :py:attr:`datadir`
        """
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


    def get_filenames(self, starttime, endtime):
        """
        Get filenames within time range.
        """
        
        self.logger.debug("Looking for data file %s" % self.fmtstr)
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


    @base.decorator_assert_integer_quotient_of_wins_per_fileunit
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
            
        self.files = self.get_filenames(stime, etime)
        self.logger.debug("Expecting files {}".format(self.files))
        for fname in self.files:
            self.extend_from_file(fname)

        # Does this really go here?
        self._check_if_requested_times_are_available(stime, etime)
        
        self._check_shape_vs_time()
        self.timeax_psd = self._get_psd_datetimeax()

        if not self.timerange:
            self.filter_psds_for_times(starttimes)
        

    def trim(self, starttime=None, endtime=None, 
                fill_value=None):
        """
        Remove data outside of start/end time.

        Parameters
        ----------------
        starttime : :py:class:`obspy.core.UTCDateTime`
            If `None` set starttime is used
        endtime : :py:class:`obspy.core.UTCDateTime`
            If `None` set endtime is used
        fill_value : float, int [None]
            If given, fills arrays to requested start/endtime 
            if available data starts after starttime or ends
            before endtime. If `None`, resulting time range
            may be shorter than requested.
        """
        self.logger.debug("Trimming data if necessary")

        if starttime is None:
            starttime = self.startdate
        if endtime is None:
            endtime = self.enddate

        nbeg = int((starttime - self.startdate) / 
                self.winlen_seconds)
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
                self.winlen_seconds)
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
        """
        Set attributes `amplitude_frequencies`, `winlen_seconds`,
        `startdate`, `enddate`, `amplitudes`, `psds`, 
        `frequency_axis` to `None`.
        """
        for attr in ["amplitude_frequencies",
                    "winlen_seconds",
                    "startdate", "enddate", 
                    "amplitudes", "psds", "frequency_axis", 
                    ]:
            try:
                self.__setattr__(attr, None)
            except AttributeError:
                continue


    def _check_if_requested_times_are_available(self, stime, etime):
        """
        Reduce start/end time to those in self if out of available range.
        Issue errors
        """
        ## 
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
        """
        Select and set PSDs only for given times.

        Finds indices of times in list, extracts 
        respective PSDs and replaces :py:attr:`psds`.
        
        Parameters
        ------------
        timelist : list of UTCDateTimes
        

        Warning
        -----------
        Overrides existing psd array! Cannot be used multiple
        times without reloading the data.

        """

        self.logger.debug("len(input timelist): {:d}".format(
            len(timelist)))
        timelist = [t for t in timelist if 
                t >= self.startdate and t < self.enddate]
        #print(starttimes)
        self.logger.debug("len(timelist) within available time: {:d}".format(
            len(timelist)))
        timeax = np.array([np.datetime64(t) for t in timelist])
        time_idx = [int((t - self.startdate)/self.winlen_seconds )
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
            pr1 = "Seconds per window = {:g}".format(self.winlen_seconds)
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
    """
    Manage interpolation/data reduction of processed amplitudes 
    and PSDs. 
    
    Grandparent is :py:class:`BaseProcessedData`.

    The main routine to call is :py:meth:`.interpolate`. 
    The other methods are mostly helpers to perform the
    smoothing/downsampling, organize the iteration
    over the data base and manage IO-operations. They can
    be useful on their own, though. 

    Parameters
    ------------------
    datadir : :py:class:`pathlib.Path`
        directory of HDF5-files
    nslc_code : str
        code of format {network}.{station}.{location}.{channel}
    fileunit : {"year", "month", "day", "hour"}
        flag indicating the expected time range per file. 
        Determines the search pattern for files.
    kernel_size : int > 1
        Number of samples over which median is computed.
        1 sample corresponds to the window size of the
        original data.
    kernel_shift : int >= 1
        Number of samples by which the kernel is shifted to
        compute the next median. See :py:attr:`.kernel_shift`


    Example
    -----------

    .. code-block:: python

        from data_quality_control import analysis, dqclogging

        nslc_code = "GR.BFO..BHZ"
        datadir = "output/"
        outdir = "output/interpolated/"
        kernel_size = 6
        kernel_shift = 3

        # Set loglevel in console, no logfile
        dqclogging.configure_handlers(loglevel_console="INFO",
                                    loglevel_file=None)

        # Interpolation
        ## Initiate interpolator
        polly = analysis.Interpolator(datadir, nslc_code, 
                                    kernel_size=kernel_size, 
                                    kernel_shift=kernel_shift)

        ## Start interpolation over whole available time range
        polly.interpolate(outdir, force_new_file=True)

        # View results
        lyza = analysis.Analyzer(datadir, nslc_code,
                                    fileunit="year")

        startdate = UTC("2020-12-24")
        enddate = UTC("2021-01-15")
        lyza.get_data(startdate, enddate)

        fig = lyza.plot_spectrogram()


    Attributes
    -----------------
    kernel_size : int > 1
        Number of samples over which median is computed.
        1 sample corresponds to the window size of the
        original data.
    kernel_shift : int >= 1
        Number of samples by which the kernel is shifted to
        compute the next median. ``kernel_shift=1`` corresponds
        to a classic running median, 
        ``kernel_shift>1`` downsamples the original data. 
        Determines 
        :py:attr:`data_quality_control.base.BaseProcessedData.winlen_seconds`. 
        For example if the seismic data were processed using 
        ``winlen_seconds=3600`` the new 
        :py:attr:`data_quality_control.base.BaseProcessedData.winlen_seconds`
        is 6*3600s = 21600s.
    WINLEN_SECONDS : int
        window size in seconds of the original processed data, i.e.
        the window size over which amplitudes and psds were computed
        from the seismic data. Is set based on the value in the
        first file and checked against 
        :py:attr:`data_quality_control.base.BaseProcessedData.winlen_seconds` 
        of every subsequent file. Raises :py:exc:`RuntimeError` if it changes.
    winlen_seconds : int
        Window size in seconds of the current data. Changes with every 
        iteration, i.e. every newly read file, between value 
        of the original processed data read from file and 
        ``kernel_shift*winlen_seconds`` after the median operation 
        is applied.
    """

    def __init__(self, datadir, nslc_code,
                fileunit="year",
                kernel_size=None,
                kernel_shift=None):
        super().__init__(datadir, nslc_code, fileunit)
        self.logger = logging.getLogger(module_logger.name+
                            '.'+"Interpolator")
        self.logger.setLevel(logging.DEBUG)
        self.kernel_size = None
        self.kernel_shift = None
        self.set_kernel(kernel_size, kernel_shift)
    

    def set_kernel(self, kernel_size, kernel_shift):
        """
        Set kernel parameters. Forces integers.

        Preserves existing attribute if 
        ``kernel_size=None`` or ``kernel_shift=None``.

        Parameters
        -------------
        kernel_size : num > 1
            Number of samples over which median is computed.
            1 sample corresponds to the window size of the
            original data.
        kernel_shift : num > 1
            Number of samples by which the kernel is shifted to
            compute the next median. See :py:attr:`.kernel_shift`
        """
        if kernel_size:
            kernel_size = int(kernel_size)
        elif self.kernel_size:
            kernel_size = self.kernel_size

        if kernel_shift:
            kernel_shift = int(kernel_shift) 
        elif self.kernel_shift:
            kernel_shift = self.kernel_shift

        self.kernel_size = kernel_size  
        self.kernel_shift = kernel_shift
        self.logger.info("Set kernel_size={}, kernel_shift={}".format(
            str(self.kernel_size), str(self.kernel_shift)
        ))
        # Add check assert_integer_quotient with self.WINLEN_SECONDS?
        # If so WINLEN_SECONDS would need to be a permanent attribute.


    def _get_WINLEN_SECONDS(self, TSTA, TEND):
        """
        Read first file in list to get window size in seconds.

        Parameters
        -------------
        TSTA : UTCDateTime
            beginning of time range in which we look for data.
        TEND : UTCDateTime
            end of time range in which we look for data.

        Example
        -----------
        If you want to know the ``winlen_seconds`` in your targeted
        data base, adjust the following code snippet:

        .. code-block:: python

            nslc_code = "GR.BFO..BHZ"
            datadir = "output/"
            polly = analysis.Interpolator(datadir, nslc_code)                             
            polly._get_WINLEN_SECONDS(*polly.get_available_timerange())
            print("Winlen in processed data", polly.WINLEN_SECONDS)
        
        """
        self.logger.debug("\n\nLooking for window size")
        for tsta, tend in self.iter_time(TSTA, TEND):
            self.get_data(tsta, tend)
            self.logger.debug("Time range to get window size: {} - {}".format(tsta, tend))
            self.logger.info("Expecting window size is {:g}s".format(self.winlen_seconds))
            self.WINLEN_SECONDS = self.winlen_seconds
            break


    def _set_check_WINLEN_SECONDS(self):
        """
        Sets :py:attr:`.WINLEN_SECONDS` if not set,
        otherwise raises error if different from 
        py:attr:`.winlen_seconds`. 
        Used to monitor if window size changes between files.
        """
        if not hasattr(self, "WINLEN_SECONDS"):
            self.logger.info("Expecting window size = {:g}s".format(
                self.winlen_seconds))
            self.WINLEN_SECONDS = self.winlen_seconds
        elif self.WINLEN_SECONDS != self.winlen_seconds:
            msg = "Window size changed"
            self.logger.error(msg)
            raise RuntimeError(msg)

            
    def _check_framed_shape(self, x, X, label=""):
        """
        Used in :py:meth:`._interpolate` to check if all data
        went into frames. 

        Parameters
        ------------
        x : ndarray (1D)
            data as time series
        X : ndarray (2D)
            data as (overlapping) frames
        kernel_shift : int
            distance between two frames in samples
        label : str
            name of x for more meaningful error msg.


        Raises
        ---------
        AssertionError
            if there are samples left.

        """
        nk, ks = X.shape
        ns = (nk-1)*self.kernel_shift+ks
        assert ns == x.size, \
            "{:d} of {} timeseries remain".format(x.size-ns, label)
    
    
    def _interpolate(self):
        """
        Transform timeseries data into frames and apply
        median operation.
        """
        self.logger.debug("Running self._interpolate()")
        x = self.amplitudes
        X = util.get_overlapping_frames(x, 
                                    self.kernel_size, self.kernel_shift)
        
        #print(x.size, X.shape)
        self._check_framed_shape(x, X, "amplitude")
        amplitudes_ = np.nanmedian(X, axis=1)

        x = self.psds[:,0]
        X = util.get_overlapping_frames(
            x, self.kernel_size, self.kernel_shift)
        #print(x.size, X.shape)
        self._check_framed_shape(x, X, "psd")
        PSD_ = np.array([np.nanmedian(
                util.get_overlapping_frames(
                    x, self.kernel_size, self.kernel_shift),axis=1) 
                         for x in self.psds.T]).T
        
        return amplitudes_, PSD_
    
    
    def iter_times_kernel(self, tsta, tend):
        """
        Generator providing start and endtime of data to read
        for interpolation.

        It's a wrapper around 
        :py:meth:`data_quality_control.analysis.Analyzer.iter_time` 
        which adjusts the yielded start and endtimes by 
        ``kernel_shift``, ``kernel_size`` and ``WINLEN_SECONDS`` to
        accommodate all data to interpolate over the entire ``fileunit``.


        Yields
        -------------
        new_tsta : :py:class:`UTCDateTime`
            starttime of data to read for interpolation. Determined
            by ``new_tend`` except for first one.
        new_tend : :py:class:`UTCDateTime`
            endtime of data to read for interpolation

        """
        new_tsta = None
        new_tend = None
        for _tsta, _tend in self.iter_time(tsta, tend):
            _tend = _tend + 24*3600
            if not new_tsta:
                new_tsta = _tsta

            new_tend = _tend + (self.kernel_size-self.kernel_shift)*self.WINLEN_SECONDS
            self.logger.debug("Times adjusted to kernel: {} - {}".format(
                new_tsta, new_tend))
            
            yield new_tsta, new_tend
            
            new_tsta = new_tend + (self.kernel_shift-self.kernel_size)*self.WINLEN_SECONDS 
            
      
    def _check_kernelshiftsize(self):
        """
        Check if :py:attr:`.winlen_seconds` and
        :py:attr:`kernel_shift` if resulting increment
        does not yield an integer quotient of possible
        total durations defined by :py:attr:`.fileunit`.

        Runs 
        :py:func:`data_quality_control.util.assert_integer_quotient_of_wins_per_fileunit`

        Raises
        ----------
        UserWarning
            if ``kernel_shift`` yields inappropriate increment
        """
        try:
            util.assert_integer_quotient_of_wins_per_fileunit(
                self.winlen_seconds*self.kernel_shift, self.fileunit
            )
        except UserWarning:
            self.logger.warning(
                "{:g} is bad choice for `kernel_shift`! ".format(self.kernel_shift))
            raise UserWarning(
                "{:g} is bad choice for `kernel_shift`! ".format(self.kernel_shift) +
                "Yields new `winlen_seconds` of {:g}s ".format(self.winlen_seconds) + 
                "which does not yield integer quotient when dividing total "+
                "duration in file, i.e. effectively 1 hour or 24 hours.")


    def _get_start_endtime(self, starttime=None, endtime=None):
        """
        Replace ``None`` input parameters with 
        start/end time from all available data.
        """
        if starttime is None or endtime is None:
            self.logger.info("Looking up available timerange "+
                "because start or endtime is None.\n")
            TSTA, TEND = self.get_available_timerange()
            TSTA = UTC(TSTA.date)
            TEND = UTC(TEND.date)
            if starttime is None:
                starttime = TSTA
            if endtime is None:
                endtime = TEND
        return starttime, endtime


    def interpolate(self, outdir=".",
            starttime=None, endtime=None,
            kernel_size=None, kernel_shift=None, 
            force_new_file=False):
        """
        Performs smoothing/downsampling of database for 
        requested time range.

        Parameters
        ---------------
        outdir : str
            directory where new results are stored. Should
            be different from :py:attr:`.datadir` because
            file names will be identical and would override
            input data.
        starttime : UTCDateTime
            Begin of time range to process. If ``None`` earliest
            available time in database is used.
        endtime : UTCDateTime
            End of time range to process. If ``None`` lastest
            available time in database is used.
        kernel_size : int, None
            See :py:attr:`.kernel_size`. Overrides existing value.
            Must be set here if not set on initialization. 
        kernel_shift : int, None
            See :py:attr:`.kernel_shift`. Overrides existing value.
            Must be set here if not set on initialization.
        force_new_file : bool [False] 
            If ``True`` overwrites existing output data.
        """
        
        self.logger.info("\n\nStarting interpolate()")

        if Path(outdir) == Path(self.datadir):
            msg = ("datadir is same as outdir. " +
                    "May override input data base!")
            self.logger.warn(msg)
            raise UserWarning(msg)

        starttime, endtime = self._get_start_endtime(starttime, endtime)
        self._get_WINLEN_SECONDS(starttime, endtime)
        self.set_kernel(kernel_size, kernel_shift)
        self._check_kernelshiftsize()

        ofilemanager = base.ProcessedDataFileManager(outdir, 
                                fileunit=self.fileunit)

        self.logger.info("\n\nIterating files over time range:\n")
        for tsta, tend in self.iter_times_kernel(starttime, endtime):
            
            #self.logger.debug("Yielded {} - {}".format(tsta, tend))
            #tsta = tsta
            #tend = tend + 24*3600# + (kernel_size-kernel_shift)*self.WINLEN_SECONDS 
            self.logger.info("Interpolating {:} - {}".format(tsta, tend))
            #self.logger.debug("Getting data...")
            
            self.get_data(tsta, tend)
            self.trim(tsta, tend, fill_value=np.nan)
            self._set_check_WINLEN_SECONDS()
            
            amplitudes_, psds_ = self._interpolate()
            self.set_data(amplitudes_, psds_, self.frequency_axis)
            #print(tsta, tend, tend+(kernel_shift-kernel_size)*self.WINLEN_SECONDS)
            self.logger.info("Setting time for output:")
            self.set_time(
                tsta, tend+(self.kernel_shift-self.kernel_size)*self.WINLEN_SECONDS)
            self.winlen_seconds = self.kernel_shift*self.WINLEN_SECONDS
            self._check_shape_vs_time()
            
            ofilemanager.set_data(self)
            ofilemanager.write_data(force_new_file, force_fileunit=False)
            # self.to_file(outdir)

            self.logger.debug("\n")

        self.logger.info("Interpolation finished.")

       

np_td2datetime_td_keywords = {'timedelta64[{}]'.format(v[0]) : v.lower() for 
                              v in ["minutes", "hours", "Days", "Months", "Years"]}



