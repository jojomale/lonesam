#import configparser
from datetime import timedelta, time
from glob import glob
import os.path
import numpy as np

#from scipy.signal import welch, get_window

#from obspy.signal.filter import bandpass
#from obspy.clients.filesystem.sds import Client
#from obspy.clients.fdsn import RoutingClient
from obspy.core import UTCDateTime as UTC
# from obspy.signal import util

import plotly.graph_objects as go

import h5py

from . import processing, base, util

import logging
logger = logging.getLogger('analysis')
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)  # set level
cformatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                            datefmt='%y-%m-%d %H:%M:%S')
ch.setFormatter(cformatter)
if not logger.hasHandlers():
    logger.addHandler(ch)


class Analyzer():
    def __init__(self, sdate, edate, 
                 datadir, stationcode, fileunit="year",
                stime="00:00", etime="23:59:59:999999"):
        self.sdate = UTC(sdate).date
        self.edate = UTC(edate).date
        self._update_time(stime, etime)
        self._update_datetime()
        self.datadir = datadir
        self.stationcode = stationcode
        self.fileunit  = fileunit
        self.iter_time = util.TIME_ITERATORS[self.fileunit]
        
        # Get fmtstr of data files
        fmtstr_base, sep, fmtstr_time = util.FNAME_FMTS[self.fileunit].rpartition("_")
        self.fmtstr = (fmtstr_base.format(
                        outdir=self.datadir, **self.nslc_as_dict()) + 
                        sep + fmtstr_time)
        
        
        self.files = self.get_filenames()
        
        
        #self.get_data()
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
    

    def get_all_data(self, sdate=None, edate=None, 
                 datadir=None, stationcode=None):
        if sdate is not None:
            self.sdate = UTC(sdate)
        if edate is not None:
            self.edate = UTC(sdate)
        if datadir is not None:
            self.datadir = datadir
        if stationcode is not None:
            self.stationcode = stationcode
        self._update_datetime()
            
        files = sorted(self.get_filenames())
        if len(files) == 0:
            logger.warn("No files for %s in %s between %s and %s" %
                        (self.stationcode, self.datadir, 
                        self.sdate, self.edate))
            return
        
        # If we found files, a
        data = processing.BaseProcessedData()
        for file in files:
            data.extend_from_file(file)
        self.data = data
            
        
            
    def get_filenames(self):
        
        
        logger.info("Looking for data file %s" % self.fmtstr)
        files = []
        
        
        for starttime, endtime in self.iter_time(self.starttime, self.endtime):
            files.append(self.fmtstr.format(year=starttime.year, 
                                        month=starttime.month, 
                                        day=starttime.day,
                                        hour=starttime.hour))
   
        return sorted(files)



    # def select_longest(self, fnames):
    #     logger.debug("Found %s files for year." % 
    #                  str(len(fnames)))
    #     f, ext = os.path.splitext(fnames[0])
    #     print(f.split('_')[-1])
    #     edate = UTC(f.split('_')[-1])
    #     for _f in fnames[1:]:
    #         _f = os.path.split(
    #                 os.path.splitext(
    #                     _f)[0])[-1]
    #         _edate = UTC(_f.split('_')[-1])
    #         if _edate >= edate:
    #             edate = _edate
    #         if edate >= self.endtime:
    #             break
    #     print(f+ext)
    #     return f+ext
    
    
    def iter_files(self):
        """
        Generator that returns open h5py.File object for
        each filename in self.files.
        """
        for fname in self.files:
            logger.debug("Opening file %s" % fname)
            try:
                val = h5py.File(fname, 'r')
                # Return file object
                yield val
                # Close before proceding
                val.close()
            # Always close file before we 
            # present the error
            except:
                val.close()
                logger.error("Error while opening file %s" % fname)
                raise
                
            
    def get_data(self, 
                 stime=None, etime=None):

        """
        Currently, stime, etime doesn't seem to do anything
        usefull. In particular, it does not select a
        subset of hours!!!
        """
        
        DATA = base.BaseProcessedData()
        for fname in self.files:
            logger.debug("Loading %s" % fname)
            DATA.extend_from_file(fname)
        
        # Cut out desired time range
        self._update_time(stime, etime)
        
        
        i = int((self.starttime - DATA.startdate) / 
                    DATA.proclen_seconds)
        j = int((self.endtime + DATA.proclen_seconds - DATA.startdate) / 
                        DATA.proclen_seconds)
        
        #print(i, j)
        self.amps = DATA.amplitudes[i:j,:]

        self.psds = DATA.psds[i:j,:]
        self.freqax = DATA.frequency_axis
        self.proclen_seconds = DATA.proclen_seconds
        self.winlen_seconds = DATA.seconds_per_window
        self.nwin = self.amps.shape[1]
        return DATA
                

    def infostr(self):
        t = (self.stationcode + "<br>" +
            "{} - {}<br>".format(self.sdate, self.edate) +
            "{} - {}".format(self.stime, self.etime))
        return t


    def plot(self):
        return self.plot_amplitudes(), self.plot_psds()


    def plot_amplitudes(self, func=None):

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

        

    def plot_psds(self, func=None):

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
        # sdate = np.datetime64(self.sdate, 'h')
        # edate = np.datetime64(self.edate, 'h') + np.timedelta64(1, "D")
        # dateax = np.arange(sdate, edate, np.timedelta64(1, "D"), 
        #             dtype='datetime64')

        # dur = self.etime.hour - self.stime.hour
        # if dur <= 0:
        #     dur = dur + 24

        # timeax = np.arange(dur+1, dtype=np.timedelta64) + self.stime.hour
        
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







class Analyzer_old():
    def __init__(self, sdate, edate, 
                 datadir, stationcode, 
                stime="00:00", etime="23:59:59:999999"):
        self.sdate = UTC(sdate).date
        self.edate = UTC(edate).date
        self.stime = time(*[int(s) for s in stime.split(':')])
        self.etime = time(*[int(s) for s in etime.split(':')])
        self.datadir = datadir
        self.stationcode = stationcode
        self._update_datetime()
        self.files = self.get_filenames()
        
        #self.get_data()
        

    def _tstr2t(self, t):
        return time(*[int(s) for s in t.split(':')])
        
            
    def _update_datetime(self):
        self.starttime = UTC("{}T{}".format(self.sdate, self.stime))
        self.endtime = UTC("{}T{}".format(self.edate, self.etime))


    def _update_time(self, stime, etime):
        if stime:
            self.stime = self._tstr2t(stime)
        if etime:
            self.etime = self._tstr2t(etime)
        self._update_datetime()
    

    def get_all_data(self, sdate=None, edate=None, 
                 datadir=None, stationcode=None):
        if sdate is not None:
            self.sdate = UTC(sdate)
        if edate is not None:
            self.edate = UTC(sdate)
        if datadir is not None:
            self.datadir = datadir
        if stationcode is not None:
            self.stationcode = stationcode
        self._update_datetime()
            
        files = sorted(self.get_filenames())
        if len(files) == 0:
            logger.warn("No files for %s in %s between %s and %s" %
                        (self.stationcode, self.datadir, 
                        self.sdate, self.edate))
            return
        
        # If we found files, a
        data = processing.BaseProcessedData()
        for file in files:
            data.extend_from_file(file)
        self.data = data
            
        
            
    def get_filenames(self):
        
        filehead = os.path.join(self.datadir, self.stationcode)
        fmtstr = filehead + "_{:04d}*.hdf5"
        logger.info("Looking for data file %s" % fmtstr)
        _year = self.sdate.year
        files = []
        while _year <= self.edate.year:
            searchstr = fmtstr.format(_year)
            fnames = glob(searchstr)

            if len(fnames) > 1:
                files.append(self.select_longest(fnames))
            elif len(fnames) == 0:
                _year = _year +1
                continue
            else:
                files.append(fnames[0])

            # Get end year of latest file
            ## Remove file-ext and path
            f, ext = os.path.splitext(files[-1])
            _endtime = UTC(f.split('_')[-1])
            if _endtime.year >= self.endtime.year:
                break
            _year = _year + 1
        return files
        


    def select_longest(self, fnames):
        logger.debug("Found %s files for year." % 
                     str(len(fnames)))
        f, ext = os.path.splitext(fnames[0])
        print(f.split('_')[-1])
        edate = UTC(f.split('_')[-1])
        for _f in fnames[1:]:
            _f = os.path.split(
                    os.path.splitext(
                        _f)[0])[-1]
            _edate = UTC(_f.split('_')[-1])
            if _edate >= edate:
                edate = _edate
            if edate >= self.endtime:
                break
        print(f+ext)
        return f+ext
    
    
    def iter_files(self):
        """
        Generator that returns open h5py.File object for
        each filename in self.files.
        """
        for fname in self.files:
            logger.debug("Opening file %s" % fname)
            try:
                val = h5py.File(fname, 'r')
                # Return file object
                yield val
                # Close before proceding
                val.close()
            # Always close file before we 
            # present the error
            except:
                val.close()
                logger.error("Error while opening file %s" % fname)
                raise
                
            
    def get_data(self, datasets=['amplitudes', 'psds'], 
                 stime=None, etime=None):
        
        DATA = {k : [] for k in datasets}
        for f in self.iter_files():
            
            t0 = UTC(*f.attrs["starttime"])
            t1 = UTC(*f.attrs["endtime"])
            # We need to add 1 day because t1 is midnight.
            ndays = timedelta(seconds=t1-t0).days+1
            logger.debug("Days in file: %s" % (ndays))

            # Get number of days this year
            year = UTC(*f.attrs["starttime"]).year
            _t0 = UTC(year, 1, 1, 0, 0, 0)
            _t1 = UTC(year+1, 1, 1, 0, 0, 0)
            dt = timedelta(seconds=_t1-_t0)
            days_this_year = dt.days
            
            
            freqax = f["frequency_axis"][:]

            for dataset in datasets:
                d = f[dataset]
                shp = d.shape
                logger.debug("Shape of %s is %s" % (dataset, shp))

                # If year is incomplete, we fill missing days with nans
                if ndays < 365:
                    tmp = np.ones((days_this_year, *shp[1:]))*np.nan
                    i = t0.julday-1  # Juldays start at 1, -1 to get indices
                    j = i + shp[0]
                    tmp[i:j,:] = d[:]
                # Otherwise we can use data directly
                else:
                    tmp = d[:]
                logger.debug("Target shape of %s is %s " % (dataset, tmp.shape))
                
                
                # Extract part of the year that we need from this file
                if self.starttime > _t0:
                    i = self.starttime.julday-1
                else:
                    i = 0
                if self.endtime < _t1:
                    j = self.endtime.julday
                else:
                    j = None
                logger.debug("Day index %s and %s, Shape = %s" % (i, j, tmp[i:j].shape))
                
                # If user wants only some hours, this is temporarily
                # strains memory more than necessary. However if time
                # range crosses midnight concatenation becomes quite
                # complex. For now, we keep it simple.
                DATA[dataset].append(tmp[i:j])
        
        # Cut out desired time range
        self._update_time(stime, etime)
        h0 = self.stime.hour
        h1 = self.etime.hour+1
        if h1 == 0 or h1>23:
            h1 = 24
        

        for dataset, data in DATA.items():
            data = np.vstack(data)
            logger.debug("Final shape of %s is %s" % (dataset, data.shape))
            # If starttime < endtime indices correspond to 
            # hours and 1 time frame lies within 1 row 
            if h0 < h1:
                DATA[dataset] = data[:,h0:h1]
            # if timerange crosses midnight, data of one
            # frame lies in 2 consecutive rows.
            else:
                DATA[dataset] = np.hstack([data[0:-1,h0:],
                                data[1:, :h1]])
        if "amplitudes" in DATA:
            self.amps = DATA["amplitudes"]
        if "psds" in DATA:
            self.psds = DATA["psds"]
            self.freqax = freqax
            DATA["frequency_axis"] = freqax
        return DATA
                
    def infostr(self):
        t = (self.stationcode + "<br>" +
            "{} - {}<br>".format(self.sdate, self.edate) +
            "{} - {}".format(self.stime, self.etime))
        return t


    def plot(self):
        return self.plot_amplitudes(), self.plot_psds()


    def plot_amplitudes(self, func=None):

        title = ("Hourly 75%-amplitude<br>" + 
                    self.infostr()
            )

        if func:
            z = func(self.amps)
        else:
            z = self.amps
        dateax, timeax = self._get_time_axis()
        fig = self._plotly_3dsurface(timeax, dateax, z,
                        name="amplitudes")

        xticks = ["{:02d}:00".format(int(s) % 24) for s in timeax]
        fig.update_layout(title=title,
            scene=dict(
                xaxis=dict(title='Time', ticktext=xticks, tickvals=timeax),
                yaxis=dict(title='Date'),
                zaxis=dict(title="m/s")
            )
        )
        return fig

        

    def plot_psds(self, func=None):

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
        sdate = np.datetime64(self.sdate, 'h')
        edate = np.datetime64(self.edate, 'h') + np.timedelta64(1, "D")
        dateax = np.arange(sdate, edate, np.timedelta64(1, "D"), 
                    dtype='datetime64')

        dur = self.etime.hour - self.stime.hour
        if dur <= 0:
            dur = dur + 24

        timeax = np.arange(dur+1, dtype=np.timedelta64) + self.stime.hour
        return dateax, timeax
        






class Analyzer_old_adapted():
    def __init__(self, sdate, edate, 
                 datadir, stationcode, 
                stime="00:00", etime="23:59:59:999999"):
        self.sdate = UTC(sdate).date
        self.edate = UTC(edate).date
        self.stime = time(*[int(s) for s in stime.split(':')])
        self.etime = time(*[int(s) for s in etime.split(':')])
        self.datadir = datadir
        self.stationcode = stationcode
        self._update_datetime()
        self.files = self.get_filenames()
        
        #self.get_data()
        

    def _tstr2t(self, t):
        return time(*[int(s) for s in t.split(':')])
        
            
    def _update_datetime(self):
        self.starttime = UTC("{}T{}".format(self.sdate, self.stime))
        self.endtime = UTC("{}T{}".format(self.edate, self.etime))


    def _update_time(self, stime, etime):
        if stime:
            self.stime = self._tstr2t(stime)
        if etime:
            self.etime = self._tstr2t(etime)
        self._update_datetime()
    

    def get_all_data(self, sdate=None, edate=None, 
                 datadir=None, stationcode=None):
        if sdate is not None:
            self.sdate = UTC(sdate)
        if edate is not None:
            self.edate = UTC(sdate)
        if datadir is not None:
            self.datadir = datadir
        if stationcode is not None:
            self.stationcode = stationcode
        self._update_datetime()
            
        files = sorted(self.get_filenames())
        if len(files) == 0:
            logger.warn("No files for %s in %s between %s and %s" %
                        (self.stationcode, self.datadir, 
                        self.sdate, self.edate))
            return
        
        # If we found files, a
        data = processing.BaseProcessedData()
        for file in files:
            data.extend_from_file(file)
        self.data = data
            
        
            
    def get_filenames(self):
        
        filehead = os.path.join(self.datadir, self.stationcode)
        fmtstr = filehead + "_{:04d}.hdf5"
        logger.info("Looking for data file %s" % fmtstr)
        _year = self.sdate.year
        files = []
        while _year <= self.edate.year:
            searchstr = fmtstr.format(_year)
            fnames = sorted(glob(searchstr))
            print(fnames)
            if len(fnames) >= 1:
                files.append(fnames[0])
            elif len(fnames) == 0:
                _year = _year +1
                continue
            
            # Get end year of latest file
            ## Remove file-ext and path
            #f, ext = os.path.splitext(files[-1])
            _year = _year + 1
        return files
        
    
    
    def iter_files(self):
        """
        Generator that returns open h5py.File object for
        each filename in self.files.
        """
        for fname in self.files:
            logger.debug("Opening file %s" % fname)
            try:
                val = h5py.File(fname, 'r')
                # Return file object
                yield val
                # Close before proceding
                val.close()
            # Always close file before we 
            # present the error
            except:
                val.close()
                logger.error("Error while opening file %s" % fname)
                raise
                
            
    def get_data(self, datasets=['amplitudes', 'psds'], 
                 stime=None, etime=None):
        
        DATA = {k : [] for k in datasets}
        for f in self.iter_files():
            
            t0 = UTC(*f.attrs["starttime"])
            t1 = UTC(*f.attrs["endtime"])
            # We need to add 1 day because t1 is midnight.
            ndays = timedelta(seconds=t1-t0).days+1
            logger.debug("Days in file: %s" % (ndays))

            # Get number of days this year
            year = UTC(*f.attrs["starttime"]).year
            _t0 = UTC(year, 1, 1, 0, 0, 0)
            _t1 = UTC(year+1, 1, 1, 0, 0, 0)
            dt = timedelta(seconds=_t1-_t0)
            days_this_year = dt.days
            
            
            freqax = f["frequency_axis"][:]

            for dataset in datasets:
                d = f[dataset]
                shp = d.shape
                logger.debug("Shape of %s is %s" % (dataset, shp))

                # If year is incomplete, we fill missing days with nans
                if ndays < 365:
                    tmp = np.ones((days_this_year, *shp[1:]))*np.nan
                    i = t0.julday-1  # Juldays start at 1, -1 to get indices
                    j = i + shp[0]
                    tmp[i:j,:] = d[:]
                # Otherwise we can use data directly
                else:
                    tmp = d[:]
                logger.debug("Target shape of %s is %s " % (dataset, tmp.shape))
                
                
                # Extract part of the year that we need from this file
                if self.starttime > _t0:
                    i = self.starttime.julday-1
                else:
                    i = 0
                if self.endtime < _t1:
                    j = self.endtime.julday
                else:
                    j = None
                logger.debug("Day index %s and %s, Shape = %s" % (i, j, tmp[i:j].shape))
                
                # If user wants only some hours, this is temporarily
                # strains memory more than necessary. However if time
                # range crosses midnight concatenation becomes quite
                # complex. For now, we keep it simple.
                DATA[dataset].append(tmp[i:j])
        
        # Cut out desired time range
        self._update_time(stime, etime)
        h0 = self.stime.hour
        h1 = self.etime.hour+1
        if h1 == 0 or h1>23:
            h1 = 24
        

        for dataset, data in DATA.items():
            data = np.vstack(data)
            logger.debug("Final shape of %s is %s" % (dataset, data.shape))
            # If starttime < endtime indices correspond to 
            # hours and 1 time frame lies within 1 row 
            if h0 < h1:
                DATA[dataset] = data[:,h0:h1]
            # if timerange crosses midnight, data of one
            # frame lies in 2 consecutive rows.
            else:
                DATA[dataset] = np.hstack([data[0:-1,h0:],
                                data[1:, :h1]])
        if "amplitudes" in DATA:
            self.amps = DATA["amplitudes"]
        if "psds" in DATA:
            self.psds = DATA["psds"]
            self.freqax = freqax
            DATA["frequency_axis"] = freqax
        return DATA
                
    def infostr(self):
        t = (self.stationcode + "<br>" +
            "{} - {}<br>".format(self.sdate, self.edate) +
            "{} - {}".format(self.stime, self.etime))
        return t


    def plot(self):
        return self.plot_amplitudes(), self.plot_psds()


    def plot_amplitudes(self, func=None):

        title = ("Hourly 75%-amplitude<br>" + 
                    self.infostr()
            )

        if func:
            z = func(self.amps)
        else:
            z = self.amps
        dateax, timeax = self._get_time_axis()
        fig = self._plotly_3dsurface(timeax, dateax, z,
                        name="amplitudes")

        xticks = ["{:02d}:00".format(int(s) % 24) for s in timeax]
        fig.update_layout(title=title,
            scene=dict(
                xaxis=dict(title='Time', ticktext=xticks, tickvals=timeax),
                yaxis=dict(title='Date'),
                zaxis=dict(title="m/s")
            )
        )
        return fig

        

    def plot_psds(self, func=None):

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
        sdate = np.datetime64(self.sdate, 'h')
        edate = np.datetime64(self.edate, 'h') + np.timedelta64(1, "D")
        dateax = np.arange(sdate, edate, np.timedelta64(1, "D"), 
                    dtype='datetime64')

        dur = self.etime.hour - self.stime.hour
        if dur <= 0:
            dur = dur + 24

        timeax = np.arange(dur+1, dtype=np.timedelta64) + self.stime.hour
        return dateax, timeax
        

