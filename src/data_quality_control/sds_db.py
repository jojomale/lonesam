from glob import glob
import os
#from datetime import timedelta, time

import numpy as np

#from obspy.core import UTCDateTime as UTC
from obspy.clients.filesystem.sds import Client
from obspy.clients.fdsn import RoutingClient, Client as FDSNClient
from obspy import Inventory

from . import base, dqclogging, util
#from .analysis import Analyzer

import logging

# Create the global logger
logger = dqclogging.create_logger()
module_logger = logging.getLogger(logger.name+'.sds_db')

class SDSProcessor(base.GenericProcessor):
    def __init__(self, nslc_code,  
            inventory_or_routing_type,
            sds_root, sds_dict={}, 
            outdir='.', preprocessing=None, 
            fileunit="year", **procparams):

        #location = ""
        dataclient = Client(sds_root, **sds_dict)
        invclient = self._set_inventory(inventory_or_routing_type)

        super().__init__(nslc_code,
                dataclient, invclient, 
                outdir=outdir, preprocessing=preprocessing, 
                fileunit=fileunit, **procparams)
        self.logger = logging.getLogger(module_logger.name+
                            '.'+"SDSProcessor")
        self.logger.setLevel(logging.DEBUG)


    def _set_inventory(self, inventory_or_routing_type):
        if isinstance(inventory_or_routing_type, Inventory):
            return inventory_or_routing_type
        elif isinstance(inventory_or_routing_type, str):
            return util.get_fdsn_or_routing_client(
                inventory_or_routing_type)



    def expand_nslc(self, starttime, endtime):
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
            logger.debug("Looking for data in year %s" % str(y))
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
        self._networks = list(networks)
        self._stations = list(stations)
        self._channels = list(channels)
        self._locations = [""]
        #return networks, stations, channels



# class SDSDataBaseAnalyzer(Analyzer):
#     def __init__(self, *args, **kwargs):
#         super().__init__(*args, **kwargs)


#     # def _get_filenames(self):
        
#     #     filehead = os.path.join(self.datadir, self.stationcode)
#     #     fmtstr = filehead + "_{:04d}.hdf5"
#     #     logger.info("Looking for data file %s" % fmtstr)
#     #     _year = self.sdate.year
#     #     files = []
#     #     while _year <= self.edate.year:
#     #         searchstr = fmtstr.format(_year)
#     #         fnames = glob(searchstr)

#     #         if len(fnames) == 0:
#     #             _year = _year +1
#     #             continue
#     #         else:
#     #             files.append(fnames[0])
#     #         # if len(fnames) > 1:
#     #         #     files.append(self.select_longest(fnames))
#     #         # elif len(fnames) == 0:
#     #         #     _year = _year +1
#     #         #     continue
#     #         # else:
#     #         #     files.append(fnames[0])

#     #         # Get end year of latest file
#     #         ## Remove file-ext and path
#     #         # f, ext = os.path.splitext(files[-1])
#     #         # _endtime = UTC(f.split('_')[-1])
#     #         # if _endtime.year >= self.endtime.year:
#     #         #     break
#     #         _year = _year + 1
#     #     return files


#     def get_selected_timerange(self, 
#                 datasets=['amplitudes', 'psds'], 
#                  stime=None, etime=None):
#         """
#         Get data only between specified hours.

#         Only works on data which was processed with
#         proclen=24*3600 and winlen=3600.
#         """
        
#         DATA = {k : [] for k in datasets}
#         for f in self.iter_files():
            
#             t0 = UTC(*f.attrs["starttime"])
#             t1 = UTC(*f.attrs["endtime"])
#             # We need to add 1 day because t1 is midnight.
#             ndays = timedelta(seconds=t1-t0).days+1
#             logger.debug("Days in file: %s" % (ndays))

#             # Get number of days this year
#             year = UTC(*f.attrs["starttime"]).year
#             _t0 = UTC(year, 1, 1, 0, 0, 0)
#             _t1 = UTC(year+1, 1, 1, 0, 0, 0)
#             dt = timedelta(seconds=_t1-_t0)
#             days_this_year = dt.days
            
            
#             freqax = f["frequency_axis"][:]

#             for dataset in datasets:
#                 d = f[dataset]
#                 shp = d.shape
#                 logger.debug("Shape of %s is %s" % (dataset, shp))

#                 # If year is incomplete, we fill missing days with nans
#                 if ndays < 365:
#                     tmp = np.ones((days_this_year, *shp[1:]))*np.nan
#                     i = t0.julday-1  # Juldays start at 1, -1 to get indices
#                     j = i + shp[0]
#                     tmp[i:j,:] = d[:]
#                 # Otherwise we can use data directly
#                 else:
#                     tmp = d[:]
#                 logger.debug("Target shape of %s is %s " % (dataset, tmp.shape))
                
                
#                 # Extract part of the year that we need from this file
#                 if self.starttime > _t0:
#                     i = self.starttime.julday-1
#                 else:
#                     i = 0
#                 if self.endtime < _t1:
#                     j = self.endtime.julday
#                 else:
#                     j = None
#                 logger.debug("Day index %s and %s, Shape = %s" % (i, j, tmp[i:j].shape))
                
#                 # If user wants only some hours, this is temporarily
#                 # strains memory more than necessary. However if time
#                 # range crosses midnight concatenation becomes quite
#                 # complex. For now, we keep it simple.
#                 DATA[dataset].append(tmp[i:j])
        
#         # Cut out desired time range
#         self._update_time(stime, etime)
#         h0 = self.stime.hour
#         h1 = self.etime.hour+1
#         if h1 == 0 or h1>23:
#             h1 = 24
        

#         for dataset, data in DATA.items():
#             data = np.vstack(data)
#             logger.debug("Final shape of %s is %s" % (dataset, data.shape))
#             # If starttime < endtime indices correspond to 
#             # hours and 1 time frame lies within 1 row 
#             if h0 < h1:
#                 DATA[dataset] = data[:,h0:h1]
#             # if timerange crosses midnight, data of one
#             # frame lies in 2 consecutive rows.
#             else:
#                 DATA[dataset] = np.hstack([data[0:-1,h0:],
#                                 data[1:, :h1]])
#         if "amplitudes" in DATA:
#             self.amps = DATA["amplitudes"]
#         if "psds" in DATA:
#             self.psds = DATA["psds"]
#             self.freqax = freqax
#             DATA["frequency_axis"] = freqax
#         return DATA