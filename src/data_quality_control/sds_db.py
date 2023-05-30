from glob import glob
import os

import numpy as np

from obspy.clients.filesystem.sds import Client
from obspy import Inventory

from . import base, dqclogging, util

import logging

# Create the global logger
logger = dqclogging.create_logger()
module_logger = logging.getLogger(logger.name+'.sds_db')


class SDSProcessor(base.GenericProcessor):
    """
    Processor tailored to SDS data bases. Mainly provides
    a customized method
    """
    def __init__(self, nslc_code,  
            inventory_or_routing_type,
            sds_root, sds_dict={}, 
            outdir='.', preprocessing=None, 
            fileunit="year", 
            FMTSTR = None,
            **procparams):

        #location = ""
        dataclient = Client(sds_root, **sds_dict)
        invclient = self._set_inventory(inventory_or_routing_type)

        if FMTSTR is not None:
            dataclient.FMTSTR = FMTSTR

        super().__init__(nslc_code,
                dataclient, invclient, 
                outdir=outdir, preprocessing=preprocessing, 
                fileunit=fileunit, **procparams)
        
        self.logger = logging.getLogger(module_logger.name+
                            '.'+"SDSProcessor")
        self.logger.setLevel(logging.DEBUG)


    def _set_inventory(self, inventory_or_routing_type):
        """
        Decide if ``inventory_or_routing_type`` is an
        :py:class:`obspy.Inventory` or a filename or
        a string indicating the routing client. 
        """
        if isinstance(inventory_or_routing_type, Inventory):
            return inventory_or_routing_type
        elif isinstance(inventory_or_routing_type, str):
            return util.get_fdsn_or_routing_client(
                inventory_or_routing_type)


    def expand_nslc(self, starttime, endtime):
        """
        Expand network, station, channel attributes (which can contain
        wildcards) to determine all unique networks, stations and channels 
        available within requested time range.

        Parameters
        --------------
        starttime : UTCDatetime
            begin of time range, we only use the year
        endtime : UTCDateTime
            end of time range, we only use the year


        Sets new attributes :py:attr:`._networks`, 
        :py:attr:`._stations`, :py:attr:`._locations`, 
        :py:attr:`._channels`

        """

        networks = set()
        stations = set()
        channels = set()
        years = np.arange(starttime.year, endtime.year+1)
        for y in years:
            self.logger.debug("Looking for data in year %s" % str(y))
            for n in glob(os.path.join(self.client.sds_root, 
                                        str(y), self.network)):
                nw = os.path.split(n)[-1]
                networks.add(nw)
                self.logger.debug("%s" % nw)
                for s in glob(os.path.join(n, self.station)):
                    sn = os.path.split(s)[-1]
                    stations.add(sn)
                    self.logger.debug("%s" % sn)
                    for c in glob(os.path.join(s, self.channel)+
                                    '.{}'.format(self.client.sds_type)):
                        ch = os.path.split(c)[-1].split('.')[0]
                        channels.add(ch)
                        self.logger.debug("%s" % ch)
        self._networks = list(networks)
        self._stations = list(stations)
        self._channels = list(channels)
        self._locations = [""]
       