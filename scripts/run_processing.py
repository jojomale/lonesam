#!/usr/bin/env python
# coding: utf-8
"""
Run processing of raw seismic data. We extract
75-percentile amplitude and spectra. 
"""

import os
from obspy.core import UTCDateTime as UTC

from data_quality_control import sds_db

network = 'GR'
station = 'BFO'
location = ''
channel = 'HHZ'

overlap = 60 #3600
fmin, fmax = (4, 14)
nperseg = 2048
winlen_in_s = 3600
proclen = 24*3600

outdir = '/home/lehr/Projects/data_quality_control/processed/'

sds_root = os.path.abspath('/sds/')
inventory_routing_type = "eida-routing"


startdate = UTC("2020-01-01")
enddate = UTC("2021-12-31")

# Set verbosity: "ERROR" < "WARNING" < "INFO" < "DEBUG"
sds_db.logger.setLevel("INFO")
sds_db.base.logger.setLevel("INFO")


def main():
        processor = sds_db.SDSProcessor(
                network, station, channel, 
                inventory_routing_type, 
                sds_root, outdir=outdir, 
                overlap=overlap, nperseg=nperseg, 
                winlen_seconds=winlen_in_s, 
                proclen_seconds=proclen,
                amplitude_frequencies=(fmin, fmax))

        print(processor)

        processor.process(startdate, enddate, force_new_file=True)


if __name__ == "__main__":
        main()



