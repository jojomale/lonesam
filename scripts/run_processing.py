#!/usr/bin/env python
# coding: utf-8
"""
Run processing of raw seismic data. We extract
75-percentile amplitude and spectra. 
"""
from pathlib import Path
#import os
from obspy.core import UTCDateTime as UTC

from data_quality_control import sds_db, dqclogging

# network = 'GR'
# station = 'BFO'
# location = ''
# channel = 'BHZ'
nscl_code = "GR.BFO..BHZ"

overlap = 60 #3600
fmin, fmax = (4, 14)
nperseg = 2048
winlen_in_s = 3600
proclen = 24*3600
sampling_rate = 100

outdir = Path('../sample_output/run_processing')

sds_root = Path('../sample_sds')
inventory_routing_type = "eida-routing"


startdate = UTC("2020-12-25")
enddate = UTC("2021-01-10")

logfilename = "log/processing.log"

# Set verbosity: "ERROR" < "WARNING" < "INFO" < "DEBUG"
dqclogging.configure_handlers(sds_db.logger, "INFO", "DEBUG", 
    logfilename, use_new_file=True)



def main():
        processor = sds_db.SDSProcessor(
                nscl_code,
                inventory_routing_type, 
                sds_root, outdir=outdir, 
                overlap=overlap, nperseg=nperseg, 
                winlen_seconds=winlen_in_s, 
                proclen_seconds=proclen,
                amplitude_frequencies=(fmin, fmax),
                sampling_rate=sampling_rate)

        print(processor)

        try:
                processor.process(startdate, enddate, force_new_file=True)
        except Exception as e:
                processor.logger.exception(e)

        processor.logger.info("Processing finished. Use `run_analysis.py` to view results.")

if __name__ == "__main__":
        main()



