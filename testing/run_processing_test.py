#!/usr/bin/env python
# coding: utf-8
"""
Run processing of raw seismic data. We extract
75-percentile amplitude and spectra. 
"""
from pathlib import Path
import os
from obspy.core import UTCDateTime as UTC

from data_quality_control import sds_db, dqclogging

network = 'GR'
station = 'BFO'
location = ''
channel = 'HHZ'

nscl_code = "GR.BFO..HHZ"

overlap = 60 #3600
fmin, fmax = (4, 14)
nperseg = 2048
winlen_in_s = 3600
proclen = 24*3600
sampling_rate = 100

outdir = Path('output')

sds_root = Path('../sample_sds')
inventory_routing_type = "eida-routing"


startdate = UTC("2020-12-20")
enddate = UTC("2021-01-10")

logfilename = "log/test_processing.log"

# Set verbosity: "ERROR" < "WARNING" < "INFO" < "DEBUG"
dqclogging.configure_handlers(sds_db.logger, "INFO", "DEBUG", 
    logfilename, use_new_file=True)



def main():
    # Clean output directory
    if outdir.exists():
        for f in outdir.glob("*"):
            os.remove(f)
        outdir.rmdir()

    assert sds_root.exists(), "sample sds database {} not found".format(str(sds_root))


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
    
    assert outdir.exists(), "output directory {} was not created".format(str(outdir))

    for t0, t1 in processor.iter_time(startdate, enddate):
        fname = outdir.joinpath("{}.{}.{}.{}_{:04d}.hdf5".format(
                network, station, location, channel, 
                t0.year))
        assert fname.is_file(), \
                "Output file {} was not created".format(fname)
        assert fname.stat().st_size > 100000, \
                "file {} seems to small".format(fname)

    print("Processing finished. Use `run_analysis_test.py` to view results.")

if __name__ == "__main__":
    main()



# processor = sds_db.SDSProcessor(
#         network=network, 
#         station=station,
#         channel=channel,
#         sds_root=sds_root,
#         inventory_routing_type=inventory_routing_type,
#         outdir=outdir, 
#         fileunit="month",
#         # Processing parameters
#         overlap=overlap, nperseg=nperseg, 
#         winlen_seconds=winlen_in_s, 
#         proclen_seconds=proclen,
#         amplitude_frequencies=(fmin, fmax))

# print(processor)


# # In[9]:


# startdate = UTC("2020-12-20")
# enddate = UTC("2021-01-15")


# # In[10]:

# processor.process(startdate, enddate, force_new_file=True)
# assert outdir.exists(), "output directory {} was not created".format(str(outdir))

# for t0, t1 in processor.iter_time(startdate, enddate):
#         fname = outdir.joinpath("{}.{}.{}.{}_{:04d}-{:02d}.hdf5".format(
#                 network, station, location, channel, 
#                 t0.year, t0.month))
#         assert fname.is_file(), \
#                 "Output file {} was not created".format(fname)
#         assert fname.stat().st_size > 10000, \
#                 "file {} seems to small".format(fname)

# # `filunit="month"` produces output files with ending `YYYY-MM.hdf5`. Note that 
# # these files are only about 1/12 of the size of the yearly files, indicating
# # that they cover only one month rather than 1 year of data.



