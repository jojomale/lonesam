#!/usr/bin/env python
# coding: utf-8

"""
Plot processed data using analysis module.
"""

import numpy as np
from obspy.core import UTCDateTime as UTC

from data_quality_control import sds_db as analysis

# Set verbosity: "ERROR" < "WARNING" < "INFO" < "DEBUG"
analysis.logger.setLevel("INFO")

# Station id
network = 'GR'
station = 'BFO'
location = ''
channel = 'HHZ'

# Data source
datadir = '/home/lehr/sds/processed/'
#datadir = "../data/"


# Date range that you want to inspect
startdate = UTC("2020-12-28")
enddate = UTC("2021-01-05")


# Choose time range 
## Full time range (= all 24h)
stime, etime = "00:00", "23:00"

## Time range crossing midnight
#stime, etime = "19:00", "05:00"



def main():
    stationcode = "{}.{}.{}.{}".format(network, station, 
                                   location, channel)
    analyzer = analysis.SDSDataBaseAnalyzer(
        startdate, enddate, datadir, stationcode)

    data = analyzer.get_data(#["amplitudes", "psds"], 
                        # stime="00:00", etime="23:00"
                        stime=stime, etime=etime)

    #fig1, fig2 = analyzer.plot()
    fig1 = analyzer.plot_amplitudes()
    fig2 = analyzer.plot_psds(np.log)

    fig1.show()
    fig2.show()

if __name__=="__main__":
    main()