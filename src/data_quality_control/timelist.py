"""
This module provides routines to read time data
and convert them into a list of `UTCDateTime`s
which can be passed to `analysis.Analyzer.get_data()`
"""

import numpy as np
from obspy.core import UTCDateTime as UTC

def read_interp_winddata(fname, stime, etime, delta):
    awind = read_winddata(fname)
    return interpolate_windarray(awind, stime, etime, delta)
    
def read_winddata(fname):
    winddata = []
    with open(fname, "r") as f:
        f.readline()
        for l in f.readlines():
            l = l.strip()
            date, time, speed = l.split()
            #datetime = UTC("{}T{}".format(date, time))
            datetime = "{}T{}".format(date, time)
            winddata.append([UTC(datetime), float(speed)])

    return np.array(winddata)

def interpolate_windarray(awind, stime, etime, delta):
    xp = np.asarray([t.timestamp for 
                    t in awind[:,0]])#, dtype=np.float_)
    fp = np.asarray(awind[:,1], dtype=np.float_)

    ## Add??
    # idx = np.where(np.logical_and(xp>=stime.timestamp, 
    #                               xp<=etime.timestamp))[0]
    # if len(idx) == 0:
    #     raise RuntimeError("No wind data for requested time range.")
    # xp = xp[idx]
    # fp = fp[idx]

    x = np.arange(stime.timestamp,
                 etime.timestamp+1,
                 delta)

    f = np.interp(x, xp, fp, )
    return x, f
    