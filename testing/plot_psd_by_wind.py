from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from obspy.core import UTCDateTime as UTC

from data_quality_control import analysis, dqclogging

analysis.logger.setLevel("INFO")

# Station id
network = 'GR'
station = 'BFO'
location = ''
channel = 'HHZ'

# Data source
datadir = Path('output')
winddata = "../../testing/winddata_bfo.txt"

# Date range that you want to inspect
starttime = UTC("2000-01-01")
endtime = UTC("2022-12-31")

dqclogging.configure_handlers(analysis.logger, "INFO", "DEBUG", "dqc_analysis_test.log")

def read_interp_winddata(fname, stime, etime, delta):
    winddata = []
    with open(fname, "r") as f:
        f.readline()
        for l in f.readlines():
            l = l.strip()
            date, time, speed = l.split()
            #datetime = UTC("{}T{}".format(date, time))
            datetime = "{}T{}".format(date, time)
            if UTC(datetime) < UTC("2020-01-01"): continue
            #winddata.append([datetime, float(speed)])
            winddata.append([UTC(datetime), float(speed)])

    awind = np.array(winddata)

    #awind_tstamps = awind.copy()
    xp = np.asarray([t.timestamp for 
                    t in awind[:,0]])#, dtype=np.float_)
    fp = np.asarray(awind[:,1], dtype=np.float_)

    #stime = starttime #UTC("2020-01-01")
    #etime = endtime #UTC("2021-12-31")
    x = np.arange(stime.timestamp,
                 etime.timestamp+1,
                 delta, #lyza.winlen_seconds
                 )
    #print(x.shape)

    f = np.interp(x, xp, fp, )
    return f

def main():
    stationcode = "{}.{}.{}.{}".format(network, station, 
                                   location, channel)
    lyza = analysis.Analyzer(datadir, stationcode,
                            fileunit="year")
    
    files = lyza.get_available_datafiles()
    print("\nAll available files in {}:\n".format(lyza.datadir),
             files, "\n")
    
    av_stime, av_etime = lyza.get_available_timerange()
    print("Available time range in {}\n{} - {}\n".format(
            datadir, av_stime, av_etime)
    )

    x, f = read_interp_winddata(winddata, starttime, endtime, lyza.winlen_seconds)

    # Extract times with certain wind speed
    cmin = 3
    cmax = 6
    xsel = x[np.logical_and(f >= cmin, f <= cmax)]
    tlist = [UTC(xi) for xi in xsel]
    print(len(tlist))
    
    DATA = lyza.get_data(tlist)
    print(lyza.starttime, lyza.endtime)
    print(lyza.amps.shape)
    print("PSD shape:", lyza.psds.shape)
    print("len(tlist):", len(tlist))
    
    fig_cont = lyza.plot_spectrogram(cmap=plt.cm.gist_heat, vmin=-0.7, vmax=2.4)
    fig_cont.savefig(datadir.joinpath("spectrogram_wind_{:g}-{:g}mps.png"))


if __name__=="__main__":
    main()