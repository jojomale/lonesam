#!/usr/bin/env python
# coding: utf-8

"""
Plot processed data using analysis module.
"""
from os import remove
from pathlib import Path
import numpy as np
from obspy.core import UTCDateTime as UTC

from matplotlib import pyplot as plt

from data_quality_control import analysis, dqclogging

# Set verbosity: "ERROR" < "WARNING" < "INFO" < "DEBUG"
#analysis.logger.setLevel("INFO")

# Station id
network = 'GR'
station = 'BFO'
location = ''
channel = 'HHZ'

# Data source
datadir = Path('output')

# Date range that you want to inspect
starttime = UTC("2020-12-20")
endtime = UTC("2021-01-10")

figdir = Path("figures/")

logfilename = "log/dqc_analysis_test.log"


if Path(logfilename).is_file():
        remove(logfilename)
dqclogging.configure_handlers(analysis.logger, "INFO", "DEBUG", logfilename)


def create_random_timelist(starttime, endtime, N):
    times = np.arange(str(starttime.date), str(endtime.date),
             dtype="datetime64[h]")
    tlist = np.sort(np.random.choice(times, N, replace=False))
    return [UTC(str(t)) for t in tlist]


def main():

    assert datadir.exists(), \
        "data directory {} does not exist".format(str(datadir))

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

    # Continuous time range
    DATA = lyza.get_data(starttime, endtime)
    print("\nLoaded data for time range {} - {}:\n".format(
        starttime, endtime), DATA, "\n")
    print(lyza.timeax_psd.shape)
    

    fig_cont = lyza.plot_spectrogram()
    fig_cont.suptitle("time range")
    figfiles = [figdir.joinpath("spectrogram_timerange.png")]
    fig_cont.savefig(figfiles[-1])


    fig_amp, fig_psd = lyza.plot3d()
    
    # Uncomment to show during runtime. Opens browser
    #fig_psd.show()


    # Save 3d-Figures as html-files. Can be opened in browser.
    for flabel, fig in zip(["amp", "psd"], [fig_amp, fig_psd]):
        html = fig.to_html(include_mathjax="cdn")
        figfiles.append(figdir.joinpath("fig3d_{}.html".format(flabel)))
        with open(figfiles[-1], "w") as f:
            f.write(html)


    # Time list
    tlist = create_random_timelist(av_stime, av_etime, 20)
    DATA = lyza.get_data(tlist)
    print("Loaded data for time list", DATA, "\n")
    #print(lyza.psds.shape)
    
    fig_tlist = lyza.plot_spectrogram()
    fig_tlist.suptitle("time list")
    figfiles.append(figdir.joinpath("spectrogram_timelist.png"))
    fig_tlist.savefig(figfiles[-1])

    # fig_cont.show()
    # fig_tlist.show()
    plt.show()
    for f in figfiles:
        assert f.exists(), "Figure {} does not exist".format(str(f))

if __name__=="__main__":
    main()