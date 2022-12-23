#!/usr/bin/env python
# coding: utf-8

"""
Plot processed data using analysis module.
"""
from pathlib import Path
import numpy as np
from obspy.core import UTCDateTime as UTC

from data_quality_control import analysis, dqclogging



# Station id
network = 'GR'
station = 'BFO'
location = ''
channel = 'HHZ'


# Date range that you want to inspect
starttime = UTC("2020-12-25")
endtime = UTC("2021-01-10")


# Data source
datadir = Path('../sample_output/run_processing')
figdir = Path("../sample_output/run_processing/figures/")

logfilename = "../sample_output/run_processing/log/processing.log"

# Set verbosity: "ERROR" < "WARNING" < "INFO" < "DEBUG"
dqclogging.configure_handlers(analysis.logger, "INFO", "DEBUG", 
    logfilename, use_new_file=True)


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

    # Continuous time range
    DATA = lyza.get_data(starttime, endtime)
    print("\nLoaded data for time range {} - {}:\n".format(
        starttime, endtime), DATA, "\n")
    print(lyza.timeax_psd.shape)
    

    fig_cont = lyza.plot_spectrogram()
    fig_cont.savefig(figdir.joinpath("spectrogram_timerange.png"))


    fig_amp, fig_psd = lyza.plot3d()
    
    # Uncomment to show during runtime. Opens browser
    #fig_psd.show()


    # Save 3d-Figures as html-files. Can be opened in browser.
    for flabel, fig in zip(["amp", "psd"], [fig_amp, fig_psd]):
        html = fig.to_html(include_mathjax="cdn")
        with open(figdir.joinpath("fig3d_{}.html".format(flabel)), "w") as f:
            f.write(html)

if __name__=="__main__":
    main()