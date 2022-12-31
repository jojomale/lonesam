from pathlib import Path
import numpy as np
from obspy.core import UTCDateTime as UTC

from data_quality_control import analysis, dqclogging

analysis.logger.setLevel("INFO")

# Station id
network = 'GR'
station = 'BFO'
location = ''
channel = 'BHZ'

# Data source
datadir = Path('output')

# Date range that you want to inspect
starttime = UTC("2000-01-01")
endtime = UTC("2022-12-31")

dqclogging.configure_handlers(analysis.logger, "INFO", "DEBUG", "dqc_analysis_test.log")


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
    fig_cont.savefig(datadir.joinpath("spectrogram_timerange.png"))


if __name__=="__main__":
    main()