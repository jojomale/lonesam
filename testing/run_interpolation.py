from os import remove
from pathlib import Path
import numpy as np
from obspy.core import UTCDateTime as UTC

from matplotlib import pyplot as plt

from data_quality_control import analysis, dqclogging, base


nslc_code = "GR.BFO..BHZ"

starttime = UTC("2020-12-20")
endtime = UTC("2021-01-10")

kernel_size = 6
kernel_shift = 3

figdir = Path("figures/")
datadir = Path("output")
outdir = Path("output/interpolated")

logfilename = "log/dqc_interpolation_test.log"

startdate = UTC("2020-12-24")
enddate = UTC("2021-01-15")

#if Path(logfilename).is_file():
#        remove(logfilename)
dqclogging.configure_handlers("INFO", "DEBUG", logfilename, 
                                use_new_file=True)


def main():
    outdir.mkdir(exist_ok=True)

    # Interpolation
    ## Initiate interpolator
    polly = analysis.Interpolator(datadir, nslc_code, 
                                kernel_size=kernel_size, 
                                kernel_shift=kernel_shift)


    ## Start interpolation over whole available time range
    polly.interpolate(outdir, force_new_file=True)


    fig, axs = plt.subplots(2,1, figsize=(8, 10))

    titles = ["raw", "interpolated"]
    datadirs = [datadir, outdir,]

    for i, _dir in enumerate(datadirs):
        ax = axs[i]
        # View results
        lyza = analysis.Analyzer(_dir, nslc_code,
                                    fileunit="year")

        lyza.get_data(startdate, enddate)
        lyza.plot_spectrogram(ax=ax)
        ax.set_title(titles[i])

    plt.show()
    # lyza = analysis.Analyzer(outdir, nslc_code, )

if __name__=="__main__":
    main()