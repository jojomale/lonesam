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
kernel_shift = 5

figdir = Path("figures/")
datadir = Path("output")
outdir = Path("output/interpolated")

logfilename = "log/dqc_interpolation_test.log"


if Path(logfilename).is_file():
        remove(logfilename)
dqclogging.configure_handlers(analysis.logger, "INFO", "DEBUG", logfilename)


def main():
    outdir.mkdir(exist_ok=True)

    polly = analysis.Interpolator(datadir, nslc_code)
    polly.interpolate(kernel_size, kernel_shift, outdir)

    fig, axs = plt.subplots(2,1, figsize=(8, 10))

    titles = ["raw", "interpolated"]
    fpatterns = [ "*202*hdf5", "interpolated/*hdf5",]

    for i, fpattern in enumerate(fpatterns):
        try: 
            res = base.BaseProcessedData()
            for fname in Path("output/").glob(fpattern):
                print(fname)
                res.extend_from_file(fname)
        except RuntimeWarning:
            pass
        
        ax = axs[i]
        res.plot_psds(np.log10, ax=ax)
        ax.set_title(titles[i])
    plt.show()
    # lyza = analysis.Analyzer(outdir, nslc_code, )

if __name__=="__main__":
    main()