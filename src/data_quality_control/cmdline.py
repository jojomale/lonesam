"""
Module defining the command line interface to dataqc suite.

If you "just" want to use the commands please go to :doc:`cmdline`.

Here, we document the functions defining the command line tool.
The CLI is implemented using 
`argparse <https://docs.python.org/3/library/argparse.html#>`_.

General usage is:

.. code-block:: console

    dataqc [-h] {process,plot,available,avail,windfilter,wind} args
    

- subcommands are implemented using argparse.subparsers, but to 
  improve legibility of the code subparsers are defined in 
  functions marked by decorator `@subcommand`.
- for each subparser, there is a function `run_subcommand` 
  containing the actual program.


Warning
---------------
Functions are not intended for use outside of this module.

"""

import sys
import argparse
from pathlib import Path
import time
from datetime import timedelta
import numpy as np
from matplotlib.pyplot import show
from obspy.core import UTCDateTime as UTC

from . import base, sds_db, dqclogging, analysis, timelist
from .util import FNAME_FMTS

logger = dqclogging.create_logger()
module_logger = dqclogging.logging.getLogger(logger.name+'.cmdline')

fileunits = list(FNAME_FMTS.keys())

commons_parser = argparse.ArgumentParser(add_help=False)
commons_parser.add_argument("--loglevel", type=str,
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="set logging level",
        default="INFO")
commons_parser.add_argument("--logfile", type=Path,
        help="Give name for logfile",
        default="dataqc_process.log")
commons_parser.add_argument("--append_logfile", #type=bool,
        action='store_false', 
        #const=True, default=False, nargs="?"
        )


def run_processing(args):
    """
    Run SDS processing job based on cli-arguments.

    Based on script `scripts/run_processing.py`. 
    Computes mean amplitudes and power spectral densities
    of seismic data for given station code using 
    `dataqc.sds_db`. Results are stored in HDF5-files.

    Parameters
    -----------
    args : argparse.Namespace
        parsed cli arguments
    """
    #
    init_args = vars(args)

    proc_args = {k: init_args.pop(k) for k in 
        ['starttime', 'endtime']}

    loglevel =  init_args.pop("loglevel")
    dqclogging.configure_handlers(base.logger, loglevel, 
                loglevel, init_args.pop("logfile"), 
                use_new_file=init_args.pop("append_logfile") )

    processor = sds_db.SDSProcessor(
            **init_args
            )

    try:
        processor.process(**proc_args)
    except Exception as e:
        processor.logger.exception(e)


def run_available(args):
    """
    Requests available HDF5 data based on cli-arguments.
    
    Uses `dataqc.analysis.Analyzer.get_all_available_data()`.
    Expects station code and data directory. Station code can
    contain wildcards.
    
    Parameters
    -----------
    args : argparse.Namespace
        parsed cli arguments
    """
    
    args = vars(args)
    init_args = {k: args[k] for k in 
        ["datadir", "nslc_code", "fileunit"]}
    
    lyza = analysis.Analyzer(**init_args)
    
    files = lyza.get_available_datafiles()
    print("\nAll available files in {}:\n".format(lyza.datadir),
             files, "\n")
    
    av_stime, av_etime = lyza.get_available_timerange()
    print("Available time range in {}\n{} - {}\n".format(
            lyza.datadir, av_stime, av_etime)
    )


def run_plot(args):
    """
    Create spectrograms and plotly-figures from cli-arguments.

    Parameters
    -----------
    args : argparse.Namespace
        parsed cli arguments
    """

    figdir = args.figdir
    init_args = {k: args.__getattribute__(k) for k in ["datadir", "nslc_code", "fileunit"]}
    
    lyza = analysis.Analyzer(**init_args)
    dqclogging.configure_handlers(base.logger, args.loglevel, 
                args.loglevel, args.logfile, args.append_logfile )

    if args.timerange:
        print("timerange")
        starttime, endtime = args.timerange #args_dict.pop("timerange")
    elif args.timelist:
        print("timelist")
        starttime = read_file_as_list_of_utcdatetimes(args.timelist)
        endtime = None
        print(starttime.__class__)
    else:
        print("Using full available timerange")
        starttime, endtime = lyza.get_available_timerange()

    DATA = lyza.get_data(starttime, endtime)
    module_logger.debug("DATA contains {}".format(DATA))

    figname = "{}_{}-{}".format(lyza.stationcode, 
                        lyza.starttime.datetime, 
                        lyza.endtime.datetime)
    module_logger.debug("Figure name base: {}".format(figname))
    fig_cont = lyza.plot_spectrogram()
    fig_cont.savefig(figdir.joinpath("{}_spectrogram.png".format(figname)))
    fig_amp, fig_psd = lyza.plot3d()
    for flabel, fig in zip(["amp", "psd"], [fig_amp, fig_psd]):
        html = fig.to_html(include_mathjax="cdn")
        with open(figdir.joinpath(
                "{}_3d_{}.html".format(figname, flabel)), "w") as f:
            f.write(html)
    if args.show:
        show()


def run_windfilter(args):
    """
    Interpolate and extract list of datetimes 
    based on observation value provided by file.

    Values are first interpolated to indicated time spacing (`delta`)
    which should correspond to the window length used for processing  
    of the seimic data. The resulting list is then filtered for times
    for which minspeed <= value <= maxspeed.
    Designed usecase is a list of wind speed measurements.
    """
    args = vars(args)
    out = args.pop("out")
    func = args.pop("func")
    speed = {k : args.pop(k) for k in 
            ["minspeed", "maxspeed"]}
    if speed["maxspeed"] is None:
        speed["maxspeed"] = 99999.
    x, f = timelist.read_interp_winddata(**args)
    x = x[np.logical_and(f>=speed["minspeed"], f<=speed["maxspeed"])]
    out.write("\n".join([str(UTC(xi)) for xi in x]))


def read_file_as_list_of_utcdatetimes(f):
    """
    Read UTCDateTime-compatible lines in file handler.

    Parameters
    ------------
    f : filehandler
        handle to file containing datetimes per line.
        Other lines are ignored (e.g. if f=stdin)
    """
    # If f is stdin, we need to skip the lines
    # that are not datetimes.
    datetimes = []
    for line in f.readlines():
        line = line.strip()
        try:
            datetimes.append(UTC(line))
        except TypeError:
            continue
    return datetimes


def main():
    """
    Executed if user calls a `dataqc` cli-command.

    - Measures execution time
    - shows help 
    - sets log level of module logger
    - parses arguments
    - calls respective subcommands
    """
    t = time.time()
    parser = main_subparser()
    # If User enters only 'dataqc' we show help of 
    # main parser which lists the subprograms
    if len(sys.argv) < 2:
        parser.parse_args(["-h"])
    
    # Otherwise we call the respective subroutine
    args = parser.parse_args()
    try:
        module_logger.setLevel(args.loglevel)
    except AttributeError:
        module_logger.setLevel("INFO")

    module_logger.info(
        "CLI-Arguments for {}:\n".format(sys.argv[1]) + 
        "{}".format(args))
    
    args.func(args)
    
    runtime = timedelta(seconds=time.time()-t) 
    module_logger.info("Finished. Took {} h".format(runtime))



def main_subparser():
    """
    Returns main parser of cli `dataqc`. Collects subcommands as
    subparsers.
    """
    
    parser = argparse.ArgumentParser(
        prog="dataqc",
        description="Command line " + 
        "interface to dataqc package.\n" + 
        "Dataqc computes average amplitudes and "+
        "power spectral densities of seismic data. Results are "+
        "stored in HDF5 files to allow quick plotting.",
        epilog="Use `dataqc subcommand -h` for details and options on each command.")
    
    subparsers = parser.add_subparsers(title="subcommands", 
        help="one of the subprogram in dataqc",
        )

    # Call all functions registered as subcommands
    for cmd in SUBCOMMANDS.values():
        cmd(subparsers)
    
    return parser


SUBCOMMANDS = dict()
def subcommand(func):
    """
    Decorator that registers functions as subcommands.
    """
    SUBCOMMANDS[func.__name__] = func
    return func


@subcommand
def process(subparsers):
    process = subparsers.add_parser("process",
        parents=[commons_parser],
        description="Compute mean amplitude and spectra of seismic data " +
            "and store results as HDF5. Requires internet access.",
        )
    process.set_defaults(func=run_processing)
    process.add_argument("nslc_code", type=str, 
            help=("station code {network}.{station}.{location}.{channel}," +
            "may contain glob-style wildcards"))
    process.add_argument("inventory_routing_type", type=str,
            help="routing client for inventory",
            choices=["eida-routing", "iris-federator"]
            )
    process.add_argument("sds_root", type=Path,
            help="root-directory of sds-filesystem")
    process.add_argument("starttime", type=UTC, 
            help=("beginning of time range you want to analyze" + 
                    "Give as YYYY-MM-DDThh:mm:ss"))
    process.add_argument("endtime", type=UTC, 
            help="end of time range you want to analyze")

    process.add_argument("-o", "--outdir", type=Path, 
            help="where to put the processed data",
            default=".")
    process.add_argument("--fileunit", type=str, 
            choices=fileunits,
            help="Time span per HDF5-file. ",
            default="year")

    process.add_argument("--overlap", type=int,
            help="seconds by which the data is extended beyond time range "+
                    "to accomodate filter effects etc.",
            default=base.default_processing_params["overlap"])
    process.add_argument("--proclen", type=int,
            help="seconds to process at once, ideally duration of " +
                    "the data file",
            default=base.default_processing_params["proclen_seconds"])
    process.add_argument("--winlen-in-s", type=int,
            help="time over which amplitude and spectra are computed," +
            " in seconds",
            default=base.default_processing_params["winlen_seconds"])
    process.add_argument("--nperseg", type=int,
            help="length of segment for spectral estimation "+
            "(scipy.signal.welch), in samples ",
            default=base.default_processing_params["nperseg"])
    process.add_argument("--amplitude-frequencies", type=float, nargs=2,
            help="min and max frequency of bandpass before "+
            "amplitude analysis.",
            default=base.default_processing_params["amplitude_frequencies"] )
    process.add_argument("--sampling_rate", "--sr", type=float,
            help="Sampling rate at which data are processed. "+
            "Data are resampled if original SR is different.",
            default=base.default_processing_params["sampling_rate"] )
    

@subcommand
def plot(subparsers):
    plot = subparsers.add_parser("plot",
        parents=[commons_parser],
        description="Make 3 different figures of amplitudes and " + 
        "power spectral densities. "+
        "Plotly is used to create interactive 3D-views of amplitudes " +
        "and PSDs with time. They are saved as html-files and can be " +
        "viewed in a browser. " +
        "PSDs are also plotted as classic spectrogram, saved as png and " +
        "optionally shown as interactive matplotlib figure. "
        )
    plot.set_defaults(func=run_plot)
    plot.add_argument("nslc_code", type=str, 
            help=("station code {network}.{station}.{location}.{channel}," +
                "May *not* contain wildcards here!"))
    plot.add_argument("datadir", type=Path, 
            help="where to look for processed data",
            default=".")
    
    plot.add_argument("--fileunit", type=str, 
            choices=fileunits,
            help="Time span per HDF5-file. ",
            default="year")
    plot.add_argument("-o", "--figdir", type=Path,
            help="Where to store figures.",
            default=".")
    plot.add_argument("-s", "--show",
        action="store_true",
        help="If given spectrogram plot is opened as matplotlib figure.")

    group = plot.add_mutually_exclusive_group()
    group.add_argument("-l", "--timelist", type=argparse.FileType("r"), 
            help=("Plot spectrograms using time list." + 
                    "Can be used as flag to read times from stdin or" + 
                    "given a file with datetimes."),
            nargs="?",
            const=sys.stdin)
    group.add_argument("-r", "--timerange", type=UTC, nargs=2,
            help=("Start and end of time range you want to plot. "+ 
                "Give as YYYY-MM-DDThh:mm:ss, " + 
                "endtime can be None to use current time."),
                    )
    

@subcommand
def avail(subparsers):
    avail = subparsers.add_parser("available",
        aliases=["avail"],
        parents=[commons_parser],
        description="Print available HDF5 files and "+ 
            "covered time range for given code in datadir.",
        )
    avail.set_defaults(func=run_available)
    avail.add_argument("nslc_code", type=str, 
        help=("station code {network}.{station}.{location}.{channel}," +
            "may contain glob-style wildcards"))
    avail.add_argument("datadir", type=Path, 
            help="where to look for processed data",
            default=".")
    avail.add_argument("--fileunit", type=str, 
            help="Time span per HDF5-file. ",
            default="year")
    

@subcommand
def windfilter(subparsers):
    windfilter = subparsers.add_parser("windfilter",
        aliases=["wind"],
        #parents=[commons_parser],
        # "Extract list of datetimes based on observation value.
        # Designed usecase is a list of wind speed measurements.
        # "
        description="Interpolate and extract list of datetimes "+
        "based on observation value provided by file. " +
        "Values are first interpolated to indicated time spacing (`delta`) " + 
        "which should correspond to the window length used for processing " + 
        "of the seimic data. The resulting list is then filtered for times "+
        "for which minspeed <= value <= maxspeed. " +
        "Designed usecase is a list of wind speed measurements."
        )
    windfilter.set_defaults(func=run_windfilter)
    windfilter.add_argument("fname", type=Path, 
        help=("Name of wind data file." +
            "File must contain 3 columns: " + 
            "date, time and speed."))
    windfilter.add_argument("stime", type=UTC, 
            help="Start of time range as YYYY-MM-DDThh:mm:ss",
            )
    windfilter.add_argument("etime", type=UTC, 
            help="End of time range as YYYY-MM-DDThh:mm:ss",
            )
    windfilter.add_argument("delta", type=float, 
            help="increment of time axis to which wind data will" +
                    "be interpolated. In seconds." + 
                    "Should be same as window length over which "+
                    "amplitudes and psds were computed.",
            )
    windfilter.add_argument("minspeed", type=float,
            help="Minimum windspeed to select.")
    windfilter.add_argument("maxspeed", type=float,
            help="Maximum windspeed to select.",
            nargs="?")
    windfilter.add_argument("out", type=argparse.FileType("w"),
            nargs="?", default=sys.stdout)


if __name__ == "__main__":
    main()
    