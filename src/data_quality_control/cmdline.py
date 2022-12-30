"""
Command line interface for accessing the old processing module.

Start a processing job via CLI. 

Adaptation to newer version should be possible.

Usage
-------
.. code-block:: console

        dataqc [-h] [-n NETWORK] [-s STATION] [-c CHANNEL] 
                [--sds-root SDS_ROOT] 
                [--inventory-routing-type {eida-routing,iris-federator}]
                [--overlap OVERLAP] [--proclen PROCLEN] 
                [--winlen-in-s WINLEN_IN_S] [--nperseg NPERSEG]
                [--amplitude-frequencies AMPLITUDE_FREQUENCIES AMPLITUDE_FREQUENCIES] 
                [--configfile CONFIGFILE] [-v {DEBUG,INFO,WARNING,ERROR,CRITICAL}]
                startdate enddate outdir


Optional arguments override settings in configfile. If a parameter
is not set via opt. arg. or is defined in config-file, we use
default values.
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

logger = dqclogging.create_logger()
module_logger = dqclogging.logging.getLogger(logger.name+'.cmdline')

default_settings = dict(network = '*',
    station = '*',
    channel = '*',
    overlap = 60,
    amplitude_frequencies = (4,14),
    nperseg = 2048,
    winlen_in_s = 3600,
    proclen = 24*3600,
    sds_root = Path('.'),
    inventory_routing_type = "eida-routing",
    sds_client_dict = {})


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


# def parse_argument():
#     """
#     Get arguments from commandline and prepare for
#     use in processing.RawDataProcessor.
#     """
#     parser = argparse.ArgumentParser(description="Compute "+
#         "amplitude levels and power spectral densities from "+
#         "raw seismic data in an sds-filesystem and store in "+
#         "HDF5-files.")
    


#     parser.add_argument("network", type=str, 
#             help="network(s), may contain glob-style wildcards")
#     parser.add_argument("station", type=str, 
#             help="station(s), may contain glob-style wildcards")
#     parser.add_argument("channel", type=str, 
#             help="channel(s), may contain glob-style wildcards")
#     parser.add_argument("inventory_routing_type", type=str,
#             help="routing client for inventory",
#             choices=["eida-routing", "iris-federator"])
#     parser.add_argument("sds_root", type=Path,
#             help="root-directory of sds-filesystem")
    
#     parser.add_argument("starttime", type=UTC, 
#             help="beginning of time range you want to analyze")
#     parser.add_argument("endtime", type=UTC, 
#             help="end of time range you want to analyze")
    
    
#     parser.add_argument("-o", "--outdir", type=Path, 
#             help="where to put the processed data",
#             default=".")
#     parser.add_argument("--overlap", type=int,
#             help="seconds by which the data is extended beyond time range "+
#                     "to accomodate filter effects etc.",
#             default=base.default_processing_params["overlap"])
#     parser.add_argument("--proclen", type=int,
#             help="seconds to process at once, ideally duration of " +
#                     "the data file",
#             default=base.default_processing_params["proclen_seconds"])
#     parser.add_argument("--winlen-in-s", type=int,
#             help="time over which amplitude and spectra are computed," +
#             " in seconds",
#             default=base.default_processing_params["winlen_seconds"])
#     parser.add_argument("--nperseg", type=int,
#             help="length of segment for spectral estimation "+
#             "(scipy.signal.welch), in samples ",
#             default=base.default_processing_params["nperseg"])
#     parser.add_argument("--amplitude-frequencies", type=float, nargs=2,
#             help="min and max frequency of bandpass before "+
#             "amplitude analysis.",
#             default=base.default_processing_params["amplitude_frequencies"] )
#     parser.add_argument("--configfile", type=Path,
#             help="file with parameters. additionally given parameters "+
#             "override those from file")
#     parser.add_argument("-f", "--force-new-file",
#             action="store_true",
#             help="overrides existing files if given",)
    

#     parser.add_argument("-v", "--verbosity", type=str,
#                     choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
#                     help="set logging level",
#                     default="DEBUG")

#     print(parser.parse_args())
#     arg_dict= vars(parser.parse_args())
#     proc_args = {k: arg_dict.pop(k) for k in ['starttime', 'endtime']}
#     loglevel = arg_dict.pop("verbosity")
#     module_logger.setLevel(loglevel)
#     #print(args)
#     return arg_dict, proc_args


class DataqcMain():
    def __init__(self) -> None:
        parser = argparse.ArgumentParser(
                prog="dataqc",
                description="Command line " + 
                        "interface to dataqc package",
                usage="dataqc command options",
                epilog="Use `dataqc subcommand -h` for details and options on each command.")

        parser.add_argument("command", help="commands")
        parser.add_argument("-v", "--verbosity", type=str,
                    choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                    help="set logging level",
                    default="INFO")
        # parser.add_argument("-v","--version", help="show version and   exit", action="version", version='1.0')
        print("args",parser.parse_args())
        args = parser.parse_args(sys.argv[1:2])
        print("args",args)
        if not hasattr(self, args.command):
            print('Unrecognized subcommand')
            parser.print_help()
            exit(1)
       
        # use dispatch pattern to invoke method with same name
        getattr(self, args.command)()


    def process(self):
        parser = argparse.ArgumentParser(
                description="Compute mean amplitude and spectra of seismic data"
        )
        parser.add_argument("nslc_code", type=str, 
            help=("station code {network}.{station}.{location}.{channel}," +
                "may contain glob-style wildcards"))
        parser.add_argument("inventory_routing_type", type=str,
                help="routing client for inventory",
                choices=["eida-routing", "iris-federator"]
                )
        parser.add_argument("sds_root", type=Path,
                help="root-directory of sds-filesystem")
        parser.add_argument("starttime", type=UTC, 
                help=("beginning of time range you want to analyze" + 
                        "Give as YYYY-MM-DDThh:mm:ss"))
        parser.add_argument("endtime", type=UTC, 
                help="end of time range you want to analyze")

        parser.add_argument("-o", "--outdir", type=Path, 
                help="where to put the parsered data",
                default=".")
        parser.add_argument("--fileunit", type=str, 
                help="where to put the processed data",
                default="year")

        parser.add_argument("--overlap", type=int,
                help="seconds by which the data is extended beyond time range "+
                        "to accomodate filter effects etc.",
                default=base.default_processing_params["overlap"])
        parser.add_argument("--proclen", type=int,
                help="seconds to process at once, ideally duration of " +
                        "the data file",
                default=base.default_processing_params["proclen_seconds"])
        parser.add_argument("--winlen-in-s", type=int,
                help="time over which amplitude and spectra are computed," +
                " in seconds",
                default=base.default_processing_params["winlen_seconds"])
        parser.add_argument("--nperseg", type=int,
                help="length of segment for spectral estimation "+
                "(scipy.signal.welch), in samples ",
                default=base.default_processing_params["nperseg"])
        parser.add_argument("--amplitude-frequencies", type=float, nargs=2,
                help="min and max frequency of bandpass before "+
                "amplitude analysis.",
                default=base.default_processing_params["amplitude_frequencies"] )

        # parser.add_argument("-f", "--force-new-file",
        #         action="store_true",
        #         help="overrides existing files if given",)
        
        parser.add_argument("--loglevel", type=str,
                choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                help="set logging level",
                default="INFO")
        parser.add_argument("--logfile", type=Path,
                help="Give name for logfile",
                default="dataqc_plot.log")
        parser.add_argument("--append_logfile", #type=bool,
                action='store_false', 
                #const=True, default=False, nargs="?"
                )
        args = parser.parse_args(sys.argv[2:])
        
        run_processing(args)


    def available(self):
        parser = argparse.ArgumentParser(
                description="Print available HDF5 files and covered time range for given code in datadir"
        )
        parser.add_argument("nslc_code", type=str, 
            help=("station code {network}.{station}.{location}.{channel}," +
                "may contain glob-style wildcards"))
        parser.add_argument("datadir", type=Path, 
                help="where to look for processed data",
                default=".")
        parser.add_argument("--fileunit", type=str, 
                help="where to put the processed data",
                default="year")
        args = parser.parse_args(sys.argv[2:])
        print(args)
        run_available(args)


    def plot(self):
        parser = argparse.ArgumentParser(
                description="Plot spectrogram and/or amplitude"
        )
        parser.add_argument("nslc_code", type=str, 
            help=("station code {network}.{station}.{location}.{channel}," +
                "may contain glob-style wildcards"))
        parser.add_argument("datadir", type=Path, 
                help="where to look for processed data",
                default=".")
        
        parser.add_argument("--fileunit", type=str, 
                help="where to put the processed data",
                default="year")

        parser.add_argument("--figdir", type=Path,
                help="where to store figures",
                default=".")

        parser.add_argument("--loglevel", type=str,
                    choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                    help="set logging level",
                    default="INFO")
        parser.add_argument("--logfile", type=Path,
                help="Give name for logfile",
                default="dataqc_plot.log")
        parser.add_argument("--append_logfile", #type=bool,
                action='store_false', 
                #const=True, default=False, nargs="?"
                )
        group = parser.add_mutually_exclusive_group()
        group.add_argument("-l", "--timelist", type=argparse.FileType("r"), 
                help=("Plot spectrograms using time list." + 
                        "Can be used as flag to read times from stdin or" + 
                        "given a file with datetimes."),
                nargs="?",
                const=sys.stdin)
        group.add_argument("-r", "--timerange", type=UTC, nargs=2,
                help=("start and end of time range you want to analyze"+ 
                        "Give as YYYY-MM-DDThh:mm:ss, endtime can be None to use current time."),
                        )
        args = parser.parse_args(sys.argv[2:])
        print(args)
        if args.timerange:
            print("timerange")
        elif args.timelist:
            print("timelist")
        else:
            print("Nothing")
        run_plot(args)


    def windfilter(self):
        parser = argparse.ArgumentParser(
                description="Interpolate and extract wind"
        )
        parser.add_argument("fname", type=Path, 
            help=("name of wind data file" +
                "."))
        parser.add_argument("stime", type=UTC, 
                help="starttime",
                )
        parser.add_argument("etime", type=UTC, 
                help="endtime",
                )
        parser.add_argument("delta", type=float, 
                help="increment of time axis to which wind data will" +
                        "be interpolated. In seconds." + 
                        "Should be same as window length.",
                )
        parser.add_argument("minspeed", type=float,
                help="Minimum windspeed")
        parser.add_argument("maxspeed", type=float,
                help="Maximum windspeed",
                nargs="?")
        parser.add_argument("out", type=argparse.FileType("w"),
                nargs="?", default=sys.stdout)
        args = parser.parse_args(sys.argv[2:])
        run_windfilter(args)

    def print(self):
        parser = argparse.ArgumentParser(description="Test piping")
        parser.add_argument("timerange", nargs="?", type=str)
        parser.add_argument("timelist", 
                type=argparse.FileType("r"),
                help="list of items to print",
                nargs="?",
                default=sys.stdin)
        # parser.add_argument("timerange", type=list,
        #         help="time range", nargs=)
        args = parser.parse_args(sys.argv[2:])
        print(args)
        print(args.timerange)
        with args.timelist as input:
            newlist = [line for line in input.readlines()]

        print(newlist)

    def experiment2(self):
        parser = argparse.ArgumentParser(description="test variable input")
        parser.add_argument("time", nargs="*", default=sys.stdin)
        args = parser.parse_args(sys.argv[2:])

        if len(args.time) == 1:
            time = args.time[0]
            if isinstance(time, argparse.FileType):
                print("FileHandler")
            else:
                try:
                    starttime = UTC(time)
                    endtime = None
                except TypeError:
                    try:
                        f = open(time, "r")
                    except FileNotFoundError:
                        raise RuntimeError("Wrong input. Assuming input is file named {}".format(time))
        elif len(args.time) == 2:
            try:
                starttime = UTC(args.time[0])
                endtime = UTC(args.time[1])
            except TypeError:
                raise RuntimeError("Wrong input format")
        else:
            raise RuntimeError("Illegal number of time arguments")

                    
        print(args)




def run_processing(args):
    t = time.time()
    #print(args)
    init_args = vars(args)

    proc_args = {k: init_args.pop(k) for k in 
        ['starttime', 'endtime']}
    #print(init_args)

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

    runtime = timedelta(seconds=time.time()-t) 
    processor.logger.info("Finished. Took {} h".format(runtime))


def run_available(args):
#     print("Available")
#     print(args)
#     print(vars(args))
#     print(vars(args)["fileunit"])
    init_args = {k: args.__getattribute__(k) for k in 
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
    #args_dict = vars(args)
    
    #args = vars(args)
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
    print(DATA)

    figname = "{}_{}-{}".format(lyza.stationcode, 
                        lyza.starttime.datetime, 
                        lyza.endtime.datetime)
    print("FIGNAME", figname)
    fig_cont = lyza.plot_spectrogram()
    print(fig_cont)
    fig_cont.savefig(figdir.joinpath("{}_spectrogram.png".format(figname)))
    fig_amp, fig_psd = lyza.plot3d()
    #print(fig_amp.to_html())
    for flabel, fig in zip(["amp", "psd"], [fig_amp, fig_psd]):
        html = fig.to_html(include_mathjax="cdn")
        with open(figdir.joinpath(
                "{}_3d_{}.html".format(figname, flabel)), "w") as f:
            f.write(html)
    if args.show:
        show()


def run_windfilter(args):
    args = vars(args)
    out = args.pop("out")
    print(out)
    func = args.pop("func")
    speed = {k : args.pop(k) for k in 
            ["minspeed", "maxspeed"]}
    if speed["maxspeed"] is None:
        speed["maxspeed"] = 99999.
    x, f = timelist.read_interp_winddata(**args)
    x = x[np.logical_and(f>=speed["minspeed"], f<=speed["maxspeed"])]
    print(x)
    print(out)
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
    #DataqcMain()
    main_subparser()


def main_subparser():
    #print("Hello")
    # Main parser
    parser = argparse.ArgumentParser(
        prog="dataqc",
        description="Command line " + 
        "interface to dataqc package",
        epilog="Use `dataqc subcommand -h` for details and options on each command.")
    

    subparsers = parser.add_subparsers(title="subcommands", 
        help="one of the subprogram in dataqc",
        )

    process(subparsers)
    plot(subparsers)
    avail(subparsers)
    windfilter(subparsers)

    # If User enters only 'dataqc' we show help of 
    # main parser which lists the subprograms
    if len(sys.argv) < 2:
        parser.parse_args(["-h"])
    
    # Otherwise we call the respective subroutine
    args = parser.parse_args()
    #print(args)
    try:
        module_logger.setLevel(args.loglevel)
    except AttributeError:
        module_logger.setLevel("INFO")

    module_logger.info("CLI-Arguments:\n" + 
                "{}".format(args))
    args.func(args)
    
    #print('Finish')

def process(subparsers):
    process = subparsers.add_parser("process",
        parents=[commons_parser],
        description="Compute mean amplitude and spectra of seismic data",
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
            help="where to put the processed data",
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
    #return process

def plot(subparsers):
    plot = subparsers.add_parser("plot",
        parents=[commons_parser],
        description="Create plots",
        )
    plot.set_defaults(func=run_plot)
    plot.add_argument("nslc_code", type=str, 
            help=("station code {network}.{station}.{location}.{channel}," +
                "May *not* contain wildcards!"))
    plot.add_argument("datadir", type=Path, 
            help="where to look for processed data",
            default=".")
    
    plot.add_argument("--fileunit", type=str, 
            help="where to put the processed data",
            default="year")
    plot.add_argument("-o", "--figdir", type=Path,
            help="where to store figures",
            default=".")
    plot.add_argument("-s", "--show",
        action="store_true",
        help="if given spectrogram plot is opened.")

    group = plot.add_mutually_exclusive_group()
    group.add_argument("-l", "--timelist", type=argparse.FileType("r"), 
            help=("Plot spectrograms using time list." + 
                    "Can be used as flag to read times from stdin or" + 
                    "given a file with datetimes."),
            nargs="?",
            const=sys.stdin)
    group.add_argument("-r", "--timerange", type=UTC, nargs=2,
            help=("start and end of time range you want to analyze"+ 
                    "Give as YYYY-MM-DDThh:mm:ss, endtime can be None to use current time."),
                    )
    #return subparsers

def avail(subparsers):
    avail = subparsers.add_parser("available",
        aliases=["avail"],
        parents=[commons_parser],
        description="Print available HDF5 files and "+ 
            "covered time range for given code in datadir",
        )
    avail.set_defaults(func=run_available)
    avail.add_argument("nslc_code", type=str, 
        help=("station code {network}.{station}.{location}.{channel}," +
            "may contain glob-style wildcards"))
    avail.add_argument("datadir", type=Path, 
            help="where to look for processed data",
            default=".")
    avail.add_argument("--fileunit", type=str, 
            help="where to put the processed data",
            default="year")
    

def windfilter(subparsers):
    windfilter = subparsers.add_parser("windfilter",
        aliases=["wind"],
        #parents=[commons_parser],
        description="Interpolate and extract wind"
        )
    windfilter.set_defaults(func=run_windfilter)
    windfilter.add_argument("fname", type=Path, 
        help=("name of wind data file" +
            "."))
    windfilter.add_argument("stime", type=UTC, 
            help="starttime",
            )
    windfilter.add_argument("etime", type=UTC, 
            help="endtime",
            )
    windfilter.add_argument("delta", type=float, 
            help="increment of time axis to which wind data will" +
                    "be interpolated. In seconds." + 
                    "Should be same as window length.",
            )
    windfilter.add_argument("minspeed", type=float,
            help="Minimum windspeed")
    windfilter.add_argument("maxspeed", type=float,
            help="Maximum windspeed",
            nargs="?")
    windfilter.add_argument("out", type=argparse.FileType("w"),
            nargs="?", default=sys.stdout)



    


# def main_subparser():
#     #print("Hello")
#     # Main parser
#     parser = argparse.ArgumentParser(
#         prog="dataqc",
#         description="Command line " + 
#         "interface to dataqc package",
#         epilog="Use `dataqc subcommand -h` for details and options on each command.")
    

#     subparsers = parser.add_subparsers(title="subcommands", 
#         help="one of the subprogram in dataqc",
#         )


#     process = subparsers.add_parser("process",
#         parents=[commons_parser],
#         description="Compute mean amplitude and spectra of seismic data",
#         )
#     process.set_defaults(func=run_processing)
#     process.add_argument("nslc_code", type=str, 
#             help=("station code {network}.{station}.{location}.{channel}," +
#             "may contain glob-style wildcards"))
#     process.add_argument("inventory_routing_type", type=str,
#             help="routing client for inventory",
#             choices=["eida-routing", "iris-federator"]
#             )
#     process.add_argument("sds_root", type=Path,
#             help="root-directory of sds-filesystem")
#     process.add_argument("starttime", type=UTC, 
#             help=("beginning of time range you want to analyze" + 
#                     "Give as YYYY-MM-DDThh:mm:ss"))
#     process.add_argument("endtime", type=UTC, 
#             help="end of time range you want to analyze")

#     process.add_argument("-o", "--outdir", type=Path, 
#             help="where to put the processed data",
#             default=".")
#     process.add_argument("--fileunit", type=str, 
#             help="where to put the processed data",
#             default="year")

#     process.add_argument("--overlap", type=int,
#             help="seconds by which the data is extended beyond time range "+
#                     "to accomodate filter effects etc.",
#             default=base.default_processing_params["overlap"])
#     process.add_argument("--proclen", type=int,
#             help="seconds to process at once, ideally duration of " +
#                     "the data file",
#             default=base.default_processing_params["proclen_seconds"])
#     process.add_argument("--winlen-in-s", type=int,
#             help="time over which amplitude and spectra are computed," +
#             " in seconds",
#             default=base.default_processing_params["winlen_seconds"])
#     process.add_argument("--nperseg", type=int,
#             help="length of segment for spectral estimation "+
#             "(scipy.signal.welch), in samples ",
#             default=base.default_processing_params["nperseg"])
#     process.add_argument("--amplitude-frequencies", type=float, nargs=2,
#             help="min and max frequency of bandpass before "+
#             "amplitude analysis.",
#             default=base.default_processing_params["amplitude_frequencies"] )



#     plot = subparsers.add_parser("plot",
#         parents=[commons_parser],
#         description="Create plots",
#         )
#     plot.set_defaults(func=run_plot)
#     plot.add_argument("nslc_code", type=str, 
#             help=("station code {network}.{station}.{location}.{channel}," +
#                 "May *not* contain wildcards!"))
#     plot.add_argument("datadir", type=Path, 
#             help="where to look for processed data",
#             default=".")
    
#     plot.add_argument("--fileunit", type=str, 
#             help="where to put the processed data",
#             default="year")
#     plot.add_argument("-o", "--figdir", type=Path,
#             help="where to store figures",
#             default=".")
#     plot.add_argument("-s", "--show",
#         action="store_true",
#         help="if given spectrogram plot is opened.")

#     group = plot.add_mutually_exclusive_group()
#     group.add_argument("-l", "--timelist", type=argparse.FileType("r"), 
#             help=("Plot spectrograms using time list." + 
#                     "Can be used as flag to read times from stdin or" + 
#                     "given a file with datetimes."),
#             nargs="?",
#             const=sys.stdin)
#     group.add_argument("-r", "--timerange", type=UTC, nargs=2,
#             help=("start and end of time range you want to analyze"+ 
#                     "Give as YYYY-MM-DDThh:mm:ss, endtime can be None to use current time."),
#                     )



#     avail = subparsers.add_parser("available",
#         aliases=["avail"],
#         parents=[commons_parser],
#         description="Print available HDF5 files and "+ 
#             "covered time range for given code in datadir",
#         )
#     avail.set_defaults(func=run_available)
#     avail.add_argument("nslc_code", type=str, 
#         help=("station code {network}.{station}.{location}.{channel}," +
#             "may contain glob-style wildcards"))
#     avail.add_argument("datadir", type=Path, 
#             help="where to look for processed data",
#             default=".")
#     avail.add_argument("--fileunit", type=str, 
#             help="where to put the processed data",
#             default="year")
    


#     windfilter = subparsers.add_parser("windfilter",
#         aliases=["wind"],
#         #parents=[commons_parser],
#         description="Interpolate and extract wind"
#         )
#     windfilter.set_defaults(func=run_windfilter)
#     windfilter.add_argument("fname", type=Path, 
#         help=("name of wind data file" +
#             "."))
#     windfilter.add_argument("stime", type=UTC, 
#             help="starttime",
#             )
#     windfilter.add_argument("etime", type=UTC, 
#             help="endtime",
#             )
#     windfilter.add_argument("delta", type=float, 
#             help="increment of time axis to which wind data will" +
#                     "be interpolated. In seconds." + 
#                     "Should be same as window length.",
#             )
#     windfilter.add_argument("minspeed", type=float,
#             help="Minimum windspeed")
#     windfilter.add_argument("maxspeed", type=float,
#             help="Maximum windspeed",
#             nargs="?")
#     windfilter.add_argument("out", type=argparse.FileType("w"),
#             nargs="?", default=sys.stdout)



#     # If User enters only 'dataqc' we show help of 
#     # main parser which lists the subprograms
#     if len(sys.argv) < 2:
#         parser.parse_args(["-h"])
    
#     # Otherwise we call the respective subroutine
#     args = parser.parse_args()
#     #print(args)
#     try:
#         module_logger.setLevel(args.loglevel)
#     except AttributeError:
#         module_logger.setLevel("INFO")

#     module_logger.info("CLI-Arguments:\n" + 
#                 "{}".format(args))
#     args.func(args)
    
#     #print('Finish')

if __name__ == "__main__":
    #parse_argument()
    #print(len(sys.argv))
    #main_subparser()
    main()
    