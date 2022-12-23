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
from obspy.core import UTCDateTime as UTC

from . import base, sds_db, dqclogging

logger = dqclogging.create_logger()
module_logger = dqclogging.logging.getLogger(logger.name+'.base')

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

def parse_argument():
    """
    Get arguments from commandline and prepare for
    use in processing.RawDataProcessor.
    """
    parser = argparse.ArgumentParser(description="Compute "+
        "amplitude levels and power spectral densities from "+
        "raw seismic data in an sds-filesystem and store in "+
        "HDF5-files.")
    


    parser.add_argument("network", type=str, 
            help="network(s), may contain glob-style wildcards")
    parser.add_argument("station", type=str, 
            help="station(s), may contain glob-style wildcards")
    parser.add_argument("channel", type=str, 
            help="channel(s), may contain glob-style wildcards")
    parser.add_argument("inventory_routing_type", type=str,
            help="routing client for inventory",
            choices=["eida-routing", "iris-federator"])
    parser.add_argument("sds_root", type=Path,
            help="root-directory of sds-filesystem")
    
    parser.add_argument("starttime", type=UTC, 
            help="beginning of time range you want to analyze")
    parser.add_argument("endtime", type=UTC, 
            help="end of time range you want to analyze")
    
    
    parser.add_argument("-o", "--outdir", type=Path, 
            help="where to put the processed data",
            default=".")
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
    parser.add_argument("--configfile", type=Path,
            help="file with parameters. additionally given parameters "+
            "override those from file")
    parser.add_argument("-f", "--force-new-file",
            action="store_true",
            help="overrides existing files if given",)
    

    parser.add_argument("-v", "--verbosity", type=str,
                    choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                    help="set logging level",
                    default="DEBUG")

    print(parser.parse_args())
    arg_dict= vars(parser.parse_args())
    proc_args = {k: arg_dict.pop(k) for k in ['starttime', 'endtime']}
    loglevel = arg_dict.pop("verbosity")
    module_logger.setLevel(loglevel)
    #print(args)
    return arg_dict, proc_args


def run_processing(args1, args2):
    t = time.time()
    processor = sds_db.SDSProcessor(
            **args1
            )
    #processor.print()
    #processor.process(**args2)
    try:
        processor.process(**args2)
    except Exception as e:
        processor.logger.exception(e)

    processor.logger.info("Processing finished. Use `run_analysis.py` to view results.")

    runtime = timedelta(seconds=time.time()-t) 
    print("Finished. Took {} h".format(runtime))



def parse_argument_old():
    """
    Get arguments from commandline and prepare for
    use in processing.RawDataProcessor.
    """
    parser = argparse.ArgumentParser(description="Compute "+
        "amplitude levels and power spectral densities from "+
        "raw seismic data in an sds-filesystem and store in "+
        "HDF5-files.")
    
    parser.add_argument("startdate", type=UTC, 
            help="beginning of time range you want to analyze")
    parser.add_argument("enddate", type=UTC, 
            help="end of time range you want to analyze")
    parser.add_argument("outdir", type=Path, 
            help="where to put the processed data")

    parser.add_argument("-n", "--network", type=str, 
            help="network(s), may contain glob-style wildcards",
            default=default_settings["network"])
    parser.add_argument("-s", "--station", type=str, 
            help="station(s), may contain glob-style wildcards",
            default=default_settings["station"])
    parser.add_argument("-c", "--channel", type=str, 
            help="channel(s), may contain glob-style wildcards",
            default=default_settings["channel"])
    parser.add_argument("--sds-root", type=Path,
            help="root-directory of sds-filesystem",
            default=default_settings["sds_root"])
    parser.add_argument("--inventory-routing-type", type=str,
            help="routing client for inventory",
            choices=["eida-routing", "iris-federator"],
            default=default_settings["inventory_routing_type"])
    parser.add_argument("--overlap", type=int,
            help="seconds by which the data is extended beyond time range "+
                    "to accomodate filter effects etc.",
            default=default_settings["overlap"])
    parser.add_argument("--proclen", type=int,
            help="seconds to process at once, ideally duration of " +
                    "the data file",
            default=default_settings["proclen"])
    parser.add_argument("--winlen-in-s", type=int,
            help="time over which amplitude and spectra are computed," +
            " in seconds",
            default=default_settings["winlen_in_s"])
    parser.add_argument("--nperseg", type=int,
            help="length of segment for spectral estimation "+
            "(scipy.signal.welch), in samples ",
            default=default_settings["nperseg"])
    parser.add_argument("--amplitude-frequencies", type=float, nargs=2,
            help="min and max frequency of bandpass before "+
            "amplitude analysis.",
            default=default_settings["amplitude_frequencies"] )
    parser.add_argument("--configfile", type=Path,
            help="file with parameters. additionally given parameters "+
            "override those from file")

    parser.add_argument("-v", "--verbosity", type=str,
                    choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                    help="set logging level",
                    default="DEBUG")

    print(parser.parse_args())
    arg_dict= vars(parser.parse_args())
    proc_args = {k: arg_dict.pop(k) for k in ['startdate', 'enddate', "outdir"]}
    loglevel = arg_dict.pop("verbosity")
    module_logger.setLevel(loglevel)
    #print(args)
    return arg_dict, proc_args


def run_raw_processing_old(args1, args2):
    
    processor = base.RawDataProcessor(
            **args1
            )
    #processor.print()
    processor.process(**args2)
    raise DeprecationWarning("This is an old part of the code. " + 
        "Try to switch to module `base`")


def main():
    
    # Main parser
    parser = argparse.ArgumentParser(description="Command line " + 
        "interface to dataqc package",
        epilog="Use `dataqc subcommand -h` for details and options on each command.")
    subparsers = parser.add_subparsers(title="subcommands", 
        help="one of the subprogram in dataqc",
        )

    # Define arguments for subparsers
    process = subparsers.add_parser("process_sds", 
        description="Compute avg amplitude and power spectral densities and store in HDF5",
        help="Compute files",
            aliases=["processing", "process"]
        )
#     process.add_argument("bla", type=str, 
#             help="network(s), may contain glob-style wildcards")
    
#     process.add_argument("network", type=str, 
#             help="network(s), may contain glob-style wildcards")
#     process.add_argument("station", type=str, 
#             help="station(s), may contain glob-style wildcards")
#     process.add_argument("channel", type=str, 
#             help="channel(s), may contain glob-style wildcards")
    
    process.add_argument("nslc_code", type=str, 
            help=("station code {network}.{station}.{location}.{channel}," +
                "may contain glob-style wildcards"))
    
    process.add_argument("inventory_routing_type", type=str,
            help="routing client for inventory",
            choices=["eida-routing", "iris-federator"])
    process.add_argument("sds_root", type=Path,
            help="root-directory of sds-filesystem")

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
#     process.add_argument("--configfile", type=Path,
#             help="file with parameters. additionally given parameters "+
#             "override those from file")
    
    
    process.add_argument("starttime", type=UTC, 
            help="beginning of time range you want to analyze")
    process.add_argument("endtime", type=UTC, 
            help="end of time range you want to analyze")
    process.add_argument("-f", "--force-new-file",
            action="store_true",
            help="overrides existing files if given",)
    process.set_defaults(func=run_processing)
    

    parser.add_argument("-v", "--verbosity", type=str,
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="set logging level",
        default="INFO")

    args = parser.parse_args()

    print(parser.parse_args())
    arg_dict= vars(parser.parse_args())
    print(arg_dict)
    proc_args = {k: arg_dict.pop(k) for k in ['starttime', 'endtime']}
    #loglevel = arg_dict.pop("verbosity")



    # If User enters only 'eida' we show help of 
    # main parser which lists the subprograms
    if len(sys.argv) < 2:
        parser.parse_args(["-h"])
    
    # Otherwise we call the respective subroutine
    args.func(arg_dict, proc_args)
    
    print('Finish')

if __name__ == "__main__":
    #parse_argument()

    main()
    