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


import argparse
import pathlib
import time
from datetime import timedelta
from obspy.core import UTCDateTime as UTC

from . import processing, base, sds_db

# default_settings = dict(network = '*',
#     station = '*',
#     channel = '*',
#     overlap = 60,
#     amplitude_frequencies = (4,14),
#     nperseg = 2048,
#     winlen_in_s = 3600,
#     proclen = 24*3600,
#     sds_root = os.path.abspath('.'),
#     inventory_routing_type = "eida-routing",
#     sds_client_dict = {})

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
    parser.add_argument("sds_root", type=pathlib.Path,
            help="root-directory of sds-filesystem")
    
    parser.add_argument("starttime", type=UTC, 
            help="beginning of time range you want to analyze")
    parser.add_argument("endtime", type=UTC, 
            help="end of time range you want to analyze")
    
    
    parser.add_argument("-o", "--outdir", type=pathlib.Path, 
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
    parser.add_argument("--configfile", type=pathlib.Path,
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
    processing.logger.setLevel(loglevel)
    #print(args)
    return arg_dict, proc_args


def run_raw_processing(args1, args2):
    
    processor = sds_db.SDSProcessor(
            **args1
            )
    #processor.print()
    processor.process(**args2)



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
    parser.add_argument("outdir", type=pathlib.Path, 
            help="where to put the processed data")

    parser.add_argument("-n", "--network", type=str, 
            help="network(s), may contain glob-style wildcards",
            default=processing.default_settings["network"])
    parser.add_argument("-s", "--station", type=str, 
            help="station(s), may contain glob-style wildcards",
            default=processing.default_settings["station"])
    parser.add_argument("-c", "--channel", type=str, 
            help="channel(s), may contain glob-style wildcards",
            default=processing.default_settings["channel"])
    parser.add_argument("--sds-root", type=pathlib.Path,
            help="root-directory of sds-filesystem",
            default=processing.default_settings["sds_root"])
    parser.add_argument("--inventory-routing-type", type=str,
            help="routing client for inventory",
            choices=["eida-routing", "iris-federator"],
            default=processing.default_settings["inventory_routing_type"])
    parser.add_argument("--overlap", type=int,
            help="seconds by which the data is extended beyond time range "+
                    "to accomodate filter effects etc.",
            default=processing.default_settings["overlap"])
    parser.add_argument("--proclen", type=int,
            help="seconds to process at once, ideally duration of " +
                    "the data file",
            default=processing.default_settings["proclen"])
    parser.add_argument("--winlen-in-s", type=int,
            help="time over which amplitude and spectra are computed," +
            " in seconds",
            default=processing.default_settings["winlen_in_s"])
    parser.add_argument("--nperseg", type=int,
            help="length of segment for spectral estimation "+
            "(scipy.signal.welch), in samples ",
            default=processing.default_settings["nperseg"])
    parser.add_argument("--amplitude-frequencies", type=float, nargs=2,
            help="min and max frequency of bandpass before "+
            "amplitude analysis.",
            default=processing.default_settings["amplitude_frequencies"] )
    parser.add_argument("--configfile", type=pathlib.Path,
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
    processing.logger.setLevel(loglevel)
    #print(args)
    return arg_dict, proc_args


def run_raw_processing_old(args1, args2):
    
    processor = processing.RawDataProcessor(
            **args1
            )
    #processor.print()
    processor.process(**args2)
    raise DeprecationWarning("This is an old part of the code. " + 
        "Try to switch to module `base`")


def main():
    t = time.time()
    args1, args2 = parse_argument()
    run_raw_processing(args1, args2)
    runtime = timedelta(seconds=time.time()-t) 
    print("Finished. Took {} h".format(runtime))
    #print(args)
    


if __name__ == "__main__":
    #parse_argument()
    
    main()
    