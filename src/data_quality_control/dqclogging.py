"""
Output messages of the package are implemented using the
`logging <https://docs.python.org/3/library/logging.html#>`_ module.

To manage logging behavior in scripts and notebooks use 
:py:func:`configure_handlers`.


For general information in logging in Python see:

    - https://docs.python.org/3/howto/logging.html#logging-basic-tutorial
    - https://docs.python.org/3/howto/logging-cookbook.html#using-logging-in-multiple-modules

"""

#import os
import logging
import logging.handlers

def create_logger():
    """
    Set logger for the package.

    Should not be used directly.
    """
    # Try to get the package name, may not work for python <3.9 versions
    try: 
        if __package__ is None and __name__ != "__main__":
            loggername = __name__.split('.')[0]
        elif __package__ == "":
            loggername = "dataqc"
        else:
            loggername = __package__
    except UnboundLocalError:
        print("Error, using ", __name__.split('.')[0])
        loggername = __name__.split('.')[0]
    
    # print("LOGGERNAME", loggername)
     # create main logger
    logger = logging.getLogger(loggername)
    logger.setLevel(logging.DEBUG)

    # Set handler for console if no handler is present
    if not logger.hasHandlers():
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)  # set level
        cformatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
             datefmt='%y-%m-%d %H:%M:%S')
        ch.setFormatter(cformatter)
        logger.addHandler(ch)
    return logger


def configure_handlers(loglevel_console, loglevel_file=None, 
                        logfilename=None, use_new_file=False
                        ):
    """
    Configure logging behavior.

    This is the function you want to call in scripts or notebooks.
    It sets logging levels and handlers, i.e. where to the information
    is displayed. By default, standard output (console) is used. Additionally,
    a log file can be defined. Information levels can be different.


    Parameters
    -----------------
    loglevel_console : str {"DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"}
        level of detail for console output
    loglevel_file : None or str {"DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"}
        Level of detail for file logging. If ``None`` logging to file is 
        omitted
    logfilename : str
        Path & name of logfile. If ``None`` logging to file is 
        omitted
    use_new_file : bool [False]
        If ``True`` and logfile exists, it will be overwritten. Otherwise
        new messages are appended.

    
    Note
    --------
    For logging to file both paramters ``loglevel_file`` and ``logfilename``
    must not be ``None``.

    """
    logger = create_logger()
    #print(logger)
    
    # Remove any existing handlers
    for hdl in logger.handlers:
        logger.removeHandler(hdl)

    # Create handlers
    ## console handler
    ch = logging.StreamHandler()
    ch.setLevel(loglevel_console)  # set level
    cformatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%y-%m-%d %H:%M:%S')
    ch.setFormatter(cformatter)
    logger.addHandler(ch)

    ## file handler
    if logfilename and loglevel_file:
        if use_new_file:
            filemode = "w"
        else:
            filemode = "a"
        fh = logging.FileHandler(logfilename, mode=filemode)
        fh.setLevel(loglevel_file)
        hformatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        fh.setFormatter(hformatter)
        logger.addHandler(fh)
        logger.info("Find log file at %s" % logfilename)

