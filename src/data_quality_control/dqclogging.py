
import os
import logging
import logging.handlers

def create_logger():
    """
    Manage logging behavior

    Warning
    ---------
    Filenames (and paths probably too) must not contain "." otherwise
    the automatic deletion of backup files does not work

    Notes
    -----------
    - See https://docs.python.org/3/howto/logging.html#logging-basic-tutorial
    - https://docs.python.org/3/howto/logging-cookbook.html#using-logging-in-multiple-modules
    - TimedRotatingFileHandler splits the path and filename base at "." and then
        uses regular expressions
     
    """
    # print("eida_logger NAME", __name__)
    # print("eida_logger PCKG", __package__)
    #print("MODULE", __name__.__module__)

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
        cformatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                                    datefmt='%y-%m-%d %H:%M:%S')
        ch.setFormatter(cformatter)
        logger.addHandler(ch)
    return logger


def configure_handlers(logger, loglevel_console, loglevel_file, 
                        logfilename, use_new_file=False
                      #  log_timeunit, log_backupcount, log_interval
                        ):
    
    if use_new_file:
        filemode = "w"
    else:
        filemode = "a"
    
    # Remove any existing handlers
    for hdl in logger.handlers:
        logger.removeHandler(hdl)

    # Create handlers
    ## console handler
    ch = logging.StreamHandler()
    ch.setLevel(loglevel_console)  # set level

    ## file handler
    #fh = logging.handlers.TimedRotatingFileHandler(
    #        os.path.join(logfilename, 'eida_availability_log'), 
    #        when=log_timeunit, backupCount=log_backupcount, interval=log_interval)
    fh = logging.FileHandler(logfilename, mode=filemode)
    fh.setLevel(loglevel_file)
    
    ## create formatter
    cformatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                                    datefmt='%y-%m-%d %H:%M:%S')
    hformatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    ### add formatter to ch
    ch.setFormatter(cformatter)
    fh.setFormatter(hformatter)

    # add handlers to main logger, if it doesn't has some already.
    # This happens if different modules are used by a script because
    # each one calls create_logger once.
    # We don't need to add the handlers again, messages are propagated
    # up to main logger. Otherwise messages are duplicated.
    #if not logger.hasHandlers():
    logger.addHandler(ch)
    logger.addHandler(fh)

    
    logger.info("Find log file at %s" % logfilename)
        #os.path.join(logfilename, 'dataqc.log'))