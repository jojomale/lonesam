
import os
import logging
import logging.handlers

def create_logger():
    """
    Manage logging behavior


    Notes
    -----------
    - See https://docs.python.org/3/howto/logging.html#logging-basic-tutorial
    - https://docs.python.org/3/howto/logging-cookbook.html#using-logging-in-multiple-modules
    - TimedRotatingFileHandler splits the path and filename base at "." and then
        uses regular expressions
     
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



# def create_console_handler(loglevel):
#     pass




# def create_file_handler(loglevel, logfilename, 
#                         use_new_file=False):
    
#     if not logfilename:
#         logfilename = "dataqc.log"
#     ## file handler
#     if use_new_file:
#         filemode = "w"
#     else:
#         filemode = "a"
#     fh = logging.FileHandler(logfilename, mode=filemode)
#     fh.setLevel(loglevel)
#     hformatter = logging.Formatter(
#         '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
#     fh.setFormatter(hformatter)
#     return fh