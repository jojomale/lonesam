How to make svn ignore jupyter checkpoints:
```bash
svn propset svn:ignore "*.pyc*
*.ipynb_checkpoints*
*.ipynb_checkpoints/*
*checkpoint*" notebooks --recursive
```

Why install did not work properly at first:
in setup.py I used at first the same settings as in eidaqc:
`packages = find_packages(where="src")`
In dataqc however pip did not find the source files which is
why after installation in other environments the modules could
not be found.
Instead, this worked:

```python
package_dir={'': 'src'}, 
packages=["data_quality_control"],  #
```

I suspect, find_packages did not work because the directory
with the source files `src/data_quality_control` has a different
name than the name of the package in the setup (`dataqc`).


# Rediscovery
- module `processing` is deprecated
- possibly, the `Analyzer` class from module `analysis` should go to `base` because
    because it is parent class to `sds_db.SDSDataBaseAnalyzer`. Then `analysis can be deprecated as well
- Analyzer in sds_db can extract time ranges, including those passing midnight


# Development

- reading psds from timelist: interpolation of times should occur outside of Analyzer class
    Usecase: Winddata, available at 6h-intervals. We want to extract times with certain ranges
    of velocities, e.g. all datetimes when 2 < v < 3 m/s.
    PSDs are available in 1h intervals. So we should interpolate the wind time series first to
    the higher interval, *then* select the times from the interpolated speeds.

    The selection of the wind data has nothing to do with the properties of the PSDs except that
    we need to know the target interval `winlen_seconds`. But how we select the times from the 
    wind data is specific to that data set. So for now, I suggest, that the Analyzer only needs 
    to receive the final timelist. And we can provide a separate functionality that inrpolates the
    wind (or whatever) data to the right time grid and then allows customized selection. Maybe we
    could create a basic function, which can take the Analyzer as input and a tuple of times and 
    values to interpolate. But extracting time and value from the raw data most likely needs a
    customized function because you never know how tims are provided, which columns to use etc.
    For selection, we also might require an intermediate, manual step, although we could probably
    build something standardized for simpler selection rules.

- Treat error due to merging of different sampling rates (GR.BFO..HHZ 30-Nov-2011)


# CLI
- dataqc
    - process
        Options:
            - code
            - starttime, endtime
            - processing paramters
            - inventory_routing_type
            - sds_root
    - available
        Options:
            - code
    - plot 
        needs to be able to differentiate between timelist and time 
        range as input

        Options:
            - code
            - starttime, endtime or timelist
            - -a, --a :plot all
            - -s, --spectrogram
            - --amplitude
            - --spectrogram3d

        Calls:
            - dataqc plot code 2020-12-01
            - dataqc plot code timelist.txt
            - dataqc plot code 2020-12-01 2020-12-02 2020-12-03 ...
            - dataqc windfilter wind.txt 2020-12-01 2021-01-01 3600 0 5| dataqc plot code 

# To dos
- Logger managment
    - Make logfile optional
    - Make logfile destination optional
    - Set logger from commandline

- pass timelists or time range to Analyzer.get_data() via command line: Would be cool if 
the CLI would determine on its own whether the user gave 
    - 1 starttime, assuming time range until today
    - 1 starttime, 1 endtime, defining a time range
    - a file with a time list
    - take time list from pipe
However mutually exclusive mandatory arguments are not exactly planned in argparse.
As of now, I think, we would need to build a relatively complicated if-else and try-except
sequence to distinguish all the above cases.
Maybe it is easier to pass them all as optional arguments and if no optional ars are given,
we plot the entire available time range
