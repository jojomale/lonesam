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