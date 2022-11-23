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
