README

Purpose
===============
This package provides routines and workflows zu extract
amplitude and spectral information from seismic data
on a regular basis for quality control.

Two pieces of information are computed:
- the 75%-percentile of the amplitude over a given time window 
  (we use 1 hour)
- the power spectral density (PSD), computed by Welch's method
  over the same time window

Raw seismic data are fed into the workflow via obspy-clients.
The processed data are stored in HDF5-Files, 1 file covering
a selected period, e.g. 1 year.


Relevant time intervals
--------------------------
The module *base* provides a low-level interface for processing
the data.
We distinguish **three relevant time intervals**:
- The time window (parameter ``winlen``) over which amplitude 
  and psd are computed. This is the smallest window. We use 1 hour.
- The processing length (parameter ``proclen``) which is the
  time range that is processed at once. We use 1 day because files
  in our database cover 1 day. This is the amount of data
  that is loaded at once for computations, thus should fit into RAM
- The time period that is written into a single output file
  (`fileunit``). We use 1 year.

These three intervals determine the shape of the resulting
data arrays.

### Amplitude

The amplitude information yields 1 data point per time window.
For 1 processing length, this results in a 1d-array with length
= ``n_winlen``, thus the number of time windows per processing
length.

In the output file, data from each processing length are stacked
vertically, resulting in a **2d-array** with 
shape (``n_proclen``, ``n_winlen``). ``n_proclen`` is the number
of processed time ranges per period covered in 1 file.

For example, with our values of ``winlen_seconds = 3600``,
``proclen_seconds = 24*3600`` and ``fileunit="year"``, we
get an amplitude arrays with shape ``(365, 24)`` or ``(366, 24)``
for leap years.

### PSD

In contrast to the amplitude, the PSD results already in a 1d-array.
Its length is determined by the settings for the FFT 
(parameter ``nperseg``), i.e. the number of frequencies (``n_freqs``).
Hence, the output of the PSD-computation has 1 additional dimension.
So per processing length, shape is ``(n_winlen, n_freqs)`` and per
file unit, shape is ``(n_proclen``, ``n_winlen, n_freqs)``.

In our case, we used ``nperseg=2048``, which results in 
``n_freqs=1025``, so the shapes are per processing length 
``(24, 1025)`` and per file unit ``(365, 24, 1025)`` or 
``(366, 24, 1025)``.




Installation
=====================
Dependencies
-------------
- obspy
- h5py >= 3.3
- plotly

We also use 
- numpy
- scipy
but they are "included" in obspy.


Install
--------
Depending on how you manage your Python, 
you may want to install the 
dependencies first, 
activate your environments, etc.

To install dataqc from BGR's SVN-server:

    ```bash
    $ pip install svn+svn://svn.hannover.bgr.de/station_quality_control/trunk/data_quality_control#egg=dataqc
    ```

This can take some time because the repo is 
quite large. If you only want the package, 
you can download an archived version:

 ```bash 
  # Make directory for source archive
  mkdir dataqc
  cd dataqc
  
  # Download package as archive
  svn export svn://svn.hannover.bgr.de/station_quality_control/trunk/data_quality_control/dist/dataqc-1.0.0.tar.gz

  # Install 
  pip install dataqc-1.0.0.tar.gz
  ```


To obtain the repo including all sample data,
scripts and tests use:

  ```bash
    $ mkdir my_repos
    $ cd my_repos
    $ svn checkout svn://svn.hannover.bgr.de/station_quality_control/trunk/data_quality_control#egg=dataqc
    $ cd station_quality_control/trunk/data_quality_control
  ```

Then use either for simple installation:
  
  ```bash 
    $ pip install .
  ```

or for an editable installation:

  ```bash 
    $ pip install -e .
  ```


Suggested editable installation with 
conda including download of repo:

  ```bash
  conda create -n dataqc -c conda-forge \
  pip obspy ipykernel h5py>=3.3 plotly
  conda activate dataqc

  mkdir my_svn_repos
  cd my_svn_repos
  svn checkout svn://svn.hannover.bgr.de/station_quality_control 
  cd station_quality_control/trunk/data_quality_control
  pip install -e . 
  ```



Documentation
=================
To view the HTML-Documentation open thi file 
(in SVN-repository) in a browser:

station_quality_control/trunk/data_quality_control/docs/build/html/index.html





Usage
============
Core functionalities can be a accessed via CLI.
For more flexiblity, we recommend using the
API to build customized processing and plotting
scipts. 


API
--------
For starting points, take a look at:
- station_quality_control/trunk/data_quality_control/scripts
- station_quality_control/trunk/data_quality_control/notebooks/usage_demo.ipynb


CLI
-----------
General use:

  ```bash
  dataqc [-h] {process,plot,available,avail,windfilter,wind} ...
  ```

Use `-h` option on subcommands for details on arguments. E.g.
`dataqc avail -h`.


Compute mean amplitudes and power spectral
densities from raw seismic data and store as 
HDF5-files:

  ```bash
  dataqc process [-h] [--loglevel {DEBUG,INFO,WARNING,ERROR,CRITICAL}] [--logfile LOGFILE] [--append_logfile] [-o OUTDIR]
                        [--fileunit {year,month,day,hour}] [--overlap OVERLAP] [--proclen PROCLEN] [--winlen-in-s WINLEN_IN_S]
                        [--nperseg NPERSEG] [--amplitude-frequencies AMPLITUDE_FREQUENCIES AMPLITUDE_FREQUENCIES]
                        [--sampling_rate SAMPLING_RATE]
                        nslc_code {eida-routing,iris-federator} sds_root starttime endtime
  ```

View available HDF5-data in a directory:
  ```bash
  dataqc available [-h] [--loglevel {DEBUG,INFO,WARNING,ERROR,CRITICAL}] [--logfile LOGFILE] [--append_logfile] [--fileunit FILEUNIT] nslc_code datadir
  ```

Plot results:
  ```bash
  dataqc plot [-h] [--loglevel {DEBUG,INFO,WARNING,ERROR,CRITICAL}] [--logfile LOGFILE] [--append_logfile] [--fileunit {year,month,day,hour}] [-o FIGDIR] [-s] [-l [TIMELIST]
                   | -r TIMERANGE TIMERANGE]
                   nslc_code datadir
  ```

Filter a list of observables for times with specific values.
Observations are interpolated to given time increment (should
be set to the time window used for spectral computations)
```bash
dataqc windfilter [-h] fname stime etime delta minspeed [maxspeed] [out]
```