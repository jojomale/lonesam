"""
Prepare seismic data for data quality control.

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

Amplitude
^^^^^^^^^^^^^^
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

PSD
^^^^^^^^
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
- h5py
- plotly

We also use 
- numpy
- scipy
but they are "included" in obspy.


Install
--------
Depending on how you manage your Python, you may want to install the 
dependencies first.

To install dataqc from BGR's SVN-server:

    ```bash
    $ pip install svn+svn://svn.hannover.bgr.de/station_quality_control/trunk/data_quality_control#egg=dataqc
    ```

From source distribution:

    ```bash
    $ pip install dataqc-1.0.0.tar.gz
    ```

Usage
============
Process raw data from commandline:

```bash
dataqc -n GR -s BFO -c BHZ --sds-root sample_sds/ -v INFO 2020-360 2021-005 tmp_data/
```
of from script in scripts/run_processing.py

View processed data using scripts/analysis.py


Documentation
=================
HTML-Documentation in SVN-repository:

station_quality_control/trunk/data_quality_control/docs/build/html/index.html


Known issues during processing
==================================
Change of sampling rate
-------------------------
The sampling rate of a station may change over the years. 
If that happens the mseed-file in the
SDS contains multiple traces with different sampling rates. 
These can not be merged without
resampling to a common rate. We choose to resample the stream to the highest rate. 
This happens in the routing `util.process_streams()` which calls 
`util.merge_different_samplingrates`.

"""