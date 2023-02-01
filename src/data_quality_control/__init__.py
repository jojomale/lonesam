"""
A Python Package to compute and view longterm spectrograms and 
amplitude levels of seismic data.

The package is build from established Python libraries like obspy, 
scipy and h5py. Matplotlib and plotly are used for data visualization. 
The core functionalities are accessible via command line 
interface while the underlying API allows for more customized workflows.

Purpose
===============
Continuous monitoring of data quality is a major issue in 
seismology because the achievement of robust scientific results 
depends on the reliability of the underlying data resources.
This package provides means to perform a systematic analysis of 
noise data in the time and frequency domain.

Two pieces of information are computed:

- the 75%-percentile of the amplitude over a given time window 
  (we use 1 hour)
- the power spectral density (PSD), computed by Welch's method
  over the same time window

The tool is designed to process large amounts of channels 
and years of data.
Raw seismic data are fed into the workflow via obspy-clients.
The processed data are stored in HDF5-Files, 1 file covering
a selected period, e.g. 1 year. Depending on the size of the 
data set, this processing takes minutes to hours. The results 
are stored in rapidly accessible HDF5-files.
They can be visualized using color-coded matrix displays (spectrograms) 
and interactive 3D-figures. The resulting figures give insight 
to characteristic noise patterns at the station and possible 
noise sources, like various forms of anthropogenic noise or 
wind generated noise. Furthermore, changes in noise levels 
or noise patterns are easily detectable.
Furthermore, the processed data can easily be restricted to 
selected times, e.g. to investigate the influence of day/night 
cycles or to obtain wind-speed dependent spectrograms.

Using a median filter, the spectral and amplitude data can be 
smoothed and/or downsampled. This is especially useful and 
recommended when plotting extremly long time series.




Installation
=============
Depending on how you manage your Python, 
you may want to install the 
dependencies first, 
activate your environments, etc.

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
After successful installation the command `dataqc` should be available.

To install dataqc from GitHub:

.. code-block:: console 

  $ pip install git+https://github.com/jojomale/lonesam.git


Editable installation with 
conda including download of repo:

.. code-block:: console

  conda create -n lonesam -c conda-forge \
  pip obspy h5py>=3.3 plotly
  conda activate lonesam

  git clone git+https://github.com/jojomale/lonesam.git
  cd lonesam

  pip install -e . 


If you don't use conda use only the last 3 commands.




Documentation
====================
https://lonesam.readthedocs.io/en/latest/index.html


In order to re-create the documentation from the repository
you need:

- `sphinx <https://www.sphinx-doc.org>`_
- `nbsphinx </https://nbsphinx.readthedocs.io>`_
- `nbsphinx-link <https://nbsphinx-link.readthedocs.io>`_
- `jupyter-lab <https://jupyterlab.readthedocs.io>`_ or `jupyter <https://jupyter-notebook.readthedocs.io>`_
- nbconvert (usually comes along with jupyter)

Assuming your repo sits in directory ``lonesam``, do:

.. code-block:: console

  $ cd lonesam
  $ sphinx-build docs/source/ docs/build/html

Then open ``lonesam/docs/build/html/index.html`` in your browser.


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

.. code-block:: console 

  dataqc [-h] {process,plot,available,avail,windfilter,wind} ...


Use `-h` option on subcommands for details on arguments. E.g.
`dataqc avail -h`.


Compute mean amplitudes and power spectral
densities from raw seismic data and store as 
HDF5-files:

.. code-block:: console 

  dataqc process [-h] [--loglevel {DEBUG,INFO,WARNING,ERROR,CRITICAL}] [--logfile LOGFILE] [--append_logfile] [-o OUTDIR]
                      [--fileunit {year,month,day,hour}] [--overlap OVERLAP] [--proclen PROCLEN] [--winlen-in-s WINLEN_IN_S]
                      [--nperseg NPERSEG] [--amplitude-frequencies AMPLITUDE_FREQUENCIES AMPLITUDE_FREQUENCIES]
                      [--sampling_rate SAMPLING_RATE]
                      nslc_code {eida-routing,iris-federator} sds_root starttime endtime


View available HDF5-data in a directory:
.. code-block:: console 

  dataqc available [-h] [--loglevel {DEBUG,INFO,WARNING,ERROR,CRITICAL}] [--logfile LOGFILE] [--append_logfile] [--fileunit FILEUNIT] nslc_code datadir


Plot results:

.. code-block:: console 

  dataqc plot [-h] [--loglevel {DEBUG,INFO,WARNING,ERROR,CRITICAL}] [--logfile LOGFILE] [--append_logfile] [--fileunit {year,month,day,hour}] [-o FIGDIR] [-s] [-l [TIMELIST]
                  | -r TIMERANGE TIMERANGE]
                  nslc_code datadir


Filter a list of observables for times with specific values.
Observations are interpolated to given time increment (should
be set to the time window used for spectral computations)

.. code-block:: console 

  dataqc windfilter [-h] fname stime etime delta minspeed [maxspeed] [out]





Technicalities
==========================
Time intervals
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


Array dimensions
--------------------
Amplitude
^^^^^^^^^^^^^^
The amplitude information yields 1 data point per time window.
Hence the resulting array is a 1d-array of shape (N,), N being
the total number of time windows in the data range.

For example, with our values of ``winlen_seconds = 3600`` and 
``fileunit="year"``, we
get an amplitude array with shape ``(365*24,)`` or ``(366*24,)``
for leap years.

PSD
^^^^^^^^
In contrast to the amplitude, the PSD results already in a 1d-array.
Its length is determined by the settings for the FFT 
(parameter ``nperseg``), i.e. the number of frequencies (``n_freqs``).
``n_freqs = nperseg // 2 + 1``.
Hence, the output of the PSD-computation has 1 additional dimension, 
thus ``(n_wins, n_freqs)``.

In our case, we used ``nperseg=2048``, which results in 
``n_freqs=1025``, so the shapes are per processing length 
``(24, 1025)`` and per file unit ``'year'`` are ``(365*24, 1025)`` or 
``(366*24, 1025)``.





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