#!/usr/bin/env python
# coding: utf-8

# In[1]:
from importlib import reload
import os
from datetime import timedelta

import numpy as np

import h5py

from obspy.clients.filesystem.sds import Client
from obspy.clients.fdsn import RoutingClient
from obspy.core import UTCDateTime as UTC
from obspy.signal import util

import matplotlib.pyplot as plt
plt.style.use('tableau-colorblind10')

from data_quality_control import processing
from data_quality_control import base


# In[4]:


overlap = 60 #3600
fmin, fmax = (4, 14)
nperseg = 2048
winlen_in_s = 3600
proclen = 24*3600


# In[8]:


reload(base)
pp = base.ProcessingParameters(overlap=overlap, 
                             amplitude_frequencies=(fmin, fmax),
                             nperseg=nperseg,
                             winlen_seconds=winlen_in_s,
                             proclen_seconds=proclen)
print(pp)


# In[9]:


network = 'GR'
station = 'BFO'
location = ''
channel = 'BHZ'

outdir = '.'

sds_root = os.path.abspath('../sample_sds/')
inventory_routing_type = "eida-routing"


# In[10]:


sdsclient = Client(sds_root)
invclient = RoutingClient(inventory_routing_type)

sdsclient.get_all_nslc()


# In[33]:


reload(base)
proc = base.GenericProcessor(network, station, location, channel, 
                             sdsclient, invclient, fileunit="year",
                             procparams=pp)

proc


# In[34]:


startdate = UTC("2020-12-20")
#enddate = UTC("2020-12-31")
enddate = UTC("2021-01-15")


# In[35]:


proc.process(startdate, enddate)


# In[ ]:


for s, e in base.TIME_ITERATORS["year"](startdate, enddate):
    print(s, e)


# In[ ]:



