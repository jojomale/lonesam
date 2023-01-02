
envname="tdi2" #test_dataqc_install

conda create -n $envname

source activate $envname
conda install -c conda-forge pip obspy ipykernel h5py<=3.3 plotly

pip install svn+svn://svn.hannover.bgr.de/station_quality_control/trunk/data_quality_control#egg=dataqc
# pip install ~/svn/station_quality_control/trunk/data_quality_control/dist/dataqc-1.0.0.tar.gz


## Download repo, install from there in editable mode
# mkdir my_svn_repos
# cd my_svn_repos
# svn checkout svn://svn.hannover.bgr.de/station_quality_control 
# pip install -e . station_quality_control/trunk/data_quality_control