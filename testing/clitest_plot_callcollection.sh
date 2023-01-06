conda activate dataqc

figdir="figures/"
datadir=output/
code="GR.BFO..BHZ"
# Time range
##################

# Simplest form: plot full available time range
dataqc plot $code $datadir --figdir $figdir

# Plot given time range
dataqc plot $code $datadir --figdir $figdir \
-r 2020-12-28 2021-01-05 \


# Try to plot time range outside available rage
dataqc plot $code $datadir --figdir $figdir \
-r 2021-02-01 2021-03-01 \


# Logging managment #
#####################

dataqc plot $code $datadir \
--figdir $figdir --loglevel INFO --logfile log/happy.log


# Time list #
#############

# Make time list from winddata
dataqc windfilter winddata_bfo.txt 2020-12-10 2021-01-11 \
3600 0 2 wind_0-2.txt

# Plot from file
dataqc plot $code $datadir --figdir $figdir \
-l wind_0-2.txt 

# Pipe list 
dataqc windfilter winddata_bfo.txt 2020-12-10 2021-01-11 3600 0 2 | \
dataqc plot $code $ --figdir $figdir -l


# Expects timelist from stdin but doesn't receive one. 
# Must be terminated by Ctrl+C
dataqc plot $code $datadir --figdir $figdir -l 


