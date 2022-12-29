conda activate dataqc

# Time range
##################

# Simplest form: plot full available time range
dataqc plot GR.BFO..HHZ ../sample_output/run_processing/ \
--figdir cli_figures

# Plot given time range
dataqc plot GR.BFO..HHZ ../sample_output/run_processing/ \
-r 2020-12-28 2021-01-05 \
--figdir cli_figures

# Try to plot time range outside available rage
dataqc plot GR.BFO..HHZ ../sample_output/run_processing/ \
-r 2021-02-01 2021-03-01 \
--figdir cli_figures


# Time list #
#############

# Make time list from winddata
dataqc windfilter winddata_bfo.txt 2020-12-10 2021-01-11 \
3600 0 2 wind_0-2.txt

# Plot from file
dataqc plot -l wind_0-2.txt GR.BFO..HHZ ../sample_output/run_processing/\
--figdir cli_figures

# Pipe list 
dataqc windfilter winddata_bfo.txt 2020-12-10 2021-01-11 3600 0 2 | \
dataqc plot GR.BFO..HHZ ../sample_output/run_processing/ -l \
--figdir cli_figures

# Expects timelist from stdin but doesn't receive one. 
# Must be terminated by Ctrl+C
dataqc plot GR.BFO..HHZ ../sample_output/run_processing/ -l \
--figdir cli_figures


# Logging managment #
#####################

dataqc plot GR.BFO..HHZ ../sample_output/run_processing/ \
--figdir cli_figures --verbosity INFO --logfile happy.log
