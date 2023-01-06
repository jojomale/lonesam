conda activate dataqc

figdir="figures/"
datadir=output/
code="GR.BFO..BHZ"
# Time range
##################

# Simplest form: plot full available time range
dataqc plot $code $datadir --figdir $figdir -s


