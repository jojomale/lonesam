datadir=output/
nslc_code="GR.BFO..BHZ"
outdir='output/smoothed'
kernelsize="6"
kernelshift="3"

logfilename="log/test_processing.log"
loglevel="INFO"

conda activate dataqc
dataqc smooth $nslc_code $datadir $outdir \
   ${kernelsize} ${kernelshift}\
   -f \
   --loglevel $loglevel \
   --logfile $logfilename

dataqc plot-amplitudes $nslc_code $outdir --figdir $figdir -s

