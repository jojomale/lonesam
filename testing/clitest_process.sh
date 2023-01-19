nscl_code=*.*..BHZ
# overlap=60 #3600
# fmin=4
# fmax=(4, 14)
# nperseg=2048
# winlen_in_s=3600
# proclen=24*3600
sampling_rate=20

outdir='output'
sds_root='../sample_sds'
inventory_routing_type='eida-routing'

startdate="2020-12-20"
enddate="2021-01-10"

logfilename="log/test_processing.log"
loglevel="INFO"

conda activate dataqc
dataqc process -o $outdir \
 $nscl_code \
 $inventory_routing_type $sds_root \
 $startdate $enddate \
    --loglevel $loglevel \
 --logfile $logfilename \
 --sampling-rate $sampling_rate