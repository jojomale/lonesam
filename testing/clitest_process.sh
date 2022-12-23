nscl_code=*.*..HHZ
overlap=60 #3600
fmin=4
fmax=(4, 14)
nperseg=2048
winlen_in_s=3600
proclen=24*3600
sampling_rate=100

outdir='output'
sds_root='../sample_sds'
inventory_routing_type='eida-routing'


startdate="2020-12-20"
enddate="2021-01-10"

logfilename="log/test_processing.log"

conda activate dataqc
dataqc process $nscl_code \
 $inventory_routing_type $sds_root \
 $startdate $enddate