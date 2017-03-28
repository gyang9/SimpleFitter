#!/bin/bash

# The following line ensures that any failed command causes the script to exit. 
set -e

# Pass in: (1) a' min., (2) a' step size, (3) number of a' pts, (4) b' min., 
# (5) b' step size, (6) number of b' pts, (7) c' min., (8) c' step size, 
# (9) number of c' pts, (10) s2t13 min., (11) s2t13 step size, 
# (12) number of s2t13 pts, (13) job name prefix, (14) path to executable,
# (15) configuration file  

nargs=15
if [ "$#" -ne $nargs ]; then
    echo "This script expects $nargs arguments. Exiting!"
    exit
fi

amin=$1
astep=$2
apts=$3
bmin=$4
bstep=$5
bpts=$6
cmin=$7
cstep=$8
cpts=$9
s2t13min=${10}
s2t13step=${11}
s2t13pts=${12}
jobstr=${13}
execpath=${14}
configfile=${15}

if [ ! -f "$execpath" ]; then
    echo "No file exists at $execpath. Exiting!"
    exit
fi

# CHANGE THIS FOR EACH SET OF JOBS!
logdir="/sps/dchooz/guang/binningPaper/${jobstr}/"

if [ ! -d "$logdir" ]; then
    echo "Log directory $logdir does not exist. Exiting!"
    exit
fi

for (( si=0; si<$s2t13pts; si++ )); do
    s2t13=`echo $s2t13min + $si\*$s2t13step | bc`
    echo "s2t13: $s2t13"
    for (( ai=0; ai<$apts; ai++ )); do
	aval=`echo $amin + $astep\*$ai | bc`
	
	for (( bi=0; bi<$bpts; bi++ )); do
	    bval=`echo $bmin + $bstep\*$bi | bc`

	    jobname="${jobstr}_s${si}_a${ai}_b${bi}"
	    qsub -P P_dchooz -o "$logdir/${jobname}.txt" -j y -l sps=1,fsize=4G,vmem=4G -N $jobname submit_grid_job.sh $configfile $s2t13 $aval $bval $cmin $cstep $cpts $jobstr $execpath
 	    
	done
    done
done
