#!/bin/bash

# Arguments passed in: (1) configfile (2) s2t13, (3) a', (4) b', (5) c' minimum, (6) c' step
# size, (7) number of c' pts, (8) outfile prefix, (9) path to executable.

# The following line ensures that any failed command causes the script to exit.
set -e

nargs=9
if [ "$#" -ne $nargs ]; then
    echo "This script expects $nargs arguments. Exiting!"
    exit
fi

if [ -z $FFITS_DATA ]; then
    echo 'FFITS_DATA environment variable not set. Exiting!'
    exit
fi

configfile=$1
s2t13=$2
aval=$3
bval=$4
cmin=$5
cstep=$6
cpts=$7
prefix=$8
execpath=$9

if [ ! -f "$execpath" ]; then
    echo "No file exists at $execpath. Exiting!"
    exit
fi

# Copy input files to worker machine. 
rm -rf HIII thirdPub twoDetector
cp -r /sps/dchooz/guang/binningPaper/*hh .
cp -r /sps/dchooz/guang/binningPaper/*cc .
cp -r /sps/dchooz/guang/binningPaper/*C .
cp -r /sps/dchooz/guang/binningPaper/*txt .
cp -r /sps/dchooz/guang/binningPaper/*root .

# Ensure input files are read from worker machine.
FFITS_DATA=`pwd`

for (( ci=0; ci<$cpts; ci++ )); do
    cval=`echo $cmin + $cstep\*$ci | bc`
# Construct output file name from inputs to executable. Negative numbers
# have their negative stripped and replaced with the letter 'n.'
    if [ `echo "$aval >= 0" | bc` -eq 1 ]; then
	if [ `echo "$aval == 0" | bc` -eq 1 ]; then
	    astring='0'
	else
	    astring=$aval
	fi
    else
	astring='n'${aval:1}
    fi
    if [ `echo "$bval >= 0" | bc` -eq 1 ]; then
	if [ `echo "$bval == 0" | bc` -eq 1 ]; then
	    bstring='0'
	else
	    bstring=$bval
	fi
    else
	bstring='n'${bval:1}
    fi
    if [ `echo "$cval >= 0" | bc` -eq 1 ]; then
	if [ `echo "$cval == 0" | bc` -eq 1 ]; then
	    cstring='0'
	else
	    cstring=$cval
	fi
    else
	cstring='n'${cval:1}
    fi
    
	outfilename=$prefix'_s'$s2t13'_'$astring'_'$bstring'_'$cstring'.txt'
# or, for a grid with a Li dimension, use something like:
#	outfilename=$prefix'_li'$s2t13'_'$astring'_'$bstring'_'$cstring'.txt'
    
	# Call executable
    eval $execpath $configfile $s2t13 $aval $bval $cval $outfilename

	# You added a public key from IN2P3 to your Nevis authorized_keys, so 
	# you can now scp to Nevis without a password.
  #  if ! scp $outfilename houston.nevis.columbia.edu:/a/data/riverside/rcarr/DC3grid; then
	    # Try to copy to Nevis twice before writing locally (scp sometimes fails).
    #	if ! scp $outfilename houston.nevis.columbia.edu:/a/data/riverside/rcarr/DC3grid/; then
	  #  cp $outfilename /sps/dchooz/guang/grid_fit/test/
        # cp ./nFitBins*.txt /sps/dchooz/guang/grid_fit/HIII_nFitBins/
          cp ./s*.txt /sps/dchooz/guang/binningPaper/test315
#	fi
 #   fi
done
