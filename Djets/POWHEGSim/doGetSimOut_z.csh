#!/bin/bash
# ---------------------
# A.Mohanty
# Utrecht University
# auro.mohanty@cern.ch
# ---------------------

count=0
while [ $count -lt $1 ]
do
	root -l -b -q getSimSpectra_z.C'("'$2'",'$count','$3','$4','$5','$6', "'$7'","'$8'","'$9'")'
	((count++))
done

#for count in `seq 0 1 $1`; do
#	root -l -b -q getSimSpectra.C'("'$2'",'$count','$3','$4','$5','$6', "'$7'","'$8'","'$9'")'
#done
