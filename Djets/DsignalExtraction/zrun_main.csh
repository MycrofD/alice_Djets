#!/bin/bash

#--------------------------------------------------
#  Author A.Mohanty
#  Utrecht University
#  auro.mohanty@cern.ch
#--------------------------------------------------
#-- script to run D-jets z_||,ch analysis


zrun_file=zrun4.csh
unoriginal=unoriginal_bash
cat $unoriginal > zrun.csh
cat $zrun_file >> zrun.csh

#EOS_local=/eos/user/a/amohanty/
EOS_local=

# Running for all required bins of jet pT
#----------------------------------------- #{5.0, 7.0, 10.0, 15.0, 36.0, 5.0, 15.0, 30.0, 15.0, 50.0, 10.0, 16.0, 36.0};
for R in 3 4 6; do
#for R in 4; do
	for thing in 1 2 3 9 6; do
	#for thing in 6; do
		bash zrun.csh $thing $R ${EOS_local} 
	done
	#bash zrun4D.csh $R ${EOS_local} 
done
