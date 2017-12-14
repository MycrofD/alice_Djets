#!/bin/bash


for count in `seq 0 1 8`; do
	#jet pT spectra, no D cut
	#root -l -b -q getSimSpectra.C'('$count',1,1,0,0,0)'
	#jet pT spectra (with D cut)
	#root -l -b -q getSimSpectra.C'('$count',1,1,0,1,0)'
	#jet pT spectra (with D cut) efficiency scaled
	root -l -b -q getSimSpectra.C'('$count',1,1,0,1,1)'
	#get D pT spectra
	#root -l -b -q getSimSpectra.C'('$count',1,0,0,1,0)'
	#get D pT spectra with jet pt cut 4-40
	#root -l -b -q getSimSpectra.C'('$count',1,0,1,1,0)'
done







