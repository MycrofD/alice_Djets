#!/bin/bash

#outPlotDir= plots
#mkdir -p $outPlotDir

	#charm, jet pT spectrum, jet pT cut (0), D pt Cut, eff
	#root -l -b -q plotSimSpectra.C'(0,1,0,0,0,"'$outPlotDir'")'
	

	#beauty, jet pT spectrum, jet pT cut (0), D pt Cut, eff
	root -l -b -q plotSimSpectra.C'(1,1,0,0,0)'
	root -l -b -q plotSimSpectra.C'(1,1,0,1,0)'	
	root -l -b -q plotSimSpectra.C'(1,1,0,1,1)'

	#beauty, D pT spectrum, jet pT cut (0), D pt Cut, eff
	root -l -b -q plotSimSpectra.C'(1,0,0,0,0)'
	root -l -b -q plotSimSpectra.C'(1,0,1,0,0)'

	#charm, jet pT spectrum, jet pT cut (0), D pt Cut, eff
	root -l -b -q plotSimSpectra.C'(0,1,0,0,0)'
	root -l -b -q plotSimSpectra.C'(0,1,0,1,0)'	

	#charm, D pT spectrum, jet pT cut (0), D pt Cut, eff
	root -l -b -q plotSimSpectra.C'(0,0,0,0,0)'
	root -l -b -q plotSimSpectra.C'(0,0,1,0,0)'
