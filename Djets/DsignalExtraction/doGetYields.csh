#!/bin/bash

# B. Trzeciak -- script to run D-jets analysis chain for the efficiency-corrected signal extraction

## get efficiency corrected yields, the 7th argument is if eff corrected, the 8th argument is if tree has a postfix, and the 9th the postfix


./getYields.csh FASTwoSDD cuts806_preliminary/AnalysisResults_FASTwoSDD.root  testOutSignal AnalysisResults_fast_MCEffRMHijing_806_cutSys.root testOutMC 1 0 0 1
