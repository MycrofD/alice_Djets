
- Prepare root files with an output of your analysis: data + MC for efficiencies and response matrix 
- If a separete bkg.fluctuation matrix is needed, prepare it first
- Prepare configs:

	-  config *.h file, examples are: Dzero.h or Dstar.h 
		- adjust variables to your needs: set up the system, D meson specie, D bins that you want to use, used R parameter, signal and SB ranges, if you want reflections (for D0), sigma in, etc ...
		- IMPORTANT -- set up also names (and description) of your non-prompt (needed for the B feed-down subtraction) and prompt (needed only at the last step when comparing the measured x-section to a model, otherwise can be skipped) D-jet simulation files (path to the files needs to be given in diffrerent place)
		- in case of running with reflections for D0, you must already have the reflection files prepared


- running scripts: here you need to define your configuration
	- run.sh: the main script running the whole analysis -- steps to be run -  settings are passed via a script that runs this one, usually this script doesn't have to be modified, unless you want to skip one of the analysis steps
		-- !! you may want to change the jet pT ranges for the efficiency evaluation: jetpteffmin, jetpteffmax; and how the denominator is define.
	- run_analysis.csh: this script runs run.sh, define here paths to data and MC output files, file names, the *.h configuration file name, if reflections are used, if external bkg. fluctuation matrix is used ...
		-- !!!! setup number of POWHEG sim files !
	- run_main.csh: this script runs run_analysis.sh, define here unfolding algorithm, reg. parameter for the unfolding, prior, prior type and bkg. 
Fluctuation matrix type (this you need to set once you use an external bkg. fluctuation matrix, then based on this type name of a root file with the corresponding bkg. fluctuation matrix will be set in the run_analysis.csh script - you should set path to your file). Set up also flags and RELEVANT SCRIPTS if you want to run analysis for systematic unc. evaluation.


#### All analyses happen inside `Djets`. Enter `Djets`. 
$ `bash run_main.csh`

##### Edit `run_analysis.csh` for toggling between default runs, and systematics-runs.
##### Change `outputdirectorybase` accordingly.
Example: for multi-trial systematics, 

0. 	Change jet radius at two places: `run_analysis.csh` and `configDzero_pp.h`.
1.  Use `outputdirectorybase=$HOME/Work/alice/analysis/pp5TeV/D0jet/results_APW/FinalSys/RawSysFinal_DzeroR0${R}_`
    Comment other `outputdirectorybase` variables.
2.  `$ bash run_main.csh` # This should succeed up to D-signal extraction, and fail for B-feed down steps.
3.  Check `run_settings.csh`. It contains flags for toggling between different systematics types.
    For multi-trial, use `flagMulti=1`.
4. 	Enter `DsignalExtraction`.
5. 	Run `getYields.csh`.

