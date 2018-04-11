# alice_Djets  

Where to store your files.  
The paths:  
1. _AnalysisResults.root_ at $HOME/.../out/...(QM2018).../  
2. _Efficiencies_ at $HOME/...(ALICE_HeavyFlavour/work/jets2).../alice_Djets/Djets/efficiency/  
3. _Reflections_ (reflections_fitted_DoubleGaus.root) at $HOME/...(ALICE_HeavyFlavour).../out/...(QM2018).../  



-----------------------------
For step-2:  
D-jet signal extraction file:  

all path variables:  
1.

------------------------------


### step-1  
----------  
The data file _AnalysisResults.root_ is analysed by [**SideBand subtraction**](Djets/DsignalExtraction/signalExtraction_SB.C) method.  
Depending on options provided, **efficiency correction** and **reflections** are taken care of.  

_Djets/DsignalExtraction/signalExtraction_SB.C_ creates  
1. _Djets/signalExtraction/plots/_,
2. _Djets/signalExtraction/plotsNoEff/_, 
3. _Djets/signalExtraction/JetPtSpectra_SB_eff.root_, 
4. _Djets/signalExtraction/JetPtSpectra_SB_noEff.root_  

### step-2  
----------  
**B-feed down**, **unfolding**   
To separate contribution from beauty quarks, _b-feed down_ is subtracted. 
The resulting jet-pt spectrum is _unfolded_ to correct it for _detector inefficiencies_.  

Following script chain is employed:  
_./run_pPbMain.csh -> ./run_pPbDzeroAnalysis.csh ->  ./run.csh_ 

What do the scripts do?  
1. _./run_pPbMain.csh_:  
runs _./run_pPbDzeroAnalysis.csh_  

2. _./run_pPbDzeroAnalysis.csh_:  
creates _outputdirectorybase_ where ... are stored.   
... are then used in ...  

3. _./run.csh_:   
mkdir -p $outdirBase/$outdir  


### step-3  
----------  

**systematic uncertainty:unfolding**, **systematic uncertainty:cut variations**, **systematic uncertainty:raw yield extraction**
