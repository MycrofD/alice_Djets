## Recap:  
file paths were updated in the _run.csh_ file.  

--------------------------

### what was done  
beauty simulations from grid were downloaded to out/QM2018/beauty.  
./run_pPbMain.csh  
### what happened  
ran ./run_pPbDzeroAnalysis.csh  
which in turn ran ./run.csh  
### what we have now  
E-LoadBackgroundMatrix: Background matrix file $HOME/.../alice_Djets/Djets/ResponseMatrix/BkgRM03/RandCones_BkgM_Djet5Excl.root not found.  

Solution:  
copy paste BkgRM03/* to the place it is looking for, i.e. alice_Djets/Djets/ResponseMatrix/BkgRM03/*   

### what else is being done...   
Issue: RooUnfold path needs to be changed.  
How solution was found: checked for the _alien_ file path for RooUnfold in _run.csh_.   
Solution: Go to _subtractFD.C_ in _FDsubtraction_  
look at ~line:583. 
> gSystem->Load("/home/.../RooUnfold-1.1.1/libRooUnfold.so");  

Issue: RooUnfold path needs to be changed, again, in _unfold_Bayes.C_.     
Solution: Go to *unfold_Bayes.C* in _unfolding_  
look at ~line:45. 
> gSystem->Load("/home/.../RooUnfold-1.1.1/libRooUnfold.so");  

### new issue  
lack of charm simulations led to ``Illegal pointer to class object... in finalJetSpectra.C''  

### more problems  
