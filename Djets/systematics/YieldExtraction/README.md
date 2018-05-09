# General Description
For each *p*~T~ bin defined by the user, the code evaluates a series of variations of the fit condition of the D-meson invariant mass distribution, retrieving the mean and the RMS of the yields from each variations.
In case of the efficiency scaling approach, the RMS of the yields for each (efficiency scaled) *p*~T,jet~ bin is already considered as the absolute yield uncertainty on that bin of the jet spectrum. The relative uncertainty is simply the absolute uncertainty divided by the mean value of the yield.
In case of the side-band scaling, for each fit variation in a given p~T,D~ bin, the value of Gaussian mean, sigma and background in the signal region are stored. Then, to build the jet *p*~T~ spectrum, a fixed, user-defined amount (*N*) of variations are chosen, independently and randomly in each *p*~T,D~ bin, and for each case the corresponding bkg-subtracted jet spectrum (using the signal region, sideband region and bkg normalization factor stored for that particular variation) as in the standard analysis. The RMS of the *N* variations of the jet yield in each *p*~T~ bin is considered as the absolute yield uncertainty on that bin. Again, the relative uncertainty is simply the absolute uncertainty divided by the mean value of the yield.

## Execution
- Before the execution, the steering macro needs to be configured, by setting the proper values in the `SetInputParametersDzero` / `SetInputParametersDstar` methods. In particular:
  - Set the names of input files, `TTree`/`THnSparse`, etc.
  - Set the number and edges of *p*~T,D~ bins *p*~T,jet~ bins. For the D^*^, the latter are used only for the efficiency scale method, since in the SB subtraction method the *p*~T,jet~ bins are already determined by the `THnSparse` projections (unless you call `SetRebinSpectrumIfSBApproach(kTRUE)`, see later on).
  - Set the D efficiency values for each *p*~T,D~ bin
  - For the D^0^: set the value of min, max and bin width for the invariant mass plots (for the D^*^, these come automatically from the `THnSparse` projections).
  - For the D^*^, SB approach: choose whether to rebin the jet p~T~ spectrum with the user-defined binning (calling `SetRebinSpectrumIfSBApproach(kTRUE)`) or keep the `THnSparse` binning
  - Define the configuration of the trials (look at the next section for details).
- To execute the code, the class has to be compiled first, and the steering macro loaded:

~~~~   C++
gROOT->LoadMacro("AliDJetRawYieldUncertainty.cxx++")
.L ExtractDJetRawYieldUncertainty.C
~~~~

- Then, independently of the approach chosen, the method `EvaluateBinPerBinUncertainty` has to be run for each p~T~ bin (i.e. *p*~T,jet~ for eff.scale method, *p*~T,D~ for SB subtraction method). `ptmin` and `ptmax` are the *p*~T~ limits of the bin to be looked at; for the D^*^, the *z* range for the projections can be also varied by setting `zmin` and `zmax`.
- Finally, the yield uncertainty is obtained by executing (once) the method `ExtractDJetRawYieldUncertainty`. In case of the SB subtraction approach, `nTrials` is the number of random variations to be picked from the total list of variations performed in the *p*~T,D~ bins, in order to get the final uncertainty.
- The code produces a .root file with the absolute value of the yield uncertainty plus another .root file with the mean of the yield from the variations plus its relative uncertainty.
- NOTE: In case a reflection template has to be added to the invariant mass fit function (D^0^ analysis), refer to the reflection section below.

## Configuration of the trials
There's no precise recipe (even in D2H on how to define the variations, and on their number). In the following, some information on the different possibilities, which can be enabled via the steering macro, are given.

* `sigmafixed_DPtBins` $\to$ array with the values of the sigma of the invariant mass fits in each *p*~T,D~ bin, taken from MC with only real D mesons
* `sigmafixed_JetPtBins` $\to$ array with the values of the sigma of the invariant mass fits in each *p*~T,ch~ ~jet~ bin, taken from MC with only real D mesons
* `chi2cut` $\to$ upper threshold of the chi^2/ndf for the fit. If the fit is worse than this threshold, is excluded from further operations.
* `meansigmaVar[6] = {kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE}` $\to$ Variations of the mean and sigma configuration, can be activated or excluded. From left to right: fixed sigma, fixed sigma (increased by 15%), fixed sigma (decreased by 15%), free mean and sigma, free sigma and fixed mean, fixed mean and sigma


* `bkgVar[8] = {kTRUE,kFALSE,kTRUE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE}` $\to$ Variations of the mean and sigma configuration, can be activated or excluded.
  From left to right: exponential, straight line, pol2, pol3, pol4, pol5, power law, power law*exponential

* `nRebinSteps=2`  $\to$ Variations of the binning of the mass spectrum. Set the number of cases and the rebin values in the array`rebinStep[2]={1,2}`

* `nMinMassSteps=2` $\to$ Variations of the lower edge of the fit range. Set the number of cases and the values in the array `minMassStep[2]={1.72,1.74}`

* `nMaxMassSteps=2`  $\to$ Variations of the lower edge of the fit range. Set the number of cases and the values in the array `maxMassStep[2]={2.00,2.03}`

* `nStepsBC=2`  $\to$ If `nStepsBC=0`, independently of what you set in thearray, skips evaluation of yield via bin counting. Otherwise, enables it, in the number of sigma from the peak center set in the array.
  `nSigmasBC[2]={3.5,4.0}`

* The following is a technicality for the `MultiTrial` code, just set `nmask` to the product of active variations of `meansigmaVar` and `bkgVar` (12 in this example) and fill the array mask with an equivalent amount of `1`'s.

  ~~~C++
  Int_t nmask = 12;
  Bool_t mask[12] =    {
       1,1,   // fixed sigma (Expo, Lin, Pol2, Pol3, Pol4, Pol5, PowLaw, PowLaw*Exp)
       1,1,   // fixed sigma+15%
       1,1,   // fixed sigma-15%
       1,1,   // free sigma, free mean
       1,1,   // free sigma, fixed mean
       1,1};  // fixed mean, fixed sigma
  ~~~

## Notes

* After executing the `EvaluateBinPerBinUncertainty` method, for each bin, always check the control canvas which is printed by the code, considering that:
  - If there are fits with $\chi^2$ over the value set in the macro, these will be chopped, and you'll get empty bins. This will often happen in case a not suitable bkg function is chosen
  - If the `sigmafixed` value is not properly set, not only you will bias the variations done with the fixed sigma, but you could get empty bins since the cases with a fit $\sigma >2$, or $\sigma<0.5$ from the `sigmafixed` value will be chopped.
* For the SB subtraction approach with the method `ExtractDJetRawYieldUncertainty`, in building the final uncertatinty never set a value of random variations greater than the minimum number of succesful variations for the varios *p*~T,D~ bins analysed. This, unless you allow repetition of randomly chosen trials, via `allowRepet=kTRUE`, which I won't anyway suggest.

## Reflections (D^0^ only)

The code implements the possibility of setting a template for the fitting of the reflection component of the invariant mass distribution in a given _p_~T~ (D or jet) bin.

This feature enters in play only during the bin-by-bin fitting phase (when calling the `EvaluateBinPerBinUncertainty` function), and can be enabled by setting as `kTRUE` the flag `refl` in the function arguments.

The ingredients needed are:

- A 1D histogram with the template distribution in the _p_~T~ bin under analysis
- The reflection/true signal ratio in the mass fit range OR the true signal mass distribution in the _p_~T~ bin under analysis (the ratio can be indeed fixed by hand, or evaluated by the reflection and true signal histograms).

An example of these plots (in several _p_~T~ bins) is given as an attachment in the JIRA ticket (file `reflections_fitted_DoubleGaus.root`). Note that these histograms are (strongly) _p_~T~-dependent. It's not mandatory that these histograms have the same binning and mass range of the data mass distribution.

The input file names and the histogram names are set in the `SetInputParameterDzero` section, together with:

- The value of the reflection/true signal ratio, if fixed from the user (first argument of `SetValueOfReflOverSignal`)
- The mass fit range of the bin under study, to automatically evaluate that ratio from the input histograms (in this case the first argument has to be set to -1)

### After these settings

It is possible to smoothen the reflection templates retrived from the MC analysis using the `FitReflDistr` in the steering macro: it needs an input file with the templates for each *p*~T~ bin inside (named `histRfl_N` and `histSgn_N`, with N the *p*~T~ bin number from 0 to `nbins`), the number of _p_~T~ bins (`nbins`) and a distribution to be used as guideline for the template smoothing (choose among `DoubleGaus`, `pol3`, `pol6`, `gaus` with the first as default).

This function can be used also to produce template variations (using different guide functions), which, together with variations of the refl/true signal ratio (by fixing it externally to different values) can help in evaluating the impact of the reflection template on the yield extraction, and its stability against template variations.