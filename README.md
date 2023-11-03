# monalisa

Code, data, and model objects used to generate the results in *"That's not the Mona Lisa! How to interpret spatial capture-recapture density surface estimates"* by Ian Durbach, Rishika Chopara, David Borchers, Rachel Phillip, Koustubh Sharma, and Ben Stevenson (to appear in Biometrics).

The main analysis in the paper is an illustration of the different SCR density surfaces that uses the Mona Lisa as the intensity surface of the point process that generates animal activity centers (the "true" surface, if you like). These illustrations are done in both a frequentist and Bayesian context, so show the the misleading inferences that can be drawn from these surfaces do not depend on which context you work within.

**To reproduce the frequentist Mona Lisa analysis**

Scripts are found in the root directory.

- run *mona-analysis.R*. This creates `secr` inputs from the image *hires_mona.jpg* and saves them in *output/mona_results.RData* and fits SCR models, which are saved in *output/mona_results.RData*. Note these outputs are already in the output folder, and running the code again will overwrite these outputs. 
- run *mona-plots.R*. This produces plots found in the paper.
- the helper file *run_secr.R* is called by *mona-analysis.R*. 

**To reproduce the Bayesian Mona Lisa analysis**

Scripts are found in the 'bayesian_code' directory.

- run *Plots_Code.R* to generate Bayesian versions of the main illustrative examples and plots. This produces the output and plots in Web Appendix B of the paper. Note that the working directory should be the 'bayesian_code' folder.
- run *UncertaintyPlots_Code.R* to generate plots of the uncertainty associated with some of the density surfaces. This produces the output and plots in Web Appendix C.

A script for the introductory example shown in Figure 1 of the paper is in *intro_example/linear-poisson-plots.R*