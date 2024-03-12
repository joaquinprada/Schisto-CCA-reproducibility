# Schisto-CCA-reproducibility
Data and Scripts for the analysis of the reproducibility of CCA in two settings in Uganda. The analysis is presented in the manuscript by E. Kabbas-Piñango et al., "Reproducibility matters: intra- and inter-sample variation of the point-of-care circulating cathodic antigen test in two Schistosoma mansoni endemic areas in Uganda", Philos Trans R Soc Lond B Biol Sci. 2023 Oct 9;378(1887):20220275. doi: 10.1098/rstb.2022.0275.  

Description of the Files:

-----
**DataMayuge.csv** contains the raw data from Mayuge.

Kato-Katz data is coded as follows: KK_X_Y, where *X* codes the day (either 1, 2 or 3) and *Y* codes the repeat (either A or B). 

Point of Care Circulating Cathodic Antigen (POC-CCA) G-score is coded as follows POC_CCA_X, where *X* codes the day (either 1, 2 or 3). A subscript _dup2 and _dup3 indicates the second and third duplicate for each day.

**DataTororo.csv** contains the raw data from Tororo. 
It is coded as Mayuge above.

-----
**ScriptMayuge.R** contains the annotated code to run the latent class analysis model for Mayuge data. The script also contains the code to generate the top-half of Figure 1.

**ScriptTororo.R** contains the annotated code to run the latent class analysis model for Tororo data. The script also contains the code to generate the bottom-half of Figure 1.

**SimulationError.R** contains the annotated code to do the simulations and the analysis required to generate Figures 2 and 3.

**SimulationErrorFunctions.R** contains some additional functions to run the script above (SimulationError.R).

-----
**kints.csv** is an additional file with posterior draws from the model run in Clark et al., 2022, Translating From Egg- to Antigen-Based Indicators for Schistosoma mansoni Elimination Targets: A Bayesian Latent Class Analysis Study.
