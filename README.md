# Spectrally Deconfounded Random Forests

This repository contains all the code used for the paper "[*Spectrally Deconfounded Random Forests*](https://doi.org/10.1080/10618600.2025.2569602)" by Markus Ulmer, Cyrill Scheidegger, and Peter Bühlmann (2025).

All the provided R-code uses the R-package [*SDModels*](https://markusul.github.io/SDModels/) v1.0.13, 
providing functionality to estimate SDForests and analyze the estimated functions. 
You can download this version from CRAN with `devtools::install_version("SDModels", version = "1.0.13", repos = "http://cran.us.r-project.org")`.

All the simulations were run on the [*Euler*](https://scicomp.ethz.ch/wiki/Euler) using the batch jobs in the 
`slurm` folder. See [*here*](https://scicomp.ethz.ch/wiki/Euler_applications_and_libraries_ubuntu) and sessionInfo.txt for details about the used libraries.

The experiments will need the folders `simulation_study/results` and `cBench/semiSimResults` to save their output.

-   `simulation_study` contains all the code to produce and visualize the simulations

-   `cBench` contains all the code to produce and visualize the single-cell gene expression experiments

-   The cBench data can be found here: <https://github.com/causalbench/causalbench>

Reproducing the different figures:
Most figures are created using simulation_study/visualization. R. The figures in Section 6 are generated using cBench/R/PlainSDF_plots.R and cBench/R/semiSimPlots.R. For the individual figures, you need to run the following experiments:

- Figure 3-6 run simulation_study/default_szenario.R
- Figure 7 run slurm/dimensions.sh
- Figure 8 and 9 run cBench/R/PlainSDF.R
- Figure 10 run slurm/semiSim.sh
- Appendix A Approximation of splitting criteria run simulation_study/fast.R
- Appendix C Non-linear confounding run simulation_study/nonlinearConfounding_2.R
