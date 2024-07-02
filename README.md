# Composite likelihood and Stochastic approximations

This repo stores the code to reproduce the results described in the paper "When composite likelihood meets stochastic approximation".

The code relies on two custom packages, whose source code is stored in the `helper_packages/` folder.
The code for the synthetic experiments is stored in the folders `isi_sims/` and `gf_sims/`, while the real data analysis can be found in the `mental_health_data/` folder.

The code for the synthetic experiments is divided in three `R` scripts, called `*_sims.R` , `raw_output.R` and `paper_plots.R` .

-   `*_sims.R` creates the sub-folders used to store the results according to the simulation setting label provided. It then computes point estimates for the stochastic algorithm under different sampling schemes and for the numerical optimiser. Finally, it computes the asymptotic variances at different running lengths of the algorithm as described in the paper.

-   `raw_output.R` reads th objects created by `*_sims.R` and computes quantities to track estimators performance both in terms of pointwise convergence and inference, as described in the paper. It also saves some raw visualisations.

-   `paper_plots.R` reads the objects created by `raw_output.R` and refines the plots as they appear both in the main paper and in the supplementary material.

In addition, the folder `gf_sims/` contains also the file `rpl.R` which provides the code to reproduce the supplementary simulation comparison with the randomized pairwise likelihood estimator (Mazo et al. 2024)

The code for the real data analysis is collected in one single script called `analysis.R` . It reads the data stored in the `mental_health_data/data/` folder and carries out the analysis described in the paper. The original data were publicly available at [https://catalog.data.gov/dataset](https://catalog.data.gov/dataset/national-epidemiologic-survey-on-alcohol-and-related-conditions-nesarcwave-1-20012002-and-) but are now accessible only via request. Because of this, we provide a mock dataset in the `mental_health_data/data/` folder to carry out the analysis. Such dataset is generated from an Ising model using the numerical estimates computed on the real data as true parameters.

