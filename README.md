# extraSuperpower
R package intended as a workflow for simulation based two-way factorial design sample size calculation.

We provide a function to create a matrix of mean values that can be used as input for the Superpower ANOVA_design function for independent of repeated measures two-factor designs. The input to create this matrix is a reference mean, number of levels of factor A, number of levels of factor B and expected effect magnitudes for each. Interaction can also be easily modeled.

It depends on packages MASS, afex, reshape2 and ggplot2. Installation of package Superpower from CRAN is highly recommended for dowstream analysis.

To install.

devtools::install_github("luisrmacias/extraSuperpower")

See ?mean_sd_matrix for examples.
