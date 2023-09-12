# extraSuperpower
R package that prepares input for two-way factorial design sample size calculation with a separate workflow for independent and repeated measures experiments. For independent measures, effect size and the power attained with a sample size *n* is calculated by exact methods. For repeated measures, cell means can be easily used as input for simulation based sample size calculation with package ``Superpower``.

We provide a function to create a matrix of mean values that can be used as input for the ``Superpower ANOVA_design`` function for repeated measures two-factor designs. The input to create this matrix is a reference mean, number of levels of factor A, number of levels of factor B and expected effect magnitudes for each. Interaction can also be easily modeled.

``extraSuperpower`` depends on packages ``MASS``, ``afex``, ``reshape2`` and ``ggplot2``. Installation of package ``Superpower`` from CRAN is highly recommended for dowstream analysis.

To install, install dependencies and ``devtools``, if not done already.

devtools::install_github("luisrmacias/extraSuperpower")

library(extraSuperpower)

?mean_sd_matrix    ## example of one between, one within (repeated measure) factorial design sample size calculation with ``extraSuperpower`` and ``Superpower``

?plot_powercurve    ## example of independent measurements power curve plotting
