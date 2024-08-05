# extraSuperpower
R package for two-way factorial design sample size calculation. This is performed in three steps.
1. Calculate expected outcomes into a cell mean model.
2. Simulates the data
3. Estimates the power for a given sample size
These steps allow for independent and repeated measures experiments. 

For the first step we provide a function to create mean values and standard deviation matrices. For repeated measures designs correlation and covariance matrices are also generated. For the second step separate functions are used to simulate independent and repeated measures experiments. 

``extraSuperpower`` depends on packages ``MASS``, ``afex``, ``fGarch``, ``truncnorm``, ``lmPerm``, ``ez``, ``nparLD``, ``Rfit``, ``stringr``, ``reshape2``, ``ggplot2`` and ``ggpubr``. Installation of package ``Superpower`` from CRAN is highly recommended for dowstream analysis.

To install, install dependencies and ``devtools``, if not done already.

``devtools::install_github("luisrmacias/extraSuperpower")``

``library(extraSuperpower)``

``?calculate_mean_matrix    ## example of one between, one within (repeated measure) factorial design simulation``

``?test_power_overkn    ## example of independent measurements sample size calculation with plot``
