# extraSuperpower
R package for two-way factorial design sample size calculation. This is performed in three steps. The package:
1. Calculates expected outcomes into a cell mean model.
2. Simulates the data
3. Estimates the power for a given sample size

These steps allow for independent and repeated measures experiments with balanced or unbalanced design. 

For the first step we provide a function to create mean values and standard deviation matrices. For repeated measures designs correlation and covariance matrices are also generated. For the second step separate functions are used to simulate independent and repeated measures experiments. Once the two-way factorial study is simulated under different sample sizes, the power under different statistical tests for these sample sizes can be estimated.

``extraSuperpower`` depends on packages ``MASS``, ``Matrix``, ``rlist``, ``afex``, ``fGarch``, ``truncnorm``, ``lmPerm``, ``ez``, ``nparLD``, ``Rfit``, ``stringr``, ``reshape2``, ``ggplot2`` and ``ggpubr``. 

To install, install dependencies and ``devtools``, if not done already.

``devtools::install_github("luisrmacias/extraSuperpower")``

Further installation and user guide in the ``vignette`` (vignettes folder).
Installation with vignettes takes a bit longer, nothing unbearable.

``devtools::install_github("luisrmacias/extraSuperpower", build_vignettes = TRUE)``
``browseVignettes(package = "extraSuperpower")``

First steps:

``library(extraSuperpower)``

``?calculate_mean_matrix    ## example of one between, one within (repeated measure) factorial design simulation``

``?test_power_overkn    ## example of independent measurements sample size calculation with plot``
