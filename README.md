# How Simulation Parameters Impact Polygenic Scoring Accuracy Across Various Methods

## Group 18
- Keng-Chi Chang
- Po-Chun Wu
## Introduction
Polygenic Risk Scores (PRS) hold promise for personalized medicine but face challenges in prediction accuracy. This study investigates how simulation parameters affect PRS accuracy using Clumping + Thresholding, BASIL (Lasso-based method in snpnet package), and BayesR (Bayesian Multiple Regression-based method in GCTB package). By manipulating parameters such as effect size distribution and heritability, we aim to understand optimal settings for each method. Our findings will inform the selection of PRS methods for precise disease risk prediction in clinical and research settings, advancing personalized healthcare interventions.
## Methods
- Clumping + Thresholding (C+T)
- [Batch Screening Iterative Lasso (BASIL)](https://github.com/rivas-lab/snpnet)
- [Bayesian Multiple Regression](https://cnsgenomics.com/software/gctb/#Overview)
## Data
- [1000 Genomes](https://www.internationalgenome.org)
## Tools 
- [PLINK 2](https://www.cog-genomics.org/plink/2.0/)
## Implementation
## Challenges
## Remaining Work 
## Results
## Discussion
## Related Work

## References
- Qian J, Tanigawa Y, Du W, Aguirre M, Chang C, Tibshirani R, et al. (2020) A fast and scalable framework for large-scale and ultrahigh-dimensional sparse regression with application to the UK Biobank. PLoS Genet 16(10): e1009141. https://doi.org/10.1371/journal.pgen.1009141
- Lloyd-Jones, L.R., Zeng, J., Sidorenko, J. et al. (2019) Improved polygenic prediction by Bayesian multiple regression on summary statistics. Nat Commun 10, 5086. https://doi.org/10.1038/s41467-019-12653-0
- The 1000 Genomes Project Consortium. (2015) A global reference for human genetic variation. Nature 526, 68–74 https://doi.org/10.1038/nature15393
