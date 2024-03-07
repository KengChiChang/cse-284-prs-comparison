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
## Data Preprocessing
### Genotypes (1000 Genomes)
- Goal: get the dataset of all chromosomes from 1000 Genomes with quality control
- Tool: PLINK
- Data Quality Control
    - `--maf 0.01 --snps-only just-acgt --max-alleles 2 --rm-dup exclude-all`
- Create new binary binary fileset
    - `--make-bed`
    - `.bed` `.bim` `.fam`
- Dataset
    - Raw data: `1000g_raw_combined/`
    - Cleaned data: `1000g_raw_combined/`
    - Separated by Population: `1000g_by_population/`
    - Separated by Super-population: `1000g_by_superpopulation`
- Reference: https://dougspeed.com/1000-genomes-project/
### Phenotypes (Simulation)
- Goal: simulate phenotypes with different settings
- Tool: [LDAK](https://dougspeed.com/simulations/)
- Options
    - `--power`: power indicates the probability that we can detect an association at the desired significance level, given that there is actually an association.
    - `--her`: heritability describes how much variation in the phenotype is described by genetics.
    - `--num-phenos`: the number of phenotypes to generate.
    - `--num-causals`: the number of predictors contributing to each phenotype 
    - Notes: LDAK will pick causal predictors at random. Also, LDAK will by default sample effect sizes from a standard normal distribution.
- Our simulation settings
    - Power: {-0.25, -1}
        - Notes: GCTA Model with power of -1 vs. Human Default Model with power of -0.25
    - Heritability: {0.1, 0.3, 0.5, 0.7, 0.9}
    - Number of causal SNPs = {1, 10, 100, 1000, 10000} 
    - Total combinations: 2 * 5 * 5 = 50
- Reference: https://dougspeed.com/simulations/
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
