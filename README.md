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
    - 1000 Genomes Project phase three call set on GRCh37
    - `data/igsr_samples.tsv`
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
    - Raw data: `data/1000g_raw_combined/`
    - Cleaned data: `data/1000g_raw_combined/`
    - Separated by Population: `data/1000g_by_population/`
    - Separated by Super-population: `data/1000g_by_superpopulation/`
- Reference: https://dougspeed.com/1000-genomes-project/
### Phenotypes (Simulation)
- Goal: simulate phenotypes with different settings
- Tool: [LDAK](https://dougspeed.com/simulations/) 
- Options
    - `--make-phenos`: perform simulation
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
    - File naming: `power=x_her=y_num-causals=z.phen`
    - Columns: FID IID PHEN1 PHEN2 PHEN3 ...
- Dataset
    - `data/1000g_pheno/` 
- Reference: https://dougspeed.com/simulations/
## Methods
### C+T
### BASIL
### BayesR
## Evaluation
### Data Splitting
- For our study on PRS, we adopted two distinct approaches to split our data into training, validation, and testing sets, ensuring a comprehensive evaluation of our methods. 
    - In one approach, we drew inspiration from Problem Set 4 (PS4), where data was segmented by population. This method allowed us to examine the influence of genetic diversity on PRS performance by closely mirroring real-world population structures. 
    - Alternatively, we explored a more conventional split, randomly dividing the data into 80% training, 10% validation, and 10% testing segments, with a keen focus on maintaining population and sex balance across these sets. This approach offers insights into the effects of random sampling on the robustness and generalizability of PRS methods. 
- Dataset
    - `data/split_ids/`
        - `by_population/`
        - `by_superpopulation/`
## Results 
## Discussion
## Related Work
- Frank Dudbridge (2013) provided a comprehensive overview of the power and predictive accuracy of PRS, highlighting the statistical principles underpinning these methods and their implications for disease risk prediction. 
- The 1000 Genomes Project Consortium (2015) provided a landmark achievement in human genetics, offering a global reference for human genetic variation. This comprehensive dataset has been instrumental in PRS research, allowing for the assessment of methodological performance across diverse genetic backgrounds.
- Vilhjálmsson et al. (2015) introduced a significant advancement with LDPred, a Bayesian method considering linkage disequilibrium and polygenic traits. This approach marked a step forward in PRS accuracy, emphasizing the complexity of genetic architecture in disease prediction.
- Khera et al. (2018) showcased the practical application of PRS, particularly in predicting coronary artery disease, underscoring the importance of large reference datasets for improving prediction accuracy. Their work emphasized the clinical utility of PRS, driving further research into methodological improvements.
- Lloyd-Jones et al. (2019) refined the Bayesian approach to PRS with ==BayesR==, demonstrating the method's capability to enhance predictive power by incorporating SNP effect sizes more accurately. This study contributed to the growing preference for Bayesian models in complex disease prediction.
- Ge et al. (2019) evaluated various PRS methods, revealing performance variability across diseases and traits. Their findings stressed the necessity for flexible PRS methodologies adaptable to different genetic architectures.
- Qian et al. (2020) presented a scalable framework for Lasso-based sparse regression with ==BASIL== in the snpnet package, efficiently handling ultra-high-dimensional genetic data. This approach represents a leap forward in managing large-scale datasets for PRS calculation.
- Yi Ding (2023) brought critical attention to the variability of polygenic scoring accuracy across different genetic ancestries. Their results highlight the need to move away from discrete genetic ancestry clusters towards the continuum of genetic ancestries when considering PGSs.


## References
- Dudbridge F. (2013) Power and Predictive Accuracy of Polygenic Risk Scores. PLOS Genetics 9(3): e1003348. https://doi.org/10.1371/journal.pgen.1003348
- The 1000 Genomes Project Consortium. (2015) A global reference for human genetic variation. Nature 526, 68–74 https://doi.org/10.1038/nature15393
- Vilhjálmsson, B. J., Yang, J., Finucane, H. K. et al. (2015) Modeling Linkage Disequilibrium Increases Accuracy of Polygenic Risk Scores. American journal of human genetics. 97(4), 576–592. https://doi.org/10.1016/j.ajhg.2015.09.001
- Khera, A.V., Chaffin, M., Aragam, K.G. et al. (2018) Genome-wide polygenic scores for common diseases identify individuals with risk equivalent to monogenic mutations. Nat Genet 50, 1219–1224. https://doi.org/10.1038/s41588-018-0183-z
- Lloyd-Jones, L.R., Zeng, J., Sidorenko, J. et al. (2019) Improved polygenic prediction by Bayesian multiple regression on summary statistics. Nat Commun 10, 5086. https://doi.org/10.1038/s41467-019-12653-0
- Ge, T., Chen, CY., Ni, Y. et al. (2019) Polygenic prediction via Bayesian regression and continuous shrinkage priors. Nat Commun 10, 1776. https://doi.org/10.1038/s41467-019-09718-5
- Qian J, Tanigawa Y, Du W, Aguirre M, Chang C, Tibshirani R, et al. (2020) A fast and scalable framework for large-scale and ultrahigh-dimensional sparse regression with application to the UK Biobank. PLoS Genet 16(10): e1009141. https://doi.org/10.1371/journal.pgen.1009141
- Ding, Y., Hou, K., Xu, Z. et al. (2023) Polygenic scoring accuracy varies across the genetic ancestry continuum. Nature 618, 774–781. https://doi.org/10.1038/s41586-023-06079-4
---
## Challenges
- In tackling our study on PRS, we ran into a common difficulty: the massive size of the 1000 Genomes dataset. With limited time and computing power, analyzing the entire dataset isn't feasible. So, we're taking a practical approach. We'll divide the data by super-population and population categories, much like we did in Problem Set 4 (PS4). This way, we can still conduct our GWAS effectively. By focusing on specific groups, we aim to strike a balance between thorough analysis and resource constraints. 
- Also, another challenge lies in determining the optimal settings for phenotype simulation. While we recognize the importance of varying simulation parameters to comprehensively evaluate PRS performance, we face uncertainties regarding the number of settings required to yield sufficient insights. 
## Remaining Work 
- We have completed the essential stages of data preprocessing and splitting for training, validation, and testing dataset. Additionally, our initial testing encompassed the evaluation of three key methods: Clumping + Thresholding (C+T), BASIL, and BayesR, yielding preliminary insights. Moving forward, our focus will shift towards extending our GWAS across a broader dataset. This expansion aims to deepen our understanding of how different phenotype simulations and the application of various PRS methods influence the resulting $R^2$ score. Furthermore, we are committed to enhancing the accessibility of our findings by visualizing results through tables and graphs. This approach will facilitate easier interpretation and comprehension, ensuring that our research contributes effectively to the broader discourse within the field.
