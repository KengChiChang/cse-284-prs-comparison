# How Simulation Parameters Impact Polygenic Scoring Accuracy Across Three Methods

## Group 18
- Keng-Chi Chang
- Po-Chun Wu
## Introduction
Polygenic Risk Scores (PRS) hold promise for personalized medicine but face challenges in prediction accuracy. This study investigates how simulation parameters affect PRS accuracy using Clumping + Thresholding, BASIL (Lasso-based method in `snpnet` package), and BayesR (Bayesian Multiple Regression-based method in `GCTB` package). By manipulating parameters such as heritability, number of causal SNPs, and how minor allel frequency affects heritability, we aim to understand whether certain PRS method might have advantage over others, and under what conditions. Our findings can inform the selection of PRS methods for precise disease risk prediction in clinical and research settings, advancing personalized healthcare interventions.


## Methods
- Clumping + Thresholding (C+T) using Plink
- [Batch Screening Iterative Lasso (BASIL)](https://github.com/rivas-lab/snpnet)
- [Bayesian Multiple Regression (BayesR)](https://cnsgenomics.com/software/gctb/#Overview)


## Hypothesis
1. Since Clumping + Thresholding does not directly model the additive effects of the contribution from mutliple SNPs, we expect that C+T would perform the worst when the number of causal SNPs increases.
2. Both BASIL and BayesR should be better at capturing additive effects from many SNPs, while shrinkage-based method might incorrectly shrink some effects to zero.

Consider heritability model of the form [(Speed, Holmes, and Balding 2023)](https://www.nature.com/articles/s41588-020-0600-y)
$$\mathbb{E}\left[h^2_j\right]=\sum_j w_j\left[f_j(1-f_j)\right]^{(1+\alpha)},$$
where $w_j$ is the weight for SNP $j$ and $f_j$ is its minor allele frequency (MAF). For simplicity, we assume that weight $w_j$ is a constant.

3. When $\alpha=-1$, the expected heritability reduces to $\sum_j w_j$.
SNPs have linear effect on heritability, regardless of minor allele frequency. Since this data generating process is closer to Lasso, we expect BASIL to perform better under this regime.
4. When $\alpha=-0.25$, the expected heritability becomes $\sum_j w_j\left[f_j(1-f_j)\right]^{0.75}$. SNPs that are more common would have higher contribution to heritability. We expect BayesR would capture this nonlinear relationship better than Lasso based BASIL.

The current published papers comparing these methods have mixed findings. Most of them either only use simulated genotypes, simulate phenotypes on a narrow demography ([Lloyd-Jones et al. 2019](https://doi.org/10.1038/s41467-019-12653-0); [Qian et al. 2020](https://doi.org/10.1371/journal.pgen.1009141)), or does not systematically evaluate alternative scenarios ([Ding et al. 2023](https://doi.org/10.1038/s41586-023-06079-4)). To our knowledge, this is the first systematic comparison on real genomic data across ancestries.


## Directory for Data and Computation
We are analyzing real genome data from [1000 Genomes Project](https://www.internationalgenome.org) phase three call set on GRCh37. Since the size of the data is huge, for the ease of data storage and collaborating, we store our data on Google Drive and use Google Colab for compute (by reading/copying files directly from Google Drive).

[Here](https://drive.google.com/drive/folders/1vONriV2u1j2BinWGxUhBwLnFN94u3Mig?usp=drive_link) is the link to our Google Drive project folder. The files are organized as follows:

```
/CSE-284-Final-Project/
├── data
│   ├── 1000g_by_superpopulation ← Genome (all autosomes) split by superpopulation
│   │   ├── AFR_{train,val,test}.{bed,bim,fam,pgen,pvar.zst,psam}
│   │   ⋮    ⋮
│   │   └── SAS_{train,val,test}.{bed,bim,fam,pgen,pvar.zst,psam}
│   ├── 1000g_by_population ← Genome (all autosomes) split by perpopulation
│   │   ├── ACB_all.{bed,bim,fam,pgen,pvar.zst,psam}
│   │   ⋮    ⋮
│   │   └── YRI_all.{bed,bim,fam,pgen,pvar.zst,psam}
│   ├── 1000g_combined ← Genome (all autosomes) 
│   │   └── 1000g.{bed,bim,fam,pgen,pvar.zst,psam}
│   ├── 1000g_raw_combined ← Raw genome before cleaning and quality control
│   │   ├── all_phase3.{pgen.zst,pvar.zst}
│   │   ├── phase3_corrected.psam
│   │   ├── raw.{bed,bim,fam}
│   │   └── clean.{bed,bim,fam}
│   ├── 1000g_pheno ← Simulated phenotype by parameters
│   │   └── power={-1,-0.25}_her={0.1,0.3,0.5,0.7,0.9}_num-causals={1,10,100,1000,10000}.phen
│   ├── igsr_samples.tsv ← Reference table for population and sex information
│   ├── split_samples.csv ← Reference table for sample splitting
│   ├── split_ids ← Ids after sample splitting (for `--keep`)
│   │   ├── by_superpopulation
│   │   └── by_population
│   └── test_data ← Genome and phenotype data for benchmark (from PS4)
│       ├── CEU_chr19_normed.{bed,bim,fam,pgen,pvar.zst,psam}
│       ├── GBR_chr19_normed.{bed,bim,fam,pgen,pvar.zst,psam}
│       └── kgvcf_ldl.{phe,phen,pheno}
└── notebook ← Data processing and analysis
    ├── 0_download_1000genome_data_and_simulate_phenotypes.ipynb
    ├── 1_convert_to_pgen.ipynb
    ├── 2_convert_to_bed.ipynb
    ├── 4_split_train_val_test.ipynb
    ├── 5_simulate_phenotype.ipynb
    ├── 6_update_phenotype_format.ipynb
    ├── 11_test_snpnet.ipynb
    ├── 12_test_GCTB.ipynb
    └── 15_test_C+T.ipynb
```


## Data Preprocessing

### [Genotypes (1000 Genomes)](https://github.com/KengChiChang/cse-284-prs-comparison/blob/main/notebook/0_download_1000genome_data_and_simulate_phenotypes.ipynb) [![image](https://img.shields.io/badge/Colab-F9AB00?style=for-the-badge&logo=googlecolab&color=525252)](https://colab.research.google.com/drive/1fsrn4J2BjSbKy6PdnAUmo6Rz7Hbbclbb?usp=sharing)
- Goal: get the dataset of all autosomes from 1000 Genomes with quality control
- Tool: PLINK
- Data Quality Control
    - `--maf 0.01 --snps-only just-acgt --max-alleles 2 --rm-dup exclude-all`
- Convert to binary formats for different packages
    - For BASIL/`snpnet`: `--make-pgen` generates `.pgen` `.pvar` `.psam`
    - For BayesR/`GCTB`: `--make-bed` generates `.bed` `.bim` `.fam`
- Dataset
    - Raw data: `data/1000g_raw_combined/`
    - Cleaned data: `data/1000g_combined/`
    - Split by Population: `data/1000g_by_population/`
    - Split by Super-population: `data/1000g_by_superpopulation/`
- Reference: https://dougspeed.com/1000-genomes-project/

### [Phenotypes (Simulation)](https://github.com/KengChiChang/cse-284-prs-comparison/blob/main/notebook/5_simulate_phenotype.ipynb) [![image](https://img.shields.io/badge/Colab-F9AB00?style=for-the-badge&logo=googlecolab&color=525252)](https://colab.research.google.com/drive/1SwnUGFM03mQmzqktkIUt9pW-OBqze0jc?usp=sharing)
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

### [Sample Splitting](https://github.com/KengChiChang/cse-284-prs-comparison/blob/main/notebook/4_split_train_val_test.ipynb) [![image](https://img.shields.io/badge/Colab-F9AB00?style=for-the-badge&logo=googlecolab&color=525252)](https://colab.research.google.com/drive/16fbgP6IK38J7wIoPl_Zj-Ts7Wk_rFqX3?usp=sharing)
- For our study on PRS, we adopted two distinct approaches to split our data into training, validation, and testing sets, ensuring a comprehensive evaluation of our methods. 
    - We randomly dividing samples from **each superpopulation** into 80% training, 10% validation, and 10% test sets, while keeping each set has balanced number of samples from each population and sex. 
    - Alternatively, we drew inspiration from Problem Set 4 (PS4), where data was segmented by **population**. This method allowed us to examine the influence of genetic diversity on PRS performance by closely mirroring real-world population structures. 
- Dataset
    - `data/split_ids/`
        - `by_population/`
        - `by_superpopulation/`



## Details

$$PRS_i = \sum_{j \in S} \beta_j X_{ij}$$

### [C+T](https://github.com/KengChiChang/cse-284-prs-comparison/blob/main/notebook/13_test_C+T.ipynb) [![image](https://img.shields.io/badge/Colab-F9AB00?style=for-the-badge&logo=googlecolab&color=525252)](https://colab.research.google.com/drive/1qmOKXFYZ5QZ4ZXcYXd8Hxry66RztcvgA?usp=sharing)
- Tool: PLINK
- Reference: UCSD CSE 284 Week 6 Lecture Slides, PS3, and PS4
1. Perform a GWAS to estimate per-variant effect sizes ($\beta$'s)
    - `--linear --maf 0.05`
2. Perform LD clumping (pruning) to get an independent set of SNPs
    - `--clump gwas.assoc.linear`
    - `--clump-p1 0.0001`: significance threshold for index SNPs   
    - `--clump-r2 0.1`: LD threshold for clumping
    - `–clump-kb 250`:  physical distance threshold for clumping
3. Choose all SNPs with $p$-value less than the threshold $T$
4. Computer PRS and evalueate accuracy ($R^2$)
    - `--score`
5. Computer final PRS with optimal $T$ and evaluate on a separate testing dataset
- Repeat steps 2-4 to find out which $T$ & clumping parameters work the best (validation dataset)

### [BASIL](https://github.com/KengChiChang/cse-284-prs-comparison/blob/main/notebook/11_test_snpnet.ipynb) [![image](https://img.shields.io/badge/Colab-F9AB00?style=for-the-badge&logo=googlecolab&color=525252)](https://colab.research.google.com/drive/1OR4CZl0jsxGFqEJ0MFgDkcrLgv4YNv94?usp=sharing)

- Tool: [`snpnet`](https://github.com/junyangq/snpnet?tab=readme-ov-file)
- Reference: [Vignette of the `snpnet` `R` package](https://github.com/junyangq/snpnet/blob/master/vignettes/vignette.pdf) by Junyang Qian and Trevor Hastie
1. (Already preprocessed) convert genotypes to `.pgen`
2. Train `snpnet()` model, get `fit_snpnet` object in `R`
  - `niter = 100` controls the number of iterations
  - `use.glmnetPlus = TRUE` recommended for faster computation
3. Predict phenotype using `predict_snpnet()`, get `pred_snpnet` object
  - `fit = fit_snpnet` to use the fitted model
  - `new_genotype_file` to test on test set
4. Evaluate $R^2$ in `pred_snpnet$metric`

### [BayesR](https://github.com/KengChiChang/cse-284-prs-comparison/blob/main/notebook/12_test_GCTB.ipynb)  [![image](https://img.shields.io/badge/Colab-F9AB00?style=for-the-badge&logo=googlecolab&color=525252)](https://colab.research.google.com/drive/1WGhJaDTsY1LcJuOcysQHrEndwxssRcwD?usp=sharing)

- Tool: [`GCTB`](https://cnsgenomics.com/software/gctb/#SBayesRTutorial)
- Reference: [Tutorial: Practical 4 Bayesian methods for genomic prediction](https://cnsgenomics.com/data/teaching/GNGWS23/module5/Practical4_Bayes.html) by `GCTB` package maintainer Jian Zeng
1. (Already preprocessed) convert genotypes to `.bed`
2. Train BayesR model
  - `--chain-length 10000` controls total length of Markov Chain
  - `--burn-in 2000` controls number of burn-ins
3. Get posterior effects for each SNP in `.snpRes` format
4. Send `.snpRes` to `plink` for scoring of PRS scores
5. Calculate R^2 against ground truth


## Evaluation
- $R^2$ score on the test set in the context of GWAS and PRS provides a measure of how well a polygenic risk score can predict (true) phenotypes.


## Results on Benchmark

As a proof of concept and testing of data/input requirements for different software packages, we use/transform data from PS4: 
- Genotype `CEU_chr19_normed.{vcf,bed,pgen}` is used as training set
- Genotype `GBR_chr19_normed.{vcf,bed,pgen}` is used as test set
- Phenotype `kgvcf_ldl.{phen,pgeno}` is used as ground truth phenotype


| Method         | Training $R^2$ | Testing $R^2$ |
| -------------- | -------------- | ------------- |
| C+T            | 0.98           | 0.64          |
| BASIL/`snpnet` | 0.99           | 0.65          |
| BayesR/`GCTB`  | 0.99           | 0.89          |

Table above reports the training and test set performance across the three methods. All three methods overfit on their own dataset. C+T and BASIL achieve similar performance on test set, while BayesR performs much better on the test set.


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
- Another complication is that different packages expect different inputs and data format. For example, `snpnet` expects `.pgen` while `GCTB` expects `.bed`. Also, `plink1.9` and `plink2` also has different input expectations and modifiers. This makes data preparation taking a huge chunk of time.
- Also, another challenge lies in determining the optimal settings for phenotype simulation. While we recognize the importance of varying simulation parameters to comprehensively evaluate PRS performance, we face uncertainties regarding the number of settings required to yield sufficient insights. 


## Remaining Work 
- We have completed the essential stages of data preprocessing and splitting for training, validation, and testing dataset. Additionally, our initial testing encompassed the evaluation of three key methods: Clumping + Thresholding (C+T), BASIL, and BayesR, yielding preliminary insights. Moving forward, our focus will shift towards extending our GWAS across a broader dataset. This expansion aims to deepen our understanding of how different phenotype simulations and the application of various PRS methods influence the resulting $R^2$ score. Furthermore, we are committed to enhancing the accessibility of our findings by visualizing results through tables and graphs. This approach will facilitate easier interpretation and comprehension, ensuring that our research contributes effectively to the broader discourse within the field.