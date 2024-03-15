# How Simulation Parameters Impact Polygenic Scoring Accuracy Across Three Methods

<img src="https://github.com/KengChiChang/cse-284-prs-comparison/blob/40b83cdf9b805f753f73a6b8fc9978ecefa33691/figure/result_by_numcausal_method_1.png?raw=true" width="50%"/><img src="https://github.com/KengChiChang/cse-284-prs-comparison/blob/40b83cdf9b805f753f73a6b8fc9978ecefa33691/figure/scores_by_pop_method_1.png?raw=true" width="50%"/>


> Keywords: GWAS (Genome-wide Association Studies), Phenotype Simulation, PRS (Polygenic Risk Score), C+T (Clumping and Thresholding), BASIL (Lasso Regression), and BayesR (Bayesian Regression)

## Group 18
- Keng-Chi Chang
- Po-Chun Wu


## Introduction
Polygenic Risk Scores (PRS) hold promise for personalized medicine but face challenges in prediction accuracy. This study investigates how simulation parameters affect PRS accuracy using Clumping + Thresholding, BASIL (Lasso-based method in `snpnet` package), and BayesR (Bayesian Multiple Regression-based method in `GCTB` package). By manipulating parameters such as heritability, number of causal SNPs, and how minor allel frequency affects heritability, we aim to understand whether certain PRS method might have advantage over others, and under what conditions. Our findings can inform the selection of PRS and phenotype simulation methods for precise disease risk prediction in clinical and research settings, advancing personalized healthcare interventions.


## Methods
- Clumping + Thresholding (C+T) (PLINK 2)
- [Batch Screening Iterative Lasso (BASIL)](https://github.com/rivas-lab/snpnet)
- [Bayesian Multiple Regression (BayesR)](https://cnsgenomics.com/software/gctb/#Overview)


## Hypothesis
1. Since Clumping + Thresholding does not directly model the additive effects of the contribution from multiple SNPs, we expect that C+T would perform the worst when the number of causal SNPs increases.
2. Both BASIL and BayesR should be better at capturing additive effects from many SNPs, while shrinkage-based method might incorrectly shrink some effects to zero.

Consider heritability model of the form [(Speed, Holmes, and Balding 2020)](https://www.nature.com/articles/s41588-020-0600-y)
$$\mathbb{E}\left[h^2_j\right]=\sum_j w_j\left[f_j(1-f_j)\right]^{(1+\alpha)},$$
where $w_j$ is the weight for SNP $j$ and $f_j$ is its minor allele frequency (MAF). For simplicity, we assume that weight $w_j$ is a constant.

3. When $\alpha=-1$, the expected heritability reduces to $\sum_j w_j$.
SNPs have linear effect on heritability, regardless of minor allele frequency. Since this data generating process is closer to Lasso, we expect BASIL to perform better under this regime.
4. When $\alpha=-0.25$, the expected heritability becomes $\sum_j w_j\left[f_j(1-f_j)\right]^{0.75}$. SNPs that are more common would have higher contribution to heritability. We expect BayesR would capture this nonlinear relationship better than Lasso based BASIL.

The current published papers comparing these methods have mixed findings. Most of them either only use simulated genotypes, simulate phenotypes on a narrow demography ([Lloyd-Jones et al. 2019](https://doi.org/10.1038/s41467-019-12653-0); [Qian et al. 2020](https://doi.org/10.1371/journal.pgen.1009141)), or does not systematically evaluate alternative scenarios ([Ding et al. 2023](https://doi.org/10.1038/s41586-023-06079-4)). To our knowledge, this is the first systematic comparison on real genomic data across ancestries/populations.


## Directory for Data Storage and Computation
We are analyzing real genome data from [1000 Genomes Project](https://www.internationalgenome.org) Phase Three call set on GRCh37. Since the size of the data is huge, for the ease of data storage and collaborating, we store our data on Google Drive and use Google Colab for compute (by reading/copying files directly from Google Drive).
- Notes: We have updated the related notebooks and the data of calculation results on GitHub. However, due to the large size of other raw data, we have only placed it on Google Drive.

[Here](https://drive.google.com/drive/folders/1vONriV2u1j2BinWGxUhBwLnFN94u3Mig?usp=drive_link) is the link to our Google Drive project folder. 

GitHub Directory:
```
.
├── data
│   ├── BASIL
│   │   └── chr19_ldl_pheno ← results of BASIL
│   ├── BayesR
│   │   └── chr19_ldl_pheno ← results of BayesR
│   ├── CT
│   │   └── chr19_ldl_pheno ← results of C+T
│   ├── chr19_ldl_pheno ← Phenotype Simulation
│   └── split_ids ← Data splitting (train, validation, test)
│       ├── by_population
│       └── by_superpopulation 
├── figure ← analysis and graphs
└── notebook ← codes and scripts
```

Part of the files in Google Drive Project are organized as follows:

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
    ├── 00_download_1000genome_data_and_simulate_phenotypes.ipynb
    ├── 01_convert_to_pgen.ipynb
    ├── 02_convert_to_bed.ipynb
    ├── 04_split_train_val_test.ipynb
    ├── 05_simulate_phenotype.ipynb
    ├── 06_update_phenotype_format.ipynb
    ├── 07_convert_chr19.ipynb
    ├── 08_LDAK_chr19_pheno_simulation.ipynb
    ├── 09_test_GCTA_pheno_simulation.ipynb
    ├── 11_test_snpnet.ipynb
    ├── 12_test_GCTB.ipynb
    ├── 13_test_C+T.ipynb
    ├── 21_1000g_snpnet.ipynb
    ├── 22_1000g_GCTB.ipynb
    ├── 23_1000g_CT_population.ipynb
    ├── 24_CT_test_PS4.ipynb
    ├── 31_snpnet_superpopulation.ipynb
    ├── 32_GCTB_superpopulation.ipynb
    ├── 33_CT_superpopulation.ipynb
    ├── 34_CT_population.ipynb
    ├── 41_chr19_prs_by_superpopulation_check_progress.ipynb
    ├── 42_chr19_prs_by_superpopulation_analysis.Rmd
    ├── 42_chr19_prs_by_superpopulation_analysis.html
    ├── 42_chr19_prs_by_superpopulation_analysis.md
    └── 43_analysis.ipynb
```
- File naming of `notebook` directory
    - 0x: Data Preprocessing
    - 1x: Model testing
    - 2x: Run models on all autosomes
    - 3x: Run models on Chr 19 only (final setting)
    - 4x: Results and analysis

## Data Preprocessing

### [Genotypes (1000 Genomes)](https://github.com/KengChiChang/cse-284-prs-comparison/blob/main/notebook/00_download_1000genome_data_and_simulate_phenotypes.ipynb) [![image](https://img.shields.io/badge/Colab-F9AB00?style=for-the-badge&logo=googlecolab&color=525252)](https://drive.google.com/file/d/1XIm8UKT1MxdjKfQDr75cxLtmREFvJpvl/view?usp=share_link)
- Goal: get the dataset of all autosomes from 1000 Genomes with quality control
- Tool: [`PLINK 1.9`](https://www.cog-genomics.org/plink/) and [`PLINK 2.0`](https://www.cog-genomics.org/plink/2.0/)
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

<img src="https://github.com/KengChiChang/cse-284-prs-comparison/blob/main/figure/1000_genomes_superpopulation_pie_chart.png?raw=true" width="50%"/>

(Number of Samples from Each Superpopulation on 1000 Genomes Project Phase 3 Dataset)

### [Phenotypes (Simulation)](https://github.com/KengChiChang/cse-284-prs-comparison/blob/main/notebook/05_simulate_phenotype.ipynb) [![image](https://img.shields.io/badge/Colab-F9AB00?style=for-the-badge&logo=googlecolab&color=525252)](https://drive.google.com/file/d/1XQ9W02ikpiQipvJ_-ByaTbNWD2dXpACE/view?usp=share_link)
- Goal: simulate phenotypes with different settings
- Tool: [LDAK](https://dougspeed.com/simulations/) 
- Options
    - `--make-phenos`: perform simulation
    - `--power` ($\alpha$) : specify how predictors are scaled
        - Diffirent from the *power* mentioned in the class that indicates the probability  we can detect an association at the desired significance level, given that there is actually an association. 
    - `--her`: heritability describes how much variation in the phenotype is described by genetics.
    - `--num-phenos`: the number of phenotypes to generate.
    - `--num-causals`: the number of predictors contributing to each phenotype.
    - ` --causals <causalsfile>`: specify which predictors are causal for each phenotype.
    - ` --effects <effectsfile>`: specify the effect sizes in a file. 
- Our simulation settings
    - Number of Phenotypes: 20
    - Power: {-0.25, -1}
        - Notes: GCTA Model with power of -1 vs. Human Default Model with power of -0.25
    - Heritability: {0.1, 0.3, 0.5, 0.7, 0.9}
    - Number of causal SNPs = {1, 10, 100, 250, 512} 
    - Total combinations: 2 * 5 * 5 = 50
    - File naming: `power=x_her=y_num-causals=z.pheno`
    - Columns: FID, IID, PHEN1, PHEN2, PHEN3, ...
- Notes: At first, we didn't specify the `<causalsfile>` and `<effectsfile>`. By default, LDAK will pick causal predictors at random and sample effect sizes from a standard normal distribution. However, the results of PRS using randomly selected causal SNPs and effect sizes are horrible ($R^2$ of testing data close to $0$). Therefore, we adopted the `ldl_prs.txt` data from UCSD CSE 284 PRS Activity, which is the LDL PRS data from this [paper](https://www.nature.com/articles/s41586-021-04064-3) and its publicly available [data](https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/prs_weights/GLGC_2021_ALL_LDL_PRS_weights_PT.txt). Also, since there are 512 SNPs on Chromosome 19 listed in `ldl_prs.txt`, the maximal number of causal SNPs is 512.
- Dataset
    - `data/chr19_ldl_pheno/` 
- Reference: https://dougspeed.com/simulations/

### [Sample Splitting](https://github.com/KengChiChang/cse-284-prs-comparison/blob/main/notebook/04_split_train_val_test.ipynb) [![image](https://img.shields.io/badge/Colab-F9AB00?style=for-the-badge&logo=googlecolab&color=525252)](https://colab.research.google.com/drive/16fbgP6IK38J7wIoPl_Zj-Ts7Wk_rFqX3?usp=share_link)
<!-- - For our study on PRS, we adopted two distinct approaches to split our data into training, validation, and testing sets, ensuring a comprehensive evaluation of our methods.  -->
- We randomly select samples from **each superpopulation** into 70% training, 15% validation, and 15% test sets, while keeping each set has balanced number of samples from each population and sex. 
<!--     - Alternatively, we drew inspiration from Problem Set 4 (PS4), where data was segmented by **population**. This method allowed us to examine the influence of genetic diversity on PRS performance by closely mirroring real-world population structures.  -->
- Dataset
    - `data/split_ids/`
        - `by_population/`
        - `by_superpopulation/`

## Details

$$PRS_i = \sum_{j \in S} \beta_j X_{ij}$$

### [C+T](https://github.com/KengChiChang/cse-284-prs-comparison/blob/main/notebook/33_CT_superpopulation.ipynb) [![image](https://img.shields.io/badge/Colab-F9AB00?style=for-the-badge&logo=googlecolab&color=525252)](https://colab.research.google.com/drive/1qmOKXFYZ5QZ4ZXcYXd8Hxry66RztcvgA?usp=share_link)
- Tool: [`PLINK 2.0`](https://www.cog-genomics.org/plink/2.0/)
- Reference: UCSD CSE 284 Week 6 Lecture Slides, PS3, and PS4
1. Perform a GWAS to estimate per-variant effect sizes ($\beta$'s)
    - `--linear --maf 0.05`
    $$Y = \beta_j X_j + \epsilon$$
2. Perform LD clumping (pruning) to get an independent set of SNPs
    - `--clump gwas.assoc.linear`
    - `--clump-p1 0.0001`: significance threshold for index SNPs   
    - `--clump-r2 0.1`: LD threshold for clumping
    - `–clump-kb 250`:  physical distance threshold for clumping
3. Choose all SNPs with $p$-value less than the threshold $T$
4. Computer PRS and evalueate accuracy ($R^2$)
    - `--score`
5. Computer final PRS with optimal $T$ and evaluate on a separate testing dataset
- Repeat steps 2-4 to find out which $T$ work the best (validation dataset)
- Notes: When running C+T for AMR (superpopulation) dataset from 1000 Genomes, we encountered a problem that `PLINK 2.0` will not run the `--score` command when there are less than 50 samples. Since there are not enough samples for validation data. We choose the best p value threhold $T$ in the training data as the $T$ for testing data. Thus, the performance may be underrated.

### [BASIL](https://github.com/KengChiChang/cse-284-prs-comparison/blob/main/notebook/31_snpnet_superpopulation.ipynb) [![image](https://img.shields.io/badge/Colab-F9AB00?style=for-the-badge&logo=googlecolab&color=525252)](https://colab.research.google.com/drive/1XpSwudemAC_NB_fxiJ_w0RzWly330Aj-?usp=share_link)
$$Minimize \sum_i (y_i - \hat{y}_i)^2 + \lambda \sum_j |\hat{\beta}_j|$$
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

### [BayesR](https://github.com/KengChiChang/cse-284-prs-comparison/blob/main/notebook/32_GCTB_superpopulation.ipynb)  [![image](https://img.shields.io/badge/Colab-F9AB00?style=for-the-badge&logo=googlecolab&color=525252)](https://colab.research.google.com/drive/18Zza9rGNuba7vBDF-in7mumzRYwTbHqk?usp=share_link)

- Tool: [`GCTB`](https://cnsgenomics.com/software/gctb/#SBayesRTutorial)
- Reference: [Tutorial: Practical 4 Bayesian methods for genomic prediction](https://cnsgenomics.com/data/teaching/GNGWS23/module5/Practical4_Bayes.html) by `GCTB` package maintainer Jian Zeng
1. (Already preprocessed) convert genotypes to `.bed`
2. Train BayesR model
  - `--chain-length 10000` controls total length of Markov Chain
  - `--burn-in 2000` controls number of burn-ins
3. Get posterior effects for each SNP in `.snpRes` format
4. Send `.snpRes` to `plink` for scoring of PRS scores
5. Calculate $R^2$ against the ground truth

## Experiments
### All autosomes with Random Effect Sizes
- $R^2$ of testing data close to $0$ 🥲 
- Probably because there are not enough samples to model the complicated relaitonships between variants and phenotypes.
### Chr 19 with Random Effect Sizes
- $R^2$ of testing data close to $0$ 🥲 
- Probably because there are not enough samples to model the complicated relaitonships between variants and phenotypes.
### Chr 19 with Fixed Effect Sizes
- It works 😃
    - Superpopulation: {'AFR', 'AMR', 'EAS', 'EUR', 'SAS'}
    - Power: {-0.25, -1}
    - Heritability: {0.1, 0.3, 0.5, 0.7, 0.9}
    - Number of causal SNPs = {1, 10, 100, 250, 512} 
    - Output files
        - `combined_predict_{superpopulation}_power={power}_her={heritability}_num-causals={number of causal SNPs}_pheno={number of phenotypes}.csv`
            - Columns: IID, Predicted Phenotype, Actual Phenotype, Phenotype Index
        - `combined_result_{superpopulation}_power={power}_her={heritability}_num-causals={number of causal SNPs}_pheno={number of phenotypes}.csv`
            - Columns: Superpopulation,	Power, Heritability, Number of Causal SNPs,	Phenotype Index, Training $R^2$, Validation $R^2$, Testing $R^2$
- Execution Time (on Google Colab)
    - C + T: $\approx$ 6 hours
    - BASIL:
    - BayesR: 

## Evaluation
- $R^2$ (coefficient of determination) on the testing dataset in the context of GWAS and PRS provides a measure of how well a polygenic risk score can predict (true) phenotypes.


## Results on Benchmark

As a proof of concept and testing of data/input requirements for different software packages, we use/transform data from PS4: 
- Genotype `CEU_chr19_normed.{vcf,bed,pgen}` is used as training set
- Genotype `GBR_chr19_normed.{vcf,bed,pgen}` is used as test set
- Phenotype `kgvcf_ldl.{phen,pgeno}` is used as ground truth phenotype


| Method         | Training $R^2$ | Testing $R^2$ |
| -------------- | -------------- | ------------- |
| C+T            | 0.98           | 0.64          |
| BASIL/`snpnet` | 0.99           | 0.65          |
| BayesR/`GCTB`  | 0.99           | <mark>0.89</mark>          |

Table above reports the training and test set performance across the three methods. All three methods overfit on their own dataset. C+T and BASIL achieve similar performance on test set, while BayesR performs much better on the test set.

## Results

- $R^2$ results on training dataset

| Superpopulation | AFR  | AMR  | EAS  | EUR  | SAS  |
| --------------- | ---- | ---- | ---- | ---- | ---- |
| C+T             | <mark>0.90</mark> | <mark>0.71</mark> | <mark>0.89</mark> | <mark>0.92</mark> |<mark>0.93</mark> |
| BASIL/`snpnet`  | 0.59 | 0.43 | 0.84 | 0.51 | 0.48 |
| BayesR/`GCTB`   |<mark>0.92</mark> | <mark>0.93</mark> | <mark>0.90</mark> | <mark>0.91</mark> | <mark>0.92</mark> |

- $R^2$ results on testing dataset

| Superpopulation | AFR | AMR | EAS | EUR | SAS |
| --------------- | --- | --- | --- | --- | --- |
| C+T             |  0.24   |  0.08   |  0.20   |  0.21   |  0.15   |
| BASIL/`snpnet`  |  0.27   |  0.11   | <mark>0.45</mark>    |  0.16   |   0.13  |
| BayesR/`GCTB`   |  <mark>0.33</mark> | <mark>0.23</mark>   |  0.23   |   <mark>0.29</mark>  |  <mark>0.27</mark>   |

## Discussions 
### Challenges and limitations
- `PLINK 1.9` and `PLINK 2.0` have different input formats and modifiers. This makes data preparation and preprocessing taking a huge chunk of time. In this project, we prioritize using `PLINK 2.0` since `PLINK 2.0` runs faster than `PLINK 1.9` in most cases, and only switch to `PLINK 1.9` when `PLINK 2.0` fails to accomplish what we expect or require significant modifications.
- The number of samples of each superpopulation from 1000 Genomes Project vary a lot. Therefore, the comparison results across superpopulations may not be representative.

### Future Work
- Different distributions of effect sizes and causal SNPs from different traits.

## Related Work
- Frank Dudbridge (2013) provided a comprehensive overview of the power and predictive accuracy of PRS, highlighting the statistical principles underpinning these methods and their implications for disease risk prediction. 
- The 1000 Genomes Project Consortium (2015) provided a landmark achievement in human genetics, offering a global reference for human genetic variation. This comprehensive dataset has been instrumental in PRS research, allowing for the assessment of methodological performance across diverse genetic backgrounds.
- Vilhjálmsson et al. (2015) introduced a significant advancement with LDPred, a Bayesian method considering linkage disequilibrium and polygenic traits. This approach marked a step forward in PRS accuracy, emphasizing the complexity of genetic architecture in disease prediction.
- Khera et al. (2018) showcased the practical application of PRS, particularly in predicting coronary artery disease, underscoring the importance of large reference datasets for improving prediction accuracy. Their work emphasized the clinical utility of PRS, driving further research into methodological improvements.
- Lloyd-Jones et al. (2019) refined the Bayesian approach to PRS with BayesR, demonstrating the method's capability to enhance predictive power by incorporating SNP effect sizes more accurately. This study contributed to the growing preference for Bayesian models in complex disease prediction.
- Ge et al. (2019) evaluated various PRS methods, revealing performance variability across diseases and traits. Their findings stressed the necessity for flexible PRS methodologies adaptable to different genetic architectures.
- Qian et al. (2020) presented a scalable framework for Lasso-based sparse regression with BASIL in the snpnet package, efficiently handling ultra-high-dimensional genetic data. This approach represents a leap forward in managing large-scale datasets for PRS calculation.
- Yi Ding (2023) brought critical attention to the variability of polygenic scoring accuracy across different genetic ancestries. Their results highlight the need to move away from discrete genetic ancestry clusters towards the continuum of genetic ancestries when considering PGSs.


## References
- Dudbridge F. (2013) Power and Predictive Accuracy of Polygenic Risk Scores. PLOS Genetics 9(3): e1003348. https://doi.org/10.1371/journal.pgen.1003348
- The 1000 Genomes Project Consortium. (2015) A global reference for human genetic variation. Nature 526, 68–74 https://doi.org/10.1038/nature15393
- Vilhjálmsson, B. J., Yang, J., Finucane, H. K. et al. (2015) Modeling Linkage Disequilibrium Increases Accuracy of Polygenic Risk Scores. American journal of human genetics. 97(4), 576–592. https://doi.org/10.1016/j.ajhg.2015.09.001
- Khera, A.V., Chaffin, M., Aragam, K.G. et al. (2018) Genome-wide polygenic scores for common diseases identify individuals with risk equivalent to monogenic mutations. Nat Genet 50, 1219–1224. https://doi.org/10.1038/s41588-018-0183-z
- Lloyd-Jones, L.R., Zeng, J., Sidorenko, J. et al. (2019) Improved polygenic prediction by Bayesian multiple regression on summary statistics. Nat Commun 10, 5086. https://doi.org/10.1038/s41467-019-12653-0
- Ge, T., Chen, CY., Ni, Y. et al. (2019) Polygenic prediction via Bayesian regression and continuous shrinkage priors. Nat Commun 10, 1776. https://doi.org/10.1038/s41467-019-09718-5
- Speed, D., Holmes, J. & Balding, D.J. (2020) Evaluating and improving heritability models using summary statistics. Nat Genet 52, 458–462. https://doi.org/10.1038/s41588-020-0600-y
- Qian J, Tanigawa Y, Du W, Aguirre M, Chang C, Tibshirani R, et al. (2020) A fast and scalable framework for large-scale and ultrahigh-dimensional sparse regression with application to the UK Biobank. PLoS Genet 16(10): e1009141. https://doi.org/10.1371/journal.pgen.1009141
- Ding, Y., Hou, K., Xu, Z. et al. (2023) Polygenic scoring accuracy varies across the genetic ancestry continuum. Nature 618, 774–781. https://doi.org/10.1038/s41586-023-06079-4
- [UCSD CSE 284 Lecture Slides and Problem Sets](https://canvas.ucsd.edu/courses/53280)
---

<!-- ## For Peer Review
### Challenges
- In tackling our study on PRS, we ran into a common difficulty: the massive size of the 1000 Genomes dataset. With limited time and computing power, analyzing the entire dataset isn't feasible. So, we're taking a practical approach. We'll divide the data by super-population and population categories, much like we did in Problem Set 4 (PS4). This way, we can still conduct our GWAS effectively. By focusing on specific groups, we aim to strike a balance between thorough analysis and resource constraints. 
- Another complication is that different packages expect different inputs and data format. For example, `snpnet` expects `.pgen` while `GCTB` expects `.bed`. Also, `plink1.9` and `plink2` also has different input expectations and modifiers. This makes data preparation taking a huge chunk of time.
- Also, another challenge lies in determining the optimal settings for phenotype simulation. While we recognize the importance of varying simulation parameters to comprehensively evaluate PRS performance, we face uncertainties regarding the number of settings required to yield sufficient insights. 


### Remaining Work 
- We have completed the essential stages of data preprocessing and splitting for training, validation, and testing dataset. Additionally, our initial testing encompassed the evaluation of three key methods: Clumping + Thresholding (C+T), BASIL, and BayesR, yielding preliminary insights. Moving forward, our focus will shift towards extending our GWAS across a broader dataset. This expansion aims to deepen our understanding of how different phenotype simulations and the application of various PRS methods influence the resulting $R^2$ score. Furthermore, we are committed to enhancing the accessibility of our findings by visualizing results through tables and graphs. This approach will facilitate easier interpretation and comprehension, ensuring that our research contributes effectively to the broader discourse within the field. -->
