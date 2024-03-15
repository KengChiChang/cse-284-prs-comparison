library("tidyverse")
library("glue")
library("snpnet")

LOG_ID           = '2024-3-14-03' # nolint: assignment_linter, single_quotes_linter.
POP_LIST         = c('EUR')
POWER_LIST       = c(-1, -0.25)
HER_LIST         = c(0.3, 0.7)
NUM_CAUSALS_LIST = c(512, 100, 10, 1)
PHENO_RANGE      = '1-20'

PARENT_DIR        = "/Users/kengchichang/Library/CloudStorage/Dropbox/GitHub/cse-284-prs-comparison"
CHR19             = glue("{PARENT_DIR}/data/1000g_by_chrom/chr19")
SPLIT_DIR         = glue("{PARENT_DIR}/data/split_ids/by_superpopulation")
GENO_DIR          = glue("{PARENT_DIR}/data/chr19_by_superpopulation")
PHEN_OLD_DIR      = glue("{PARENT_DIR}/data/chr19_ldl_pheno")
PHEN_DIR          = glue("{PARENT_DIR}/data/chr19_ldl_pheno_corrected")
MODEL_DIR         = glue("{PARENT_DIR}/data/BASIL/chr19_ldl_pheno/models3")
RESULT_DIR        = glue("{PARENT_DIR}/data/BASIL/chr19_ldl_pheno")
# GDRIVE_RESULT_DIR = 'data/drive/MyDrive/CSE-284-Final-Project/data/BASIL/chr19_ldl_pheno'
pheno_col_names = paste0('PHENO', c(1:100))

configs = list(
  save = FALSE,  # save intermediate results per iteration (default FALSE)
  nCores = 2,  # number of cores available (default 1)
  niter = 30,  # max number of iterations (default 50)
  # prevIter = 15,  # if we want to start from some iteration saved in results.dir
  use.glmnetPlus = TRUE,  # recommended for faster computation
  verbose = FALSE,
  # early.stopping = FALSE,  # whether to stop based on validation performance (default TRUE)
  early.stopping = TRUE,
  # stopping.lag = 2,
  plink2.path = "/Users/kengchichang/bin/plink2",   # path to plink2 program
  zstdcat.path = "/opt/homebrew/bin/zstdcat",  # path to zstdcat program
  # meta.dir = "meta/",
  # save.dir = "save/",
  results.dir = MODEL_DIR
)

for (POP in POP_LIST) {
    print(glue("-------------- Start: {POP} --------------"))
    # combine train test split
    pop_samples = rbind(
        read_delim(glue('{SPLIT_DIR}/{POP}_train.txt'),
                    delim=' ', col_names=c('FID', 'IID'), skip=1, col_types='cc') |>
                    mutate(split='train'),
        read_delim(glue('{SPLIT_DIR}/{POP}_val.txt'),
                    delim=' ', col_names=c('FID', 'IID'), skip=1, col_types='cc') |>
                    mutate(split='val'),
        read_delim(glue('{SPLIT_DIR}/{POP}_test.txt'),
                    delim=' ', col_names=c('FID', 'IID'), skip=1, col_types='cc') |>
                    mutate(split='test')) |>
        select(-FID)
    write_csv(pop_samples['IID'], glue("{SPLIT_DIR}/{POP}_all.txt"))
    # generate pgen for superpopulation
    system(glue(
    "plink2 \\
        --make-pgen vzs 'psam-cols=+fid' \\
        --bfile {CHR19} \\
        --const-fid 0 \\
        --keep {SPLIT_DIR}/{POP}_all.txt \\
        --out {GENO_DIR}/{POP}_all \\
        --silent \\
        &> /dev/null"
        ))
    for (NUM_CAUSALS in NUM_CAUSALS_LIST) {
        for (HER in HER_LIST) {
            for (POWER in POWER_LIST) {
                params = glue("power={POWER}_her={HER}_num-causals={NUM_CAUSALS}")
                pheno_file = glue("{PHEN_OLD_DIR}/{params}.pheno")
                # transform phenotype to suit snpnet inpu format
                pheno_df = read_delim(pheno_file, delim=' ', skip=1,
                                      col_names=c('FID', 'IID', pheno_col_names), 
                                      col_types=paste0('cc', paste0(rep('n', 100), collapse=''))) |>
                                      select_if(~!(all(is.na(.)) | all(. == "")))
                pheno_df = pheno_df |>
                    inner_join(pop_samples, by=c('IID')) |>
                    select(all_of(c("FID", "IID", "split", pheno_col_names)))
                write_csv(pheno_df, glue("{PHEN_DIR}/{POP}_{params}.phe"))

                PHENO_INDEX_START        = as.integer(strsplit(PHENO_RANGE, '-')[[1]][1])
                PHENO_INDEX_END          = as.integer(strsplit(PHENO_RANGE, '-')[[1]][2])
                combined_result_file     = glue('{RESULT_DIR}/combined_result_{POP}_{params}_pheno={PHENO_INDEX_START}-{PHENO_INDEX_END}.csv')
                combined_prediction_file = glue('{RESULT_DIR}/combined_predict_{POP}_{params}_pheno={PHENO_INDEX_START}-{PHENO_INDEX_END}.csv')
                combined_result_data     = list()
                combined_prediction_data = list()
                # get models, one for each phenotype
                for (PHENO_INDEX in PHENO_INDEX_START:PHENO_INDEX_END) {
                    print(glue("---- Train: pop={POP}, power={POWER}, her={HER}, num-causals={NUM_CAUSALS}, pheno={PHENO_INDEX}"))
                    # train snpnet
                    fit_snpnet = snpnet(
                        genotype.pfile = glue('{GENO_DIR}/{POP}_all'),
                        phenotype.file = glue('{PHEN_DIR}/{POP}_{params}.phe'),
                        phenotype = glue('PHENO{PHENO_INDEX}'),
                        family = "gaussian",
                        split.col = "split", 
                        mem = 48000,  # amount of memory available (MB), recommended
                        configs = configs
                    )
                    # get predictions
                    pred_snpnet = predict_snpnet(
                        fit = fit_snpnet,
                        new_genotype_file = glue('{GENO_DIR}/{POP}_all'),
                        new_phenotype_file = glue('{PHEN_DIR}/{POP}_{params}.phe'),
                        phenotype = glue('PHENO{PHENO_INDEX}'),
                        split_col = "split", 
                        split_name = c("train", "val", "test"),
                        configs = configs
                        )
                    train_score = tail(as.vector(unlist(pred_snpnet$metric$train)), n=1)
                    print(glue("-- Train score {train_score}: pop={POP}, power={POWER}, her={HER}, num-causals={NUM_CAUSALS}, pheno={PHENO_INDEX}"))
                    val_score = tail(as.vector(unlist(pred_snpnet$metric$val)), n=1)
                    test_score = tail(as.vector(unlist(pred_snpnet$metric$test)), n=1)
                    print(glue("-- Test score {test_score}: pop={POP}, power={POWER}, her={HER}, num-causals={NUM_CAUSALS}, pheno={PHENO_INDEX}"))
                    # resulst to save
                    result_row = data.frame(
                        "population/superpopulation"=POP,
                        "power"=POWER,
                        "her"=HER,
                        "num_causals"=NUM_CAUSALS,
                        "pheno_num"=PHENO_INDEX,
                        "train"=train_score,
                        "val"=val_score,
                        "test"=test_score)
                    combined_result_data[[PHENO_INDEX]] = result_row
                    # prediction to save
                    m1 = do.call(rbind, strsplit(rownames(pred_snpnet$prediction$test), '_'))
                    IID = m1[, 2]
                    predict = as.vector(pred_snpnet$prediction$test[,ncol(pred_snpnet$prediction$test)])
                    actual = pred_snpnet$response$test
                    predict_df = data.frame(IID, predict, actual,
                                            pheno_num=PHENO_INDEX)
                    combined_prediction_data[[PHENO_INDEX]] = predict_df
                    
                }
                
                # Write combined prediction data to file
                combined_prediction_df = bind_rows(combined_prediction_data)
                write_csv(combined_prediction_df, combined_prediction_file)
                print(glue("-- Combined prediction data has been written to: {combined_prediction_file}"))

                # Write combined result data to file
                combined_result_df = bind_rows(combined_result_data)
                write_csv(combined_result_df, combined_result_file)
                print(glue("-- Combined result data has been written to: {combined_result_file}"))

                # system(glue("cp {combined_result_file} {GDRIVE_RESULT_DIR}"))
                # system(glue("cp {combined_prediction_file} {GDRIVE_RESULT_DIR}"))
                system(glue("rm -rf {MODEL_DIR}/*"))

            }
        }
    }
    system(glue("rm -rf {GENO_DIR}/{POP}*"))
}