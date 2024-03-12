#!/bin/bash
# Script for C+T Analysis

# Filter association results for training set
cat train.assoc.linear | awk '( (($7<10) && ($7!="NA")) || ($1~/CHR/))' > train.assoc.linear.filt

# Generate range list for clumping thresholds
echo "0.000001 0 0.000001" > range_list 
echo "0.00001 0 0.00001" >> range_list 
echo "0.0001 0 0.0001" >> range_list 
echo "0.001 0 0.001" >> range_list 
echo "0.05 0 0.05" >> range_list
echo "0.1 0 0.1" >> range_list

# Perform clumping for training set
plink \
    --vcf train_normed.vcf.gz \
    --clump-p1 1 \
    --clump-r2 0.1 \
    --clump-kb 250 \
    --clump train.assoc.linear.filt \
    --out train

# Extract clumped SNPs and valid SNPs for validation set
awk 'NR!=1{print $3}' train.clumped > train.valid.snp
cat train.clumped | awk '{print $3 "\t" $5}' > train.pvals

# Calculate PRS for validation set
plink \
    --vcf val_normed.vcf.gz \
    --score train.assoc.linear.filt 2 4 7 header sum \
    --q-score-range range_list train.pvals \
    --extract train.valid.snp \
    --out train_val

# Calculate PRS for test set
plink \
    --vcf test_normed.vcf.gz \
    --score train.assoc.linear.filt 2 4 7 header sum \
    --q-score-range range_list train.pvals \
    --extract train.valid.snp \
    --out train_test
