
# Real Data Example

We show Step 6 and 7 as in [README](../README.md) when applying TAGC and WAGC to real data from the UKB study. Here the bfiles used for LD matrix moments computation are unimputed genotype bfiles of the UKB study.

Both examples below use the trait of weight (bioimpedance measure) (UKB Data Field: 23098).

## Within-Ancestry GC (WAGC)
The population is white (regardless of British or non-British). The GWAS summary statistics of trait A of population I is generated using UKB white British subjects, and the individual level observations of trait B of population II is for UKB white non-British subjects. Population I is considered the same as population II, thus we do not distinguish between them in the code of WAGC.

Code below is run in R. Trait A and B are the same trait, weight (bioimpedance measure) (UKB Data Field: 23098).
```{bash}
library(TAGC)

# LD moment, output of LD_block_moment_all_WAGC
LD_mom = c(461466, 2122484, 110792909)

# heritability estimate of trait A of population
h2 = 0.42449

# sample size of GWAS of trait A of population
n = 344046

# correction factor of TAGC
tagc_bias = bias_factor(LD_mom[1], LD_mom[2], LD_mom[3], h2, h2, n)
tagc_bias
# [1] 5.324878

# TAGC point estimate
covar_path = 'data/test_output/covar/ukb_500k_unrelated_british_dec01_2019_gcta.covar.qcovar'
pheno_path = 'data/ukb_test_ukww_10k_unrelated_Aug23_2023.phenb'
prs_path = 'data/test_output/prs/fastgwa_unimputed_pheno10_ukbo49k.profile'
pheno_col_num = 11 # the 11th column in pheno_path contains phenotype weight

rv = TAGC_gc_se_CI(TAGC_bias = tagc_bias,
                   covar_path = covar_path,
                   pheno_path = pheno_path,
                   pheno_col_num = pheno_col_num,
                   prs_path = prs_path,
                   prs_col_num = 6, # the 6th column in prs_path contains PRS
                   number_bs_samples=500, # 500 bootstrap samples for bootstrap SE
                   bs_sample_size=1500, # each bootstrap sameple has size 1500
                   CI_conf_level=0.95 # 95% confidence interval 
                   )

rv
#         P_value       GC0         SE0        GC    GC_SE    GC_CI1   GC_CI2
# 1 3.034978e-134 0.1604448 0.006453209 0.8543491 0.117692 0.5886357 1.066361
#   N_obs
# 1 18486
```



## Trans-Ancestry GC (TAGC)
Population I is white British, and population II is Asian. Trait A and B are the same trait, weight (bioimpedance measure) (UKB Data Field: 23098).

```{bash}
library(TAGC)

# LD moment, output of LD_block_moment_all_TAGC
LD_mom = c(461466, 1464320, 38452232)

# heritability estimate of trait A of population I (the trait without individual level observations for popn II)
h2a = 0.42449

# heritability estimate of trait B of population II (the trait with individual level observations for popn II)
h2b = 0.634465

# sample size of GWAS of trait A of population I
n = 344046

# correction factor of TAGC
tagc_bias = bias_factor(LD_mom[1], LD_mom[2], LD_mom[3], h2a, h2b, n)
tagc_bias
# [1] 3.822642


# TAGC point estimate
covar_path = 'data/test_output/covar/ukb_500k_unrelated_british_dec01_2019_gcta.covar.qcovar'
pheno_path = 'data/test_output/pheno/ukb_test_52phenotypes_asian_9k_unrelated.phenb'
prs_path = 'data/test_output/prs/fastgwa_unimputed_pheno6_ukbo49k.profile'
pheno_col_num = 11 # the 11th column in pheno_path contains phenotype weight

rv = TAGC_gc_se_CI(TAGC_bias = tagc_bias,
                   covar_path = covar_path,
                   pheno_path = pheno_path,
                   pheno_col_num = pheno_col_num,
                   prs_path = prs_path,
                   prs_col_num = 6, # the 6th column in prs_path contains PRS
                   number_bs_samples=500, # 500 bootstrap samples for bootstrap SE
                   bs_sample_size=1500, # each bootstrap sameple has size 1500
                   CI_conf_level=0.95 # 95% confidence interval 
                   )

rv
#       P_value       GC0        SE0       GC      GC_SE    GC_CI1    GC_CI2
# 1 5.40941e-87 0.1909251 0.00953748 0.729838 0.08629865 0.5654883 0.8947309
#   N_obs
# 1  7920
```

## Interpret the output
- `P_value` is the p-value of the linear regression, where the population II trait B is regressed on the the population II trait A PRS.
- `GC0` is the original correlation between trait B observation and trait A PRS in the linear regression.
- `SE0` is the standard error of the linear regression slope.
- `GC` is the trans-ancestry genetic correlation (TAGC) estimate of heel bone mineral density.
- `GC_SE` is the standard error of the TAGC estimate based on the 500 bootstrap samples.
- `GC_CI1` and `GC_CI2` is the quantile-based confidence interval using the bootstrap samples. Here `GC_CI1` is the 2.5% quantile and `GC_CI2` is the 97.5% quantile of the 500 TAGC estimates of the 500 bootstrap samples. 
- `N_obs` is the number of complete observations in the linear regression.





