
# Simulated Data Example

We show Step 6 and 7 as in [README](../README.md) when applying TAGC and WAGC to simulated data, where the genotype file is from the UKB study and contains HapMap3 SNPs, and the phenotype is simulated using [simu](https://github.com/precimed/simu). The true heritability of both traits are 0.6, and the true genetic correlation is 0.5. No covariates are included in the phenotype generation, thus covariates are not input to `TAGC_gc_se_CI` in both WAGC and TACG examples here.


## Within-Ancestry GC (WAGC)
Both population I and II are European.

```{bash}
library(TAGC)

LD_mom = c(990376, 19765193, 921944801)
h2 = 0.6
n = 340000
tagc_bias = bias_factor(LD_mom[1], LD_mom[2], LD_mom[3], h2, h2, n)
tagc_bias
# [1] 2.073844


rv = TAGC_gc_se_CI(TAGC_bias = tagc_bias,
    pheno_path = "pheno_3_train_n480k_ldsc_h05_c01_m05.pheno_1000",
    pheno_col_num = 4, 
    prs_path = paste0("fastgwa_ldsc_c01_m05_pheno3_ukbb10k.profile_1000"),
    prs_col_num = 6, number_bs_samples=500, bs_sample_size=1500, CI_conf_level=0.95)
rv
#        P_value       GC0        SE0        GC      GC_SE    GC_CI1    GC_CI2
# 1 2.643551e-19 0.2787554 0.03039973 0.5780952 0.04776215 0.4793125 0.6730096
#   N_obs
# 1  1000
```



## Trans-Ancestry GC (TAGC)

Population I is European and population II is Asian.

```{bash}
library(TAGC)

LD_mom = c(990376, 16345610, 736704003)
h2 = 0.6 # true heritability
n = 340000 # population I GWAS sample size
tagc_bias = bias_factor(LD_mom[1], LD_mom[2], LD_mom[3], h2, h2, n)
tagc_bias
# [1] 2.245349

rv = TAGC_gc_se_CI(TAGC_bias = tagc_bias,
    pheno_path = "pheno_3_train_n480k_ldsc_h05_c01_m05.pheno_946",
    pheno_col_num = 4, 
    prs_path = "fastgwa_ldsc_c01_m05_pheno3_ukbo49k.profile_946",
    prs_col_num = 6, number_bs_samples=500, bs_sample_size=1500, CI_conf_level=0.95)

rv
#        P_value       GC0        SE0        GC      GC_SE    GC_CI1    GC_CI2
# 1 5.167309e-13 0.2318691 0.03166022 0.5206271 0.04842782 0.4327686 0.6206826
#   N_obs
# 1   946
```




## Interpret the output
- `P_value` is the p-value of the linear regression, where the population II trait B is regressed on the the population II trait A PRS.
- `GC0` is the original correlation between trait B observation and trait A PRS in the linear regression.
- `SE0` is the standard error of the linear regression slope.
- `GC` is the trans-ancestry genetic correlation (TAGC) estimate of heel bone mineral density.
- `GC_SE` is the standard error of the TAGC estimate based on the 500 bootstrap samples.
- `GC_CI1` and `GC_CI2` is the quantile-based confidence interval using the bootstrap samples. Here `GC_CI1` is the 2.5% quantile and `GC_CI2` is the 97.5% quantile of the 500 TAGC estimates of the 500 bootstrap samples. 
- `N_obs` is the number of complete observations in the linear regression.


