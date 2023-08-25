# TAGC: Estimating trans-ancestry genetic correlation with unbalanced data resources

TAGC is a novel estimation method of using genetic predicted observations to estimate trans-ancestry genetic correlations, which describes how genetic architecture of complex traits varies among populations, in genomewide association studies (GWAS).

Our new estimator corrects for prediction errors caused by high-dimensional weak GWAS signals, while addressing the heterogeneity of GWAS data across ethnicities, such as linkage disequilibrium (LD) differences, which can lead to biased findings in homogeneity-agnostic analyses. Moreover, our estimator only requires one population to have a large GWAS sample size, and the second population can only have a much smaller number of participants (for example, hundreds). It is designed to specifically address the unbalanced data resources such that the GWAS sample size for European populations is usually larger than that of non-European ancestry groups.


## Installation
```{bash}
# dependency
install.packages("BEDMatrix")
install.packages("data.table")
install.packages("DescTools")
install.packages("psych")
# install
library("devtools")
install_github("FSSKM/TAGC_review", ref="main")
```


## Output
TAGC generates a point estimate and corresponding confidence intervals for the trans-ancestry genetic correlation between trait A of population I and trait B of population II. 

It is assumed that GWAS summary statistics are available for trait A of population I, and individual level observations are available for trait B of population II. Individual level observations of trait A of population II is not needed, and the TAGC uses the PRS of trait A for population II instead of individual level observations.

Population II typically is the smaller population, and population I is typically the larger population with higher quality GWAS. 

The TAGC estimate will be a corrected estimate based on the naïve correlation between polygenic risk score (PRS) of trait A of population II and the individual observation of trait B of population II, and TAGC aymptotically unbiasedly estimates the Pearson correlation of population-specific genetic effects.


## Input
TAGC requires the following input and outputs an estimate of trans-ancestry genetic correlations, and resampling-based standard error.
- LD block matrix moments
    * User can use our pre-computed LD block matrix moments and need to match the GWAS summary statistics with our SNP ID list accordingly, or
    * user can compute the LD block matrix moments of their own bfile using our TAGC LD moment function.
- LD block boundaries (if user computes their own LD block matrix moments)
    * User can use our provided [LD block boundaries](TAGC_data/LD_boundary), or
    * user can use their own LD block boundaries.
- Polygenic risk score (PRS) generated using GWAS summary statistics (GWAS as a marginal estimator of SNP effects) of one trait (PRS of trait A of population II).
- Heritability estimate of phenotypes
    * We recommend [GCTA](https://yanglab.westlake.edu.cn/software/gcta/#GREMLanalysis) if individual level data is available, and [LDSC](https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation#estimating-heritability-genetic-correlation-and-the-ld-score-regression-intercept) if only summary-level data is available. A brief tutorial for both methods are provided below.
- Individual level observation of the other trait (trait B of population II).
- (Optional) Individual level observation of the covaraites corresponding to trait B of population II.


## Input Data Format
Both the phenotype file of trait B of population II, the covariates of population II are assumed to have the first column containing subject IDs. The covariates file should only have the first column containing subject IDs.



## Additional Data
We provide the following data for Step 1 (matching GWAS summary statistics), and Step 7 (TAGC/WAGC estimate):
- Lists of SNPs that were used in computation of our pre-computed LD matrix moments
    * [HapMap3](TAGC_data/hm3_snp.id.list)
    * [Unimputed UKB bfile](TAGC_data/UKB_unimputed_snp.id.list)
- Pre-computed LD matrix moments
    * LD matrix moment for WAGC/TAGC input:
        * HapMap3: [WAGC](TAGC_data/LD_moments/hm3/UKB_eur/WAGC_moment), [TAGC](TAGC_data/LD_moments/hm3/UKB_eur_asian/TAGC_moment)
        * Unimputed UKB bfile: [WAGC](TAGC_data/LD_moments/unimputed/UKB_eur/WAGC_moment), [TAGC](TAGC_data/LD_moments/unimputed/UKB_eur_asian/TAGC_moment)
    * LD matrix moment for each block (aggregates to the above input using our R function `LD_block_moment_all_WAGC`/`LD_block_moment_all_TAGC`):
        * HapMap3: [WAGC](TAGC_data/LD_moments/hm3/UKB_eur/each_block), [TAGC](TAGC_data/LD_moments/hm3/UKB_eur_asian/each_block)
        * Unimputed UKB bfile: [WAGC](TAGC_data/LD_moments/unimputed/UKB_eur/each_block), [TAGC](TAGC_data/LD_moments/unimputed/UKB_eur_asian/each_block)

We also provide the LD matrix used in our simulation analysis in our paper, which is a $10000 \times 10000$ block of the HapMap3 LD of Chromosome 22 of European population and that of Eastern Asian population. 
- EUR: [dropbox link](https://www.dropbox.com/scl/fi/8him3d4228sz42tamjpwi/LD-EUR-CHR22-hapmap3.csv?rlkey=gunrhrz92hhmhrkkpeh4dc2xv&dl=0)
- EAS: [dropbox link](https://www.dropbox.com/scl/fi/xjck63iguxjgad7ctptls/LD-EAS-CHR22-hapmap3.csv?rlkey=jec0rcy22nxj1g7tzs30rbnuy&dl=0)




## Tutorial
We will walk through a simulated data example and a real data example, where the former uses a simulated phenotype with known heritability, and the latter uses a real phenotype with unknown heritability.

We first demonstrate the pre-processing steps that are required for both examples.

### Step 1: Match GWAS summary statistics with our snp list (necessary if using our pre-computed LD block matrix moments)
User can skip this step if computing their own LD block matrix moments. Otherwise, use the code below to remove SNPs that are not in our provided SNP list:
```{bash}
awk 'NR==FNR { ids[$1]=1; next } FNR==1 || $2 in ids' ${TAGC_snp_list} ${input_gwas_path} > ${output_gwas_path}
```
where `${TAGC_snp_list}` is the list of SNPs included in our bfile and used in our LD block matrix moments, available with two versions (see Additional Data section above). `${input_gwas_path}` is the path to user's GWAS summary statistics file, and `${output_gwas_path}` is the desired output path.


### Step 2: Confounding check and correction using LDSC (optional)
User can check for confounding in the GWAS summary statistics using [LDSC](https://github.com/bulik/ldsc) and correct the confounding using LDSC output.
```{bash}
# run LDSC on the GWAS summary statistics
munge_sumstats.py \
--${gwas_path} \
--N ${gwas_sample_size} \
--out ${ldsc_output} \
--merge-alleles w_hm3.snplist
```
As a rule of thumb, if the `Lambda GC` in `${ldsc_output}.log` is larger than 1.1, then the user can consider correcting the confounding by dividing the effect size in the GWAS summary statistics by `Lambda GC`:
```{bash}
awk -v lambda_gc="${lambda_gc}" -v k="$k" '{$k = $k / lambda_gc} 1' ${gwas_path} > ${gwas_by_lambda_path}
```
where the `k`-th column in the GWAS summary statistics file (`${gwas_path}`) is the effect size (beta), and `{lambda_gc}` is the `Lambda GC` in `${ldsc_output}.log`. 


### Step 3: PRS generation
Polygenic risk score using effect size in GWAS summary statistics is one of the input required by TAGC. Run the PLINK command below to generate the PRS. See more details at https://www.cog-genomics.org/plink/1.9/score.
```{bash}
${plink_1.9_path} \
  --bfile ${bfile_path}
  --out ${prs_output_filename}
  --score ${gwas_path}
```
where the path to the bfile is `${bfile_path}`, and `${gwas_path}` is the path to the GWAS summary statistics that is processed in Step 1 and 2.


### Step 4: Heritability estimation for both traits
#### GCTA (individual level data)
See the [GCTA GREML tutorial](https://yanglab.westlake.edu.cn/software/gcta/#GREMLanalysis) for more details. Below is a sample command.
```{bash}
gcta64 \
    --grm ${your_grm} \
    --reml \
    --pheno ${your_phenotype} \
    --thread-num 8 \
    --mpheno ${your_pheno_column} \
    --covar ${your_categorical_covar} \
    --qcovar ${your_quantitative_covar} \
    --reml-lrt 1 \
    --out ${GCTA_output_path}
```
where `${your_pheno_column}` is the column number of the phenotype in the phenotype file `${your_phenotype}`, and the first two columns of `${your_phenotype}` should be FID and IID. Similarly for `${your_categorical_covar}` and `${your_quantitative_covar}`.

#### LDSC (GWAS summary statistics)
See the [LDSC tutorial](https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation#estimating-heritability-genetic-correlation-and-the-ld-score-regression-intercept) for more details. Below is a sample command.
```{bash}
ldsc.py \
--rg ${popn1_trait}.sumstats.gz,${popn2_trait}.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/
--out ${LDSC_output_path} \
```
where the summary statistics of the two traits of interest are `{popn1_trait}.sumstats.gz` and `${popn2_trait}.sumstats.gz`. The heritability estimation will be in `${LDSC_output_path}.log`.



We demonstrate Step 5 and 6 in detail using 1000 Genome bfiles [here](docs/LD_block_matrix_moment_example.md).

### Step 5: LD block generation (if not using our pre-computed LD block matrix moments)
It is recommended to use slurm to run the LD block generation in parallel. For example, if using our provided LD block matrix boundary (`TAGC_data/LD_boundary`), the user can submit 253 jobs as in the example command below:

```{bash}
for i in `seq 1 253`
do
sbatch LD_block_bfile.sh $i
done
```
where the content of `LD_block_bfile.sh` is
```{bash}
# SBATCH header

i=${1}
${plink_1.9_path} \
    --bfile ${bfile_path} \
    --extract TAGC_data/LD_boundary/file${i} \ # or replace with your own LD boundary file
    --make-bed \
    --out ${LD_block_bfile_output_path}/${LD_block_bfile_output_filename} \
    --range
```
User can also use their own LD block boundary.


Step 6 and 7 below depend on whether the user considers population I different from population II. We show detailed steps


### Step 6: LD block matrix moment computation (if not using our pre-computed LD block matrix moments)
If the user has a large bfile of many SNPs (e.g. multiple millions), it is recommended to use slurm to run the LD block matrix moment computation in parallel. For example, if the user used our provided LD block matrix boundary (TAGC_data/LD_boundary/file{1..253}) in Step 5, then the user can again submit 253 jobs as in the example command below:
```{bash}
for i in `seq 1 253`
do
sbatch LD_block_moment.sh $i
done
```
where `LD_block_moment.sh` calls the `TAGC_LD_block_moment` or `WAGC_LD_block_moment` function from the `TAGC` R package. If population I is different from population II (e.g. white vs African), then the user should use `TAGC_LD_block_moment` ('**T**rans-**A**ncestry **G**enetic **C**orrelation'). Otherwise, use `WAGC_LD_block_moment` ('**W**ithin-**A**ncestry **G**enetic **C**orrelation'). 

We provide a detailed [example](docs/LD_block_matrix_moment_example.md) using the publicly available 1000 Genome data. 



### Step 7: TAGC (Trans-Ancestry) and WAGC (Within-Ancestry) estimate, SE, and CI
See 'Simulated data example' and 'Real data example' below. Both simulated data example and real data example show the trans-ancestry and within-ancestry genetic correlation estimate.


## Simulated data example
We provide two detailed examples using simulated data [here](docs/simulated_data_example.md), which is aligned with the simulation analysis in our TAGC paper.


## Real data example
We provide two detailed examples using real data [here](docs/real_data_example.md), which is aligned with the real data analysis in our TAGC paper.





## FAQ
- **Q.** Why is your TAGC point estimate sometimes outside the $[-1,1]$ interval?

- **A.** In real data applications, this 'out-of-boundary' phenomenon might occur when the true genetic correlation is near the boundary of the parameter space. For example, if a trait exhibits significant genetic similarity between two populations, the underlying true genetic correlation would be near to, but less than, 1. Owing to sample variation and the inherent uncertainty in the estimation process, some point estimates may exceed 1 in real-world applications. If we enforce a strict boundary to restrict the estimates to less than 1, we risk introducing bias and potentially underestimating the true genetic correlation. We therefore integrated uncertainty quantification into our estimators through the use of resampling-based confidence intervals, and the TAGC confidence interval [`GC_CI1`, `GC_CI2`] should contain $\pm 1$, if the TAGC point estimate is beyond $[-1,1]$. Our `TAGC` R package also gives a warning when the TAGC estimate is outside $[-1,+1]$.


- **Q.** Why the confounding bias check and correction in Step 2?

- **A.** Confounding bias can yield inflated test statistics in GWAS summary statistics and can be deteced by the intercept of LDSC ([LDSC paper](https://www.nature.com/articles/ng.3211), [LDSC software](https://github.com/bulik/ldsc)), and the GWAS summary statistics can be corrected by dividing GWAS effect size by `Lambda GC` of LDSC software output. Besides, user should adjust necessary covariates such as age, sex, genetic PCs, to reduce confounding bias.


- **Q.** Does TAGC apply to binary traits?

- **A.** TAGC supports both GLM ($-\log$(odds ratio)) and LM GWAS summary statistics. Regarding the estimation of heritability, when using LDSC, please make sure to use the observed scale h2 estimate and not to use the liability scale h2 estimate. 


- **Q.** How will the TAGC result change when we have GWAS for the smaller population and individual phenotype data for the larger population? For example, use GWAS summary statistics from a small Asian population and individual observations from a large European population for genetic correlation estimate.

- **A.** TAGC applies to the setup where the GWAS sample size is small and the individual observation sample size is large. In this case, the TAGC estimate remains unbiased, though the variance of the TAGC estimate gets larger due to the small GWAS sample size and noisier summary statistics.





