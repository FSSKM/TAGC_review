


# LD Block Matrix Example

Here we show Step 5 and 6 in [README](../README.md) when computing LD matrix moments for TAGC and WAGC estimates. We use the 1000 Genome EUR and EAS genotype data. 

The 1000G EUR reference panel is available [here](https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip). The 1000G EAS reference panel is available [here](https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eas.zip).



### LD block matrix generation
Below we use PLINK 1.9 to generate LD blocks from 1000G EUR/EAS bfiles. 
```{bash}
# generate 253 LD blocks from 1000G bfile EUR
for i in `seq 1 253`
do
plink_1.9
  --bfile g1000_eur_allchr_ukb_unimputed_v2 \
  --extract TAGC_code/TAGC_data/LD_boundary/file${i} \
  --make-bed \
  --out g1000_eur_ukb_unimputed_file${i} \
  --range
done

# generate 253 LD blocks from 1000G bfile EAS
for i in `seq 1 253`
do
plink_1.9
  --bfile g1000_eas_allchr_ukb_unimputed_v2 \
  --extract TAGC_code/TAGC_data/LD_boundary/file${i} \
  --make-bed \
  --out g1000_eas_ukb_unimputed_file${i} \
  --range
done
```


## Within-Ancestry GC (WAGC) LD block matrix moment computation
Code below is run in R and can be parallelized. `WAGC_LD_block_moment` returns the path to the .txt file containing the moments of the LD block matrix. `LD_block_moment_all_WAGC` returns a vector of the three matrix moments needed to compute WAGC estimate.

```{bash}
library(TAGC)

# compute matrix moments for each LD block
for (j in 1:253) {
    # bfile_path contains both of the 253 EUR LD blocks
    bfile_path = '.'
    bfile_name = paste0('g1000_eur_ukb_unimputed_file', j)
    output_path = 'test_output_g1000_eur'

    rv = WAGC_LD_block_moment(
        bfile_path = bfile_path, 
        bfile_name = bfile_name,
        output_path = output_path)
}
rv # for block 253
# [1] "test_output_g1000_eur/WAGC_LD_block_moment/g1000_eur_ukb_unimputed_file253_traces.txt"

# compute matrix moments using LD block moments
LD_mom = LD_block_moment_all_WAGC(
    LD_block_mom_path = 'test_output_g1000_eur/WAGC_LD_block_moment', 
    n1 = 503, # 503 subjects in 1000G EUR
    output_path = 'test_output_g1000_eur', 
    output_filename = 'g1000_eur_WAGC')
LD_mom
# [1]   281430  1252953 47253520
```
User can choose to only use a subset of subjects in the bfile for LD matrix moment computation by providing the path to a list of desired subject IDs to the `WAGC_LD_block_moment` function, argument `popn_list_path`.




## Trans-Ancestry GC (TAGC) LD block matrix moment computation
Here we show the LD matrix moments computation for TAGC, where population I is EUR and population II is EAS. The R code below has a for loop and can be parallelized. 

`TAGC_LD_block_moment` returns the path to the .txt file containing the moments of the LD block matrix. `LD_block_moment_all_TAGC` returns a vector of the three matrix moments needed to compute TAGC estimate.

```{bash}
library(TAGC)

# compute matrix moments for each LD block
for (j in 1:253) {
    bfile_name1 = paste0('g1000_eur_ukb_unimputed_file', j)
    bfile_name2 = paste0('g1000_eas_ukb_unimputed_file', j)
    output_path = 'test_output_g1000_eur_eas'

    rv = TAGC_LD_block_moment(
        bfile_path_popn1 = '.', 
        bfile_path_popn2 = '.', 
        bfile_name_popn1 = bfile_name1, 
        bfile_name_popn2 = bfile_name2, 
        output_path = output_path)
}
rv # for block 253
# [1] "test_output_g1000_eur_eas/TAGC_LD_block_moment/g1000_eur_ukb_unimputed_file1_g1000_eas_ukb_unimputed_file1_TAGC_traces.txt"

# compute matrix moments using LD block moments
LD_mom2 = LD_block_moment_all_TAGC(
    LD_block_mom_path = 'test_output_g1000_eur_eas/TAGC_LD_block_moment', 
    n1 = 503, # 503 subjects in 1000G EUR
    output_path = 'test_output_g1000_eur_eas',
    output_filename = 'g1000_eur_eas_TAGC')
LD_mom2
# [1]   281430.0   668493.2 15315197.7
```
Similarly, user can choose to only use a subset of subjects in the bfile of population I (`bfile_name1`) and another subset of population II (`bfile_name_popn2`) for LD matrix moment computation by providing the path to two lists of desired subject IDs to the `TAGC_LD_block_moment` function, argument `popn1_list_path` and `popn2_list_path`.

`LD_mom2` is the three matrix moments needed for TAGC bias correction and estimate. `LD_block_moment_all_TAGC` also outputs a .txt file at `test_output_g1000_eur_eas/g1000_eur_eas_TAGC.txt` that contains these three matrix moments.





