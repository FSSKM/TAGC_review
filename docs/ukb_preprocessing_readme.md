
# UKB data preprocessing
Here we outline the steps to obtain UKB phenotypes and covariates used in our real data analysis section.


## Data downloading

Once your data application is approved by the UK Biobank, the data can be downloaded from the https://ams.ukbiobank.ac.uk/ams/. A detailed instruction for downloading various data types is available in the [UK Biobank Data Access Guide](https://biobank.ctsu.ox.ac.uk/~bbdatan/Data_Access_Guide.pdf). The UNIX program `ukbconv` can also be downloaded from this page: https://biobank.ndph.ox.ac.uk/showcase/download.cgi

```{bash}
### preparation
wget  -nd  biobank.ndph.ox.ac.uk/ukb/util/ukbmd5
wget  -nd  biobank.ndph.ox.ac.uk/ukb/util/ukbconv
wget  -nd  biobank.ndph.ox.ac.uk/ukb/util/ukbfetch
wget  -nd  biobank.ndph.ox.ac.uk/ukb/util/ukblink
wget  -nd  biobank.ndph.ox.ac.uk/ukb/ukb/utilx/encoding.dat
wget  -nd  biobank.ndph.ox.ac.uk/ukb/util/ukbunpack
```

Then, the genetic data used in this study can be downloaded following the steps in "Section 3: Download utilities for bulk & genomics data and returns" of the [UK Biobank Data Access Guide](https://biobank.ctsu.ox.ac.uk/~bbdatan/Data_Access_Guide.pdf). Similarly, the phenotype and covariate data used in this study can be downloaded following the steps in "Section 2: The main dataset" of the UK Biobank Data Access Guide.



## Reprocessing steps of genetic data
After downloading the genetic data, we perform standard QC steps with Plink (https://www.cog-genomics.org/plink/2.0/) for both genotyping and imputed data. 

Here is an example script for the imputed data. 
```{bash}
for j in `seq 1 22`
do
sbatch  job_qc_ukb.sl  $j
done
```

where the content of `job_qc_ukb.sl` is below:
```{bash}
###
#!/bin/bash

#SBATCH -N 1
#SBATCH --mem 100g
#SBATCH -n 1
#SBATCH -t 2-00:00:00
#SBATCH --mail-type=end

echo "$j"

~/plink2 --bgen ~/ukb_imp_chr${1}_v3.bgen --sample ~/ukbxxxxx_imp_chr22_v3_sxxxx.sample  --mind 0.10 --maf 0.01 --geno 0.10 --hwe 0.0000001  --make-bed --out ~/ukb_imp_chr${1}_v3_qcd
```



## Preprocessing steps of phenotype and covariate data

After downloading the data, you will have a file called ukb12345.enc_ukb, where 12345 stands for your application ID (you actual application ID will be different). Here are the steps to obtain the phenotype and covariate data. 

First, you can convert the phenotype and covariates data using ukbconv and the following command. The result will be a text file called ukb12345.csv.

```{bash}
ukbconv ukb12345.enc_ukb csv -oukb12345
mv ukb12345.enc_ukb.csv ukb12345.csv
```

Next, we will construct a file containing the phenotypes and covariates you used for the data analysis. There are multiple options. Below is an example in R:

```{bash}
library(data.table)
data0<-fread("~/ukb12345.csv")
data1<-data0[,c(1,(22+1),(23+1),(88+1),(89+1):(92+1),(93+1):(96+1),(97+1):(100+1),
                (10651+1):(10654+1),
                (10663+1):(10666+1),(11002+1),(11003+1),(11011+1):(11050+1),
                (14996+1):(14999+1),(15026+1):(15031+1),(13530+1):(13531+1))]
fwrite(data1, "~/ukb12345_variables.csv", row.names = F)
```

Here you will need to find the Data-Field ID for the variables, and then identify which columns are useful for your analysis. 
Below is the list of the UKB phenotypes are used in our real data analysis and their Data-Field IDs (if available). The UKB website of Data-Field xxxxx can be accessed at https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=xxxxx. If more variables are needed, you can search on https://biobank.ndph.ox.ac.uk/showcase/search.cgi.
- Heel bone mineral density: 78
- Hand grip strength (left): 46
- Hand grip strength (right): 47
- Waist: 48
- Hip: 49
- BMI (manual measure): 21001
- BMI (bioimpedance measure): 23104
- Weight (bioimpedance measure): 23098
- Weight: (manual measure) 21002
- Body fat percentage: 23099
- HDL: 30760
- Platelet count: 30080
- Platelet distribution width: 30110
- Red blood cell count: 30010
- Diastolic blood pressure: 4079
- Systolic blood pressure: 4080
- Pulse rate: 102
- Vascular heart problems: 6150
- Blood clot problems: 6152
- Age at first live birth: 2754
- Age when periods started: 2714
- Time to identify matches: 20023
- Age completed education: 845
- Depression sum score
- Neuroticism sum score
- Anxiety sum score

The last three "sum score" variables were generated from original data fields as follows. 

```{bash}
Depression_sum_score<-as.matrix(apply(cbind(data1$X20510.0.0,data1$X20513.0.0,data1$X20508.0.0,data1$X20507.0.0,
                                            data1$X20517.0.0,data1$X20518.0.0,data1$X20511.0.0,data1$X20519.0.0,
                                            data1$X20514.0.0),1,sum))

Neuroticism_sum_score<-as.matrix(apply(cbind(data1$X2040.0.0,data1$X2030.0.0,data1$X2020.0.0,data1$X2010.0.0,
                                             data1$X2000.0.0,data1$X1990.0.0,data1$X1980.0.0,data1$X1970.0.0,
                                             data1$X1960.0.0,data1$X1950.0.0,data1$X1940.0.0,data1$X1930.0.0,
                                             data1$X1920.0.0),1,sum))

Anxiety_sum_score<-as.matrix(apply(cbind(data1$X20520.0.0,data1$X20515.0.0,data1$X20506.0.0,
                                         data1$X20505.0.0,data1$X20509.0.0,data1$X20516.0.0,
                                         data1$X20512.0.0),1,sum))          
```
For example, `data1$X20510.0.0` corresponds to Data-Field 20510, which can be found at https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=20510.
