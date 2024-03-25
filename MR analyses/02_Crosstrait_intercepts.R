## Project: EA_MH_CBS_MR 
## Script purpose: Run genetic correlations and LDSC intercepts between traits 
## to identify sample overlap between GWAS 
## Starting date: 06-08-2021
## Author: Perline Demange & Michel Nivard

# Set up ####
library(data.table)
library(tidyverse)
#library(devtools)
#install_github("GenomicSEM/GenomicSEM")
library(GenomicSEM)

setwd("C:/Users/user/Dropbox/CBS - MR/02_data_analysis/Summary_statistics/")

# 1. Munge summary statistics #################

## 1.1 EA ###### 

### 1.1.1 EA without 23andMe ######
munge(files = "EA/GWAS_EA_excl23andMe.txt", 
      hm3 = "C:/Users/user/OneDrive - Vrije Universiteit Amsterdam/Documents/References_Data/w_hm3.snplist",
      trait.names = "EA_excl23andMe",
      N = 766345,
      info.filter = 0.9, maf.filter = 0.01)

### 1.1.2 EA meta ######
munge(files = "EA/meta_education_lee_23andMe_betas_20210813.txt", 
      hm3 = "C:/Users/user/OneDrive - Vrije Universiteit Amsterdam/Documents/References_Data/w_hm3.snplist",
      trait.names = "EA",
      N = 1131883,
      info.filter = 0.9, maf.filter = 0.01)

### 1.1.3 EA siblings  ######
#### 1.1.3.1 Reformat the summary statistics ########
#install.packages("vcfR")
library(vcfR)
#vcf <- read.vcfR("EA/withinsib_EA_ieu-b-4835.vcf") 
vcf <- read.vcfR("../Summary_statistics/EA/Withinsib/ieu-b-4836.vcf") #from https://gwas.mrcieu.ac.uk/datasets/ieu-b-4836/
vcf@meta
data_snp <- vcf@fix
data_effect <- vcf@gt
head(data_snp)
data_snp <- as.data.frame(data_snp)
data_effect <- as.data.frame(data_effect)
head(data_effect)

#get the effect from the vcr 
library(stringr)
data_effect_cl <- str_split_fixed(data_effect$`ieu-b-4836`, ":", 5)
head(data_effect_cl)
data_effect_cl <- as.data.frame(data_effect_cl)
colnames(data_effect_cl) <- c("effect", "SE", "logpval", "SS", "ID") # this depends of the format here https://github.com/MRCIEU/gwas-vcf-specification

# Check potential missing data
data_snp %>% summarise_all(~ sum(is.na(.)))
data_effect_cl %>% summarise_all(~ sum(is.na(.)))

#combine
Education_wf_ieu_4836 <- cbind(data_snp, data_effect_cl) # careful, this should not be merged based on Id, as there are multiple identical rsid
head(Education_wf_ieu_4836)

names(Education_wf_ieu_4836)[13] <- "ID2"
Education_wf_ieu_4836 %>% summarise_all(~ sum(is.na(.)))
Education_wf_ieu_4836[Education_wf_ieu_4836$ID == "Na",]
Education_wf_ieu_4836[Education_wf_ieu_4836$POS == "241863494",] #this row has ID Na but doesnt appear when looking for Na 
# it looks like when ID is na, ID2 is empty 

Education_wf_ieu_4836 <- Education_wf_ieu_4836 %>% mutate_all(na_if,"")
#this should solve readability issues, and NA should be excluded in merging 

Education_wf_ieu_4836$CHROM <- as.numeric(Education_wf_ieu_4836$CHROM)
Education_wf_ieu_4836$POS <- as.numeric(Education_wf_ieu_4836$POS)
Education_wf_ieu_4836$effect <- as.numeric(Education_wf_ieu_4836$effect)
Education_wf_ieu_4836$SE <- as.numeric(Education_wf_ieu_4836$SE)
Education_wf_ieu_4836$logpval <- as.numeric(Education_wf_ieu_4836$logpval)
Education_wf_ieu_4836$SS <- as.numeric(Education_wf_ieu_4836$SS)

Education_wf_ieu_4836$pval <- 10^-(Education_wf_ieu_4836$logpval) # log pvalue is - log pvalue
summary(Education_wf_ieu_4836$pval)


write.table(Education_wf_ieu_4836, "EA/Withinsib/Education_wf_ieu_4836_reformat.txt", row.names = F, quote = F)
write.csv2(Education_wf_ieu_4836, "EA/Withinsib/Education_wf_ieu_4836_reformat.csv", row.names = F, quote = F)

#remove unnecessary columns 
Education_wf_ieu_4836_minimal <- Education_wf_ieu_4836[c("ID", "REF", 
                                                         "ALT", "effect",
                                                         "pval")]

write.table(Education_wf_ieu_4836_minimal, "EA/Withinsib/Education_wf_ieu_4836_reformat_minimal.txt", row.names = F, quote = F)

test <- fread("../Summary_statistics/EA/Withinsib/Education_wf_ieu_4836_reformat.txt")

#### 1.1.3.1 Munging ########
munge(files = "../Summary_statistics/EA/Withinsib/Education_wf_ieu_4836_reformat_minimal.txt", 
      hm3 = "C:/Users/user/OneDrive - Vrije Universiteit Amsterdam/Documents/References_Data/w_hm3.snplist",
      trait.names = "EAwf2209",
      N = 99970, #this is the sample size given by the ieu, if the SS col is sample size than the max is higher 
      info.filter = 0.9, maf.filter = 0.01)

lookup <-fread("EAwf2209.sumstats.gz")
tail(lookup)
summary(lookup$Z) #all Nas 


## 1.2 Psychopathologies ######

files <- c("ASD_iPSYCH-PGC_ASD_Nov2017.gz",
           "Anorexia_pgcAN2.2019-07.vcf.tsv.gz",
           "Anxiety/TotAnx_effect_sumstats.gz",
           "Anxiety/META_UKBB_iPSYCH_ANGST_wNcol.sumstats",
           "ADHD/daner_adhd_meta_filtered_NA_iPSYCH23_PGC11_sigPCs_woSEX_2ell6sd_EUR_Neff_70.meta.gz"
)

files2 <- c("Bipolar/pgc-bip2021-BDI.vcf.tsv.txt",
            "Bipolar/pgc-bip2021-BDII.vcf.tsv.txt"
)

files3 <- c("MDD/MDD2018_ex23andMe.gz",
            "OCD/ocd_aug2017.gz",
            "pts_eur_freeze2_overall.results"
)

trait.names <- c("ASD_Grove",
                 "Ano_Watson",
                 "Anx_Purves_UKB",
                 "Anx_Purves_meta",
                 "ADHD_Demontis"
)
trait.names2 <- c("BIP_I_Mullins",
                  "BIP_II_Mullins"
)
trait.names3 <- c("MDD_Wray",
                  "OCD_Arnold",
                  "PTSD_Nievergelt" 
)




sample.size <- c(46351,
                 72517,
                 83566,
                 114019,
                 53293
)
sample.size2 <- c(475038,
                  370856
)
sample.size3 <- c(480359,
                  9725,
                  174659
)


munge(files = files, 
      hm3 = "C:/Users/user/OneDrive - Vrije Universiteit Amsterdam/Documents/References_Data/w_hm3.snplist",
      trait.names = trait.names,
      N = sample.size,
      info.filter = 0.9, maf.filter = 0.01)
munge(files = files2, 
      hm3 = "C:/Users/user/OneDrive - Vrije Universiteit Amsterdam/Documents/References_Data/w_hm3.snplist",
      trait.names = trait.names2,
      N = sample.size2,
      info.filter = 0.9, maf.filter = 0.01)
munge(files = files3, 
      hm3 = "C:/Users/user/OneDrive - Vrije Universiteit Amsterdam/Documents/References_Data/w_hm3.snplist",
      trait.names = trait.names3,
      N = sample.size3,
      info.filter = 0.9, maf.filter = 0.01)

munge(files = "MDD/PGC_UKB_depression_genome-wide.txt", 
      hm3 = "C:/Users/user/OneDrive - Vrije Universiteit Amsterdam/Documents/References_Data/w_hm3.snplist",
      trait.names = "MDD_Howard",
      N =  500199,
      info.filter = 0.9, maf.filter = 0.01)

munge(files = "Bipolar/pgc-bip2021-all.vcf.tsv.txt", 
      hm3 = "C:/Users/user/OneDrive - Vrije Universiteit Amsterdam/Documents/References_Data/w_hm3.snplist",
      trait.names = "BIP_all_Mullins",
      N =  413466,
      info.filter = 0.9, maf.filter = 0.01)

munge(files =  "SCZ/PGC3_SCZ_wave3_public.v2.tsv", 
      hm3 = "C:/Users/user/OneDrive - Vrije Universiteit Amsterdam/Documents/References_Data/w_hm3.snplist",
      trait.names = "SCZ_PGC3",
      N =  161405,
      info.filter = 0.9, maf.filter = 0.01)

munge(files =  "Alc/pgc_alcdep.eur_discovery.aug2018_release_rsid.txt", 
      hm3 = "C:/Users/user/OneDrive - Vrije Universiteit Amsterdam/Documents/References_Data/w_hm3.snplist",
      trait.names = "Alc_Walters",
      N =  46568,
      info.filter = 0.9, maf.filter = 0.01)

check <- fread("ADHD/Demontis 2023/ADHD2022_iPSYCH_deCODE_PGC.meta.gz")


munge(files =  "ADHD/Demontis 2023/ADHD2022_iPSYCH_deCODE_PGC.meta.gz", 
      hm3 = "C:/Users/user/OneDrive - Vrije Universiteit Amsterdam/Documents/References_Data/w_hm3.snplist",
      trait.names = "ADHD_Demontis_2023",
      N =  225534,
      info.filter = 0.9, maf.filter = 0.01)

# 2. Run LDSC to obtain matrix with LD cross trait intercept  ##################

trait.file <- c("ASD_Grove.sumstats.gz",
                "Ano_Watson.sumstats.gz",
                "Anx_Purves_UKB.sumstats.gz",
                "Anx_Purves_meta.sumstats.gz",
                "ADHD_Demontis.sumstats.gz",
                "BIP_all_Mullins.sumstats.gz",
                "BIP_I_Mullins.sumstats.gz",
                "BIP_II_Mullins.sumstats.gz",
                "MDD_Howard.sumstats.gz",
                "MDD_Wray.sumstats.gz",
                "OCD_Arnold.sumstats.gz",
                "PTSD_Nievergelt.sumstats.gz",
                "SCZ_PGC3.sumstats.gz",
                "Alc_Walters.sumstats.gz",
                "ADHD_Demontis_2023.sumstats.gz"
)
sample.prevs <- c(0.396582598,
                  0.23431747,
                  0.304585597,
                  0.280453258,
                  0.378717655,
                  0.101379557,
                  0.052753674,
                  0.018284725,
                  0.341376132,
                  0.281993259,
                  0.276401028,
                  0.132898963,
                  0.417521142,
                  0.248,
                  0.171552848
                  
)
pop.prevs <- c(0.012,
               0.04,
               0.16,
               0.16,
               0.05,
               0.02,
               0.01,
               0.01,
               0.15,
               0.15,
               0.015,
               0.3,
               0.01,
               0.159,
               0.05
)
trait.Ids <- c("ASD_Grove",
                 "Ano_Watson",
                 "Anx_Purves_UKB",
                 "Anx_Purves_meta",
                 "ADHD_Demontis",
                 "BIP_all_Mullins",
                 "BIP_I_Mullins",
                 "BIP_II_Mullins",
                 "MDD_Howard",
                 "MDD_Wray",
                 "OCD_Arnold",
                 "PTSD_Nievergelt",
                 "SCZ_PGC3",
                 "Alc_Walters",
                 "ADHD_Demontis_2023"
)

for (i in 1:length(trait.file)){ 
traits <- c(trait.file[i], "EA.sumstats.gz")
sample.prev <- c(sample.prevs[i],NA)
population.prev <- c(pop.prevs[i],NA)
ld <- "C:/Users/user/OneDrive - Vrije Universiteit Amsterdam/Documents/References_Data/eur_w_ld_chr/"
wld <- "C:/Users/user/OneDrive - Vrije Universiteit Amsterdam/Documents/References_Data/eur_w_ld_chr/"
trait.names<-c(trait.Ids[i],"EA")
ldsc(traits, sample.prev, population.prev, ld, wld, trait.names)
}

ld <- "C:/Users/user/OneDrive - Vrije Universiteit Amsterdam/Documents/References_Data/eur_w_ld_chr/"
wld <- "C:/Users/user/OneDrive - Vrije Universiteit Amsterdam/Documents/References_Data/eur_w_ld_chr/"

for (i in 1:length(trait.file)){ 
   traits <- c(trait.file[i], "EAwf2209.sumstats.gz")
   sample.prev <- c(sample.prevs[i],NA)
   population.prev <- c(pop.prevs[i],NA)
   trait.names<-c(trait.Ids[i],"EAwf2209")
   ldsc(traits, sample.prev, population.prev, ld, wld, trait.names)
}

# 3. Get pvalues of the CTI #####
# I got the CTI manually in the log file 

cti <- fread("Cross Trait Intercepts/cti.csv")
cti$Z <- cti$cti/cti$se
cti
cti$p <-2*pnorm(-abs(cti$Z))

write.table(cti, "cti_pval.csv", row.names = F)

cti <- fread("Cross Trait Intercepts/cti_wf_2209.csv")
cti$Z <- cti$cti/cti$se
cti
cti$p <-2*pnorm(-abs(cti$Z))

write.table(cti, "cti_wf_2209_pval.csv", row.names = F)


# 4. Check some additional correlations to make sure of proper strand direction #####
## 4.1 Individual checks ####
traits <- c("BIP_all_Mullins.sumstats.gz", "SCZ_PGC3.sumstats.gz")
sample.prev <- c(0.101379557, 0.417521142)
population.prev <- c(0.02,0.01)
ld <- "C:/Users/user/OneDrive - Vrije Universiteit Amsterdam/Documents/References_Data/eur_w_ld_chr/"
wld <- "C:/Users/user/OneDrive - Vrije Universiteit Amsterdam/Documents/References_Data/eur_w_ld_chr/"
trait.names<-c("BIP_all_Mullins","SCZ_PGC3")
ldsc(traits, sample.prev, population.prev, ld, wld, trait.names)
# correlation is positive 

traits <- c("Ano_Watson.sumstats.gz","Anx_Purves_meta.sumstats.gz")
sample.prev <- c(0.23431747, 0.280453258)
population.prev <- c(0.04,0.16)
ld <- "C:/Users/user/OneDrive - Vrije Universiteit Amsterdam/Documents/References_Data/eur_w_ld_chr/"
wld <- "C:/Users/user/OneDrive - Vrije Universiteit Amsterdam/Documents/References_Data/eur_w_ld_chr/"
trait.names<-c("Ano_Watson","Anx_Purves_meta")
ldsc(traits, sample.prev, population.prev, ld, wld, trait.names)
#correlation is negative!!!!

traits <- c("OCD_Arnold.sumstats.gz","Anx_Purves_meta.sumstats.gz")
sample.prev <- c(0.276401028, 0.280453258)
population.prev <- c(0.015,0.16)
ld <- "C:/Users/user/OneDrive - Vrije Universiteit Amsterdam/Documents/References_Data/eur_w_ld_chr/"
wld <- "C:/Users/user/OneDrive - Vrije Universiteit Amsterdam/Documents/References_Data/eur_w_ld_chr/"
trait.names<-c("OCD_Arnold","Anx_Purves_meta")
ldsc(traits, sample.prev, population.prev, ld, wld, trait.names)
#correlation is positive 

traits <- c("Ano_Watson.sumstats.gz","OCD_Arnold.sumstats.gz")
sample.prev <- c(0.23431747, 0.276401028)
population.prev <- c(0.04,0.015)
ld <- "C:/Users/user/OneDrive - Vrije Universiteit Amsterdam/Documents/References_Data/eur_w_ld_chr/"
wld <- "C:/Users/user/OneDrive - Vrije Universiteit Amsterdam/Documents/References_Data/eur_w_ld_chr/"
trait.names<-c("Ano_Watson","Anx_Purves_meta")
ldsc(traits, sample.prev, population.prev, ld, wld, trait.names)
# correlation is negative!!! Anorexia strand direction is flipped 

traits <- c("ASD_Grove.sumstats.gz","Anx_Purves_meta.sumstats.gz")
sample.prev <- c(0.396582598, 0.280453258)
population.prev <- c(0.012,0.16)
ld <- "C:/Users/user/OneDrive - Vrije Universiteit Amsterdam/Documents/References_Data/eur_w_ld_chr/"
wld <- "C:/Users/user/OneDrive - Vrije Universiteit Amsterdam/Documents/References_Data/eur_w_ld_chr/"
trait.names<-c("ASD_Grove","Anx_Purves_meta")
ldsc(traits, sample.prev, population.prev, ld, wld, trait.names)
#correlation is positive 

## 4.2 All disorders ######


trait.file <- c(trait.file, "EA.sumstats.gz", "EAwf2209.sumstats.gz")
sample.prevs <- c(sample.prevs, NA, NA)
pop.prevs <- c(pop.prevs, NA, NA)
trait.Ids <- c(trait.Ids, "EA", "EAwf2209")

LDSCoutput<-ldsc(traits=trait.file, 
                 sample.prev=sample.prevs, 
                 population.prev=pop.prevs, 
                 ld="C:/Users/user/OneDrive - Vrije Universiteit Amsterdam/Documents/References_Data/eur_w_ld_chr/",
                 wld="C:/Users/user/OneDrive - Vrije Universiteit Amsterdam/Documents/References_Data/eur_w_ld_chr/",
                 trait.names=trait.Ids)

#optional command to save the output as a .RData file for later use
save(LDSCoutput,file="LDSCoutput_matrix_psychdisorders_EA.RData")

load("LDSCoutput_matrix_psychdisorders_EA.RData")
# plot
corrplot::corrplot(cov2cor(LDSCoutput$S),is.corr=F)

#Anorexia GWAS is reverse coded (so previously presented results are sign-switched)
#Anx_Purves_UKb is reverse coded (as mentionned in their readme) BUT NOT Anx_Purves_meta
#As I already knew EAwf2209 is reverse coded 

# to output the standard errors of the ld-score regression in the order they are listed in the genetic covariance (i.e., S) matrix
k<-nrow(LDSCoutput$S)
SE<-matrix(0, k, k)
SE[lower.tri(SE,diag=TRUE)] <- sqrt(diag(LDSCoutput$V))
SE
