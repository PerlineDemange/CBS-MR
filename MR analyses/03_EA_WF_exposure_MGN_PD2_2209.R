## Project: EA_MH_CBS_MR 
## Script purpose: Run MR analyses with EA as exposure 
## Author: Perline Demange & Michel Nivard

# Set up #####
library(data.table)
library(tidyverse)
library(readxl)

library(devtools)
#install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR) # https://mrcieu.github.io/TwoSampleMR/index.html 

setwd("C:/Users/user/Dropbox/CBS - MR/MR analyses/")

# 1. Load data Exposure EA #############

# #* 1.1. EA Within Family ----------------

ea_wf <- fread("../Summary_statistics/EA/Withinsib/Education_wf_ieu_4836_reformat.txt")
head(ea_wf)

#** 1.1.1 Select SNPS as instruments, use Lee et al Sup table 2 as source for instruments -------
# SNP in Sup Table 2 are independent SNPs, so no need to clump 

Leeetal_lead_SNPs_EA <- read_excel("../Summary_statistics/EA/Leeetal_lead_SNPs_EA.xlsx", skip = 1)
head(Leeetal_lead_SNPs_EA)
 
ea_wf_instrument_bp <- ea_wf[paste0(ea_wf$CHROM,ea_wf$POS) %in%
                            paste0(Leeetal_lead_SNPs_EA$Chr,
                                   Leeetal_lead_SNPs_EA$Position), ] #1254

ea_wf_instrument <- ea_wf[ea_wf$ID %in% Leeetal_lead_SNPs_EA$SNP, ] #1177
# More instruments if look at chrom pos, but it is safer to use the rsid


#** 1.1.2 Check that the instruments are valid #####
# We define valid instruments for the Within-sibling gwas as p < 0.05
summary(ea_wf_instrument$pval)

nrow(ea_wf_instrument[ea_wf_instrument$pval <0.05,]) # 578
ea_wf_instrument <- ea_wf_instrument[ea_wf_instrument$pval <0.05,]

 
# #** 1.1.3 Reformat for MR ----------
ea_wf_exp_dat <- format_data( ea_wf_instrument, 
                           type="exposure",   
                           snp_col = "ID",
                           beta_col = "effect",
                           se_col = "SE",
                           effect_allele_col = "ALT", #this is odd, but the ref as effect allele doesnt make sense, also in regards to previous versions of the withinsib gwas
                           other_allele_col = "REF",
                           pval_col = "pval")


ea_wf_exp_dat$exposure <- "EAwf"
head(ea_wf_exp_dat)
 
# ** 1.1.4 rescale to the standard deviation of "years of schooling" from Lee et al. ####
ea_wf_exp_dat$beta.exposure <- 3.9*ea_wf_exp_dat$beta.exposure 
ea_wf_exp_dat$se.exposure <- 3.9*ea_wf_exp_dat$se.exposure
 

#* 1.2 EA GWAS (including 23andMe, meta-analysis from CogNonCog project) --------
# Meta-analysis was described and checked in document Meta-analysis EA.R

ea_meta <- fread("../Summary_statistics/EA/meta_education_lee_23andMe_130820211.txt") # Perline's MA AUGST 13TH
head(ea_meta)

#** 1.2.1 Select significant SNPS as instruments ----------------
ea_meta_instrument <- ea_meta[ea_meta$`P-value` <= 0.00000005, ]

#** 1.2.2 compute beta & se  ----------------

ea_meta_instrument$Beta <- ea_meta_instrument$Zscore / 
  sqrt(2*ea_meta_instrument$Freq1*
         (1 - ea_meta_instrument$Freq1)*
         (ea_meta_instrument$Weight+ ea_meta_instrument$Zscore*
            ea_meta_instrument$Zscore))  
ea_meta_instrument$SE <- 1/
  sqrt(2*ea_meta_instrument$Freq1*
         (1 - ea_meta_instrument$Freq1)*
         (ea_meta_instrument$Weight+ 
            ea_meta_instrument$Zscore*ea_meta_instrument$Zscore))  
  
#** 1.2.4 scale beta to "years of education", which has a specific standard deviation  (3.9). #################

ea_meta_instrument$Beta <- 3.9*ea_meta_instrument$Beta
ea_meta_instrument$SE <- 3.9*ea_meta_instrument$SE
  
rm(ea_meta)
rm(ea_wf)

#** 1.2.2 Reformat for MR ---------------------
ea_meta_exp_dat <- format_data(ea_meta_instrument, 
                             type="exposure",   
                             snp_col = "MarkerName",
                             beta_col = "Beta",
                             se_col = "SE",
                             effect_allele_col = "Allele1",
                             other_allele_col = "Allele2",
                             pval_col = "P-value",
                             eaf_col = "Freq1")
ea_meta_exp_dat$exposure <- "EA"


# 2. Load Outcomes ########################

#* 2.1. ASD #################
# (read me for the file: https://figshare.com/articles/dataset/asd2019/14671989?file=28169289 ) #################

asd <- fread("../Summary_statistics/ASD_iPSYCH-PGC_ASD_Nov2017.gz")
head(asd)

# convert the OR to log(OR) or the beta of a logistic regression:
asd$beta <- log(asd$OR)

#** 2.1.1. Clump EA specifically for ASD ######

ea_meta_exp_dat_asd <- ea_meta_exp_dat[ea_meta_exp_dat$SNP %in% asd$SNP,]

ea_meta_exp_clumped_asd <- clump_data(ea_meta_exp_dat_asd,
                                  clump_kb = 1000, 
                                  clump_r2 = 0.001, 
                                  clump_p1 = 1,
                                  clump_p2 = 1,
                                  pop = "EUR"
)
# 1000000 0.1 1393 SNPS
# 1000 0.1: 1402 SNPs 
# 1000 0.01: 990 SNPS
# 1000 0.001: 841 SNPS
# 10000 0.001:515 SNPs

ea_wf_exp_dat_asd <- ea_wf_exp_dat[ea_wf_exp_dat$SNP %in% asd$SNP,]

ea_wf_exp_dat_clumped_asd <- clump_data(ea_wf_exp_dat_asd,
                                      clump_kb = 1000, 
                                      clump_r2 = 0.001, 
                                      clump_p1 = 1,
                                      clump_p2 = 1,
                                      pop = "EUR"
)



#** 2.1.2. Extract SNPS for ASD outcome and format: #################
asd_dat <- format_data(asd, 
                        type="outcome",
                        snps= ea_meta_exp_clumped_asd$SNP,
                        snp_col = "SNP",
                        beta_col = "beta",
                        se_col = "SE",
                        effect_allele_col = "A1",
                        other_allele_col = "A2",
                        pval_col = "P")
asd_dat$outcome <- "ASD"


asd_dat.wf <- format_data(asd, 
                       type="outcome",
                       snps= ea_wf_exp_dat_clumped_asd$SNP,
                       snp_col = "SNP",
                       beta_col = "beta",
                       se_col = "SE",
                       effect_allele_col = "A1",
                       other_allele_col = "A2",
                       pval_col = "P")
asd_dat.wf$outcome <- "ASD"

rm(asd)


#** 2.1.3. Harmonize exposure and outcome data : #################

dat.asd <- harmonise_data(
  exposure_dat = ea_meta_exp_clumped_asd, 
  outcome_dat = asd_dat, action = 2
)
nb.snps.asd <- nrow(dat.asd)

dat.asd.wf <- harmonise_data(
  exposure_dat = ea_wf_exp_dat_clumped_asd, 
  outcome_dat = asd_dat.wf, action = 2
)
nb.snps.asd.wf <- nrow(dat.asd.wf)

#** 2.1.4. Run pre-registered MR: #################
# function in file MR_analyses)functions.R 

asd.results <- run_MR_analyses(dat.asd)
asd.results
asd.wf.results <- run_MR_analyses(dat.asd.wf)
asd.wf.results

save(asd.results, asd.wf.results, file = "asd.MR.results.2209.RData")
load("asd.MR.results.2209.RData")

#** 2.1.5. Figures #################
p1 <- mr_scatter_plot(asd.results[[1]], dat.asd)
p1[[1]]


#* 2.2. ADHD 2018  #################
# (for format of the file see: https://docs.google.com/document/d/1TWIhr8-qpCXB13WCXcU1_HDio8lC_MeWoAg2jlggrtU/edit )
adhd <- fread("../Summary_statistics/ADHD/daner_adhd_meta_filtered_NA_iPSYCH23_PGC11_sigPCs_woSEX_2ell6sd_EUR_Neff_70.meta.gz")
head(adhd)

# take log(OR) to get the beta, s.e. is already on the log(OR) scale according to description of the daner file format.
adhd$beta<- log(adhd$OR)


#** 2.2.1. Clump EA specifically for ADHD ######

ea_meta_exp_dat_adhd <- ea_meta_exp_dat[ea_meta_exp_dat$SNP %in% adhd$SNP,]


ea_meta_exp_clumped_adhd <- clump_data(ea_meta_exp_dat_adhd,
                                      clump_kb = 1000,
                                      clump_r2 = 0.001,
                                      clump_p1 = 1,
                                      clump_p2 = 1,
                                      pop = "EUR"
)

ea_wf_exp_dat_adhd <- ea_wf_exp_dat[ea_wf_exp_dat$SNP %in% adhd$SNP,]

ea_wf_exp_dat_clumped_adhd <- clump_data(ea_wf_exp_dat_adhd,
                                        clump_kb = 1000, 
                                        clump_r2 = 0.001, 
                                        clump_p1 = 1,
                                        clump_p2 = 1,
                                        pop = "EUR"
)




#** 2.2.2 Extract SNPS for outcome and format: #################

adhd_dat <- format_data(adhd, 
                        type="outcome",
                        snps= ea_meta_exp_clumped_adhd$SNP,
                        snp_col = "SNP",
                        beta_col = "beta",
                        se_col = "SE",
                        effect_allele_col = "A1",
                        other_allele_col = "A2",
                        eaf_col = "FRQ_A_19099",
                        pval_col = "P",
                        samplesize_col = "Neff")
adhd_dat$outcome <- "ADHD"

adhd_dat.wf <- format_data(adhd, 
                        type="outcome",
                        snps= ea_wf_exp_dat_clumped_adhd$SNP,
                        snp_col = "SNP",
                        beta_col = "beta",
                        se_col = "SE",
                        effect_allele_col = "A1",
                        other_allele_col = "A2",
                        eaf_col = "FRQ_A_19099",
                        pval_col = "P",
                        samplesize_col = "Neff")
adhd_dat.wf$outcome <- "ADHD"


rm(adhd)


#** 2.2.3 Harmonize exposure and outcome data #################
dat.adhd <- harmonise_data(
  exposure_dat = ea_meta_exp_clumped_adhd, 
  outcome_dat = adhd_dat, action = 2
)
head(dat.adhd)
nb.snps.adhd <- nrow(dat.adhd)

dat.adhd.wf <- harmonise_data(
  exposure_dat = ea_wf_exp_dat_clumped_adhd, 
  outcome_dat = adhd_dat.wf, action = 2
)
head(dat.adhd.wf)
nb.snps.adhd.wf <- nrow(dat.adhd.wf)

#** 2.2.4 MR analyses ADHD #############

adhd.results <- run_MR_analyses(dat.adhd)
adhd.results
adhd.wf.results <- run_MR_analyses(dat.adhd.wf)
adhd.wf.results

save(adhd.results, adhd.wf.results, file = "adhd.MR.results.2209.RData")
load("adhd.MR.results.2209.RData")

#** 2.2.5 Figures ADHD ###########
res_single <- mr_singlesnp(dat.adhd)
p4 <- mr_funnel_plot(res_single)
p4[[1]]

res  <- mr(dat.adhd, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_mode", "mr_weighted_median"))
p1 <- mr_scatter_plot(res, dat.adhd) 
p1[[1]]

res_loo <- mr_leaveoneout(dat.adhd)
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]

#* 2.2.bis ADHD 2023  #################
adhd2023 <- fread("../Summary_statistics/ADHD/Demontis 2023/ADHD2022_iPSYCH_deCODE_PGC.meta.gz")
head(adhd2023)

# take log(OR) to get the beta, s.e. is already on the log(OR) scale according to description of the daner file format.
adhd2023$beta<- log(adhd2023$OR)


#** 2.2.1. Clump EA specifically for ADHD ######

ea_meta_exp_dat_adhd2023 <- ea_meta_exp_dat[ea_meta_exp_dat$SNP %in% adhd2023$SNP,]


ea_meta_exp_clumped_adhd2023 <- clump_data(ea_meta_exp_dat_adhd2023,
                                       clump_kb = 1000,
                                       clump_r2 = 0.001,
                                       clump_p1 = 1,
                                       clump_p2 = 1,
                                       pop = "EUR"
)

ea_wf_exp_dat_adhd2023 <- ea_wf_exp_dat[ea_wf_exp_dat$SNP %in% adhd2023$SNP,]

ea_wf_exp_dat_clumped_adhd2023 <- clump_data(ea_wf_exp_dat_adhd2023,
                                         clump_kb = 1000, 
                                         clump_r2 = 0.001, 
                                         clump_p1 = 1,
                                         clump_p2 = 1,
                                         pop = "EUR"
)




#** 2.2.2 Extract SNPS for outcome and format: #################

adhd2023_dat <- format_data(adhd2023, 
                        type="outcome",
                        snps= ea_meta_exp_clumped_adhd2023$SNP,
                        snp_col = "SNP",
                        beta_col = "beta",
                        se_col = "SE",
                        effect_allele_col = "A1",
                        other_allele_col = "A2",
                        eaf_col = "FRQ_A_19099",
                        pval_col = "P",
                        samplesize_col = "Neff")
adhd2023_dat$outcome <- "ADHD2023"

adhd2023_dat.wf <- format_data(adhd2023, 
                           type="outcome",
                           snps= ea_wf_exp_dat_clumped_adhd2023$SNP,
                           snp_col = "SNP",
                           beta_col = "beta",
                           se_col = "SE",
                           effect_allele_col = "A1",
                           other_allele_col = "A2",
                           eaf_col = "FRQ_A_19099",
                           pval_col = "P",
                           samplesize_col = "Neff")
adhd2023_dat.wf$outcome <- "ADHD2023"


rm(adhd2023)


#** 2.2.3 Harmonize exposure and outcome data #################
dat.adhd2023 <- harmonise_data(
  exposure_dat = ea_meta_exp_clumped_adhd2023, 
  outcome_dat = adhd2023_dat, action = 2
)
head(dat.adhd2023)
nb.snps.adhd2023 <- nrow(dat.adhd2023)

dat.adhd2023.wf <- harmonise_data(
  exposure_dat = ea_wf_exp_dat_clumped_adhd2023, 
  outcome_dat = adhd2023_dat.wf, action = 2
)
head(dat.adhd2023.wf)
nb.snps.adhd2023.wf <- nrow(dat.adhd2023.wf)

#** 2.2.4 MR analyses ADHD #############

adhd2023.results <- run_MR_analyses(dat.adhd2023)
adhd2023.results
adhd2023.wf.results <- run_MR_analyses(dat.adhd2023.wf)
adhd2023.wf.results

save(adhd2023.results, adhd2023.wf.results, file = "adhd2023.MR.results.2209.RData")
load("adhd2023.MR.results.2209.RData")

#** 2.4.5 Figures ADHD ###########
res_single <- mr_singlesnp(dat.adhd2023)
p4 <- mr_funnel_plot(res_single)
p4[[1]]

res  <- mr(dat.adhd2023, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_mode", "mr_weighted_median"))
p1 <- mr_scatter_plot(res, dat.adhd2023) 
p1[[1]]

res_loo <- mr_leaveoneout(dat.adhd2023)
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]


#* 2.3 Anorexia Nervosa #################

an <- fread("../Summary_statistics/Anorexia_pgcAN2.2019-07.vcf.tsv.gz")
head(an)

# #Get the header which contains the readme
# an_header <- file("../Summary_statistics/pgcAN2.2019-07.vcf.tsv")
# header <- readLines(an_header, n=100)
# sink("../Summary_statistics/HeadersOf_pgcAN2.2019-07.vcf.tsv.txt")
# cat(paste0(header, collapse="\n"))
# sink()

#** 2.3.1. Clump EA specifically for Anorexia ######
ea_meta_exp_dat_an <- ea_meta_exp_dat[ea_meta_exp_dat$SNP %in% an$ID,]


ea_meta_exp_clumped_an <- clump_data(ea_meta_exp_dat_an,
                                       clump_kb = 1000,
                                       clump_r2 = 0.001,
                                       clump_p1 = 1,
                                       clump_p2 = 1,
                                       pop = "EUR"
)

ea_wf_exp_dat_an <- ea_wf_exp_dat[ea_wf_exp_dat$SNP %in% an$ID,]

ea_wf_exp_dat_clumped_an <- clump_data(ea_wf_exp_dat_an,
                                         clump_kb = 1000, 
                                         clump_r2 = 0.001, 
                                         clump_p1 = 1,
                                         clump_p2 = 1,
                                         pop = "EUR"
)


#** 2.3.2  Extract SNPS for outcome and format: #################

# Correlation matrix with other traits suggests this GWAS was reverse coded 
an_dat <- format_data(an, 
                        ncontrol_col = "NCON",
                        ncase_col = "NCAS",
                        type="outcome",
                        snps= ea_meta_exp_clumped_an$SNP,
                        snp_col = "ID",
                        beta_col = "BETA",
                        se_col = "SE",
                        effect_allele_col = "ALT",
                        other_allele_col = "REF",
                        pval_col = "PVAL",
                        )
an_dat$outcome <- "ANNO"

an_dat.wf <- format_data(an, 
                      ncontrol_col = "NCON",
                      ncase_col = "NCAS",
                      type="outcome",
                      snps= ea_wf_exp_dat_clumped_an$SNP,
                      snp_col = "ID",
                      beta_col = "BETA",
                      se_col = "SE",
                      effect_allele_col = "ALT",
                      other_allele_col = "REF",
                      pval_col = "PVAL",
)
an_dat.wf$outcome <- "ANNO"

rm(an)

#** 2.3.3 Harmonize exposure and outcome data #################

dat.an <- harmonise_data(
  exposure_dat = ea_meta_exp_clumped_an, 
  outcome_dat = an_dat, action = 2
)
head(dat.an)
nb.snps.an <- nrow(dat.an)

dat.an.wf <- harmonise_data(
  exposure_dat = ea_wf_exp_dat_clumped_an, 
  outcome_dat = an_dat.wf, action = 2
)
head(dat.an.wf)
nb.snps.an <- nrow(dat.an.wf)

#** 2.3.4 MR analyses Anorexia Nervosa #############
an.results <- run_MR_analyses(dat.an)
an.wf.results <- run_MR_analyses(dat.an.wf)

save(an.results, an.wf.results, file = "an.MR.results.2209.RData")
load("an.MR.results.2209.RData")

#* 2.4 Anxiety Disorder #################

#Purves metanalysis 

anx <- fread("../Summary_statistics/Anxiety/META_UKBB_iPSYCH_ANGST_wNcol.sumstats")
head(anx)

#** 2.4.1. Clump EA specifically for Anxiety ######

ea_meta_exp_dat_anx <- ea_meta_exp_dat[ea_meta_exp_dat$SNP %in% anx$SNP,]


ea_meta_exp_clumped_anx <- clump_data(ea_meta_exp_dat_anx,
                                       clump_kb = 1000,
                                       clump_r2 = 0.001,
                                       clump_p1 = 1,
                                       clump_p2 = 1,
                                       pop = "EUR"
)

ea_wf_exp_dat_anx <- ea_wf_exp_dat[ea_wf_exp_dat$SNP %in% anx$SNP,]

ea_wf_exp_dat_clumped_anx <- clump_data(ea_wf_exp_dat_anx,
                                       clump_kb = 1000, 
                                       clump_r2 = 0.001, 
                                       clump_p1 = 1,
                                       clump_p2 = 1,
                                       pop = "EUR"
)




#** 2.4.2 Extract SNPS for outcome and format: ############

anx_dat <- format_data(anx, 
                       samplesize_col = "TotalSampleSize",
                       type="outcome",
                       snps= ea_meta_exp_clumped_anx$SNP,
                       snp_col = "SNP",
                       beta_col = "Effect",
                       se_col = "StdErr",
                       effect_allele_col = "Allele1",
                       other_allele_col = "Allele2",
                       pval_col = "P",
)
anx_dat$outcome <- "ANX"
#warning   Duplicated SNPs present in exposure data for phenotype 'outcome. Just keeping the first instance:
#rs7040995

anx_dat.wf <- format_data(anx, 
                       samplesize_col = "TotalSampleSize",
                       type="outcome",
                       snps= ea_wf_exp_dat_clumped_anx$SNP,
                       snp_col = "SNP",
                       beta_col = "Effect",
                       se_col = "StdErr",
                       effect_allele_col = "Allele1",
                       other_allele_col = "Allele2",
                       pval_col = "P",
)
anx_dat.wf$outcome <- "ANX"

rm(anx)



#** 2.4.3 Harmonize exposure and outcome data #################
dat.anx <- harmonise_data(
  exposure_dat = ea_meta_exp_clumped_anx, 
  outcome_dat = anx_dat, action = 2
)
head(dat.anx)

dat.anx.wf <- harmonise_data(
  exposure_dat = ea_wf_exp_dat_clumped_anx, 
  outcome_dat = anx_dat.wf, action = 2
)
head(dat.anx.wf )

#** 2.4.4 MR analyses Anxiety Disorder #############
anx.results <- run_MR_analyses(dat.anx)
anx.wf.results <- run_MR_analyses(dat.anx.wf)

save(anx.results, anx.wf.results, file = "anx.MR.results.2209.RData")
load("anx.MR.results.2209.RData")

#* 2.5 Bipolar Disorder #################

# #get header where readme is
# header <- file("../Summary_statistics/Bipolar/pgc-bip2021-all.vcf.tsv.gz")
# header <- readLines(header, n=100)
# sink("../Summary_statistics/Bipolar/HeadersOf_pgc-bip2021-all.vcf.tsv.txt")
# cat(paste0(header, collapse="\n"))
# sink()


bip <- fread("../Summary_statistics/Bipolar/pgc-bip2021-all.vcf.tsv.gz")
#write.table(bip, "pgc-bip2021-all.vcf.tsv.txt") #resave data to use with genomic sem munging
head(bip)

#** 2.5.1. Clump EA specifically for BIP ######

ea_meta_exp_dat_bip <- ea_meta_exp_dat[ea_meta_exp_dat$SNP %in% bip$ID,]


ea_meta_exp_clumped_bip <- clump_data(ea_meta_exp_dat_bip,
                                      clump_kb = 1000,
                                      clump_r2 = 0.001,
                                      clump_p1 = 1,
                                      clump_p2 = 1,
                                      pop = "EUR"
)

ea_wf_exp_dat_bip <- ea_wf_exp_dat[ea_wf_exp_dat$SNP %in% bip$ID,]

ea_wf_exp_dat_clumped_bip <- clump_data(ea_wf_exp_dat_bip,
                                        clump_kb = 1000, 
                                        clump_r2 = 0.001, 
                                        clump_p1 = 1,
                                        clump_p2 = 1,
                                        pop = "EUR"
)




#** 2.5.2 Extract SNPS for outcome and format: ############

bip_dat <- format_data(bip,
                       ncase_col = "NCAS",
                       ncontrol_col = "NCOL",
                       type="outcome",
                       snps= ea_meta_exp_clumped_bip$SNP,
                       snp_col = "ID",
                       beta_col = "BETA",
                       se_col = "SE",
                       effect_allele_col = "A1",
                       other_allele_col = "A2",
                       pval_col = "PVAL"
)
bip_dat$outcome <- "BIP"

bip_dat.wf <- format_data(bip,
                       ncase_col = "NCAS",
                       ncontrol_col = "NCOL",
                       type="outcome",
                       snps= ea_wf_exp_dat_clumped_bip$SNP,
                       snp_col = "ID",
                       beta_col = "BETA",
                       se_col = "SE",
                       effect_allele_col = "A1",
                       other_allele_col = "A2",
                       pval_col = "PVAL"
)
bip_dat.wf$outcome <- "BIP"

rm(bip)

#** 2.5.3 Harmonize exposure and outcome data #################

dat.bip <- harmonise_data(
  exposure_dat = ea_meta_exp_clumped_bip, 
  outcome_dat = bip_dat, action = 2
)
head(dat.bip)

dat.bip.wf <- harmonise_data(
  exposure_dat = ea_wf_exp_dat_clumped_bip, 
  outcome_dat = bip_dat.wf, action = 2
)
head(dat.bip.wf )

#** 2.5.4 MR analyses for Bipolar disorder #############
bip.results <- run_MR_analyses(dat.bip)
bip.wf.results <- run_MR_analyses(dat.bip.wf)

save(bip.results, bip.wf.results, file = "bip.MR.results.2209.RData")
load("bip.MR.results.2209.RData")

#* 2.6 Bipolar disorder Type-1 #################

#read header to get read me
# header <- file("../Summary_statistics/Bipolar/pgc-bip2021-BDI.vcf.tsv.gz")
# header <- readLines(header, n=100)
# sink("../Summary_statistics//Bipolar/HeadersOf_pgc-bip2021-BDI.vcf.tsv.txt")
# cat(paste0(header, collapse="\n"))
# sink()

bip1 <- fread("../Summary_statistics/bipolar/pgc-bip2021-BDI.vcf.tsv.gz")
#write.table(bip1, "pgc-bip2021-BDI.vcf.tsv.txt") #resave data to use with genomic sem munging
head(bip1)

#** 2.6.1. Clump EA specifically for BIP1 ######

ea_meta_exp_dat_bip1 <- ea_meta_exp_dat[ea_meta_exp_dat$SNP %in% bip1$ID,]


ea_meta_exp_clumped_bip1 <- clump_data(ea_meta_exp_dat_bip1,
                                      clump_kb = 1000,
                                      clump_r2 = 0.001,
                                      clump_p1 = 1,
                                      clump_p2 = 1,
                                      pop = "EUR"
)

ea_wf_exp_dat_bip1 <- ea_wf_exp_dat[ea_wf_exp_dat$SNP %in% bip1$ID,]

ea_wf_exp_dat_clumped_bip1 <- clump_data(ea_wf_exp_dat_bip1,
                                        clump_kb = 1000, 
                                        clump_r2 = 0.001, 
                                        clump_p1 = 1,
                                        clump_p2 = 1,
                                        pop = "EUR"
)




#** 2.6.2 Extract SNPS for outcome and format: ############
bip1_dat <- format_data(bip1,
                       ncase_col = "NCAS",
                       ncontrol_col = "NCOL",
                       type="outcome",
                       snps= ea_meta_exp_clumped_bip1$SNP,
                       snp_col = "ID",
                       beta_col = "BETA",
                       se_col = "SE",
                       effect_allele_col = "A1",
                       other_allele_col = "A2",
                       pval_col = "PVAL"
)
bip1_dat$outcome <- "bip1"

bip1_dat.wf <- format_data(bip1,
                        ncase_col = "NCAS",
                        ncontrol_col = "NCOL",
                        type="outcome",
                        snps= ea_wf_exp_dat_clumped_bip1$SNP,
                        snp_col = "ID",
                        beta_col = "BETA",
                        se_col = "SE",
                        effect_allele_col = "A1",
                        other_allele_col = "A2",
                        pval_col = "PVAL"
)
bip1_dat.wf$outcome <- "bip1"


rm(bip1)


#** 2.6.3 Harmonize exposure and outcome data #################
dat.bip1 <- harmonise_data(
  exposure_dat = ea_meta_exp_clumped_bip1, 
  outcome_dat = bip1_dat, action = 2
)
head(dat.bip1)

dat.bip1.wf <- harmonise_data(
  exposure_dat = ea_wf_exp_dat_clumped_bip1, 
  outcome_dat = bip1_dat.wf, action = 2
)
head(dat.bip1.wf )

#** 2.6.4 MR analyses Bipolar Disorder Type 1 #############
bip1.results <- run_MR_analyses(dat.bip1)
bip1.wf.results <- run_MR_analyses(dat.bip1.wf)

save(bip1.results, bip1.wf.results, file = "bip1.MR.results.2209.RData")
load("bip1.MR.results.2209.RData")

###* 2.7 Bipolar Disorder type 2 #################

# #read header 
# header <- file("../Summary_statistics/Bipolar/pgc-bip2021-BDII.vcf.tsv.gz")
# header <- readLines(header, n=100)
# sink("../Summary_statistics//Bipolar/HeadersOf_pgc-bip2021-BDII.vcf.tsv.txt")
# cat(paste0(header, collapse="\n"))
# sink()


bip2 <- fread("../Summary_statistics/bipolar/pgc-bip2021-BDII.vcf.tsv.gz")
#write.table(bip2, "pgc-bip2021-BDII.vcf.tsv.txt") #resave data to use with genomic sem munging
head(bip2)

#** 2.7.1. Clump EA specifically for BIP2 ######

ea_meta_exp_dat_bip2 <- ea_meta_exp_dat[ea_meta_exp_dat$SNP %in% bip2$ID,]


ea_meta_exp_clumped_bip2 <- clump_data(ea_meta_exp_dat_bip2,
                                      clump_kb = 1000,
                                      clump_r2 = 0.001,
                                      clump_p1 = 1,
                                      clump_p2 = 1,
                                      pop = "EUR"
)

ea_wf_exp_dat_bip2 <- ea_wf_exp_dat[ea_wf_exp_dat$SNP %in% bip2$ID,]

ea_wf_exp_dat_clumped_bip2 <- clump_data(ea_wf_exp_dat_bip2,
                                         clump_kb = 1000, 
                                         clump_r2 = 0.001, 
                                         clump_p1 = 1,
                                         clump_p2 = 1,
                                         pop = "EUR"
)


#** 2.7.2 Extract SNPS for outcome and format: ############

bip2_dat <- format_data(bip2,
                        ncase_col = "NCAS",
                        ncontrol_col = "NCOL",
                        type="outcome",
                        snps= ea_meta_exp_clumped_bip2$SNP,
                        snp_col = "ID",
                        beta_col = "BETA",
                        se_col = "SE",
                        effect_allele_col = "A1",
                        other_allele_col = "A2",
                        pval_col = "PVAL"
)
bip2_dat$outcome <- "bip2"

bip2_dat.wf <- format_data(bip2,
                        ncase_col = "NCAS",
                        ncontrol_col = "NCOL",
                        type="outcome",
                        snps= ea_wf_exp_dat_clumped_bip2$SNP,
                        snp_col = "ID",
                        beta_col = "BETA",
                        se_col = "SE",
                        effect_allele_col = "A1",
                        other_allele_col = "A2",
                        pval_col = "PVAL"
)
bip2_dat.wf$outcome <- "bip2"

rm(bip2)
gc()

#** 2.7.3 Harmonize exposure and outcome data #################
dat.bip2 <- harmonise_data(
  exposure_dat = ea_meta_exp_clumped_bip2, 
  outcome_dat = bip2_dat, action = 2
)
head(dat.bip2)

dat.bip2.wf <- harmonise_data(
  exposure_dat = ea_wf_exp_dat_clumped_bip2, 
  outcome_dat = bip2_dat.wf, action = 2
)
head(dat.bip2.wf )

#** 2.7.4 MR analyses Bipolar Disorder Type 2 #############
bip2.results <- run_MR_analyses(dat.bip2)
bip2.wf.results <- run_MR_analyses(dat.bip2.wf)

save(bip2.results, bip2.wf.results, file = "bip2.MR.results.2209.RData")
load("bip2.MR.results.2209.RData")

#* 2.8 MDD #################

# Howard gwas
mdd <- fread("../Summary_statistics/MDD/PGC_UKB_depression_genome-wide.txt") 
head(mdd)

#** 2.8.1. Clump EA specifically for MDD ######

ea_meta_exp_dat_mdd <- ea_meta_exp_dat[ea_meta_exp_dat$SNP %in% mdd$MarkerName,]


ea_meta_exp_clumped_mdd <- clump_data(ea_meta_exp_dat_mdd,
                                      clump_kb = 1000,
                                      clump_r2 = 0.001,
                                      clump_p1 = 1,
                                      clump_p2 = 1,
                                      pop = "EUR"
)

ea_wf_exp_dat_mdd <- ea_wf_exp_dat[ea_wf_exp_dat$SNP %in% mdd$MarkerName,]
ea_wf_exp_dat_clumped_mdd <- clump_data(ea_wf_exp_dat_mdd,
                                         clump_kb = 1000, 
                                         clump_r2 = 0.001, 
                                         clump_p1 = 1,
                                         clump_p2 = 1,
                                         pop = "EUR"
)

#** 2.8.2 Extract SNPS for outcome and format: ############

mdd_dat <- format_data(mdd,
                        type="outcome",
                        snps= ea_meta_exp_clumped_mdd$SNP,
                        snp_col = "MarkerName",
                        beta_col = "LogOR",
                        se_col = "StdErrLogOR",
                        effect_allele_col = "A1",
                        other_allele_col = "A2",
                        pval_col = "P",
                       eaf= "Freq"
)
mdd_dat$outcome <- "MDD"

mdd_dat.wf <- format_data(mdd,
                       type="outcome",
                       snps= ea_wf_exp_dat_clumped_mdd$SNP,
                       snp_col = "MarkerName",
                       beta_col = "LogOR",
                       se_col = "StdErrLogOR",
                       effect_allele_col = "A1",
                       other_allele_col = "A2",
                       pval_col = "P",
                       eaf= "Freq"
)
mdd_dat.wf$outcome <- "MDD"
  
rm(mdd)
gc()

#** 2.8.3 Harmonize exposure and outcome data #################

dat.mdd <- harmonise_data(
  exposure_dat = ea_meta_exp_clumped_mdd, 
  outcome_dat = mdd_dat, action = 2
)
head(dat.mdd) 

dat.mdd.wf <- harmonise_data(
  exposure_dat = ea_wf_exp_dat_clumped_mdd, 
  outcome_dat = mdd_dat.wf, action = 2
)
head(dat.mdd.wf )

#** 2.8.4 MR analyses for MR #############
mdd.results <- run_MR_analyses(dat.mdd)
mdd.wf.results <- run_MR_analyses(dat.mdd.wf)

save(mdd.results, mdd.wf.results, file = "mdd.MR.results.2209.RData")
load("mdd.MR.results.2209.RData")

#* 2.9 OCD #################

ocd <- fread("../Summary_statistics/OCD/ocd_aug2017.gz")
head(ocd)

ocd$beta <- log(ocd$OR)

#** 2.9.1. Clump EA specifically for OCD ######

ea_meta_exp_dat_ocd <- ea_meta_exp_dat[ea_meta_exp_dat$SNP %in% ocd$SNP,]


ea_meta_exp_clumped_ocd <- clump_data(ea_meta_exp_dat_ocd,
                                      clump_kb = 1000,
                                      clump_r2 = 0.001,
                                      clump_p1 = 1,
                                      clump_p2 = 1,
                                      pop = "EUR"
)

ea_wf_exp_dat_ocd <- ea_wf_exp_dat[ea_wf_exp_dat$SNP %in% ocd$SNP,]
ea_wf_exp_dat_clumped_ocd <- clump_data(ea_wf_exp_dat_ocd,
                                        clump_kb = 1000, 
                                        clump_r2 = 0.001, 
                                        clump_p1 = 1,
                                        clump_p2 = 1,
                                        pop = "EUR"
)


#** 2.9.2 Extract SNPS for outcome and format: ############

ocd_dat <- format_data(ocd,
                       type="outcome",
                       snps= ea_meta_exp_clumped_ocd$SNP,
                       snp_col = "SNP",
                       beta_col = "beta",
                       se_col = "SE",
                       effect_allele_col = "A1",
                       other_allele_col = "A2",
                       pval_col = "P"
)
ocd_dat$outcome <- "ocd"

ocd_dat.wf <- format_data(ocd,
                       type="outcome",
                       snps= ea_wf_exp_dat_clumped_ocd$SNP,
                       snp_col = "SNP",
                       beta_col = "beta",
                       se_col = "SE",
                       effect_allele_col = "A1",
                       other_allele_col = "A2",
                       pval_col = "P"
)
ocd_dat.wf$outcome <- "ocd"

rm(ocd)

#** 2.9.3 Harmonize exposure and outcome data #################
dat.ocd <- harmonise_data(
  exposure_dat = ea_meta_exp_clumped_ocd, 
  outcome_dat = ocd_dat, action = 2
)
head(dat.ocd) ## fairly steep SNP dropout....

dat.ocd.wf <- harmonise_data(
  exposure_dat = ea_wf_exp_dat_clumped_ocd, 
  outcome_dat = ocd_dat.wf, action = 2
)
head(dat.ocd.wf )

#** 2.9.1 MR analyses for OCD #############
ocd.results <- run_MR_analyses(dat.ocd)
ocd.wf.results <- run_MR_analyses(dat.ocd.wf)

save(ocd.results, ocd.wf.results, file = "ocd.MR.results.2209.RData")
load("ocd.MR.results.2209.RData")

#* 2.10 PTSD #################

ptsd <- fread("../Summary_statistics/pts_eur_freeze2_overall.results")
head(ptsd)

ptsd$beta <- log(ptsd$OR)

#** 2.10.1. Clump EA specifically for PTSD ######

ea_meta_exp_dat_ptsd <- ea_meta_exp_dat[ea_meta_exp_dat$SNP %in% ptsd$SNP,]


ea_meta_exp_clumped_ptsd <- clump_data(ea_meta_exp_dat_ptsd,
                                      clump_kb = 1000,
                                      clump_r2 = 0.001,
                                      clump_p1 = 1,
                                      clump_p2 = 1,
                                      pop = "EUR"
)

ea_wf_exp_dat_ptsd <- ea_wf_exp_dat[ea_wf_exp_dat$SNP %in% ptsd$SNP,]
ea_wf_exp_dat_clumped_ptsd <- clump_data(ea_wf_exp_dat_ptsd,
                                        clump_kb = 1000, 
                                        clump_r2 = 0.001, 
                                        clump_p1 = 1,
                                        clump_p2 = 1,
                                        pop = "EUR"
)


#** 2.10.2 Extract SNPS for outcome and format: ############

ptsd_dat <- format_data(ptsd,
                        ncase_col = "Nca",
                        ncontrol_col = "Ncol",
                        type="outcome",
                        snps= ea_meta_exp_clumped_ptsd$SNP,
                        snp_col = "SNP",
                        beta_col = "beta",
                        se_col = "SE",
                        effect_allele_col = "A1",
                        other_allele_col = "A2",
                        pval_col = "P"
)
ptsd_dat$outcome <- "ptsd"

ptsd_dat.wf <- format_data(ptsd,
                        ncase_col = "Nca",
                        ncontrol_col = "Ncol",
                        type="outcome",
                        snps= ea_wf_exp_dat_clumped_ptsd$SNP,
                        snp_col = "SNP",
                        beta_col = "beta",
                        se_col = "SE",
                        effect_allele_col = "A1",
                        other_allele_col = "A2",
                        pval_col = "P"
)
ptsd_dat.wf$outcome <- "ptsd"

rm(ptsd)

#** 2.10.3 Harmonize exposure and outcome data #################
dat.ptsd <- harmonise_data(
  exposure_dat = ea_meta_exp_clumped_ptsd, 
  outcome_dat = ptsd_dat, action = 2
)
head(dat.ptsd) ## fairly steep SNP dropout....

dat.ptsd.wf <- harmonise_data(
  exposure_dat = ea_wf_exp_dat_clumped_ptsd, 
  outcome_dat = ptsd_dat.wf, action = 2
)
head(dat.ptsd.wf)

#** 2.10.4 MR analyses for PTSD #############
ptsd.results <- run_MR_analyses(dat.ptsd)
ptsd.wf.results <- run_MR_analyses(dat.ptsd.wf)

save(ptsd.results, ptsd.wf.results, file = "ptsd.MR.results.2209.RData")
load("ptsd.MR.results.2209.RData")

#* 2.11 Schizophrenia  #################

# scz2 <- fread("../Summary_statistics/SCZ/CLOZUK_PGC2noclo.METAL.assoc.dosage.fix") # from pgc website
# 
# scz <- fread("../Summary_statistics/SCZ/SCZ_clozuk_pgc2.meta.sumstats.txt.gz") # from walter website migth be clozuk extra
# head(scz2)
# 
# #PGC doc is probably only PGC2, while the sumstats on Walter group website contains CLOZUK 
# #PGC sumstats have 110000 less SNPS, but they have RSID. Walter sumstats have freq od A1 


#SCZ3 on the PGC website, but no readme 
scz <- fread("../Summary_statistics/SCZ/PGC3_SCZ_wave3_public.v2.tsv", fill=T) 
head(scz)

scz$beta <- log(scz$OR)

#** 2.10.1. Clump EA specifically for SCZ ######

ea_meta_exp_dat_scz <- ea_meta_exp_dat[ea_meta_exp_dat$SNP %in% scz$SNP,]


ea_meta_exp_clumped_scz <- clump_data(ea_meta_exp_dat_scz,
                                       clump_kb = 1000,
                                       clump_r2 = 0.001,
                                       clump_p1 = 1,
                                       clump_p2 = 1,
                                       pop = "EUR"
)

ea_wf_exp_dat_scz <- ea_wf_exp_dat[ea_wf_exp_dat$SNP %in% scz$SNP,]
ea_wf_exp_dat_clumped_scz <- clump_data(ea_wf_exp_dat_scz,
                                         clump_kb = 1000, 
                                         clump_r2 = 0.001, 
                                         clump_p1 = 1,
                                         clump_p2 = 1,
                                         pop = "EUR"
)



#** 2.11.2 Extract SNPS for outcome and format: ############

scz_dat <- format_data(scz,
                       type="outcome",
                       snps= ea_meta_exp_clumped_scz$SNP,
                       snp_col = "SNP",
                       beta_col = "beta",
                       se_col = "SE",
                       effect_allele_col = "A1",
                       other_allele_col = "A2",
                       pval_col = "P"
)
scz_dat$outcome <- "scz"


scz_dat.wf <- format_data(scz,
                       type="outcome",
                       snps= ea_wf_exp_dat_clumped_scz$SNP,
                       snp_col = "SNP",
                       beta_col = "beta",
                       se_col = "SE",
                       effect_allele_col = "A1",
                       other_allele_col = "A2",
                       pval_col = "P"
)
scz_dat.wf$outcome <- "scz"


rm(scz)
gc()

#** 2.11.3 Harmonize exposure and outcome data #################
dat.scz <- harmonise_data(
  exposure_dat = ea_meta_exp_clumped_scz, 
  outcome_dat = scz_dat, action = 2
)
head(dat.scz) ## fairly steep SNP dropout....

dat.scz.wf <- harmonise_data(
  exposure_dat = ea_wf_exp_dat_clumped_scz, 
  outcome_dat = scz_dat.wf, action = 2
)
head(dat.scz.wf )

#** 2.11.1 MR analyses Schizophrenia  #############
scz.results <- run_MR_analyses(dat.scz)
scz.wf.results <- run_MR_analyses(dat.scz.wf)

save(scz.results, scz.wf.results, file = "scz.MR.results.2209.RData")
load("scz.MR.results.2209.RData")

#** 2.11.2 Figures Schizophrenia ###########
res_single <- mr_singlesnp(dat.scz)
p4 <- mr_funnel_plot(res_single)
p4[[1]]

res  <- mr(dat.scz, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_mode", "mr_weighted_median"))
p1 <- mr_scatter_plot(res, dat.scz) 
p1[[1]]

res_loo <- mr_leaveoneout(dat.scz)
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]



#*2.12 Alcohol dependence ################
# alc <- fread("Alc/pgc_alcdep.eur_discovery.aug2018_release.txt.gz") 
# head(alc)
# alc <- alc %>% separate(SNP, c("rsID", NA), extra = "drop", fill = "right") # snp column contains pos and alleles, keep only rsid
#write.table(alc, "Alc/pgc_alcdep.eur_discovery.aug2018_release_rsid.txt", row.names = F)

alc <- fread("../Summary_statistics/Alc/pgc_alcdep.eur_discovery.aug2018_release_rsid.txt")

# missing frequency of the effect allele, take the one from reference 1000G
# only use it to compute beta and SE
ref <- fread("C:/Users/user/OneDrive - Vrije Universiteit Amsterdam/Documents/References_Data/reference.1000G.maf.0.005.txt.gz")
alc_ref <- merge(alc, ref, by.x="rsID", by.y="SNP", suffixes= c("alc", "ref"))
head(alc_ref)
alc_ref$MAF2 <- ifelse(alc_ref$A1alc == alc_ref$A2ref & alc_ref$A2alc == alc_ref$A1ref, 1-alc_ref$MAF, alc_ref$MAF)
alc_ref$A1ref2 <- ifelse(alc_ref$A1alc == alc_ref$A2ref & alc_ref$A2alc == alc_ref$A1ref, alc_ref$A2ref, alc_ref$A1ref)
alc_ref$A2ref2 <- ifelse(alc_ref$A1alc == alc_ref$A2ref & alc_ref$A2alc == alc_ref$A1ref, alc_ref$A1ref, alc_ref$A2ref)
alc_ref$diffferentallele <- ifelse(alc_ref$A1alc == alc_ref$A1ref2 & alc_ref$A2alc == alc_ref$A2ref2, T, F)
table(alc_ref$diffferentallele) # all are true
alc <- alc_ref[,c("rsID", "CHRalc", "BPalc", "A1alc", "A2alc", "Z", "P", "Weight", "MAF2")]
names(alc) <- c("rsID", "CHR", "BP", "A1", "A2", "Z", "P", "Weight", "MAF")

alc$beta <- alc$Z / sqrt(2*alc$MAF*(1 - alc$MAF)*(alc$Weight+ alc$Z*alc$Z))  
alc$SE <- 1/sqrt(2*alc$MAF*(1 - alc$MAF)*(alc$Weight+ alc$Z*alc$Z))  

head(alc)
rm(ref)

#write.table(alc, "../Summary_statistics/Alc/pgc_alcdep.eur_discovery.aug2018_release_rsid_betase.txt", row.names = F)

#** 2.12.1. Clump EA specifically for alcohol ######

ea_meta_exp_dat_alc <- ea_meta_exp_dat[ea_meta_exp_dat$SNP %in% alc$rsID,]


ea_meta_exp_clumped_alc <- clump_data(ea_meta_exp_dat_alc,
                                      clump_kb = 1000,
                                      clump_r2 = 0.001,
                                      clump_p1 = 1,
                                      clump_p2 = 1,
                                      pop = "EUR"
)

ea_wf_exp_dat_alc <- ea_wf_exp_dat[ea_wf_exp_dat$SNP %in% alc$rsID,]
ea_wf_exp_dat_clumped_alc <- clump_data(ea_wf_exp_dat_alc,
                                        clump_kb = 1000, 
                                        clump_r2 = 0.001, 
                                        clump_p1 = 1,
                                        clump_p2 = 1,
                                        pop = "EUR"
)



#** 2.11.2 Extract SNPS for outcome and format: ############

alc_dat <- format_data(alc,
                       type="outcome",
                       snps= ea_meta_exp_clumped_alc$SNP,
                       snp_col = "rsID",
                       beta_col = "beta",
                       se_col = "SE",
                       effect_allele_col = "A1",
                       other_allele_col = "A2",
                       pval_col = "P"
)
alc_dat$outcome <- "alc"


alc_dat.wf <- format_data(alc,
                          type="outcome",
                          snps= ea_wf_exp_dat_clumped_alc$SNP,
                          snp_col = "rsID",
                          beta_col = "beta",
                          se_col = "SE",
                          effect_allele_col = "A1",
                          other_allele_col = "A2",
                          pval_col = "P"
)
alc_dat.wf$outcome <- "alc"


#** 2.11.3 Harmonize exposure and outcome data #################
dat.alc <- harmonise_data(
  exposure_dat = ea_meta_exp_clumped_alc, 
  outcome_dat = alc_dat, action = 2
)

dat.alc.wf <- harmonise_data(
  exposure_dat = ea_wf_exp_dat_clumped_alc, 
  outcome_dat = alc_dat.wf, action = 2
)


#** 2.11.1 MR analyses alcohol #############
alc.results <- run_MR_analyses(dat.alc)
alc.wf.results <- run_MR_analyses(dat.alc.wf)

save(alc.results, alc.wf.results, file = "alc.MR.results.2209.RData")
load("alc.MR.results.2209.RData")

#*3. Combine and save all results ######

all.results <- rbind(asd.results[[1]], 
                     adhd.results[[1]], 
                     adhd2023.results[[1]],
                     an.results[[1]], 
                     anx.results[[1]], 
                     bip.results[[1]], 
                     bip1.results[[1]],
                     bip2.results[[1]], 
                     mdd.results[[1]], 
                     ocd.results[[1]],
                     ptsd.results[[1]], 
                     scz.results[[1]],
                     alc.results[[1]]
                     )
iv.all.results <- all.results[all.results$method == "Inverse variance weighted",]
ggplot(data= iv.all.results, aes(x=outcome, y=or))+ 
  geom_col()+ 
  geom_errorbar(ymin=iv.all.results$or_lci95, ymax=iv.all.results$or_uci95)


all.wf.results <- rbind(asd.wf.results[[1]], 
                     adhd.wf.results[[1]], 
                     adhd2023.wf.results[[1]], 
                     an.wf.results[[1]], 
                     anx.wf.results[[1]], 
                     bip.wf.results[[1]], 
                     bip1.wf.results[[1]],
                     bip2.wf.results[[1]], 
                     mdd.wf.results[[1]], 
                     ocd.wf.results[[1]],
                     ptsd.wf.results[[1]], 
                     scz.wf.results[[1]],
                     alc.wf.results[[1]]
)
iv.all.wf.results <- all.wf.results[all.wf.results$method == "Inverse variance weighted",]
ggplot(data= iv.all.wf.results, aes(x=outcome, y=or))+ 
  geom_col()+ 
  geom_errorbar(ymin=iv.all.wf.results$or_lci95, ymax=iv.all.wf.results$or_uci95)


head(iv.all.wf.results)
iv.all.wf.results$GWAS <- "within"
iv.all.results$GWAS <- "population"

iv.all.both.results <- rbind(iv.all.results, iv.all.wf.results)
iv.all.both.results$GWAS<- factor(iv.all.both.results$GWAS, levels = c("within","population"))

save(iv.all.both.results, file="ivw.all.both.results.RData")




dotCOLS = c("#a6d8f0", "#f9b282")
barCOLS= c("#008fd5", "#de6b35")
p <- ggplot(iv.all.both.results, aes(x=outcome, y=or, ymin=or_lci95, ymax=or_uci95, col=GWAS, fill=GWAS)) + 
  geom_linerange(size=9, position=position_dodge(width=0.5))+
  geom_hline(yintercept=1, lty=2)+
  geom_point(size=7,shape=21,colour="white", stroke=0.5, position=position_dodge(width=0.5))+
  scale_fill_manual(values=barCOLS)+
  scale_color_manual(values=dotCOLS)+
  #scale_x_discrete(limits = rev(order_traits), name="Disorders")+ #
  scale_y_continuous(name= "Odds ratio")+ 
  coord_flip()+ 
  theme_minimal(base_size = 22)+
  theme(panel.grid.minor=element_blank(), legend.title=element_blank(), axis.title.y=element_blank())
p


all.sensitivity <- rbind(asd.results[[2]], 
                         adhd.results[[2]], 
                         adhd2023.results[[2]], 
                         an.results[[2]], 
                         anx.results[[2]], 
                         bip.results[[2]], 
                         bip1.results[[2]],
                         bip2.results[[2]], 
                         mdd.results[[2]], 
                         ocd.results[[2]],
                         ptsd.results[[2]], 
                         scz.results[[2]], 
                         alc.results[[2]]
)

all.Q <- rbind(asd.results[[3]], 
               adhd.results[[3]], 
               adhd2023.results[[3]], 
               an.results[[3]], 
               anx.results[[3]], 
               bip.results[[3]], 
               bip1.results[[3]],
               bip2.results[[3]], 
               mdd.results[[3]], 
               ocd.results[[3]],
               ptsd.results[[3]], 
               scz.results[[3]], 
               alc.results[[3]]
)

all.wf.sensitivity <- rbind(asd.wf.results[[2]], 
                         adhd.wf.results[[2]], 
                         adhd2023.wf.results[[2]], 
                         an.wf.results[[2]], 
                         anx.wf.results[[2]], 
                         bip.wf.results[[2]], 
                         bip1.wf.results[[2]],
                         bip2.wf.results[[2]], 
                         mdd.wf.results[[2]], 
                         ocd.wf.results[[2]],
                         ptsd.wf.results[[2]], 
                         scz.wf.results[[2]],
                         alc.wf.results[[2]]
)

all.wf.Q <- rbind(asd.wf.results[[3]], 
               adhd.wf.results[[3]], 
               adhd2023.wf.results[[3]], 
               an.wf.results[[3]], 
               anx.wf.results[[3]], 
               bip.wf.results[[3]], 
               bip1.wf.results[[3]],
               bip2.wf.results[[3]], 
               mdd.wf.results[[3]], 
               ocd.wf.results[[3]],
               ptsd.wf.results[[3]], 
               scz.wf.results[[3]],
               alc.wf.results[[3]]
)

write.csv(all.results, "all_MR_results_221020.csv")
write.csv(all.wf.results, "all_MR_wf_results_221020.csv")
write.csv(all.sensitivity, "all_MR_sensitivity_221020.csv")
write.csv(all.Q, "all_MR_Q_221020.csv")
write.csv(all.wf.sensitivity, "all_MR_wf_sensitivity_221020.csv")
write.csv(all.wf.Q, "all_MR_wf_Q_221020.csv")


all.results.wide <- reshape(all.results, idvar = "outcome" , timevar = "method", direction = "wide")
write.csv(all.results.wide, "all_MR_results_wide_221020.csv")

all.wf.results.wide <- reshape(all.wf.results, idvar = "outcome" , timevar = "method", direction = "wide")
write.csv(all.wf.results.wide, "all_MR_wf_results_wide_221020.csv")


all.Q.wide <- reshape(all.Q, idvar = "outcome" , timevar = "method", direction = "wide")
write.csv(all.Q.wide, "all_MR_Q_wide_221020.csv")


all.wf.Q.wide <- reshape(all.wf.Q , idvar = "outcome" , timevar = "method", direction = "wide")
write.csv(all.wf.Q.wide, "all_MR_wf_Q_wide_221020.csv")

