## Project: EA_MH_CBS_MR 
## Script purpose: Run MR analyses, mental health as exposure 
## Author: Perline Demange & Michel Nivard


# Libraries #####
library(data.table)
library(tidyverse)
library(readxl)

library(devtools)
#install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR) # https://mrcieu.github.io/TwoSampleMR/index.html 


# Load Exposure ##################
#can not run everything at once due to memory issues, split into several groups

asd <- fread("../Summary_statistics/ASD_iPSYCH-PGC_ASD_Nov2017.gz")
adhd <- fread("../Summary_statistics/ADHD/daner_adhd_meta_filtered_NA_iPSYCH23_PGC11_sigPCs_woSEX_2ell6sd_EUR_Neff_70.meta.gz")
an <- fread("../Summary_statistics/Anorexia_pgcAN2.2019-07.vcf.tsv.gz")
anx <- fread("../Summary_statistics/Anxiety/META_UKBB_iPSYCH_ANGST_wNcol.sumstats")
bip <- fread("../Summary_statistics/Bipolar/pgc-bip2021-all.vcf.tsv.gz")
bip1 <- fread("../Summary_statistics/bipolar/pgc-bip2021-BDI.vcf.tsv.gz")
bip2 <- fread("../Summary_statistics/bipolar/pgc-bip2021-BDII.vcf.tsv.gz")
mdd <- fread("../Summary_statistics/MDD/PGC_UKB_depression_genome-wide.txt") 
ocd <- fread("../Summary_statistics/OCD/ocd_aug2017.gz")
ptsd <- fread("../Summary_statistics/pts_eur_freeze2_overall.results")
scz <- fread("../Summary_statistics/SCZ/PGC3_SCZ_wave3_public.v2.tsv", fill=T) 
alc <- fread("../Summary_statistics/Alc/pgc_alcdep.eur_discovery.aug2018_release_rsid_betase.txt")
adhd2023 <- fread("../Summary_statistics/ADHD/Demontis 2023/ADHD2022_iPSYCH_deCODE_PGC.meta.gz")


# convert the OR to log(OR) or the beta of a logistic regression if necessary #####
asd$beta <- log(asd$OR)
adhd$beta<- log(adhd$OR)
adhd2023$beta<- log(adhd2023$OR)
#an is beta already 
#anx is beta already 
#bip is beta already
#bip1 already
#bip2 already
#mdd already 
ocd$beta <- log(ocd$OR)
ptsd$beta <- log(ptsd$OR)
scz$beta <- log(scz$OR)
#alc beta is ready, built in script EA to MH 

# select instruments ########
# If the number of instruments after clumping is <5, we select suggestive hits instead 
asd_instrument <- asd[asd$P <= 0.00001, ] #4028
adhd_instrument <- adhd[adhd$P <= 0.00000005, ] #317
adhd2023_instrument <- adhd2023[adhd2023$P <= 0.00000005, ] #1428
an_instrument <- an[an$PVAL <= 0.00000005, ] # 327
anx_instrument <- anx[anx$P <= 0.00001, ] 
bip_instrument <- bip[bip$PVAL <= 0.00000005, ]#
bip1_instrument <- bip1[bip1$PVAL <= 0.00000005, ]
bip2_instrument <- bip2[bip2$PVAL <= 0.00001, ]
mdd_instrument <- mdd[mdd$P <= 0.00000005, ]
ocd_instrument <- ocd[ocd$P <= 0.00001, ]
ptsd_instrument <- ptsd[ptsd$P <= 0.00001, ]
scz_instrument <- scz[scz$P <= 0.00000005, ]
alc_instrument <- alc[alc$P <= 0.00001, ]
# Format ####
asd_exp_dat <- format_data(asd_instrument, 
                               type="exposure",   
                               snp_col = "SNP",
                               beta_col = "beta",
                               se_col = "SE",
                               effect_allele_col = "A1",
                               other_allele_col = "A2",
                               pval_col = "P")
adhd_exp_dat <- format_data(adhd_instrument, 
                            type="exposure",   
                            snp_col = "SNP",
                            beta_col = "beta",
                            se_col = "SE",
                            effect_allele_col = "A1",
                            other_allele_col = "A2",
                            pval_col = "P")
adhd2023_exp_dat <- format_data(adhd2023_instrument, 
                            type="exposure",   
                            snp_col = "SNP",
                            beta_col = "beta",
                            se_col = "SE",
                            effect_allele_col = "A1",
                            other_allele_col = "A2",
                            pval_col = "P")
an_exp_dat <- format_data(an_instrument, 
                      ncontrol_col = "NCON",
                      ncase_col = "NCAS",
                      type="exposure",
                      snp_col = "ID",
                      beta_col = "BETA",
                      se_col = "SE",
                      effect_allele_col = "ALT",
                      other_allele_col = "REF",
                      pval_col = "PVAL",
)
anx_exp_dat <- format_data(anx_instrument, 
                       samplesize_col = "TotalSampleSize",
                       type="exposure",
                       snp_col = "SNP",
                       beta_col = "Effect",
                       se_col = "StdErr",
                       effect_allele_col = "Allele1",
                       other_allele_col = "Allele2",
                       pval_col = "P",
)
bip_exp_dat <- format_data(bip_instrument,
                       ncase_col = "NCAS",
                       ncontrol_col = "NCOL",
                       type="exposure",
                       snp_col = "ID",
                       beta_col = "BETA",
                       se_col = "SE",
                       effect_allele_col = "A1",
                       other_allele_col = "A2",
                       pval_col = "PVAL"
)
bip1_exp_dat <- format_data(bip1_instrument,
                       ncase_col = "NCAS",
                       ncontrol_col = "NCOL",
                       type="exposure",
                       snp_col = "ID",
                       beta_col = "BETA",
                       se_col = "SE",
                       effect_allele_col = "A1",
                       other_allele_col = "A2",
                       pval_col = "PVAL"
)
bip2_exp_dat <- format_data(bip2_instrument,
                        ncase_col = "NCAS",
                        ncontrol_col = "NCOL",
                        type="exposure",
                        snp_col = "ID",
                        beta_col = "BETA",
                        se_col = "SE",
                        effect_allele_col = "A1",
                        other_allele_col = "A2",
                        pval_col = "PVAL"
)
mdd_exp_dat <- format_data(mdd_instrument,
                       type="exposure",
                       snp_col = "MarkerName",
                       beta_col = "LogOR",
                       se_col = "StdErrLogOR",
                       effect_allele_col = "A1",
                       other_allele_col = "A2",
                       pval_col = "P",
                       eaf= "Freq"
)
ocd_exp_dat <- format_data(ocd_instrument,
                       type="exposure",
                       snp_col = "SNP",
                       beta_col = "beta",
                       se_col = "SE",
                       effect_allele_col = "A1",
                       other_allele_col = "A2",
                       pval_col = "P"
)
ptsd_exp_dat <- format_data(ptsd_instrument,
                        ncase_col = "Nca",
                        ncontrol_col = "Ncol",
                        type="exposure",
                        snp_col = "SNP",
                        beta_col = "beta",
                        se_col = "SE",
                        effect_allele_col = "A1",
                        other_allele_col = "A2",
                        pval_col = "P"
)
scz_exp_dat <- format_data(scz_instrument,
                       type="exposure",
                       snp_col = "SNP",
                       beta_col = "beta",
                       se_col = "SE",
                       effect_allele_col = "A1",
                       other_allele_col = "A2",
                       pval_col = "P"
)
alc_exp_dat <- format_data(alc_instrument,
                           type="exposure",
                           snp_col = "rsID",
                           beta_col = "beta",
                           se_col = "SE",
                           effect_allele_col = "A1",
                           other_allele_col = "A2",
                           pval_col = "P"
)
asd_exp_dat$exposure <- "ASD"
adhd_exp_dat$exposure <- "ADHD"
adhd2023_exp_dat$exposure <- "ADHD2023"
an_exp_dat$exposure <- "ANNO"
anx_exp_dat$exposure <- "ANX"
bip_exp_dat$exposure <- "BIP"
bip1_exp_dat$exposure <- "BIP1"
bip2_exp_dat$exposure <- "BIP2"
mdd_exp_dat$exposure <- "MDD"
ocd_exp_dat$exposure <- "OCD"
ptsd_exp_dat$exposure <- "PTSD"
scz_exp_dat$exposure <- "SCZ"
alc_exp_dat$exposure <- "alc"
#clump #####
asd_exp_dat_clumped <- clump_data(asd_exp_dat,
                                        clump_kb = 1000, 
                                        clump_r2 = 0.001, 
                                        clump_p1 = 1,
                                        clump_p2 = 1,
                                        pop = "EUR"
)
#64
adhd_exp_dat_clumped <- clump_data(adhd_exp_dat,
                                   clump_kb = 1000, 
                                   clump_r2 = 0.001, 
                                   clump_p1 = 1,
                                   clump_p2 = 1,
                                   pop = "EUR"
) #11
adhd2023_exp_dat_clumped <- clump_data(adhd2023_exp_dat,
                                   clump_kb = 1000, 
                                   clump_r2 = 0.001, 
                                   clump_p1 = 1,
                                   clump_p2 = 1,
                                   pop = "EUR"
) #27
an_exp_dat_clumped <- clump_data(an_exp_dat,
                                   clump_kb = 1000, 
                                   clump_r2 = 0.001, 
                                   clump_p1 = 1,
                                   clump_p2 = 1,
                                   pop = "EUR"
) #9
anx_exp_dat_clumped <- clump_data(anx_exp_dat,
                                 clump_kb = 1000, 
                                 clump_r2 = 0.001, 
                                 clump_p1 = 1,
                                 clump_p2 = 1,
                                 pop = "EUR"
) 
bip_exp_dat_clumped <- clump_data(bip_exp_dat,
                                  clump_kb = 1000, 
                                  clump_r2 = 0.001, 
                                  clump_p1 = 1,
                                  clump_p2 = 1,
                                  pop = "EUR"
) 
bip1_exp_dat_clumped <- clump_data(bip1_exp_dat,
                                  clump_kb = 1000, 
                                  clump_r2 = 0.001, 
                                  clump_p1 = 1,
                                  clump_p2 = 1,
                                  pop = "EUR"
) 
bip2_exp_dat_clumped <- clump_data(bip2_exp_dat,
                                  clump_kb = 1000, 
                                  clump_r2 = 0.001, 
                                  clump_p1 = 1,
                                  clump_p2 = 1,
                                  pop = "EUR"
) 
mdd_exp_dat_clumped <- clump_data(mdd_exp_dat,
                                   clump_kb = 1000, 
                                   clump_r2 = 0.001, 
                                   clump_p1 = 1,
                                   clump_p2 = 1,
                                   pop = "EUR"
) 
ocd_exp_dat_clumped <- clump_data(ocd_exp_dat,
                                   clump_kb = 1000, 
                                   clump_r2 = 0.001, 
                                   clump_p1 = 1,
                                   clump_p2 = 1,
                                   pop = "EUR"
) 
ptsd_exp_dat_clumped <- clump_data(ptsd_exp_dat,
                                  clump_kb = 1000, 
                                  clump_r2 = 0.001, 
                                  clump_p1 = 1,
                                  clump_p2 = 1,
                                  pop = "EUR"
) 
scz_exp_dat_clumped <- clump_data(scz_exp_dat,
                                  clump_kb = 1000, 
                                  clump_r2 = 0.001, 
                                  clump_p1 = 1,
                                  clump_p2 = 1,
                                  pop = "EUR"
) 
alc_exp_dat_clumped <- clump_data(alc_exp_dat,
                                  clump_kb = 1000, 
                                  clump_r2 = 0.001, 
                                  clump_p1 = 1,
                                  clump_p2 = 1,
                                  pop = "EUR"
) 

rm(asd)
rm(adhd)
rm(adhd2023)
rm(an)
rm(anx)
rm(bip)
rm(bip1)
rm(bip2)
rm(mdd)
rm(ocd)
rm(ptsd)
rm(scz)
rm(alc)
rm(asd_instrument)
rm(adhd_instrument)
rm(an_instrument)
rm(anx_instrument)
rm(bip_instrument)
rm(bip1_instrument)
rm(bip2_instrument)
rm(mdd_instrument)
rm(ocd_instrument)
rm(ptsd_instrument)
rm(scz_instrument)
rm(asd_exp_dat)
rm(adhd_exp_dat)
rm(an_exp_dat)
rm(anx_exp_dat)
rm(bip_exp_dat)
rm(bip1_exp_dat)
rm(bip2_exp_dat)
rm(mdd_exp_dat)
rm(mdd_exp_dat)
nb_snps_clumped <- c(nrow(asd_exp_dat_clumped), 
                     nrow(adhd_exp_dat_clumped), 
                     nrow(adhd2023_exp_dat_clumped), 
                     nrow(an_exp_dat_clumped),
                     nrow(anx_exp_dat_clumped), 
                     nrow(bip_exp_dat_clumped), 
                     nrow(bip1_exp_dat_clumped), 
                     nrow(bip2_exp_dat_clumped), 
                     nrow(mdd_exp_dat_clumped), 
                     nrow(ocd_exp_dat_clumped), 
                     nrow(ptsd_exp_dat_clumped), 
                     nrow(scz_exp_dat_clumped),
                     nrow(alc_exp_dat_clumped))
#25
# Load outcome

#* 1.2 EA GWAS (including 23andMe, meta-analysis from CogNonCog project) --------
# Meta-analysis was described and checked in document Meta-analysis EA.R

ea_meta <- fread("../Summary_statistics/EA/meta_education_lee_23andMe_130820211.txt") # Perline's MA AUGST 13TH
head(ea_meta)

#** 1.2.2 compute beta & se  ----------------
ea_meta$Beta <- ea_meta$Zscore / sqrt(2*ea_meta$Freq1*(1 - ea_meta$Freq1)*(ea_meta$Weight+ ea_meta$Zscore*ea_meta$Zscore))  
ea_meta$SE <- 1/sqrt(2*ea_meta$Freq1*(1 - ea_meta$Freq1)*(ea_meta$Weight+ ea_meta$Zscore*ea_meta$Zscore))  
#** 1.2.4 scale beta to "years of education", which has a specific standard deviation  (3.9). #################
ea_meta$Beta <- 3.9*ea_meta$Beta
ea_meta$SE <- 3.9*ea_meta$SE


# EA within family 
ea_wf <- fread("../Summary_statistics/EA/Withinsib/Education_wf_ieu_4836_reformat.txt") 
head(ea_wf)
# scale beta to "years of education", which has a specific standard deviation  (3.9). #################
ea_wf$Beta <- 3.9*ea_wf$effect
ea_wf$SE <- 3.9*ea_wf$SE


# Format outcome for MR with each exposure ---------------------
ea_meta_dat_asd <- format_data(ea_meta, 
                       type="outcome",
                       snps= asd_exp_dat_clumped$SNP,
                       snp_col = "MarkerName",
                       beta_col = "Beta",
                       se_col = "SE",
                       effect_allele_col = "Allele1",
                       other_allele_col = "Allele2",
                       pval_col = "P-value",
                       eaf_col = "Freq1")
ea_wf_dat_asd <- format_data(ea_wf, 
                               type="outcome",
                               snps= asd_exp_dat_clumped$SNP,
                               snp_col = "ID",
                               beta_col = "Beta",
                               se_col = "SE",
                               effect_allele_col = "ALT", #correct ordering
                               other_allele_col = "REF",
                               pval_col = "pval")
ea_meta_dat_adhd <- format_data(ea_meta, 
                                type="outcome",
                                snps= adhd_exp_dat_clumped$SNP,
                                snp_col = "MarkerName",
                                beta_col = "Beta",
                                se_col = "SE",
                                effect_allele_col = "Allele1",
                                other_allele_col = "Allele2",
                                pval_col = "P-value",
                                eaf_col = "Freq1")
ea_wf_dat_adhd <- format_data(ea_wf, 
                              type="outcome",
                              snps= adhd_exp_dat_clumped$SNP,
                              snp_col = "ID",
                              beta_col = "Beta",
                              se_col = "SE",
                              effect_allele_col = "ALT",
                              other_allele_col = "REF",
                              pval_col = "pval")
ea_meta_dat_adhd2023 <- format_data(ea_meta, 
                                type="outcome",
                                snps= adhd2023_exp_dat_clumped$SNP,
                                snp_col = "MarkerName",
                                beta_col = "Beta",
                                se_col = "SE",
                                effect_allele_col = "Allele1",
                                other_allele_col = "Allele2",
                                pval_col = "P-value",
                                eaf_col = "Freq1")
ea_wf_dat_adhd2023 <- format_data(ea_wf, 
                              type="outcome",
                              snps= adhd2023_exp_dat_clumped$SNP,
                              snp_col = "ID",
                              beta_col = "Beta",
                              se_col = "SE",
                              effect_allele_col = "ALT",
                              other_allele_col = "REF",
                              pval_col = "pval")
ea_meta_dat_an <- format_data(ea_meta, 
                                type="outcome",
                                snps= an_exp_dat_clumped$SNP,
                                snp_col = "MarkerName",
                                beta_col = "Beta",
                                se_col = "SE",
                                effect_allele_col = "Allele1",
                                other_allele_col = "Allele2",
                                pval_col = "P-value",
                                eaf_col = "Freq1")
ea_wf_dat_an <- format_data(ea_wf, 
                              type="outcome",
                              snps= an_exp_dat_clumped$SNP,
                              snp_col = "ID",
                              beta_col = "Beta",
                              se_col = "SE",
                              effect_allele_col = "ALT",
                              other_allele_col = "REF",
                              pval_col = "pval")
ea_meta_dat_anx <- format_data(ea_meta, 
                              type="outcome",
                              snps= anx_exp_dat_clumped$SNP,
                              snp_col = "MarkerName",
                              beta_col = "Beta",
                              se_col = "SE",
                              effect_allele_col = "Allele1",
                              other_allele_col = "Allele2",
                              pval_col = "P-value",
                              eaf_col = "Freq1")
ea_wf_dat_anx <- format_data(ea_wf, 
                            type="outcome",
                            snps= anx_exp_dat_clumped$SNP,
                            snp_col = "ID",
                            beta_col = "Beta",
                            se_col = "SE",
                            effect_allele_col = "ALT",
                            other_allele_col = "REF",
                            pval_col = "pval")
ea_meta_dat_bip <- format_data(ea_meta, 
                               type="outcome",
                               snps= bip_exp_dat_clumped$SNP,
                               snp_col = "MarkerName",
                               beta_col = "Beta",
                               se_col = "SE",
                               effect_allele_col = "Allele1",
                               other_allele_col = "Allele2",
                               pval_col = "P-value",
                               eaf_col = "Freq1")
ea_wf_dat_bip <- format_data(ea_wf, 
                             type="outcome",
                             snps= bip_exp_dat_clumped$SNP,
                             snp_col = "ID",
                             beta_col = "Beta",
                             se_col = "SE",
                             effect_allele_col = "ALT",
                             other_allele_col = "REF",
                             pval_col = "pval")
ea_meta_dat_bip1 <- format_data(ea_meta, 
                               type="outcome",
                               snps= bip1_exp_dat_clumped$SNP,
                               snp_col = "MarkerName",
                               beta_col = "Beta",
                               se_col = "SE",
                               effect_allele_col = "Allele1",
                               other_allele_col = "Allele2",
                               pval_col = "P-value",
                               eaf_col = "Freq1")
ea_wf_dat_bip1 <- format_data(ea_wf, 
                             type="outcome",
                             snps= bip1_exp_dat_clumped$SNP,
                             snp_col = "ID",
                             beta_col = "Beta",
                             se_col = "SE",
                             effect_allele_col = "ALT",
                             other_allele_col = "REF",
                             pval_col = "pval")
ea_meta_dat_bip2 <- format_data(ea_meta, 
                               type="outcome",
                               snps= bip2_exp_dat_clumped$SNP,
                               snp_col = "MarkerName",
                               beta_col = "Beta",
                               se_col = "SE",
                               effect_allele_col = "Allele1",
                               other_allele_col = "Allele2",
                               pval_col = "P-value",
                               eaf_col = "Freq1")
ea_wf_dat_bip2 <- format_data(ea_wf, 
                             type="outcome",
                             snps= bip2_exp_dat_clumped$SNP,
                             snp_col = "ID",
                             beta_col = "Beta",
                             se_col = "SE",
                             effect_allele_col = "ALT",
                             other_allele_col = "REF",
                             pval_col = "pval")
ea_meta_dat_mdd <- format_data(ea_meta, 
                                type="outcome",
                                snps= mdd_exp_dat_clumped$SNP,
                                snp_col = "MarkerName",
                                beta_col = "Beta",
                                se_col = "SE",
                                effect_allele_col = "Allele1",
                                other_allele_col = "Allele2",
                                pval_col = "P-value",
                                eaf_col = "Freq1")
ea_wf_dat_mdd <- format_data(ea_wf, 
                              type="outcome",
                              snps= mdd_exp_dat_clumped$SNP,
                              snp_col = "ID",
                              beta_col = "Beta",
                              se_col = "SE",
                              effect_allele_col = "ALT",
                              other_allele_col = "REF",
                              pval_col = "pval")
ea_meta_dat_ocd <- format_data(ea_meta, 
                                type="outcome",
                                snps= ocd_exp_dat_clumped$SNP,
                                snp_col = "MarkerName",
                                beta_col = "Beta",
                                se_col = "SE",
                                effect_allele_col = "Allele1",
                                other_allele_col = "Allele2",
                                pval_col = "P-value",
                                eaf_col = "Freq1")
ea_wf_dat_ocd <- format_data(ea_wf, 
                              type="outcome",
                              snps= ocd_exp_dat_clumped$SNP,
                              snp_col = "ID",
                              beta_col = "Beta",
                              se_col = "SE",
                              effect_allele_col = "ALT",
                              other_allele_col = "REF",
                              pval_col = "pval")
ea_meta_dat_ptsd <- format_data(ea_meta, 
                               type="outcome",
                               snps= ptsd_exp_dat_clumped$SNP,
                               snp_col = "MarkerName",
                               beta_col = "Beta",
                               se_col = "SE",
                               effect_allele_col = "Allele1",
                               other_allele_col = "Allele2",
                               pval_col = "P-value",
                               eaf_col = "Freq1")
ea_wf_dat_ptsd <- format_data(ea_wf, 
                             type="outcome",
                             snps= ptsd_exp_dat_clumped$SNP,
                             snp_col = "ID",
                             beta_col = "Beta",
                             se_col = "SE",
                             effect_allele_col = "ALT",
                             other_allele_col = "REF",
                             pval_col = "pval")
ea_meta_dat_scz <- format_data(ea_meta, 
                               type="outcome",
                               snps= scz_exp_dat_clumped$SNP,
                               snp_col = "MarkerName",
                               beta_col = "Beta",
                               se_col = "SE",
                               effect_allele_col = "Allele1",
                               other_allele_col = "Allele2",
                               pval_col = "P-value",
                               eaf_col = "Freq1")
ea_wf_dat_scz <- format_data(ea_wf, 
                             type="outcome",
                             snps= scz_exp_dat_clumped$SNP,
                             snp_col = "ID",
                             beta_col = "Beta",
                             se_col = "SE",
                             effect_allele_col = "ALT",
                             other_allele_col = "REF",
                             pval_col = "pval")
ea_meta_dat_alc <- format_data(ea_meta, 
                               type="outcome",
                               snps= alc_exp_dat_clumped$SNP,
                               snp_col = "MarkerName",
                               beta_col = "Beta",
                               se_col = "SE",
                               effect_allele_col = "Allele1",
                               other_allele_col = "Allele2",
                               pval_col = "P-value",
                               eaf_col = "Freq1")
ea_wf_dat_alc <- format_data(ea_wf, 
                             type="outcome",
                             snps= alc_exp_dat_clumped$SNP,
                             snp_col = "ID",
                             beta_col = "Beta",
                             se_col = "SE",
                             effect_allele_col = "ALT",
                             other_allele_col = "REF",
                             pval_col = "pval")
ea_meta_dat_asd$outcome <- "EA"
ea_wf_dat_asd$outcome <- "EA_wf"
ea_meta_dat_adhd$outcome <- "EA"
ea_wf_dat_adhd$outcome <- "EA_wf"
ea_meta_dat_adhd2023$outcome <- "EA"
ea_wf_dat_adhd2023$outcome <- "EA_wf"
ea_meta_dat_an$outcome <- "EA"
ea_wf_dat_an$outcome <- "EA_wf"
ea_meta_dat_anx$outcome <- "EA"
ea_wf_dat_anx$outcome <- "EA_wf"
ea_meta_dat_bip$outcome <- "EA"
ea_wf_dat_bip$outcome <- "EA_wf"
ea_meta_dat_bip1$outcome <- "EA"
ea_wf_dat_bip1$outcome <- "EA_wf"
ea_meta_dat_bip2$outcome <- "EA"
ea_wf_dat_bip2$outcome <- "EA_wf"
ea_meta_dat_mdd$outcome <- "EA"
ea_wf_dat_mdd$outcome <- "EA_wf"
ea_meta_dat_ocd$outcome <- "EA"
ea_wf_dat_ocd$outcome <- "EA_wf"
ea_meta_dat_ptsd$outcome <- "EA"
ea_wf_dat_ptsd$outcome <- "EA_wf"
ea_meta_dat_scz$outcome <- "EA"
ea_wf_dat_scz$outcome <- "EA_wf"
ea_meta_dat_alc$outcome <- "EA"
ea_wf_dat_alc$outcome <- "EA_wf"
# get number of snps missing from ea 
nb_snp_shared_ea <- c(nrow(ea_meta_dat_asd), 
                      nrow(ea_meta_dat_adhd),
                      nrow(ea_meta_dat_adhd2023),
                      nrow(ea_meta_dat_an), 
                      nrow(ea_meta_dat_anx), 
                      nrow(ea_meta_dat_bip), 
                      nrow(ea_meta_dat_bip1), 
                      nrow(ea_meta_dat_bip2), 
                      nrow(ea_meta_dat_mdd), 
                      nrow(ea_meta_dat_ocd), 
                      nrow(ea_meta_dat_ptsd), 
                      nrow(ea_meta_dat_scz), 
                      nrow(ea_meta_dat_alc))
nb_snp_shared_ea_wf <- c(nrow(ea_wf_dat_asd), 
                      nrow(ea_wf_dat_adhd),
                      nrow(ea_wf_dat_adhd2023),
                      nrow(ea_wf_dat_an), 
                      nrow(ea_wf_dat_anx), 
                      nrow(ea_wf_dat_bip), 
                      nrow(ea_wf_dat_bip1), 
                      nrow(ea_wf_dat_bip2), 
                      nrow(ea_wf_dat_mdd), 
                      nrow(ea_wf_dat_ocd), 
                      nrow(ea_wf_dat_ptsd), 
                      nrow(ea_wf_dat_scz), 
                      nrow(ea_wf_dat_alc))



rm(ea_meta)
rm(ea_wf)

# Harmonise SNPS ######
dat.asd <- harmonise_data(
  exposure_dat = asd_exp_dat_clumped, 
  outcome_dat = ea_meta_dat_asd, action = 2
)

dat.asd.wf <- harmonise_data(
  exposure_dat = asd_exp_dat_clumped, 
  outcome_dat = ea_wf_dat_asd, action = 2
)

dat.adhd <- harmonise_data(
  exposure_dat = adhd_exp_dat_clumped, 
  outcome_dat = ea_meta_dat_adhd, action = 2
)
dat.adhd.wf <- harmonise_data(
  exposure_dat = adhd_exp_dat_clumped, 
  outcome_dat = ea_wf_dat_adhd, action = 2
)
dat.adhd2023 <- harmonise_data(
  exposure_dat = adhd2023_exp_dat_clumped, 
  outcome_dat = ea_meta_dat_adhd2023, action = 2
)
dat.adhd2023.wf <- harmonise_data(
  exposure_dat = adhd2023_exp_dat_clumped, 
  outcome_dat = ea_wf_dat_adhd2023, action = 2
)
dat.an <- harmonise_data(
  exposure_dat = an_exp_dat_clumped, 
  outcome_dat = ea_meta_dat_an, action = 2
)
dat.an.wf <- harmonise_data(
  exposure_dat = an_exp_dat_clumped, 
  outcome_dat = ea_wf_dat_an, action = 2
)
dat.anx <- harmonise_data(
  exposure_dat = anx_exp_dat_clumped, 
  outcome_dat = ea_meta_dat_anx, action = 2
)
dat.anx.wf <- harmonise_data(
  exposure_dat = anx_exp_dat_clumped, 
  outcome_dat = ea_wf_dat_anx, action = 2
)
dat.bip <- harmonise_data(
  exposure_dat = bip_exp_dat_clumped, 
  outcome_dat = ea_meta_dat_bip, action = 2
)
dat.bip.wf <- harmonise_data(
  exposure_dat = bip_exp_dat_clumped, 
  outcome_dat = ea_wf_dat_bip, action = 2
)
dat.bip1 <- harmonise_data(
  exposure_dat = bip1_exp_dat_clumped, 
  outcome_dat = ea_meta_dat_bip1, action = 2
)
dat.bip1.wf <- harmonise_data(
  exposure_dat = bip1_exp_dat_clumped, 
  outcome_dat = ea_wf_dat_bip1, action = 2
)
dat.bip2 <- harmonise_data(
  exposure_dat = bip2_exp_dat_clumped, 
  outcome_dat = ea_meta_dat_bip2, action = 2
)
dat.bip2.wf <- harmonise_data(
  exposure_dat = bip2_exp_dat_clumped, 
  outcome_dat = ea_wf_dat_bip2, action = 2
)
dat.mdd <- harmonise_data(
  exposure_dat = mdd_exp_dat_clumped, 
  outcome_dat = ea_meta_dat_mdd, action = 2
)
dat.mdd.wf <- harmonise_data(
  exposure_dat = mdd_exp_dat_clumped, 
  outcome_dat = ea_wf_dat_mdd, action = 2
)
dat.ocd <- harmonise_data(
  exposure_dat = ocd_exp_dat_clumped, 
  outcome_dat = ea_meta_dat_ocd, action = 2
)
dat.ocd.wf <- harmonise_data(
  exposure_dat = ocd_exp_dat_clumped, 
  outcome_dat = ea_wf_dat_ocd, action = 2
)
dat.ptsd <- harmonise_data(
  exposure_dat = ptsd_exp_dat_clumped, 
  outcome_dat = ea_meta_dat_ptsd, action = 2
)
dat.ptsd.wf <- harmonise_data(
  exposure_dat = ptsd_exp_dat_clumped, 
  outcome_dat = ea_wf_dat_ptsd, action = 2
)
dat.scz <- harmonise_data(
  exposure_dat = scz_exp_dat_clumped, 
  outcome_dat = ea_meta_dat_scz, action = 2
)
dat.scz.wf <- harmonise_data(
  exposure_dat = scz_exp_dat_clumped, 
  outcome_dat = ea_wf_dat_scz, action = 2
)
dat.alc <- harmonise_data(
  exposure_dat = alc_exp_dat_clumped, 
  outcome_dat = ea_meta_dat_alc, action = 2
)
dat.alc.wf <- harmonise_data(
  exposure_dat = alc_exp_dat_clumped, 
  outcome_dat = ea_wf_dat_alc, action = 2
)


# Analyses 

asd.results <- run_MR_analyses(dat.asd)
asd.results.wf <- run_MR_analyses(dat.asd.wf)
save(asd.results, asd.results.wf, file = "asd.MRMHEA.results.RData")

adhd.results <- run_MR_analyses(dat.adhd)
adhd.results.wf <- run_MR_analyses(dat.adhd.wf)
save(adhd.results, adhd.results.wf, file = "adhd.MRMHEA.results.RData")

adhd2023.results <- run_MR_analyses(dat.adhd2023)
adhd2023.results.wf <- run_MR_analyses(dat.adhd2023.wf)
save(adhd2023.results, adhd2023.results.wf, file = "adhd2023.MRMHEA.results.RData")

an.results <- run_MR_analyses(dat.an)
an.results.wf <- run_MR_analyses(dat.an.wf)
save(an.results, an.results.wf, file = "an.MRMHEA.results.RData")

anx.results <- run_MR_analyses(dat.anx)
anx.results.wf <- run_MR_analyses(dat.anx.wf)
save(anx.results, anx.results.wf, file = "anx.MRMHEA.results.RData")

bip.results <- run_MR_analyses(dat.bip)
bip.results.wf <- run_MR_analyses(dat.bip.wf)
save(bip.results, bip.results.wf, file = "bip.MRMHEA.results.RData")

bip1.results <- run_MR_analyses(dat.bip1)
bip1.results.wf <- run_MR_analyses(dat.bip1.wf)
save(bip1.results, bip1.results.wf, file = "bip1.MRMHEA.results.RData")

bip2.results <- run_MR_analyses(dat.bip2)
bip2.results.wf <- run_MR_analyses(dat.bip2.wf)
save(bip2.results, bip2.results.wf, file = "bip2.MRMHEA.results.RData")

mdd.results <- run_MR_analyses(dat.mdd)
mdd.results.wf <- run_MR_analyses(dat.mdd.wf)
save(mdd.results, mdd.results.wf, file = "mdd.MRMHEA.results.RData")

ocd.results <- run_MR_analyses(dat.ocd)
ocd.results.wf <- run_MR_analyses(dat.ocd.wf)
save(ocd.results, ocd.results.wf, file = "ocd.MRMHEA.results.RData")

ptsd.results <- run_MR_analyses(dat.ptsd)
ptsd.results.wf <- run_MR_analyses(dat.ptsd.wf)
save(ptsd.results, ptsd.results.wf, file = "ptsd.MRMHEA.results.RData")

scz.results <- run_MR_analyses(dat.scz)
scz.results.wf <- run_MR_analyses(dat.scz.wf)
save(scz.results, scz.results.wf, file = "scz.MRMHEA.results.RData")

alc.results <- run_MR_analyses(dat.alc)
alc.results.wf <- run_MR_analyses(dat.alc.wf)
save(alc.results, alc.results.wf, file = "alc.MRMHEA.results.RData")


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
all.wf.results <- rbind(asd.results.wf[[1]], 
                        adhd.results.wf[[1]], 
                        adhd2023.results.wf[[1]], 
                        an.results.wf[[1]], 
                        anx.results.wf[[1]], 
                        bip.results.wf[[1]], 
                        bip1.results.wf[[1]],
                        bip2.results.wf[[1]], 
                        mdd.results.wf[[1]], 
                        ocd.results.wf[[1]],
                        ptsd.results.wf[[1]], 
                        scz.results.wf[[1]],
                        alc.results.wf[[1]]
                        
)
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
all.wf.sensitivity <- rbind(asd.results.wf[[2]], 
                            adhd.results.wf[[2]], 
                            adhd2023.results.wf[[2]], 
                            an.results.wf[[2]], 
                            anx.results.wf[[2]], 
                            bip.results.wf[[2]], 
                            bip1.results.wf[[2]],
                            bip2.results.wf[[2]], 
                            mdd.results.wf[[2]], 
                            ocd.results.wf[[2]],
                            ptsd.results.wf[[2]], 
                            scz.results.wf[[2]],
                            alc.results.wf[[2]]
)

all.wf.Q <- rbind(asd.results.wf[[3]], 
                  adhd.results.wf[[3]], 
                  adhd2023.results.wf[[3]],
                  an.results.wf[[3]], 
                  anx.results.wf[[3]], 
                  bip.results.wf[[3]], 
                  bip1.results.wf[[3]],
                  bip2.results.wf[[3]], 
                  mdd.results.wf[[3]], 
                  ocd.results.wf[[3]],
                  ptsd.results.wf[[3]], 
                  scz.results.wf[[3]],
                  alc.results.wf[[3]]
                  
)

write.csv(all.results, "all_MR_results_rev_221020.csv")
write.csv(all.wf.results, "all_MR_wf_results_rev_221020.csv")
write.csv(all.sensitivity, "all_MR_sensitivity_rev_221022.csv")
write.csv(all.Q, "all_MR_Q_rev_221020.csv")
write.csv(all.wf.sensitivity, "all_MR_wf_sensitivity_rev_221020.csv")
write.csv(all.wf.Q, "all_MR_Q_wf_rev_221020.csv")

all.results.wide <- reshape(all.results, idvar = "exposure" , timevar = "method", direction = "wide")
write.csv(all.results.wide, "all_MR_results_wide_rev_221020.csv")

all.wf.results.wide <- reshape(all.wf.results, idvar = "exposure" , timevar = "method", direction = "wide")
write.csv(all.wf.results.wide, "all_MR_wf_results_wide_rev_221020.csv")


all.Q.wide <- reshape(all.Q, idvar = "exposure" , timevar = "method", direction = "wide")
write.csv(all.Q.wide, "all_MR_Q_wide_rev_221020.csv")


all.wf.Q.wide <- reshape(all.wf.Q , idvar = "exposure" , timevar = "method", direction = "wide")
write.csv(all.wf.Q.wide, "all_MR_wf_Q_wide_rev_221020.csv")
