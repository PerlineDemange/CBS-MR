## Project: EA_MH_CBS_MR 
## Script purpose: Run meta-analysis EA Lee et al + 23andMe 
## Author: Perline Demange 

# 1. Meta-analysis done for previous paper Demange et al. 2021 ####

ea_meta <- fread("Summary_statistics/EA/meta_education_lee_23andMe_adjust_rsid.txt", sep =" ", header=F)
head(ea_meta)
# This file was saved with both tabular and space as separator 
# Need some reformating, done below 

split2 <- str_split_fixed(ea_meta$V1, "\t", 2) #separate first column which was badly saved with tab as sep
colnames(split2) <- c("chr:pos", "rsid")
head(split2)
split2 <- print.data.frame(as.data.frame(split2), quote=FALSE) # remove the quotes 

ea_meta_final <- cbind(split2, ea_meta[, 2:ncol(ea_meta)])
names(ea_meta_final) <- c("CHR:POS", "SNP", "CHR", "POS", "Allele1", "Allele2", "Freq1", "FreqSE", "MinFreq", "MaxFreq", "Weight", "Zscore", "P", "Direction")

write.csv(ea_meta_final, "meta_education_lee_23andMe_adjust_rsid_reformat.csv", row.names=F)

# this file only contains Zscores and P, no beta and SE 
# I need to get beta and SE 


# 2. Get Beta and SE #######
# I use the equation answered here: https://www.biostars.org/p/319584/
# from this paper https://www.nature.com/articles/ng.3538
data <- fread("Summary_statistics/meta_education_lee_23andMe_adjust_rsid_reformat.csv", header=T)
head(data)

data$Beta <- data$Zscore / sqrt(2*data$Freq1*(1 - data$Freq1)*
                                  (data$Weight + data$Zscore*data$Zscore))
data$SE <- 1/sqrt(2*data$Freq1*(1 - data$Freq1)*
                    (data$Weight + data$Zscore*data$Zscore))

# Check method to get beta and SE used in Genomic SEM
data$beta_gsem <- (sign(data$Zscore)*sqrt(qchisq(data$P,1, lower = F)))/
  sqrt(data$Weight*2*data$Freq1*(1-data$Freq1))
data$SE_gsem <- abs(data$beta_gsem/(sign(data$Zscore)*
                                      sqrt(qchisq(data$P,1, lower = F))))

sign <- data[data$P < 5e-8, ]
plot(sign$Zscore, sign$Beta)
plot(sign$Zscore, sign$beta_gsem)
plot(sign$Beta, sign$beta_gsem)

# Both methods are consistent
# I use the first one 

datatosave <- data[, 1:16]
write.table(datatosave, "meta_education_lee_23andMe_betas_20210813.txt", row.names=F)


# 3. Re-do meta-analysis to double check ###########
# in LISA 
# module load 2019
# module load Metal/2011-03-25-foss-2019b
# metal < metal #file from previous project, see other Github /non-cognitive

newdata <- fread("Summary_statistics/EA/meta_education_lee_23andMe_130820211.txt", header=T)

newsign <- newdata[newdata$`P-value` < 5e-8, ]
new <- newsign[order(newsign$MarkerName),] 
plot(new$Zscore, sign$Zscore)
head(new)
head(sign)
# the old meta-analysis is identical to the new one, so it ws correct and I will use the old doc as it contained extra colums. 
