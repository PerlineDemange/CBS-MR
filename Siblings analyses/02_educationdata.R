# Script to add education variables to our population data 
# And create subset with same sex siblings 
# Author: Perline Demange
# Project: CBS-MR EA-MH 
# 


# 1. Load Educational Attainment data ####
edu <- fread("../Data/HOOGSTEOPL2019TABV1.csv")
head(edu) #11925865  300000 more than in 2018
summary(edu)

# Check data quality 
table(edu$OPLNIVSOI2016AGG4HBMETNIRWO) #no missing data (coded as 9999) 

# 2. Recode education data ####
# *2.1 Investigation of Praktijkonderwijs code 1211 ----
# Children with highest diploma being Praktijkonderwijs are mostly children with disabilities. 
# It is difficult to convert this type of higher diploma in a number of years studied 
# as children might stay there until they are 18

table(edu$OPLNIVSOI2016AGG4HBMETNIRWO)
# 6898 individuals have 1211 as highest diploma 
6869*100/11925865 # This represents 0.05 % of the population 

# I decide to code this diploma as missing value, and will exclude them from further analyses 

# *2.2 Code in number of years ------
edu$yrs <- edu$OPLNIVSOI2016AGG4HBMETNIRWO
edu$yrs[edu$yrs == 1111] <- 2
edu$yrs[edu$yrs == 1112] <- 8
edu$yrs[edu$yrs == 1211] <- NA
edu$yrs[edu$yrs == 1212] <- 12
edu$yrs[edu$yrs == 1213] <- 12
edu$yrs[edu$yrs == 1221] <- 12
edu$yrs[edu$yrs == 1222] <- 11
edu$yrs[edu$yrs == 2111] <- 13
edu$yrs[edu$yrs == 2112] <- 14
edu$yrs[edu$yrs == 2121] <- 15
edu$yrs[edu$yrs == 2131] <- 13
edu$yrs[edu$yrs == 2132] <- 14
edu$yrs[edu$yrs == 3111] <- 15
edu$yrs[edu$yrs == 3112] <- 17
edu$yrs[edu$yrs == 3113] <- 17
edu$yrs[edu$yrs == 3211] <- 18
edu$yrs[edu$yrs == 3212] <- 18
edu$yrs[edu$yrs == 3213] <- 22

barplot(table(edu$OPLNIVSOI2016AGG4HBMETNIRWO))
barplot(table(edu$yrs))


# 3. Population data with EA ######
# *3.1 Merge population data with educational data -----
gba <- fread("gba_6585_20201015.csv") #6539767
head(gba)
edu_pop_1211 <- merge(gba, edu, by=c('RINPERSOONS', "RINPERSOON")) #3306956

# *3.2 Missing rate  -----
100 - nrow(edu_pop_1211)*100/nrow(gba) #49.43%

# * 3.3 Exclude EA 1211 
edu_pop <- edu_pop_1211[!which(edu_pop_1211$OPLNIVSOI2016AGG4HBMETNIRWO == 1211),] #3305733

# 4. Sibling data with EA ##### 
# *4.1 Merge sibling data and education data ------
sib <- fread("gba_6585_FID_fullsibship_final_20210215.csv")
head(sib) #3234923

edu_sib_full <- merge(sib, edu, by=c('RINPERSOONS', "RINPERSOON")) #2182763

# *4.2 Missing rate full sibling population  ----
100 - nrow(edu_sib_full)*100/nrow(sib) #32.53%

# *4.3 Remove incomplete families ----
# ** 4.3.1 Data with 1211 -----
# get new family size 
count <- as.data.frame(table(edu_sib_full$FID)) #get number of occurences of FID
edu_sib_full$Fsize_edu <- count[[2]][match(edu_sib_full$FID, count[[1]])]  # match data fid with value in count, if true takes value of Freq 
table(edu_sib_full$Fsize_edu) 
edu_sib_1211 <- edu_sib_full[edu_sib_full$Fsize_edu > 1] #1743826

# ** 4.3.2 Data without 1211 ------
edu_sib <- edu_sib_full[!which(edu_sib_full$OPLNIVSOI2016AGG4HBMETNIRWO == 1211),] #2182137 #626 indiv with 1211 in edu_sib_full
count <- as.data.frame(table(edu_sib$FID)) #get number of occurences of FID
edu_sib$Fsize_edu <- count[[2]][match(edu_sib$FID, count[[1]])]  # match data fid with value in count, if true takes value of Freq 
table(edu_sib$Fsize_edu) 
edu_sib <- edu_sib[edu_sib$Fsize_edu > 1] #1743032

# *4.4 number of unique sibships ----
numberfamilies <- unique(edu_sib$FID) 
length(numberfamilies) #766 514 families 
nrow(edu_sib)/ length(unique(edu_sib$FID)) #Average of 2.27 sibling per sibships

# number of multiples 
summary(edu_sib$multiple) #45134


# 5. Descriptive of the data #######
# *5.1 PratijktOnderwijs 1211 ------
table(edu_pop_1211$OPLNIVSOI2016AGG4HBMETNIRWO) #1223
1223* 100/ nrow(edu_pop) #0.037 %

table(edu_sib_1211$OPLNIVSOI2016AGG4HBMETNIRWO) #476
476 *100/ nrow(edu_sib_1211) #0.027 % 

# * 5.2 Compare sibling subset with population ----
# **5.2.1 Descriptives and distributions ######
barplot(table(edu_pop$OPLNIVSOI2016AGG4HBMETNIRWO))
barplot(table(edu_sib$OPLNIVSOI2016AGG4HBMETNIRWO))
count_table_diploma <- as.data.frame(rbind(table(edu_pop$OPLNIVSOI2016AGG4HBMETNIRWO), 
                     table(edu_sib$OPLNIVSOI2016AGG4HBMETNIRWO)))

count_table_diploma$population <- c("gba_EA", "sib_EA")

barplot(table(edu_pop$yrs))
barplot(table(edu_sib$yrs))
count_table_yrs <- as.data.frame(rbind(table(edu_pop$yrs), 
                                           table(edu_sib$yrs)))

count_table_yrs$population <- c("gba_EA", "sib_EA")

descriptive_edu_pop <- as.data.frame(t(unclass(summary(edu_pop$yrs))))
descriptive_edu_pop$SD <- sd(edu_pop$yrs)
descriptive_edu_pop$sample_size <- nrow(edu_pop)
descriptive_edu_pop$population <- "gba_EA"
descriptive_edu_pop

descriptive_edu_sib <- as.data.frame(t(unclass(summary(edu_sib$yrs))))
descriptive_edu_sib$SD <- sd(edu_sib$yrs)
descriptive_edu_sib$sample_size <- nrow(edu_sib)
descriptive_edu_sib$population <- "sib_EA"
descriptive_edu_sib

descriptive_edu <- rbind(descriptive_edu_pop, descriptive_edu_sib)


# **5.2.2 F-test, t-test, and cohen'd #####
ftest <- var.test(edu_pop$yrs, edu_sib$yrs)
ttest <- t.test(edu_pop$yrs, edu_sib$yrs, var.equal=T)
#  min p-value < 2.2e-16

cohend <- (mean(edu_sib$yrs) - mean(edu_pop$yrs))/ sqrt((var(edu_sib$yrs)+var(edu_pop$yrs)/2))

descriptive_edu$ttest_p <-  ttest$p.value
descriptive_edu$ftest_p <-  ftest$p.value
descriptive_edu$cohend <- cohend

# 6. Save population data and sibling subset merged with educational data, 1211 excluded #####
#write.csv(edu_sib, "sib_edu_subset_20210219.csv", row.names = F)
edu_sib <- fread("sib_edu_subset_20210219.csv")
#write.csv(edu_pop, "gba6585_edu_subset_20210219.csv", row.names = F)
edu_pop <- fread("gba6585_edu_subset_20210219.csv")

#write.csv(descriptive_edu, "descriptive_table_edu_allpop_20220729.csv", row.names = F)
#write.csv(count_table_yrs, "descriptive_table_count_yrs_20210223.csv", row.names=F)
#write.csv(count_table_diploma, "descriptive_table_count_diploma_20210223.csv", row.names=F)
descriptive_edu <- fread("descriptive_table_edu_allpop_20220729.csv")
descriptive_edu

descriptive_edu_output <- descriptive_edu[,!c(1,6)]
write.xlsx(descriptive_edu, "descriptive_table_edu_allpop_20220729.xlsx", row.names = F)

# 7. Create subset with same-sex siblings ############
edu_sib <- fread("sib_edu_subset_20210219.csv")
head(edu_sib)
edu_sib_men <- edu_sib[edu_sib$GBAGESLACHT ==1, ] #875381
count <- as.data.frame(table(edu_sib_men$FID)) #get number of occurences of FID
edu_sib_men$Fsize_sex <- count[[2]][match(edu_sib_men$FID, count[[1]])]  # match data fid with value in count, if true takes value of Freq 
edu_sib_men <- edu_sib_men[edu_sib_men$Fsize_sex > 1] #516307

edu_sib_women <- edu_sib[edu_sib$GBAGESLACHT ==2, ] #867651
count <- as.data.frame(table(edu_sib_women$FID)) #get number of occurences of FID
edu_sib_women$Fsize_sex <- count[[2]][match(edu_sib_women$FID, count[[1]])]  # match data fid with value in count, if true takes value of Freq 
edu_sib_women <- edu_sib_women[edu_sib_women$Fsize_sex > 1] #508062


length(unique(edu_sib_men$FID) ) #240 416 families 
nrow(edu_sib_men)/ length(unique(edu_sib_men$FID)) #Average of 2.15 sibling per sibships
length(unique(edu_sib_women$FID) ) #236 138
nrow(edu_sib_women)/ length(unique(edu_sib_women$FID)) #Average of 2.15 sibling per sibships

save(edu_sib_men, file="edu_sib_men_220310.Rda")
save(edu_sib_women, file="edu_sib_women_220310.Rda")
