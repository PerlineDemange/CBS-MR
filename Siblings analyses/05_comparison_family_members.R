# Script to compare siblings/parents of patients with specific disorders.
# Author: Perline Demange
# Project: CBS-MR EA-MH 
#

# Compare education of siblings of patients ####
sib_EA_alldata <- fread("sib_EA_diagnoses_20210221.csv")
head(sib_EA_alldata)

# 1. Identify no diagnoses families ####
Any <- sib_EA_alldata[sib_EA_alldata$Any == 1,] 
family_of_any_FID <- Any$FID
family_of_any <- sib_EA_alldata[sib_EA_alldata$FID %in% family_of_any_FID,] 
no_diagnoses <- sib_EA_alldata[which(!sib_EA_alldata$FID %in% family_of_any_FID),]

# 2. Create a loop to run analyses over all diagnoses ##########
## 2.1 Set up the loop ##############
codefile <- as.data.frame(fread("H:/Data/Diagnoses/Diagnoses_code_preregistered.csv", header=T))
trait_list <- codefile$trait

sample_size_total <- c("No_disorder", nrow(no_diagnoses), NA, NA,NA)
descriptives_total <- NULL
summary_ttest <- NULL

## 2.2 Loop the analysis over all disorders ###################
for (i in 2:(length(trait_list)+1)){
  disorder_name <- colnames(sib_EA_alldata)[i]
  
  ### 2.2.1 Create datasets ############
  # get family with at least one member who has the specific disorder 
  patient <- sib_EA_alldata[sib_EA_alldata[[i]] == 1,] # 5920 people with bipolar 
  number_of_patient <- nrow(patient)
  family_of_patient_FID <- patient$FID
  family_of_patient <- sib_EA_alldata[sib_EA_alldata$FID %in% family_of_patient_FID,] # 14059 members
  family_of_patient_safe <- family_of_patient[family_of_patient$Any == 0,] #sibling of bipolar without diagnoses 6863
  number_sib_patient <- nrow(family_of_patient_safe)
  
  # Identify families where there is not specific disorder but some other diagnoses 
  no_patient <- sib_EA_alldata[which(!sib_EA_alldata$FID %in% family_of_patient_FID),] # 1728973
  no_patient_but_diagnoses <- no_patient[no_patient$Any == 1,] # disorder but not bipolar 177078
  number_of_patient_but_diagnoses <- nrow(no_patient_but_diagnoses)
  family_of_no_patient_but_diagnoses_FID <- no_patient_but_diagnoses$FID
  family_of_no_patient_but_diagnoses <- no_patient[no_patient$FID %in% family_of_no_patient_but_diagnoses_FID,] # 377089
  family_of_no_patient_but_diagnoses_safe <- family_of_no_patient_but_diagnoses[family_of_no_patient_but_diagnoses$Any == 0, ] #200011
  number_sib_no_patient_but_diagnoses <- nrow(family_of_no_patient_but_diagnoses_safe)
  
  ### 2.2.2. F and T tests #######
  #### 2.2.2.1 Siblings of patients compared with patients #####
  # F test to compare two variances 
  Ftest_patient <- var.test(family_of_patient_safe$yrs, patient$yrs)
  
  # Ttest to compare the means 
  if (Ftest_patient$p.value <0.05){ 
    Ttest_patient <- t.test(family_of_patient_safe$yrs, patient$yrs, var.equal = F)
  } else {
    Ttest_patient <- t.test(family_of_patient_safe$yrs, patient$yrs, var.equal = T)
  }
  
  #### 2.2.2.2 with families with another diagnoses #######
  # F test to compare two variances 
  Ftest_diagnoses <- var.test(family_of_patient_safe$yrs, family_of_no_patient_but_diagnoses_safe$yrs)

  # Ttest to compare the means 
  if (Ftest_diagnoses$p.value <0.05){ 
    Ttest_diagnoses <- t.test(family_of_patient_safe$yrs, family_of_no_patient_but_diagnoses_safe$yrs, var.equal = F)
  } else {
    Ttest_diagnoses <- t.test(family_of_patient_safe$yrs, family_of_no_patient_but_diagnoses_safe$yrs, var.equal = T)
  }
  
  #### 2.2.2.3 with families with no diagnoses ######
  # F test to compare two variances 
  Ftest <- var.test(family_of_patient_safe$yrs, no_diagnoses$yrs)
  
  # Ttest to compare the means 
  if (Ftest$p.value <0.05){ 
    Ttest <- t.test(family_of_patient_safe$yrs, no_diagnoses$yrs, var.equal = F)
  } else {
    Ttest <- t.test(family_of_patient_safe$yrs, no_diagnoses$yrs, var.equal = T)
  }
  
  ### 2.2.3. Descriptive measures #######
  sample_size <- c(number_of_patient, number_sib_patient, 
                   number_of_patient_but_diagnoses, number_sib_no_patient_but_diagnoses)
  # histogrames
  # hist_patient <- ggplot(patient, aes(x=yrs)) +
  #   geom_histogram(binwidth = 1) 
  # hist_sib_patient <- ggplot(family_of_patient_safe, aes(x=yrs)) +
  #   geom_histogram(binwidth = 1) 
  # hist_sib_no_patient_but_diagnoses <- ggplot(family_of_no_patient_but_diagnoses_safe, aes(x=yrs)) +
  #   geom_histogram(binwidth = 1)
  
  patient_desc <- c(mean(patient$yrs), sd(patient$yrs))
  sib_patient_desc <- c(mean(family_of_patient_safe$yrs), sd(family_of_patient_safe$yrs))
  sib_no_patient_but_diagnoses_desc <- c(mean(family_of_no_patient_but_diagnoses_safe$yrs), sd(family_of_no_patient_but_diagnoses_safe$yrs))
  
  
  ### 2.2.4 Save the data ######
  sample_size_total <- rbind(sample_size_total, c(disorder_name, sample_size))
  descriptives_total <- rbind(descriptives_total, 
                              c(disorder_name, 
                                c(patient_desc,
                                sib_patient_desc, 
                                sib_no_patient_but_diagnoses_desc, 
                                mean(no_diagnoses$yrs), 
                                sd(no_diagnoses$yrs))))
  summary_ttest <- rbind(summary_ttest,c(disorder_name, 
                                         Ftest_patient$p.value, Ttest_patient$p.value, 
                                         Ftest_diagnoses$p.value, Ttest_diagnoses$p.value, 
                                         Ftest$p.value, Ttest$p.value
                                         ))
  #full_results <- c(hist_sib_patient, hist_sib_no_patient_but_diagnoses, Ftest_diagnoses, Ftest, Ttest_diagnoses, Ttest)
} 


## 2.3 Reformat the data and save ########
# this is not the cleanest/fastest way.. 
descriptives_total <- rbind(c("No_disorder", NA, NA, NA, NA, NA, NA,NA, NA), descriptives_total)
descriptives_total <- as.data.frame(descriptives_total)
descriptives_total[,1]  <-  as.character(descriptives_total[,1])
descriptives_total[,2]  <-  as.numeric(as.character(descriptives_total[,2]))
descriptives_total[,3]  <-  as.numeric(as.character(descriptives_total[,3]))
descriptives_total[,4]  <-  as.numeric(as.character(descriptives_total[,4]))
descriptives_total[,5]  <-  as.numeric(as.character(descriptives_total[,5]))
descriptives_total[,6]  <-  as.numeric(as.character(descriptives_total[,6]))
descriptives_total[,7]  <-  as.numeric(as.character(descriptives_total[,7]))
descriptives_total[,8]  <-  as.numeric(as.character(descriptives_total[,8]))
descriptives_total[,9]  <-  as.numeric(as.character(descriptives_total[,9]))

colnames(descriptives_total) <- c("Diagnoses", 
                                  "Mean_patients", "SD_patients",
                                  "Mean_sib_patients", "SD_sib_patients",
                                  "Mean_sib_other_diagnoses", "SD_sib_other_diagnoses",
                                  "Mean_no_diagnoses", "SD_no_diagnoses")

descriptives_total$diff_nodiagnoses <- as.numeric(descriptives_total$Mean_no_diagnoses) - 
  as.numeric(descriptives_total$Mean_sib_patients)
descriptives_total$diff_otherdiagnoses <- as.numeric(descriptives_total$Mean_sib_other_diagnoses) - 
  as.numeric(descriptives_total$Mean_sib_patients)
descriptives_total$diff_patient <- as.numeric(descriptives_total$Mean_patient) - 
  as.numeric(descriptives_total$Mean_sib_patients)

summary_ttest <- rbind(c("No_disorder", NA, NA, NA, NA, NA, NA), summary_ttest)
summary_ttest <- as.data.frame(summary_ttest)
summary_ttest[,1] <- as.character(summary_ttest[,1])
summary_ttest[,2] <- as.numeric(as.character(summary_ttest[,2]))
summary_ttest[,3] <- as.numeric(as.character(summary_ttest[,3]))
summary_ttest[,4] <- as.numeric(as.character(summary_ttest[,4]))
summary_ttest[,5] <- as.numeric(as.character(summary_ttest[,5]))
summary_ttest[,6] <- as.numeric(as.character(summary_ttest[,6]))
summary_ttest[,7] <- as.numeric(as.character(summary_ttest[,7]))

colnames(summary_ttest) <- c("Diagnoses", 
                             "Ftest_patient", "Ttest_patient",
                             "Ftest_other_diagnoses", "Ttest_other_diagnoses",
                             " Ftest_no_diagnoses", "Ttest_no_diagnoses")
summary_ttest

sample_size_total <- as.data.frame(sample_size_total)
sample_size_total[,1] <- as.character(sample_size_total[,1])
sample_size_total[,2] <- as.numeric(as.character(sample_size_total[,2]))
sample_size_total[,3] <- as.numeric(as.character(sample_size_total[,3]))
sample_size_total[,4] <- as.numeric(as.character(sample_size_total[,4]))
sample_size_total[,5] <- as.numeric(as.character(sample_size_total[,5]))
rownames(sample_size_total) <- NULL
sample_size_total
colnames(sample_size_total) <- c("Diagnoses", 
                             "size_patients", "size_sib_patients","size_sib_other_diagnoses",
                                  "size_no_diagnoses")

full_results <- cbind(descriptives_total, summary_ttest[,-1], sample_size_total[,-1]) 

# save
write.csv(full_results, "comparison_family_members_results.csv", row.names=F)


## 2.4  Plot mean of each groups ########

full_long <- full_results %>% 
  pivot_longer(cols= starts_with("Mean_"), 
               names_to = "groups", 
               names_prefix = "Mean_", 
               values_to = "mean_EA")

ggplot(full_long, aes(x=Diagnoses, y=mean_EA, fill=groups))+
  geom_bar(stat="identity", position = "dodge")

full_long <- full_results %>% 
  pivot_longer(cols= starts_with("diff_"), 
               names_to = "groups", 
               names_prefix = "diff_", 
               values_to = "diff_EA")

ggplot(full_long, aes(x=Diagnoses, y=diff_EA, fill=groups))+
  geom_bar(stat="identity", position = "dodge")



# Compare parents. NOT DONE: sample size too low ##################
# SPOILER: looked at BIP, sample size of reported parental education is too low (N=0-3), did not continue the analysis 
# # 2.1 all siblings with EA data and diagnoses data ###
# sib_EA_alldata <- fread("sib_EA_diagnoses_20210221.csv")
# 
# #everybody with EA data
# edu_pop <- fread("gba6585_edu_subset_20210219.csv")
# head(edu_pop)
# edu_pop_easy <- select(edu_pop, c("RINPERSOON", "yrs"))
# 
# #identify sib with BIP 
# Bipolar <- sib_EA_alldata[sib_EA_alldata$Bipolar == 1,] # 5920 people with bipolar 
# sum(duplicated(Bipolar$FID)) #135, so need to select unique 
# 
# Bipolar <- as.data.frame(Bipolar)
# unique_fam_bipolar <- Bipolar[!duplicated(Bipolar$FID),] #5785
# 
# parents_of_bip <- select(unique_fam_bipolar, c("RINPERSOONpa", "RINPERSOONMa"))
# 
# 
# yrs_of_father <- merge(parents_of_bip, edu_pop, by.y = "RINPERSOON",  by.x = "RINPERSOONpa") # 0
# yrs_of_mother <- merge(parents_of_bip, edu_pop, by.y = "RINPERSOON",  by.x = "RINPERSOONMa") #2
# 
# 
# # 2.2 all siblings with diagnoses data 
# diagnoses_totalyears_per_person_sib <- fread("siblings_diagnoses_totalyear_per_person_20210219.csv")
# 
# head(diagnoses_totalyears_per_person_sib)
# sib <- fread("gba_6585_FID_fullsibship_final_20210215.csv")
# head(sib)
# sib_diagnoses <- merge(sib, diagnoses_totalyears_per_person_sib, by = "RINPERSOON")
# 
# #everybody with EA data
# edu_pop <- fread("gba6585_edu_subset_20210219.csv")
# head(edu_pop)
# edu_pop_easy <- select(edu_pop, c("RINPERSOON", "yrs"))
# 
# #identify sib with BIP 
# Bipolar <- sib_diagnoses[sib_diagnoses$Bipolar == 1,] # 9342 people with bipolar 
# sum(duplicated(Bipolar$FID)) 
# 
# Bipolar <- as.data.frame(Bipolar)
# unique_fam_bipolar <- Bipolar[!duplicated(Bipolar$FID),] #9118
# 
# parents_of_bip <- select(unique_fam_bipolar, c("RINPERSOONpa", "RINPERSOONMa"))
# 
# 
# yrs_of_father <- merge(parents_of_bip, edu_pop, by.y = "RINPERSOON",  by.x = "RINPERSOONpa") # 0
# yrs_of_mother <- merge(parents_of_bip, edu_pop, by.y = "RINPERSOON",  by.x = "RINPERSOONMa") #3
# 
# 
# # The sample size of parents who have reported EA is too low 
# 
# 
