# Health data 
# Diagnoses 
# Author: Perline Demange
# Project: CBS-MR EA-MH 
# 


# Diagnoses #############
# Load gba6585 to get ID of individuals 
gba <- fread("gba_6585_20201015.csv")
gba <- gba[,2 ]

# 1. Identify Disorders ####
## 1.1 Load summary tables for disorders datasets information -----
codefile <- as.data.frame(fread("H:/Data/Diagnoses/Diagnoses_code_preregistered_autism.csv", header=T))
data_list <- as.data.frame(fread("H:/Data/Diagnoses/Diagnoses_file_list.csv", header=T))


## 1.2 Run function to get disorders status (main and secondary diagnoses), per year, ####
# and save intermediate results 
gba_allyears <- get_diagnoses(codefile, data_list, gba)
# Documentation of the function get_diagnoses
# The function returns a table with one row per individuals. For every year, 
#     for every diagnoses (both main and secondary, and total), it states True or False 
#     if the individual was diagnosed. 
# Warning: Only works if the all disorders code containing the string code provided are 
# part of the disorder of interest
# Arguments: 
#   - codefile: data frame with two columns, $trait containing the name of the disorder, 
#       $code containing the DSM IV code for the disorder (if several codes they need to be separated by | )
#   - data_list: data frame with three columns: $year: year of the data, 
#       $main: data with path to the main diagnoses data,
#       $second: data with path to the secondary diagnoses data
#   - gba: column of individual ID 
# Output: data table with one column for the year, and three columns for each disorder
#     (main and secondary diagnoses, and both)


gba_allyears <- as.data.frame(gba_allyears)
rm(gba)

write.csv(gba_allyears, "gba6585_diagnoses_allyears_detailed_20230703.csv", row.names = F) # doesnt work becasue it takes forever... 
#gba_allyears <- fread("gba6585_diagnoses_allyears_detailed_20210219.csv")
#gba_allyears <- as.data.frame(gba_allyears)

##  1.3 Get full disorder status (main and secondary pulled together), per year ######
trait_list <- codefile$trait 
my_vars <- c( "RINPERSOON", "year",  trait_list)
gba_allyears_full_diagnosis <- gba_allyears[my_vars]
summary(gba_allyears_full_diagnosis$ASD) # do not summary all, dataset is too big

## 1.4 Add variable "any diagnoses" and "nb of disorders" ######
gba_allyears_full_diagnosis <- gba_allyears_full_diagnosis %>%
  mutate(nb_disorders = rowSums(across(where(is.logical))))

summary(gba_allyears_full_diagnosis$nb_disorders) 
table(gba_allyears_full_diagnosis$nb_disorders)


gba_allyears_full_diagnosis$Any <- ifelse(gba_allyears_full_diagnosis$nb_disorders > 0, T,F)
summary(gba_allyears_full_diagnosis$Any)
head(gba_allyears_full_diagnosis[gba_allyears_full_diagnosis$Any ==T,])
head(gba_allyears_full_diagnosis[gba_allyears_full_diagnosis$nb_disorders >1,])

write.csv(gba_allyears_full_diagnosis, "gba6585_diagnoses_allyears_full_diagnosis_20230703.csv", row.names = F)  # not sure it worked
#gba_allyears_full_diagnosis <- fread("gba6585_diagnoses_allyears_full_diagnosis_20210219.csv")
rm(gba_allyears)

# 2.Descriptive ###########

# trait_list needs to be changed to include "Any"
trait_list <- c(trait_list, "Any")

## 2.1. In total population gba6585 ----

# In total how many times one person had the diagnosis over the years 
#count_diagnoses_allyears_per_person <- gba_allyears_full_diagnosis %>%
#  group_by(RINPERSOON) %>%
#  summarise_if(is.logical, sum)


### 2.1.1 Get overall disorders status per person (across years) --------
diagnoses_totalyears_per_person_gba <- get_diagnoses_totalyear_per_person(gba_allyears_full_diagnosis, "gba")
# - write.csv(diagnoses_totalyears_per_person_gba, "gba_diagnoses_totalyear_per_person_20210219.csv", row.names = F)
# - diagnoses_totalyears_per_person_gba <- fread("gba_diagnoses_totalyear_per_person_20210219.csv")

### 2.2.3 Get descriptive table of diagnoses --------
descriptive_table_gba <- get_descriptive_table(gba_allyears_full_diagnosis,
                                               diagnoses_totalyears_per_person_gba, 
                                               "gba", trait_list)
descriptive_table_gba


## 2.2 In population gba6585 with EA ------
### 2.2.1 Load data and merge --------
gba_EA <- fread("gba6585_edu_subset_20210219.csv")
gba_EA_list <- gba_EA$RINPERSOON #3305733

gba_EA_allyears_full_diagnosis <- gba_allyears_full_diagnosis[which(gba_allyears_full_diagnosis$RINPERSOON %in% gba_EA_list),]
#nrow should be 3305733 * number of years
nrow(gba_EA_allyears_full_diagnosis)/6 #good

### 2.2.2 Get total diagnoses per person --------
diagnoses_totalyears_per_person_gba_EA <- get_diagnoses_totalyear_per_person(gba_EA_allyears_full_diagnosis, "gba_EA")
write.csv(diagnoses_totalyears_per_person_gba_EA, "gba_EA_diagnoses_totalyear_per_person_20230921.csv", row.names = F)
# - diagnoses_totalyears_per_person_gba_EA <- fread("gba_EA_diagnoses_totalyear_per_person_20210219.csv")

### 2.2.3 Get descriptive table of diagnoses --------
descriptive_table_gba_EA <- get_descriptive_table(gba_EA_allyears_full_diagnosis, diagnoses_totalyears_per_person_gba_EA, "gba_EA", trait_list)
descriptive_table_gba_EA


## 2.3 In total sibling pop -------
### 2.3.1 Load data and merge --------
sib <- fread("gba_6585_FID_fullsibship_final_20210215.csv")
head(gba_allyears_full_diagnosis)
sib_list <- sib$RINPERSOON #3234923

sib_allyears_full_diagnosis <- gba_allyears_full_diagnosis[which(gba_allyears_full_diagnosis$RINPERSOON %in% sib_list),]
#nrow should be 3234923 * number of years
nrow(sib_allyears_full_diagnosis)/6 #good 

### 2.3.2 Get total diagnoses per person --------
diagnoses_totalyears_per_person_sib <- get_diagnoses_totalyear_per_person(sib_allyears_full_diagnosis, "siblings")
write.csv(diagnoses_totalyears_per_person_sib, "siblings_diagnoses_totalyear_per_person_20230921.csv", row.names = F)
# - diagnoses_totalyears_per_person_sib <- fread("siblings_diagnoses_totalyear_per_person_20210219.csv")

### 2.3.3 Get descriptive table of diagnoses --------
descriptive_table_sib <- get_descriptive_table(sib_allyears_full_diagnosis, diagnoses_totalyears_per_person_sib, "siblings", trait_list)
descriptive_table_sib

## 2.4 prevalence in sibling pop with EA ---------
sib_EA <- fread("sib_edu_subset_20210219.csv")
sib_EA_list <- sib_EA$RINPERSOON #1743032

sib_EA_allyears_full_diagnosis <- gba_allyears_full_diagnosis[which(gba_allyears_full_diagnosis$RINPERSOON %in% sib_EA_list),]
#nrow should be  1743032 * number of years
nrow(sib_EA_allyears_full_diagnosis)/6 #good 

### 2.4.1 Get total diagnoses per person --------
diagnoses_totalyears_per_person_sib_EA <- get_diagnoses_totalyear_per_person(sib_EA_allyears_full_diagnosis, "siblings_EA")
write.csv(diagnoses_totalyears_per_person_sib_EA, "siblings_EA_diagnoses_totalyear_per_person_20230920.csv", row.names = F)
# - diagnoses_totalyears_per_person_sib_EA <-fread("siblings_EA_diagnoses_totalyear_per_person_20210219.csv")

### 2.4.2 Get descriptive table of diagnoses --------
descriptive_table_sib_EA <- get_descriptive_table(sib_EA_allyears_full_diagnosis, diagnoses_totalyears_per_person_sib_EA, "siblings_EA", trait_list)
descriptive_table_sib_EA

## 2.5 All populations together -------
descriptive_table <- rbind(descriptive_table_gba, descriptive_table_gba_EA, descriptive_table_sib, descriptive_table_sib_EA)
descriptive_table

write.csv(descriptive_table, "descriptive_table_diagnoses_allpop_20230707.csv", row.names = F)
write.xlsx(descriptive_table, "descriptive_table_diagnoses_allpop_20230707.xlsx", row.names = F)
#descriptive_table <- fread("descriptive_table_diagnoses_allpop_2021025.csv")

#remove ODD and Conduct for output
descriptive_table_output <- subset(descriptive_table, select= -c(Conduct, ODD))
write.csv(descriptive_table_output, "descriptive_table_diagnoses_allpop_20230707_output.csv", row.names = F)
write.xlsx(descriptive_table_output, "descriptive_table_diagnoses_allpop_20230707_output.xlsx", row.names = F)

# 3. Descriptives per individual characteristics in gba ea  subset####

#diagnoses_totalyears_per_person_gba_EA <- fread("gba_EA_diagnoses_totalyear_per_person_20210219.csv")
total_data <- merge(gba_EA, diagnoses_totalyears_per_person_gba_EA, by = "RINPERSOON") #3305733

col_subset <- c("RINPERSOON", "GBAGESLACHT", "yrs", "OPLNIVSOI2016AGG4HBMETNIRWO",  trait_list)
sum_data <- subset(total_data, select=col_subset)
sum_data

#add total eating disorders
sum_data$Eating <- ifelse(sum_data$Anorexia == 1 | sum_data$Bulimia == 1, 1, 0)
table(sum_data$Eating)
table(sum_data$Anorexia)
table(sum_data$Bulimia)

nrow(sum_data[sum_data$Anorexia ==1 & sum_data$Bulimia ==1, ]) # 91 people with both anorexia and bulimia, fits with sum Eating 

trait_list <- c(trait_list, "Eating")

## 3.1 Per sex -----
total_count <- sum_data %>%
  group_by(GBAGESLACHT) %>%
  summarise(number = n())

diagnoses_count <- sum_data %>%
  group_by(GBAGESLACHT) %>%
  summarise_all(sum)

preval <- diagnoses_count/ total_count$number
prevalence_CI <- 1.96*sqrt((preval*(1-preval))/total_count$number)
preval$measure <- "prevalence"
prevalence_CI$measure <- "prevalence_CI"
diagnoses_count$measure <- "count"
diagnoses_by_sex <- rbind(preval,prevalence_CI)
diagnoses_by_sex <- rbind(diagnoses_by_sex, as.data.frame(diagnoses_count))
diagnoses_by_sex$sample_size <- total_count$number
diagnoses_by_sex$sex <- diagnoses_count$GBAGESLACHT
diagnoses_by_sex <- subset(diagnoses_by_sex, select= c("sex", "measure", "sample_size", trait_list))
diagnoses_by_sex

write.csv(diagnoses_by_sex, "diagnoses_by_sex_gba_20230707.csv", row.names=F)
write.xlsx(diagnoses_by_sex, "diagnoses_by_sex_gba_20230707.xlsx", row.names=F)


## 3.2 Per EA in years ------
total_count <- sum_data %>%
  group_by(yrs) %>%
  summarise(number = n())

diagnoses_count <- sum_data %>%
  group_by(yrs) %>%
  summarise_all(sum)

preval <- diagnoses_count/ total_count$number
prevalence_CI <- 1.96*sqrt((preval*(1-preval))/total_count$number)
preval$measure <- "prevalence"
prevalence_CI$measure <- "prevalence_CI"
diagnoses_count$measure <- "count"
diagnoses_by_yrs <- rbind(preval,prevalence_CI)
diagnoses_by_yrs <- rbind(diagnoses_by_yrs, as.data.frame(diagnoses_count))
diagnoses_by_yrs$sample_size <- total_count$number
diagnoses_by_yrs$yrs <- diagnoses_count$yrs
diagnoses_by_yrs <- subset(diagnoses_by_yrs, select= c("yrs", "measure", "sample_size", trait_list))
diagnoses_by_yrs

write.csv(diagnoses_by_yrs, "diagnoses_by_yrs_gba_20230707.csv", row.names=F)
write.xlsx(diagnoses_by_yrs, "diagnoses_by_yrs_gba_20230707.xlsx", row.names=F)

# make it output safe 
# remove CDD and Conduct and anorexia and bulimia (keep eating disorder sum)
diagnoses_by_yrs_output <- subset(diagnoses_by_yrs, select= -c(Conduct, ODD, Anorexia, Bulimia))
diagnoses_by_yrs_output[diagnoses_by_yrs_output$measure== 'count',]
# remove make ClusterA yrs 22 NA 
diagnoses_by_yrs_output[diagnoses_by_yrs_output$yrs == 22,]$ClusterA <- NA
diagnoses_by_yrs_output

write.csv(diagnoses_by_yrs_output, "diagnoses_by_yrs_gba_20230707_output.csv", row.names=F)
write.xlsx(diagnoses_by_yrs_output, "diagnoses_by_yrs_gba_20230707_output.xlsx", row.names=F)

## 3.3 Per EA in years and per sex  ------
total_count <- sum_data %>%
  group_by(yrs, GBAGESLACHT) %>%
  summarise(number = n())

diagnoses_count <- sum_data %>%
  group_by(yrs, GBAGESLACHT) %>%
  summarise_all(sum)

preval <- diagnoses_count / total_count$number
prevalence_CI <- 1.96*sqrt((preval*(1-preval))/total_count$number)
preval$measure <- "prevalence"
prevalence_CI$measure <- "prevalence_CI"
diagnoses_count$measure <- "count"
diagnoses_by_yrs_sex <- rbind(preval,prevalence_CI)
diagnoses_by_yrs_sex <- rbind(diagnoses_by_yrs_sex, as.data.frame(diagnoses_count))
diagnoses_by_yrs_sex$sample_size <- total_count$number
diagnoses_by_yrs_sex$yrs <- diagnoses_count$yrs
diagnoses_by_yrs_sex$sex <- diagnoses_count$GBAGESLACHT
diagnoses_by_yrs_sex <- subset(diagnoses_by_yrs_sex, select= c("yrs","sex", "measure", "sample_size", trait_list))

diagnoses_by_yrs_sex
diagnoses_by_yrs_sex[diagnoses_by_yrs_sex$measure == "count",]

write.csv(diagnoses_by_yrs_sex, "diagnoses_by_yrs_sex_gba_20230707.csv", row.names=F)
write.xlsx(diagnoses_by_yrs_sex, "diagnoses_by_yrs_sex_gba_20230707.xlsx", row.names=F)

#For output 
diagnoses_by_yrs_sex_output <- subset(diagnoses_by_yrs_sex, select= -c(Conduct, ODD, Anorexia, Bulimia))
diagnoses_by_yrs_sex_output[diagnoses_by_yrs_sex_output$measure == "count",]

trait_list_sub <- trait_list
trait_list_sub <- trait_list_sub[trait_list_sub != "Conduct"]
trait_list_sub <- trait_list_sub[trait_list_sub != "ODD"]  
trait_list_sub <- trait_list_sub[trait_list_sub != "Anorexia"] 
trait_list_sub <- trait_list_sub[trait_list_sub != "Bulimia"]  

for (trait in trait_list_sub){
  print(trait)
  for (num in c(2,8,11,12,13,14,15,17,18,22)){
    print(num)
    for (sex in c(2,1)){ 
      print(sex)
      if (diagnoses_by_yrs_sex_output[diagnoses_by_yrs_sex_output$measure == "count" &
                                     diagnoses_by_yrs_sex_output$sex == sex &
                                     diagnoses_by_yrs_sex_output$yrs == num,][[trait]] <= 10){ 
         diagnoses_by_yrs_sex_output[diagnoses_by_yrs_sex_output$sex == sex & diagnoses_by_yrs_sex_output$yrs == num,][[trait]] <- NA} 
    } 
  }
}

diagnoses_by_yrs_sex_output
write.csv(diagnoses_by_yrs_sex_output, "diagnoses_by_yrs_sex_gba_20230707_output.csv", row.names=F)
write.xlsx(diagnoses_by_yrs_sex_output, "diagnoses_by_yrs_sex_gba_20230707_output.xlsx", row.names=F)

#diagnoses_by_yrs_sext <- fread("diagnoses_by_yrs_sex_gba_20210225_output.csv")
mdd_test <- subset(diagnoses_by_yrs_sext, select=c("yrs","sex", "measure", "sample_size","MDD"))
mdd_test_prev <- mdd_test[mdd_test$measure == "prevalence", ]
mdd_test_prevCI <- mdd_test[mdd_test$measure == "prevalence_CI", ]
mdd_test_prev$CI <- mdd_test_prevCI$MDD
mdd_test_prev

ggplot(data=mdd_test_prev, aes(x=yrs, y=MDD, color=sex))+
  geom_point() + 
  geom_errorbar(aes(x=yrs, ymin=MDD-CI, ymax=MDD+CI), position="dodge")


## 3.4 Per EA in diploma ------
total_count <- sum_data %>%
  group_by(OPLNIVSOI2016AGG4HBMETNIRWO) %>%
  summarise(number = n())

diagnoses_count <- sum_data %>%
  group_by(OPLNIVSOI2016AGG4HBMETNIRWO) %>%
  summarise_all(sum)

preval <- diagnoses_count / total_count$number
prevalence_CI <- 1.96*sqrt((preval*(1-preval))/total_count$number)
preval$measure <- "prevalence"
prevalence_CI$measure <- "prevalence_CI"
diagnoses_count$measure <- "count"
diagnoses_by_diploma <- rbind(preval,prevalence_CI)
diagnoses_by_diploma <- rbind(diagnoses_by_diploma, as.data.frame(diagnoses_count))
diagnoses_by_diploma$sample_size <- total_count$number
diagnoses_by_diploma$diploma <- diagnoses_count$OPLNIVSOI2016AGG4HBMETNIRWO
diagnoses_by_diploma <- subset(diagnoses_by_diploma, select= c("diploma", "measure",
                                                               "sample_size", trait_list))
diagnoses_by_diploma

diagnoses_by_diploma[diagnoses_by_diploma$measure == "count",]
write.csv(diagnoses_by_diploma, "diagnoses_by_diploma_gba_20230707.csv", row.names=F)
write.xlsx(diagnoses_by_diploma, "diagnoses_by_diploma_gba_20230707.xlsx", row.names=F)

# For output
diagnoses_by_diploma_output <- subset(diagnoses_by_diploma, select= -c(Conduct, ODD, Anorexia, Bulimia))
diagnoses_by_diploma_output[diagnoses_by_diploma_output$measure == "count",]

diagnoses_by_diploma_output[diagnoses_by_diploma_output$diploma == 3213,]$ClusterA <- NA
diagnoses_by_diploma_output
write.csv(diagnoses_by_diploma_output, "diagnoses_by_diploma_gba_20230707_output.csv", row.names=F)
write.xlsx(diagnoses_by_diploma_output, "diagnoses_by_diploma_gba_20230707_output.xlsx", row.names=F)

# 4. Descriptives per individual characteristics in siblings subset####
#diagnoses_totalyears_per_person_sib_EA <- fread("siblings_EA_diagnoses_totalyear_per_person_20210219.csv")
total_data <- merge(sib_EA, diagnoses_totalyears_per_person_sib_EA, by = "RINPERSOON")

col_subset <- c("RINPERSOON", "GBAGESLACHT", "yrs", "OPLNIVSOI2016AGG4HBMETNIRWO",  "birth_order",  trait_list)
sum_data <- subset(total_data, select=col_subset)
sum_data

#add total eating disorders
sum_data$Eating <- ifelse(sum_data$Anorexia == 1 | sum_data$Bulimia == 1, 1, 0)
table(sum_data$Eating)
table(sum_data$Anorexia)
table(sum_data$Bulimia)
table(sum_data$ASD)

nrow(sum_data[sum_data$Anorexia ==1 & sum_data$Bulimia ==1, ]) # 61 people with both anorexia and bulimia, fits with sum Eating 

trait_list <- c(trait_list, "Eating")


## 4.1 Per sex -----
total_count <- sum_data %>%
  group_by(GBAGESLACHT) %>%
  summarise(number = n())

diagnoses_count <- sum_data %>%
  group_by(GBAGESLACHT) %>%
  summarise_all(sum)

preval <- diagnoses_count*100 / total_count$number
preval$measure <- "percent"
diagnoses_count$measure <- "count"
diagnoses_by_sex <- rbind(preval,diagnoses_count)
diagnoses_by_sex$sample_size <- total_count$number
diagnoses_by_sex$sex <- c(1,2)
diagnoses_by_sex <- subset(diagnoses_by_sex, select= c("sex", "measure", "sample_size", trait_list))
diagnoses_by_sex

write.csv(diagnoses_by_sex, "diagnoses_by_sex_20230707.csv", row.names=F)
write.xlsx(diagnoses_by_sex, "diagnoses_by_sex_20230707.xlsx", row.names=F)
#diagnoses_by_sex <- fread("diagnoses_by_sex_20210219.csv")

## 4.2 Per EA in years ------
total_count <- sum_data %>%
  group_by(yrs) %>%
  summarise(number = n())

diagnoses_count <- sum_data %>%
  group_by(yrs) %>%
  summarise_all(sum)

preval <- diagnoses_count*100 / total_count$number
preval$measure <- "percent"
diagnoses_count$measure <- "count"
diagnoses_by_yrs <- rbind(preval,diagnoses_count)
diagnoses_by_yrs$sample_size <- total_count$number
diagnoses_by_yrs$yrs <- diagnoses_count$yrs
diagnoses_by_yrs <- subset(diagnoses_by_yrs, select= c("yrs", "measure", "sample_size", trait_list))
diagnoses_by_yrs

write.csv(diagnoses_by_yrs, "diagnoses_by_yrs_20230707.csv", row.names=F)
write.xlsx(diagnoses_by_yrs, "diagnoses_by_yrs_20230707.xlsx", row.names=F)
#diagnoses_by_yrs <- fread("diagnoses_by_yrs_20220729.csv")


# make it output safe 
# remove CDD and Conduct and anorexia and bulimia (keep eating disorder sum)
diagnoses_by_yrs_output <- subset(diagnoses_by_yrs, select= -c(Conduct, ODD, Anorexia, Bulimia))
diagnoses_by_yrs_output[diagnoses_by_yrs_output$measure== 'count',]
# remove make ClusterA, BIP2 yrs 22 NA 
diagnoses_by_yrs_output[diagnoses_by_yrs_output$yrs == 22,]$ClusterA <- NA
diagnoses_by_yrs_output[diagnoses_by_yrs_output$yrs == 22,]$Bipolar_2 <- NA
diagnoses_by_yrs_output

write.csv(diagnoses_by_yrs_output, "diagnoses_by_yrs_20230707_output.csv", row.names=F)
write.xlsx(diagnoses_by_yrs_output, "diagnoses_by_yrs_20230707_output.xlsx", row.names=F)

## 4.3 Per EA in years and per sex  ------
total_count <- sum_data %>%
  group_by(yrs, GBAGESLACHT) %>%
  summarise(number = n())

diagnoses_count <- sum_data %>%
  group_by(yrs, GBAGESLACHT) %>%
  summarise_all(sum)

preval <- diagnoses_count / total_count$number
prevalence_CI <- 1.96*sqrt((preval*(1-preval))/total_count$number)
preval$measure <- "prevalence"
prevalence_CI$measure <- "prevalence_CI"
diagnoses_count$measure <- "count"
diagnoses_by_yrs_sex <- rbind(preval,prevalence_CI)
diagnoses_by_yrs_sex <- rbind(diagnoses_by_yrs_sex, as.data.frame(diagnoses_count))
diagnoses_by_yrs_sex$sample_size <- total_count$number
diagnoses_by_yrs_sex$yrs <- diagnoses_count$yrs
diagnoses_by_yrs_sex$sex <- diagnoses_count$GBAGESLACHT
diagnoses_by_yrs_sex <- subset(diagnoses_by_yrs_sex, select= c("yrs","sex", "measure", "sample_size", trait_list))
diagnoses_by_yrs_sex[diagnoses_by_yrs_sex$measure == "count",]

write.csv(diagnoses_by_yrs_sex, "diagnoses_by_yrs_sex_20230707.csv", row.names=F)
write.xlsx(diagnoses_by_yrs_sex, "diagnoses_by_yrs_sex_20230707.xlsx", row.names=F)

#For output 
diagnoses_by_yrs_sex_output <- subset(diagnoses_by_yrs_sex, select= -c(Conduct, ODD, Anorexia, Bulimia))
diagnoses_by_yrs_sex_output[diagnoses_by_yrs_sex_output$measure == "count",]

trait_list_sub <- trait_list
trait_list_sub <- trait_list_sub[trait_list_sub != "Conduct"]
trait_list_sub <- trait_list_sub[trait_list_sub != "ODD"]  
trait_list_sub <- trait_list_sub[trait_list_sub != "Anorexia"] 
trait_list_sub <- trait_list_sub[trait_list_sub != "Bulimia"]  

for (trait in trait_list_sub){
  print(trait)
  for (num in c(2,8,11,12,13,14,15,17,18,22)){
    print(num)
    for (sex in c(2,1)){ 
      print(sex)
      if (diagnoses_by_yrs_sex_output[diagnoses_by_yrs_sex_output$measure == "count" &
                                      diagnoses_by_yrs_sex_output$sex == sex &
                                      diagnoses_by_yrs_sex_output$yrs == num,][[trait]] <= 10){ 
        diagnoses_by_yrs_sex_output[diagnoses_by_yrs_sex_output$sex == sex & diagnoses_by_yrs_sex_output$yrs == num,][[trait]] <- NA} 
    } 
  }
}

tail(diagnoses_by_yrs_sex_output)
write.csv(diagnoses_by_yrs_sex_output, "diagnoses_by_yrs_sex_20230707_output.csv", row.names=F)
write.xlsx(diagnoses_by_yrs_sex_output, "diagnoses_by_yrs_sex_20230707_output.xlsx", row.names=F)

## 4.4 Per EA in diploma ------
total_count <- sum_data %>%
  group_by(OPLNIVSOI2016AGG4HBMETNIRWO) %>%
  summarise(number = n())

diagnoses_count <- sum_data %>%
  group_by(OPLNIVSOI2016AGG4HBMETNIRWO) %>%
  summarise_all(sum)

preval <- diagnoses_count*100 / total_count$number
preval$measure <- "percent"
diagnoses_count$measure <- "count"
diagnoses_by_diploma <- rbind(preval,diagnoses_count)
diagnoses_by_diploma$sample_size <- total_count$number
diagnoses_by_diploma$diploma <- diagnoses_count$OPLNIVSOI2016AGG4HBMETNIRWO
diagnoses_by_diploma <- subset(diagnoses_by_diploma, select= c("diploma", "measure", "sample_size", trait_list))
diagnoses_by_diploma

write.csv(diagnoses_by_diploma, "diagnoses_by_diploma_20230707.csv", row.names=F)
write.xlsx(diagnoses_by_diploma, "diagnoses_by_diploma_20230707.xlsx", row.names=F)

# For output
diagnoses_by_diploma_output <- subset(diagnoses_by_diploma, select= -c(Conduct, ODD, Anorexia, Bulimia))
diagnoses_by_diploma_output[diagnoses_by_diploma_output$measure == "count",]

diagnoses_by_diploma_output[diagnoses_by_diploma_output$diploma == 3213,]$ClusterA <- NA
diagnoses_by_diploma_output[diagnoses_by_diploma_output$diploma == 3213,]$Bipolar_2 <- NA
diagnoses_by_diploma_output
write.csv(diagnoses_by_diploma_output, "diagnoses_by_diploma_20230707_output.csv", row.names=F)
write.xlsx(diagnoses_by_diploma_output, "diagnoses_by_diploma_20230707_output.xlsx", row.names=F)

## 4.5 Per birth order --------
sum_data$birth_order_cat <- sum_data$birth_order
sum_data$birth_order_cat[sum_data$birth_order >= 3] <- "3+"
table(sum_data$birth_order_cat)

total_count <- sum_data %>%
  group_by(birth_order_cat) %>%
  summarise(number = n())

diagnoses_count <- sum_data %>%
  group_by(birth_order_cat) %>%
  summarise_all(sum)

preval <- diagnoses_count[-1]*100 / total_count[-1]$number
preval$measure <- "percent"
preval$birth_order_cat <- diagnoses_count$birth_order_cat
diagnoses_count$measure <- "count"
diagnoses_by_birth_order <- rbind(preval,diagnoses_count)
diagnoses_by_birth_order$sample_size <- total_count$number
diagnoses_by_birth_order$birth_order <- diagnoses_count$birth_order_cat
diagnoses_by_birth_order <- subset(diagnoses_by_birth_order, select= c("birth_order", "measure", "sample_size", trait_list))
diagnoses_by_birth_order

write.csv(diagnoses_by_birth_order, "diagnoses_by_birth_order_20230707.csv", row.names=F)
write.xlsx(diagnoses_by_birth_order, "diagnoses_by_birth_order_20230707.xlsx", row.names=F)

# For output
diagnoses_by_birth_order_output <- subset(diagnoses_by_birth_order, select= -c(Conduct, ODD))
diagnoses_by_birth_order_output[diagnoses_by_birth_order_output$measure == "count",]

write.csv(diagnoses_by_birth_order_output, "diagnoses_by_birth_order_output_20230707.csv", row.names=F)
write.xlsx(diagnoses_by_birth_order_output, "diagnoses_by_birth_order_output_20230707.xlsx", row.names=F)

## 4.6 Per year of birth ######
head(sib_EA)
birthdate <- subset(sib_EA, select= c(RINPERSOON, GBAGEBOORTEJAAR))
head(birthdate)

birthdate_sum_data <- merge(sum_data, birthdate, by = "RINPERSOON")
head(birthdate_sum_data)

birthdate_sum_data$birthyear <- birthdate_sum_data$GBAGEBOORTEJAAR
birthdate_sum_data$birthyear[birthdate_sum_data$GBAGEBOORTEJAAR >= 1981] <- "81-85"
birthdate_sum_data$birthyear[birthdate_sum_data$birthyear >= 1976 & 
                               birthdate_sum_data$birthyear < 1981] <- "76-80"
birthdate_sum_data$birthyear[birthdate_sum_data$birthyear >= 1971& 
                               birthdate_sum_data$birthyear < 1976] <- "71-75"
birthdate_sum_data$birthyear[birthdate_sum_data$birthyear >= 1965& 
                               birthdate_sum_data$birthyear < 1971] <- "65-70"
table(birthdate_sum_data$birthyear)

birthdate_sum_data <- subset(birthdate_sum_data, select = -c(birth_order_cat))
total_count <- birthdate_sum_data %>%
  group_by(birthyear) %>%
  summarise(number = n())

diagnoses_count <- birthdate_sum_data %>%
  group_by(birthyear) %>%
  summarise_all(sum)


preval <- diagnoses_count[-1]*100 / total_count[-1]$number
preval$measure <- "percent"
preval$birthyear <- diagnoses_count$birthyear
diagnoses_count$measure <- "count"
diagnoses_by_birthyear <- rbind(preval,diagnoses_count)
diagnoses_by_birthyear$sample_size <- total_count$number
diagnoses_by_birthyear$birthyear <- diagnoses_count$birthyear
diagnoses_by_birthyear <- subset(diagnoses_by_birthyear, select= c("birthyear", "measure", "sample_size", trait_list))
diagnoses_by_birthyear

write.csv(diagnoses_by_birthyear, "diagnoses_by_birthyear_20230707.csv", row.names=F)
write.xlsx(diagnoses_by_birthyear, "diagnoses_by_birthyear_20230707.xlsx", row.names=F)

# For output
diagnoses_by_birthyear_output <- subset(diagnoses_by_birthyear, select= -c(Conduct, ODD))
diagnoses_by_birthyear_output[diagnoses_by_birthyear_output$measure == "count",]

write.csv(diagnoses_by_birthyear_output, "diagnoses_by_birthyear_output_20230707.csv", row.names=F)
write.xlsx(diagnoses_by_birthyear_output, "diagnoses_by_birthyear_output_20230707.xlsx", row.names=F)

# 5. Co-occurence table ------- 
## 5.1 In gba6585 #############
### 5.1.1 Percentage matrix -------
#diagnoses_totalyears_per_person_gba <- fread("gba_diagnoses_totalyear_per_person_20210219.csv")

nrow(diagnoses_totalyears_per_person_gba[diagnoses_totalyears_per_person_gba$Any ==1,])

head(diagnoses_totalyears_per_person_gba)
cooc <- as.data.frame(diagnoses_totalyears_per_person_gba[,2:24])
cooc <- cooc[, -c(3:4)] # Remove Conduct and ODD 
head(cooc)
cooc <- as.matrix(cooc)
out <- crossprod(cooc) #symetric matrix with number of individuals diagnosed with both x and y 
diag(out) # number of individuals diagnoses with each disorder 

out # 
write.table(out, "matrix_cooccurrence_count_20230707.txt")
write.xlsx(out, "matrix_cooccurrence_count_20230707.xlsx")

test <- t(apply(out,1, function(x)x/diag(out))) # columns are ratio people with the column disorder who also have the row disorder 
test_melt <- melt(test) # var 1 is row variable, so value is ratio people with var 2 disorder who also have var 1 
ggplot(data=test_melt, aes(x=Var1, y=Var2, fill=value))+ 
  geom_tile(color="white") +
  scale_fill_gradient2(low='blue', high="red", mid="white")

write.table(test, "matrix_cooccurrence_freq_20230707.txt")
write.xlsx(test, "matrix_cooccurrence_freq_20230707.xlsx")

test <- read.table("matrix_cooccurrence_freq_20210223.txt")
test


### 5.1.2 Polychoric correlations #####
library(psych)
polycor <- polychoric(diagnoses_totalyears_per_person_gba[, 2:25][, -c(3:4)]) #remove ocd and conduct
polycor
#Warning message:
#In cor.smooth(mat) : Matrix was not positive definite, smoothing was done
matrix_plycor <- polycor$rho
write.table(matrix_plycor, "matrix_polychoriccorrelation_20230707.csv")
write.xlsx(matrix_plycor, "matrix_polychoriccorrelation_20230707.xlsx")

fig <- melt(matrix_plycor)
ggplot(data=fig, aes(x=Var1, y=Var2, fill=value))+ 
  geom_tile(color="white") +
  scale_fill_gradient2(high="red", mid="white")

## 5.2 In siblings EA subset #####
### 5.2.1 Percentage matrix -------
#diagnoses_totalyears_per_person_sib_EA <- fread("siblings_EA_diagnoses_totalyear_per_person_20210219.csv")

diagnoses_totalyears_per_person_sib_EA$Eating <- ifelse(
  diagnoses_totalyears_per_person_sib_EA$Anorexia == 1 |
    diagnoses_totalyears_per_person_sib_EA$Bulimia == 1, 1, 0)


cooc <- diagnoses_totalyears_per_person_sib_EA[,2:26]
cooc <- cooc[, -c(3:4)] # Remove Conduct and ODD 
head(cooc)
cooc <- cooc[, -22]#remove any
cooc <- cooc[, -c(16:17)]#remove anorexia and bulimia
coocm <- as.matrix(cooc)
out <- crossprod(coocm) #symetric matrix with number of individuals diagnosed with both x and y 
diag(out) # number of individuals diagnoses with each disorder 

out #not clean for output, need to remove schizophrenia * eating 
write.table(out, "matrix_cooccurrence_count_sibEA_20230711.txt")
write.xlsx(out, "matrix_cooccurrence_count_sibEA_20230711.xlsx")

test <- t(apply(out,1, function(x)x/diag(out))) # columns are ratio people with the column disorder who also have the row disorder 
test_melt <- melt(test) # var 1 is row variable, so value is ratio people with var 2 disorder who also have var 1 
ggplot(data=test_melt, aes(x=Var1, y=Var2, fill=value))+ 
  geom_tile(color="white") +
  scale_fill_gradient2(low='blue', high="red", mid="white")

write.table(test, "matrix_cooccurrence_freq_sibEA_20230711.txt")
write.xlsx(test, "matrix_cooccurrence_freq_sibEA_20230711.xlsx")

### 5.1.2 Polychoric correlations #####
library(psych)
polycor <- polychoric(cooc)
polycor
matrix_plycor <- polycor$rho
write.table(matrix_plycor, 
            "matrix_polychoriccorrelation_sibEA_20230711.csv")
write.xlsx(matrix_plycor, 
            "matrix_polychoriccorrelation_sibEA_20230711.xlsx")

fig <- melt(matrix_plycor)
ggplot(data=fig, aes(x=Var1, y=Var2, fill=value))+ 
  geom_tile(color="white") +
  scale_fill_gradient2(low='blue', high="red", mid="white")

# for all these sibea results, not clean for output, need to remove schizophrenia * eating 

#6. Investigate primary vs secondary diagnoses #############

# Just to have an idea 
# Individuals can have several primary diagnoses in the timeframe (one different per year) 
# I therefore look at the total number of time a disorder was given as primary or secondary 
# This is only a crude measure, as some disorders might call for a yearly diagnosis while others not

## 6.1 In gba #####

#Keep only primary and secondary, so remove combined 
codefile <- as.data.frame(fread("H:/Data/Diagnoses/Diagnoses_code_preregistered_autism.csv",
                                header=T))
trait_list <- codefile$trait
gba_allyears_primsec <- gba_allyears[,!(names(gba_allyears) %in% trait_list)]
rm(gba_allyears)

primsec_totalyears_per_person_gba <- get_diagnoses_totalyear_per_person(gba_allyears_primsec,
                                                                          "gba")
rm(gba_allyears_primsec)

# count number of diagnosed individuals 
count_primsec <- primsec_totalyears_per_person_gba %>%
  summarise_at(vars(-RINPERSOON), sum)


# among people who were never diagnosed with the disorder as primary diagnose,
# how many got diagnosed with it as secondary 

main_names <- names(
  primsec_totalyears_per_person_gba[,grepl("main$",
                                    names(primsec_totalyears_per_person_gba))])

second_names <- names(primsec_totalyears_per_person_gba[,grepl("second$",
                                                             names(primsec_totalyears_per_person_gba))])

no_second <- primsec_totalyears_per_person_gba

for(col in 1:length(main_names)){
    no_second[,second_names[col]] <- 
      ifelse(no_second[,main_names[col]] == 1,0, 
         no_second[,second_names[col]])     
  } 


#check it does what I what want
head(primsec_totalyears_per_person_gba[primsec_totalyears_per_person_gba$Anxiety_main==1,])
head(no_second[no_second$Anxiety_main==1,])

# now sum the second only to get number of people with secondary without having main 
only_second <- no_second[, !names(no_second) %in% main_names] 

count_only_second <- only_second %>%
  summarise_at(vars(-RINPERSOON), sum)

count_only_second

count_primsec

names(count_only_second) <- paste(names(count_only_second), "only", sep="")

count <- cbind(count_primsec, count_only_second)
write.xlsx(count, "count_primary_secondary_20231030.xlsx", row.names = F)


## 6.2 Do same in SibEA ############
#get only sibEA 

sib_EA <- fread("sib_edu_subset_20210219.csv")
sib_EA_list <- sib_EA$RINPERSOON #1743032

primsec_totalyears_per_person_sib_EA <- primsec_totalyears_per_person_gba[which(primsec_totalyears_per_person_gba$RINPERSOON %in% sib_EA_list),]


# count number of diagnosed individuals 
count_primsec <- primsec_totalyears_per_person_sib_EA %>%
  summarise_at(vars(-RINPERSOON), sum)


# among people who were never diagnosed with the disorder as primary diagnose,
# how many got diagnosed with it as secondary 

main_names <- names(
  primsec_totalyears_per_person_sib_EA[,grepl("main$",
                                           names(primsec_totalyears_per_person_sib_EA
                                                ))])

second_names <- names(primsec_totalyears_per_person_sib_EA[,grepl("second$",
                                                               names(primsec_totalyears_per_person_sib_EA
                                                                    ))])

no_second <- primsec_totalyears_per_person_sib_EA

for(col in 1:length(main_names)){
  no_second[,second_names[col]] <- 
    ifelse(no_second[,main_names[col]] == 1,0, 
           no_second[,second_names[col]])     
} 


#check it does what I what want
head(primsec_totalyears_per_person_sib_EA
    [primsec_totalyears_per_person_sib_EA
                                        $Anxiety_main==1,])
head(no_second[no_second$Anxiety_main==1,])

# now sum the second only to get number of people with secondary without having main 
only_second <- no_second[, !names(no_second) %in% main_names] 

count_only_second <- only_second %>%
  summarise_at(vars(-RINPERSOON), sum)

count_only_second

count_primsec

names(count_only_second) <- paste(names(count_only_second), "only", sep="")

count <- cbind(count_primsec, count_only_second)
write.xlsx(count, "count_primary_secondary_sibEA_20231030.xlsx", row.names = F)




