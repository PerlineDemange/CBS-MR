# Functions 
# Author: Perline Demange
# Project: CBS-MR EA-MH 
# 


get_diagnoses <- function(codefile, data_list, gba){
  # Documentation:
  # The function returns a table with one row per individuals. For every year, 
  #     for every diagnoses (both main and secondary, and total), it states True or False 
  #     if the individual was diagnosed. 
  # Warning: Only works if the all disorders code containing the string code provided are 
  # part of the disorder of interest
  # Arguments: 
  #   - codefile: data frame with two columns, $trait containing the name of the disorder, 
  #       $code containing the DSM IV code for the disorder
  #   - data_list: data frame with three columns: $year: year of the data, 
  #       $main: data with path to the main diagnoses data,
  #       $second: data with path to the secondary diagnoses data
  #   - gba: column of individual ID 
  # Output: data table with one column for the year, and three columns for each disorder
  #     (main and secondary diagnoses, and both)
  gba_allyears <- NULL 
  trait_list <- codefile$trait #do this outside of the loop will save time 
  code_list <- codefile$code
  year <- data_list$year
  file_list_main <- data_list$main
  file_list_second <- data_list$second
  for (y in 1:nrow(data_list)) { #for each year
    gba_year <- gba
    gba_year$year <- year[y] #get year
    print(year[y])
    main <- fread(file_list_main[y]) #load data main diagnosis for the correct year
    second <- fread(file_list_second[y]) #load data second diagnosis for the correct year
    for(x in 1:nrow(codefile)){ #for each trait
      trait <- trait_list[x]
      print(trait)
      code <- code_list[x]
      main_diagnosis <- main[grep(code, main$GGZDBCHoofddiagnoseDSMIV),]$RINPERSOON # get all ID of people diagnosed with the disorder in main diagnosis 
      second_diagnosis <- second[grep(code, second$GGZDBCnevendiagnoseDSMIV),]$RINPERSOON # get all ID of people diagnosed with the disorder in second diagnosis
      gba_year[, paste(trait, "main", sep="_")] <- ifelse(gba$RINPERSOON %in% main_diagnosis, T, F) 
      gba_year[, paste(trait, "second", sep="_")] <- ifelse(gba$RINPERSOON %in% second_diagnosis, T, F) 
      gba_year[, paste(trait)] <- ifelse(gba$RINPERSOON %in% main_diagnosis| gba$RINPERSOON %in% second_diagnosis, T, F)
      #print(gba_year)
    }
    gba_allyears <- rbind(gba_allyears, gba_year)
  }
  return(gba_allyears)
}


get_diagnoses_totalyear_per_person <- function(full_diagnosis_data, population_name){
  # Documentation
  # Returns diagnoses per person for all years pulled together, 1 if they were diagnosed at least once, 0 if not
  # Arguments: 
  #   - full_diagnosis_data is the dataframe of diagnoses (main and secondary pulled together) for every person for every year
  #   - population_name is the name of the subset of the population used 
  # Careful: When running summarise, the columns are not logical anymore but just containing 0 or 1 !
  diagnoses_totalyears_per_person <- full_diagnosis_data %>%
    group_by(RINPERSOON) %>%
    summarise_if(is.logical, max)
  diagnoses_totalyears_per_person <- as.data.frame(diagnoses_totalyears_per_person)
  return(diagnoses_totalyears_per_person)
}

# old version without prevalcen CI 
# get_descriptive_table <- function(full_diagnosis_data, totalyear_diagnosis_data, population_name, trait_list){
#   # Documentation 
#   # Function returns a table countaining the count and prevalence of the diagnoses per year and in total 
#   # Arguments
#   #   - full_diagnosis_data is the dataframe of diagnoses (main and secondary pulled together) 
#   #       for every person for every year
#   #   - totalyear_diagnosis_data is the output of get_diagnoses_totalyear_per_person, 
#   #       diagnoses per person for all years pulled together
#   #   - population_name is the name of the subset of the population used 
#   #   - trait_list is codefile$trait. Codefile is a data frame with two columns, $trait containing the name of the disorder, 
#   #       $code containing the DSM IV code for the disorder
#   # 1. Count number of persons diagnosed at least once for each disorder for all years -----
#   count_diagnoses_totalyears <- totalyear_diagnosis_data %>%
#     summarise_at(vars(-RINPERSOON), sum)
#   count_diagnoses_totalyears$year <- 'total'
#   
#   # 2. Count number of persons diagnosed for each disorder every year ----
#   count_diagnoses_everyyears <- full_diagnosis_data %>%
#     group_by(year) %>%
#     summarise_if(is.logical, sum)
#   count_diagnoses_total <- rbind(count_diagnoses_totalyears, count_diagnoses_everyyears)
#   count_diagnoses_total$sample_size <- nrow(totalyear_diagnosis_data)
#   count_diagnoses_total$population <- population_name
#   count_diagnoses_total$measure <- "count"
#   
#   # 3. Prevalence calculation -----
#   sample_size <- nrow(totalyear_diagnosis_data)
#   prevalence <- count_diagnoses_total[, 1:length(trait_list)]*100/sample_size 
#   prevalence$year <- count_diagnoses_total$year
#   prevalence$sample_size <- nrow(totalyear_diagnosis_data)
#   prevalence$population <- population_name
#   prevalence$measure <- "percent"
#   
#   descriptive_table <- rbind(count_diagnoses_total, prevalence)
#   descriptive_table <- format(descriptive_table, scientific=F)
#   return(descriptive_table)
#   
# }
# 



get_descriptive_table <- function(full_diagnosis_data, totalyear_diagnosis_data, population_name, trait_list){
  # Documentation 
  # Function returns a table countaining the count and prevalence of the diagnoses per year and in total 
  # Arguments
  #   - full_diagnosis_data is the dataframe of diagnoses (main and secondary pulled together) 
  #       for every person for every year
  #   - totalyear_diagnosis_data is the output of get_diagnoses_totalyear_per_person, 
  #       diagnoses per person for all years pulled together
  #   - population_name is the name of the subset of the population used 
  #   - trait_list is codefile$trait. Codefile is a data frame with two columns, $trait containing the name of the disorder, 
  #       $code containing the DSM IV code for the disorder
  # 1. Count number of persons diagnosed at least once for each disorder for all years -----
  count_diagnoses_totalyears <- totalyear_diagnosis_data %>%
    summarise_at(vars(-RINPERSOON), sum)
  count_diagnoses_totalyears$year <- 'total'
  
  # 2. Count number of persons diagnosed for each disorder every year ----
  count_diagnoses_everyyears <- full_diagnosis_data %>%
    group_by(year) %>%
    summarise_if(is.logical, sum)
  count_diagnoses_total <- rbind(count_diagnoses_totalyears, count_diagnoses_everyyears)
  count_diagnoses_total$sample_size <- nrow(totalyear_diagnosis_data)
  count_diagnoses_total$population <- population_name
  count_diagnoses_total$measure <- "count"
  
  # 3. Prevalence calculation -----
  sample_size <- nrow(totalyear_diagnosis_data)
  prevalence <- count_diagnoses_total[, 1:length(trait_list)]/sample_size #I removed the *100 to make it percent
  prevalence$year <- count_diagnoses_total$year
  prevalence$sample_size <- nrow(totalyear_diagnosis_data)
  prevalence$population <- population_name
  prevalence$measure <- "prevalence"
  
  #Get prevalence CI 
  prevalence_CI <- 1.96*sqrt((prevalence[, 1:length(trait_list)]*(1-prevalence[, 1:length(trait_list)]))/sample_size)
  prevalence_CI$year <- count_diagnoses_total$year
  prevalence_CI$sample_size <- nrow(totalyear_diagnosis_data)
  prevalence_CI$population <- population_name
  prevalence_CI$measure <- "prevalence_CI"
  
  descriptive_table <- rbind(count_diagnoses_total, prevalence, prevalence_CI)
  descriptive_table <- format(descriptive_table, scientific=F)
  return(descriptive_table)
  
}
