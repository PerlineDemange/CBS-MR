# Script for sensitivity analyses in same-sex siblings only 
# Same-sex siblings with EA data were identified in script 02
# Author: Perline Demange
# Project: CBS-MR EA-MH 


# *Diagnoses* #############
# This part is adapted from script 03 

codefile <- as.data.frame(fread("H:/Data/Diagnoses/Diagnoses_code_preregistered_autism.csv", header=T))
trait_list <-  c(codefile$trait, "Any")

# 1. Descriptives #####################
## Women ############
### 1.1 prevalence in female sibling pop with EA ---------
gba_allyears_full_diagnosis <- fread("gba6585_diagnoses_allyears_full_diagnosis_20230703.csv")
load("edu_sib_women_220310.Rda")
sib_EA_women <- edu_sib_women
rm(edu_sib_women)
sib_EA_women_list <- sib_EA_women$RINPERSOON #508062

sib_EA_women_allyears_full_diagnosis <- gba_allyears_full_diagnosis[
  which(gba_allyears_full_diagnosis$RINPERSOON %in%
          sib_EA_women_list),]
#nrow should be  nrow * number of years
nrow(sib_EA_women_allyears_full_diagnosis)/6 #good 

### 1.2 Get total diagnoses per person --------
diagnoses_totalyears_per_person_sib_EA_women <- get_diagnoses_totalyear_per_person(
  sib_EA_women_allyears_full_diagnosis, 
  "siblings_EA")
#write.csv(diagnoses_totalyears_per_person_sib_EA_women, "siblings_women_EA_diagnoses_totalyear_per_person_20231108.csv", row.names = F)

### 1.3 Get descriptive table of diagnoses --------
descriptive_table_sib_EA_women <- get_descriptive_table(sib_EA_women_allyears_full_diagnosis, 
                                                        diagnoses_totalyears_per_person_sib_EA_women,
                                                        "siblings_EA_women", 
                                                        trait_list)
descriptive_table_sib_EA_women

## Men ######
### 1.1 prevalence in men sibling pop with EA ---------
load("edu_sib_men_220310.Rda")
sib_EA_men <- edu_sib_men
rm(edu_sib_men)
sib_EA_men_list <- sib_EA_men$RINPERSOON #516307

sib_EA_men_allyears_full_diagnosis <- gba_allyears_full_diagnosis[
  which(gba_allyears_full_diagnosis$RINPERSOON %in%
          sib_EA_men_list),]
#nrow should be  nrow * number of years
nrow(sib_EA_men_allyears_full_diagnosis)/6 #good 

### 1.2 Get total diagnoses per person --------
diagnoses_totalyears_per_person_sib_EA_men <- get_diagnoses_totalyear_per_person(
  sib_EA_men_allyears_full_diagnosis,
  "siblings_EA")
write.csv(diagnoses_totalyears_per_person_sib_EA_men,"siblings_men_EA_diagnoses_totalyear_per_person_20231108.csv", row.names = F)
#diagnoses_totalyears_per_person_sib_EA_men <- fread("siblings_men_EA_diagnoses_totalyear_per_person_20220314.csv")

### 1.3 Get descriptive table of diagnoses --------
descriptive_table_sib_EA_men <- get_descriptive_table(sib_EA_men_allyears_full_diagnosis,
                                                      diagnoses_totalyears_per_person_sib_EA_men,
                                                      "siblings_EA_men",
                                                      trait_list)
descriptive_table_sib_EA_men

## Save descriptives #####
descriptive_table <- rbind(descriptive_table_sib_EA_women,
                           descriptive_table_sib_EA_men)
descriptive_table
write.csv2(descriptive_table, 
          "descriptive_table_diagnoses_allpop_samesex_20231108.csv", 
          row.names = F, quote=F)
descriptive_table_cl <- fread("descriptive_table_diagnoses_allpop_samesex_20231108.csv") #this solves issues with numbers being characters quickly
write.csv2(descriptive_table_cl, 
           "descriptive_table_diagnoses_allpop_samesex_20231108.csv", 
           row.names = F, quote=F)

# make it output safe 
# remove CDD and Conduct 
descriptive_table_cl <- subset(descriptive_table_cl, select= -c(Conduct, ODD))
descriptive_table_cl[descriptive_table_cl$measure== 'count',]
# Combine Anorexia and Bulimia into eating disorders
descriptive_table_cl$Eating <- descriptive_table_cl$Anorexia + descriptive_table_cl$Bulimia
# remove anorexia and bulimia for men 
descriptive_table_cl[descriptive_table_cl$population == "siblings_EA_men",]$Anorexia<- NA
descriptive_table_cl[descriptive_table_cl$population == "siblings_EA_men",]$Bulimia<- NA
# remove eating for men in 2016. 
descriptive_table_cl[descriptive_table_cl$population == "siblings_EA_men" & descriptive_table_cl$year == 2016,]$Eating<- NA

write.csv2(descriptive_table_cl, 
           "descriptive_table_diagnoses_allpop_samesex_20231108_clean.csv", 
           row.names = F, quote=F)
write.xlsx(descriptive_table_cl, "descriptive_table_diagnoses_allpop_samesex_20231108_clean.xlsx", row.names = F)


# 2. Analyses ###################

sib_EA <- fread("sib_edu_subset_20210219.csv")

sib_EA_women_alldata <- merge(diagnoses_totalyears_per_person_sib_EA_women, sib_EA, by= "RINPERSOON")

## Women #######
### 2.1 Population estimates #######

results <- NULL
for(x in 1:(nrow(codefile)+1)){ #for each trait + Any
  trait <- trait_list[x]
  print(trait)
  model <- glm(paste(trait, "~ yrs + GBAGEBOORTEJAAR + birth_order"), 
               family = binomial, 
               data = sib_EA_women_alldata) #need paste to solve "variable lengths differ" error
  results_trait <- as.data.frame(summary(model)$coefficients)
  results_trait <- rownames_to_column(results_trait, "Variables")
  results_trait <- results_trait %>%
    mutate(OR = exp(Estimate), 
           var.diag= diag(vcov(model)),
           OR.SE = sqrt(OR^2 * var.diag))
  results_trait$df <- model$df.residual
  results_trait$trait <- paste(trait)
  robust <- coeftest(model, 
                     vcov = vcovCL, # CL is clustered errors, seems to be what we need here (other options are HC)
                     type = "HC0", cluster = sib_EA_women_alldata$FID) # HC1 is used for lm (degree of freedom based correction), HC0 is used for anything else
  results_trait_robust <- as.data.frame(robust[,])
  results_trait_robust <- results_trait_robust %>%
    mutate(OR = exp(Estimate), 
           var.diag= results_trait_robust[,2]* results_trait_robust[,2],  #diag(vcov(model)) is equivalent to (SE)^2 
           OR.SE = sqrt(OR^2 * var.diag))
  results_trait <- cbind(results_trait, results_trait_robust)
  results <- rbind(results, results_trait)
} 

colnames(results) <-c("variable", 
                                  "estimate", "SE", "z_value",
                                  "p_value","OR", "var.diag", "OR.SE",
                                  "df", "trait", 
                                  "estimate_robust", "SE_robust", "z_value_robust",
                                  "p_value_robust","OR_robust", "var.diag_robust", "OR.SE_robust"
)



### 2.2 Within-sibling analyses #########

mean_of_sib_EA <- 
  group_by(sib_EA_women_alldata, FID) %>%
  summarize(EA_between = mean(yrs))

final <- merge(sib_EA_women_alldata, mean_of_sib_EA, by = "FID")

final$EA_within <- final$yrs - final$EA_between

results_within <- NULL
for(x in 1:(nrow(codefile)+1)){ #for each trait + Any
  trait <- trait_list[x]
  print(trait)
  model1 <- glm(paste(trait, "~ EA_within + EA_between + GBAGEBOORTEJAAR + birth_order"), 
                family = binomial,
                data = final) #need paste to solve "variable lengths differ" error
  results_trait <- as.data.frame(summary(model1)$coefficients)
  results_trait <- rownames_to_column(results_trait, "Variables")
  results_trait <- results_trait %>%
    mutate(OR = exp(Estimate), 
           var.diag= diag(vcov(model1)),
           OR.SE = sqrt(OR^2 * var.diag))
  results_trait$df <- model1$df.residual
  results_trait$trait <- paste(trait)
  robust <- coeftest(model1, 
                     vcov = vcovCL, # CL is clustered errors, seems to be what we need here (other options are HC)
                     type = "HC0", cluster = final$FID) # HC1 is used for lm (degree of freedom based correction), HC0 is used for anything else
  results_trait_robust <- as.data.frame(robust[,])
  results_trait_robust <- results_trait_robust %>%
    mutate(OR = exp(Estimate), 
           var.diag= results_trait_robust[,2]* results_trait_robust[,2],  #diag(vcov(model)) is equivalent to (SE)^2 
           OR.SE = sqrt(OR^2 * var.diag))
  results_trait <- cbind(results_trait, results_trait_robust)
  results_within <- rbind(results_within, results_trait)
} 

head(results_within)

colnames(results_within) <- c("variable", 
                              "estimate", "SE", "z_value",
                              "p_value","OR", "var.diag", "OR.SE",
                              "df", "trait", 
                              "estimate_robust", "SE_robust", "z_value_robust",
                              "p_value_robust","OR_robust", "var.diag_robust", "OR.SE_robust"
)


results$model <- "population"
results_within$model <- "sib_comparison"

results_all <- rbind(results, results_within)
write.csv(results_all, "glm_sib_EA_women_results_20231108.csv")

## Men ####

sib_EA_men_alldata <- merge(diagnoses_totalyears_per_person_sib_EA_men, sib_EA, by= "RINPERSOON")

### 2.1 Population estimates #######

results <- NULL
for(x in 1:(nrow(codefile)+1)){ #for each trait + Any
  trait <- trait_list[x]
  print(trait)
  model <- glm(paste(trait, "~ yrs + GBAGEBOORTEJAAR + birth_order"), 
               family = binomial, 
               data = sib_EA_men_alldata) #need paste to solve "variable lengths differ" error
  results_trait <- as.data.frame(summary(model)$coefficients)
  results_trait <- rownames_to_column(results_trait, "Variables")
  results_trait <- results_trait %>%
    mutate(OR = exp(Estimate), 
           var.diag= diag(vcov(model)),
           OR.SE = sqrt(OR^2 * var.diag))
  results_trait$df <- model$df.residual
  results_trait$trait <- paste(trait)
  robust <- coeftest(model, 
                     vcov = vcovCL, # CL is clustered errors, seems to be what we need here (other options are HC)
                     type = "HC0", cluster = sib_EA_men_alldata$FID) # HC1 is used for lm (degree of freedom based correction), HC0 is used for anything else
  results_trait_robust <- as.data.frame(robust[,])
  results_trait_robust <- results_trait_robust %>%
    mutate(OR = exp(Estimate), 
           var.diag= results_trait_robust[,2]* results_trait_robust[,2],  #diag(vcov(model)) is equivalent to (SE)^2 
           OR.SE = sqrt(OR^2 * var.diag))
  results_trait <- cbind(results_trait, results_trait_robust)
  results <- rbind(results, results_trait)
} 
colnames(results) <-c("variable", 
                      "estimate", "SE", "z_value",
                      "p_value","OR", "var.diag", "OR.SE",
                      "df", "trait", 
                      "estimate_robust", "SE_robust", "z_value_robust",
                      "p_value_robust","OR_robust", "var.diag_robust", "OR.SE_robust"
)





### 2.2 Within-sibling analyses #########

mean_of_sib_EA <- 
  group_by(sib_EA_men_alldata, FID) %>%
  summarize(EA_between = mean(yrs))

final <- merge(sib_EA_men_alldata, mean_of_sib_EA, by = "FID")

final$EA_within <- final$yrs - final$EA_between

results_within <- NULL
for(x in 1:(nrow(codefile)+1)){ #for each trait + Any
  trait <- trait_list[x]
  print(trait)
  model1 <- glm(paste(trait, "~ EA_within + EA_between + GBAGEBOORTEJAAR + birth_order"), 
                family = binomial,
                data = final) #need paste to solve "variable lengths differ" error
  results_trait <- as.data.frame(summary(model1)$coefficients)
  results_trait <- rownames_to_column(results_trait, "Variables")
  results_trait <- results_trait %>%
    mutate(OR = exp(Estimate), 
           var.diag= diag(vcov(model1)),
           OR.SE = sqrt(OR^2 * var.diag))
  results_trait$df <- model1$df.residual
  results_trait$trait <- paste(trait)
  robust <- coeftest(model1, 
                     vcov = vcovCL, # CL is clustered errors, seems to be what we need here (other options are HC)
                     type = "HC0", cluster = final$FID) # HC1 is used for lm (degree of freedom based correction), HC0 is used for anything else
  results_trait_robust <- as.data.frame(robust[,])
  results_trait_robust <- results_trait_robust %>%
    mutate(OR = exp(Estimate), 
           var.diag= results_trait_robust[,2]* results_trait_robust[,2],  #diag(vcov(model)) is equivalent to (SE)^2 
           OR.SE = sqrt(OR^2 * var.diag))
  results_trait <- cbind(results_trait, results_trait_robust)
  results_within <- rbind(results_within, results_trait)
} 

colnames(results_within) <- c("variable", 
                              "estimate", "SE", "z_value",
                              "p_value","OR", "var.diag", "OR.SE",
                              "df", "trait", 
                              "estimate_robust", "SE_robust", "z_value_robust",
                              "p_value_robust","OR_robust", "var.diag_robust", "OR.SE_robust"
)


results$model <- "population"
results_within$model <- "sib_comparison"

results_all <- rbind(results, results_within)
write.csv(results_all, "glm_sib_EA_men_results_20231108.csv")

## Save ####
# save again as excel file? row,anmes = false and "false? 
results_women<- fread("glm_sib_EA_women_results_20231108.csv")
write.xlsx(results_women, "glm_sib_EA_women_results_20231108.xlsx", 
           row.names=F)


write.xlsx(results_men, "glm_sib_EA_men_results_20231108.xlsx", 
           row.names=F)

#3. Figure diagnoses same sex ######
results_women <-fread("glm_sib_EA_women_results_20231108.csv", drop=1)
results_women$sex <- "women"
results_women$OR <- as.numeric(results_women$OR)
results_men <-fread("glm_sib_EA_men_results_20231108.csv", drop=1)
results_men$sex <- "men"
results_all <- rbind(results_women, results_men)

results_plot2 <- results_all[which(results_all$variable == "EA_within" | results_all$variable == "EA_between"),]
dotCOLS= c("#a6d8f0", "#f9b282")
barCOLS= c("#008fd5", "#de6b35")
#results_plot2 <- results_plot2[results_plot2$trait != "Bulimia",]
p <- ggplot(results_plot2, aes(x=trait, y=OR, ymin=OR-OR.SE*1.96, ymax=OR+OR.SE*1.96, col= variable, fill=variable)) +
  geom_linerange(size= 5, position=position_dodge(width=0.5))+ 
  geom_hline( yintercept=1, lty=2)+
  geom_point(size=3, shape=21, colour="white", stroke=0.5, position=position_dodge(width=0.5)) +
  scale_fill_manual(values=barCOLS)+
  scale_color_manual(values=dotCOLS) + 
  scale_x_discrete(name= "Disorders") + 
  scale_y_continuous(name= "Odds ratio") +
  coord_flip()+ 
  facet_wrap(~ sex)+
  theme_minimal()

p


# look up other effects
results_age <- results_all[which(results_all$variable == "GBAGEBOORTEJAAR"),]
ggplot(data=results_age, aes(x=trait, y=estimate, fill=variable))+
  geom_bar(stat="identity", position="dodge") +   facet_wrap(~ sex)+
  geom_errorbar(aes(x=trait, ymin=estimate-SE*1.96, ymax=estimate+SE*1.96), position="dodge")

results_birth_order <- results_all[which(results_all$variable == "birth_order"),]
ggplot(data=results_birth_order, aes(x=trait, y=estimate, fill=variable))+
  geom_bar(stat="identity", position="dodge") +  facet_wrap(~ sex)+
  geom_errorbar(aes(x=trait, ymin=estimate-SE*1.96, ymax=estimate+SE*1.96), position="dodge")



# *Costs* ###############
load("sib_allyears_costs.Rda")

## 1. Create datasets ##########
data_cost_sib_women <- merge(sib_EA_women, sib_allyears, by = "RINPERSOON")
data_cost_sib_men <- merge(sib_EA_men, sib_allyears, by = "RINPERSOON")

# Get new number of family size and remove families < 2 
count <- as.data.frame(table(data_cost_sib_women$FID)) #get number of occurences of FID
data_cost_sib_women$Fsize_cost <- count[[2]][match(data_cost_sib_women$FID, count[[1]])] # match data fid with value in count, if true takes value of Freq 
data_cost_sib_women <- data_cost_sib_women[data_cost_sib_women$Fsize_cost > 1,] #497219

count <- as.data.frame(table(data_cost_sib_men$FID)) #get number of occurences of FID
data_cost_sib_men$Fsize_cost <- count[[2]][match(data_cost_sib_men$FID, count[[1]])] # match data fid with value in count, if true takes value of Freq 
data_cost_sib_men <- data_cost_sib_men[data_cost_sib_men$Fsize_cost > 1,] #492816

## 2. Analyses #############

### 2.1 Create within and between variables ##########
mean_of_sib_EA_women <- 
  group_by(data_cost_sib_women, FID) %>%
  summarize(EA_between = mean(yrs))

data_cost_sib_women <- merge(data_cost_sib_women, mean_of_sib_EA_women, by = "FID")

data_cost_sib_women$EA_within <- data_cost_sib_women$yrs - data_cost_sib_women$EA_between

mean_of_sib_EA_men <- 
  group_by(data_cost_sib_men, FID) %>%
  summarize(EA_between = mean(yrs))

data_cost_sib_men <- merge(data_cost_sib_men, mean_of_sib_EA_men, by = "FID")

data_cost_sib_men$EA_within <- data_cost_sib_men$yrs - data_cost_sib_men$EA_between

### 2.2 Run and save analyses ###############

#### 2.2.1 Analyses linear model ################

modelwithin_women  <- lm(cost_mean_log ~ EA_within + EA_between + 
                      GBAGEBOORTEJAAR + birth_order,
                   data = data_cost_sib_women)

results_within_women <- as.data.frame(summary(modelwithin_women)$coef)
results_within_women <- rownames_to_column(results_within_women, "Variables")
results_within_women$df <- modelwithin_women$df.residual
results_within_women$model <- "sib_comparison"
results_within_women$sex <- "women"

model_women <- lm(cost_mean_log ~ yrs + 
               GBAGEBOORTEJAAR + birth_order, 
            data = data_cost_sib_women)
results_women <- as.data.frame(summary(model_women)$coef)
results_women <- rownames_to_column(results_women, "Variables")
results_women$df <- model_women$df.residual
results_women$model <- "population"
results_women$sex <- "women"

results_all_women <- rbind(results_women, results_within_women)


modelwithin_men  <- lm(cost_mean_log ~ EA_within + EA_between + 
                           GBAGEBOORTEJAAR + birth_order,
                         data = data_cost_sib_men)

results_within_men <- as.data.frame(summary(modelwithin_men)$coef)
results_within_men <- rownames_to_column(results_within_men, "Variables")
results_within_men$df <- modelwithin_men$df.residual
results_within_men$model <- "sib_comparison"
results_within_men$sex <- "men"

model_men <- lm(cost_mean_log ~ yrs + 
                    GBAGEBOORTEJAAR + birth_order, 
                  data = data_cost_sib_men)
results_men <- as.data.frame(summary(model_men)$coef)
results_men <- rownames_to_column(results_men, "Variables")
results_men$df <- model_men$df.residual
results_men$model <- "population"
results_men$sex <- "men"

results_all_men <- rbind(results_men, results_within_men)

results_all<- rbind(results_all_women, results_all_men)


write.csv(results_all, "lm_costs_results_samesex_20220830.csv", 
          row.names=F, quote=F)
write.xlsx(results_all, "lm_costs_results_samesex_20220830.xlsx", 
           row.names=F)


#### 2.2.2 Analyses Mixed models ####################

modelwithin_women <- lmer(cost_mean_log ~ EA_within + EA_between + 
                        GBAGEBOORTEJAAR + birth_order +(1|FID),
                     data = data_cost_sib_women)

results_within_women<- as.data.frame(summary(modelwithin_women)$coef)
results_within_women<- rownames_to_column(results_within_women, "Variables")
results_within_women$Z <- results_within_women$Estimate/results_within_women$`Std. Error`
results_within_women$Pvalue <- 2*pnorm(-abs(results_within_women$Z))
results_within_women$model <- "sib_comparison"
results_within_women$sex <- "women"

model_women <- lmer(cost_mean_log ~ yrs + 
                 GBAGEBOORTEJAAR + birth_order +(1|FID), 
              data = data_cost_sib_women)
results_women <- as.data.frame(summary(model_women)$coef)
results_women <- rownames_to_column(results_women, "Variables")
results_women$Z <- results_women$Estimate/results_women$`Std. Error`
results_women$Pvalue <- 2*pnorm(-abs(results_women$Z))
results_women$model <- "population"
results_women$sex <- "women"

results_all_women <- rbind(results_women, results_within_women)

modelwithin_men <- lmer(cost_mean_log ~ EA_within + EA_between + 
                            GBAGEBOORTEJAAR + birth_order +(1|FID),
                          data = data_cost_sib_men)

results_within_men<- as.data.frame(summary(modelwithin_men)$coef)
results_within_men<- rownames_to_column(results_within_men, "Variables")
results_within_men$Z <- results_within_men$Estimate/results_within_men$`Std. Error`
results_within_men$Pvalue <- 2*pnorm(-abs(results_within_men$Z))
results_within_men$model <- "sib_comparison"
results_within_men$sex <- "men"

model_men <- lmer(cost_mean_log ~ yrs + 
                      GBAGEBOORTEJAAR + birth_order +(1|FID), 
                    data = data_cost_sib_men)
results_men <- as.data.frame(summary(model_men)$coef)
results_men <- rownames_to_column(results_men, "Variables")
results_men$Z <- results_men$Estimate/results_men$`Std. Error`
results_men$Pvalue <- 2*pnorm(-abs(results_men$Z))
results_men$model <- "population"
results_men$sex <- "men"

results_all_men <- rbind(results_men, results_within_men)

results_all <- rbind(results_all_women, results_all_men)


write.csv(results_all, "lmer_costs_results_samesex_20220830.csv", row.names=F, quote=F)
write.xlsx(results_all, "lmer_costs_results_samesex_20220830.xlsx", 
           row.names=F)

