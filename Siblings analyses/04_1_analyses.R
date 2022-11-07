# Script to run the diagnoses analyses 
# Author: Perline Demange
# Project: CBS-MR EA-MH 
#

# 1. Population Estimates ------
# As our diagnosis is binary, we will run a fixed effects logistic regression, 
# including sex and year of birth and birthorder as covariates 

## 1.1 Finalise sibling with EA data subset #############
diagnoses_totalyears_per_person_sib_EA <- fread("siblings_EA_diagnoses_totalyear_per_person_20210219.csv") 
sib_EA <- fread("sib_edu_subset_20210219.csv")

head(diagnoses_totalyears_per_person_sib_EA)
head(sib_EA)

sib_EA_alldata <- merge(diagnoses_totalyears_per_person_sib_EA, sib_EA, by= "RINPERSOON")

#write.csv(sib_EA_alldata, "sib_EA_diagnoses_20210221.csv", row.names=F)
sib_EA_alldata <- fread("sib_EA_diagnoses_20210221.csv")

## 1.2 Analyses ##########

codefile <- as.data.frame(fread("H:/Data/Diagnoses/Diagnoses_code_preregistered.csv", header=T))
trait_list <- codefile$trait
trait_list <- c(trait_list, "Any")

results <- NULL
for(x in 1:(nrow(codefile)+1)){ #for each trait + Any
  trait <- trait_list[x]
  print(trait)
  model <- glm(paste(trait, "~ yrs + GBAGESLACHT + GBAGEBOORTEJAAR + birth_order"),
               family = binomial,
               data = sib_EA_alldata) #need paste to solve "variable lengths differ" error
  results_trait <- as.data.frame(summary(model)$coefficients)
  results_trait <- rownames_to_column(results_trait, "Variables")
  results_trait <- results_trait %>%
    mutate(OR = exp(Estimate), 
           var.diag= diag(vcov(model)),
           OR.SE = sqrt(OR^2 * var.diag))
  results_trait$df <- model$df.residual
  results_trait$trait <- paste(trait)
  results <- rbind(results, results_trait)
} 
colnames(results) <- c("variable", "estimate", "SE", "z_value", "p_value","OR", 
                       "var.diag", "OR.SE", "df", "trait")


#write.csv(results, "glm_pop_sib_EA_results_2021024.csv")

ggplot(data=results[results$variable == "yrs",], aes(x=trait, y=estimate))+
  geom_bar(stat="identity") + 
  geom_errorbar(aes(x=trait, ymin=estimate-SE*1.96, ymax=estimate+SE*1.96))

 
# 2. Within sibling regression ####################

# ## 2.1 Get ICC for EA ###################
# 
# ICCest <- function(model){
#   icc <- sqrt(diag(getVarCov(model)))^2/(sqrt(diag(getVarCov(model)))^2 + model$sigma^2)
#   as.vector(icc)
# }
# 
# 
# icc <- lme(yrs~1,
#            random=~1|FID,
#            method="ML",
#            na.action=na.omit,
#            data=sib_EA_alldata)
# ICCest(icc)

## 2.1 create within and between ##################
mean_of_sib_EA <- 
  group_by(sib_EA_alldata, FID) %>%
  summarize(EA_between = mean(yrs))
    
final <- merge(sib_EA_alldata, mean_of_sib_EA, by = "FID")

final$EA_within <- final$yrs - final$EA_between

#write.csv(final, "final_within_between_data_20210222.csv", row.names=F)
final <- fread("final_within_between_data_20210222.csv")


## 2.2 Loop analyses ##################
results_within <- NULL
for(x in 1:(nrow(codefile)+1)){ #for each trait + Any
  trait <- trait_list[x]
  print(trait)
  model1 <- glm(paste(trait, "~ EA_within + EA_between + GBAGESLACHT + GBAGEBOORTEJAAR + birth_order"), 
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
  results_within <- rbind(results_within, results_trait)
} 

# Warning messages:
#   1: In doTryCatch(return(expr), name, parentenv, handler) :
#   restarting interrupted promise evaluation
# 2: In doTryCatch(return(expr), name, parentenv, handler) :
#   restarting interrupted promise evaluation
# 3: In doTryCatch(return(expr), name, parentenv, handler) :
#   restarting interrupted promise evaluation

head(results_within)
colnames(results_within) <- c("variable", "estimate", "SE", "z_value", "p_value",
                              "OR", "var.diag", "OR.SE","df", "trait")

results$model <- "population"
results_within$model <- "sib_comparison"

results_all <- rbind(results, results_within)
#write.csv(results_all, "glm_sib_EA_results_20210224.csv")



# 3. Figures ##### 
results_all <- fread("glm_sib_EA_results_20210224.csv")
results_all$OR <- as.numeric(results_all$OR)
results_plot <- results_all[which(results_all$variable == "yrs" | results_all$variable == "EA_within" | results_all$variable == "EA_between"),]
ggplot(data=results_plot, aes(x=trait, y=estimate, fill=variable))+
  geom_bar(stat="identity", position="dodge") + 
  geom_errorbar(aes(x=trait, ymin=estimate-SE*1.96, ymax=estimate+SE*1.96), position="dodge")


ggplot(data=results_plot, aes(x=trait, y=OR, fill=variable))+
  geom_bar(stat="identity", position="dodge")  +
  geom_errorbar(aes(x=trait, ymin=OR-OR.SE*1.96, ymax=OR+OR.SE*1.96), position="dodge")

results_plot2 <- results_all[which(results_all$variable == "EA_within" | results_all$variable == "EA_between"),]
dotCOLS= c("#a6d8f0", "#f9b282")
barCOLS= c("#008fd5", "#de6b35")
p <- ggplot(results_plot2, aes(x=trait, y=OR, ymin=OR-OR.SE*1.96, ymax=OR+OR.SE*1.96, col= variable, fill=variable)) +
  geom_linerange(size= 5, position=position_dodge(width=0.5))+ 
  geom_hline( yintercept=1, lty=2)+
  geom_point(size=3, shape=21, colour="white", stroke=0.5, position=position_dodge(width=0.5)) +
  scale_fill_manual(values=barCOLS)+
  scale_color_manual(values=dotCOLS) + 
  scale_x_discrete(name= "Disorders") + 
  scale_y_continuous(name= "Odds ratio") +
  coord_flip()+ 
  theme_minimal()

p

# look up other effects
results_sex <- results_all[which(results_all$variable == "GBAGESLACHT"),]
ggplot(data=results_sex, aes(x=trait, y=estimate, fill=variable))+
  geom_bar(stat="identity", position="dodge") + 
  geom_errorbar(aes(x=trait, ymin=estimate-SE*1.96, ymax=estimate+SE*1.96), position="dodge")

results_age <- results_all[which(results_all$variable == "GBAGEBOORTEJAAR"),]
ggplot(data=results_age, aes(x=trait, y=estimate, fill=variable))+
  geom_bar(stat="identity", position="dodge") + 
  geom_errorbar(aes(x=trait, ymin=estimate-SE*1.96, ymax=estimate+SE*1.96), position="dodge")

results_birth_order <- results_all[which(results_all$variable == "birth_order"),]
ggplot(data=results_birth_order, aes(x=trait, y=estimate, fill=variable))+
  geom_bar(stat="identity", position="dodge") + 
  geom_errorbar(aes(x=trait, ymin=estimate-SE*1.96, ymax=estimate+SE*1.96), position="dodge")


-----------------------------------------
# 4. Robustness analyses #########################

## 4.1 Remove yrs= 11 ##########################

sib_EA_alldata_exc <- sib_EA_alldata[!which(sib_EA_alldata$yrs == 11), ] #1736756 #sib_all_data 1743032
# exclude families that are now incomplete
count <- as.data.frame(table(sib_EA_alldata_exc$FID)) #get number of occurences of FID
sib_EA_alldata_exc$Fsize_excl11 <- count[[2]][match(sib_EA_alldata_exc$FID, count[[1]])]  # match data fid with value in count, if true takes value of Freq 
sib_EA_alldata_exc <- sib_EA_alldata_exc[sib_EA_alldata_exc$Fsize_excl11  > 1]  #1732499

#get between and within
mean_of_sib_EA_exc <- 
  group_by(sib_EA_alldata_exc, FID) %>%
  summarize(EA_between = mean(yrs))

final_exc <- merge(sib_EA_alldata_exc, mean_of_sib_EA_exc, by = "FID")

final_exc$EA_within <- final_exc$yrs - final_exc$EA_between

write.csv(final_exc, "final_within_between_data_excl11_20220715.csv", row.names=F)


results_within_exc <- NULL
for(x in 1:(nrow(codefile)+1)){ #for each trait + Any
  trait <- trait_list[x]
  print(trait)
  model1 <- glm(paste(trait, "~ EA_within + EA_between + GBAGESLACHT + GBAGEBOORTEJAAR + birth_order"),
                family = binomial,
                data = final_exc) #need paste to solve "variable lengths differ" error
  results_trait <- as.data.frame(summary(model1)$coefficients)
  results_trait <- rownames_to_column(results_trait, "Variables")
  results_trait <- results_trait %>%
    mutate(OR = exp(Estimate), 
           var.diag= diag(vcov(model1)),
           OR.SE = sqrt(OR^2 * var.diag))
  results_trait$df <- model1$df.residual
  results_trait$trait <- paste(trait)
  results_within_exc <- rbind(results_within_exc, results_trait)
} 

results_within_exc
colnames(results_within_exc) <- c("variable", "estimate", "SE", "z_value", "p_value","OR", "var.diag", "OR.SE","df", "trait")



results <- NULL
for(x in 1:(nrow(codefile)+1)){ #for each trait + Any
  trait <- trait_list[x]
  print(trait)
  model <- glm(paste(trait, "~ yrs + GBAGESLACHT + GBAGEBOORTEJAAR + birth_order"),
               family = binomial,
               data = final_exc) #need paste to solve "variable lengths differ" error
  results_trait <- as.data.frame(summary(model)$coefficients)
  results_trait <- rownames_to_column(results_trait, "Variables")
  results_trait <- results_trait %>%
    mutate(OR = exp(Estimate), 
           var.diag= diag(vcov(model)),
           OR.SE = sqrt(OR^2 * var.diag))
  results_trait$df <- model$df.residual
  results_trait$trait <- paste(trait)
  results <- rbind(results, results_trait)
} 
colnames(results) <- c("variable", "estimate", "SE", "z_value", "p_value","OR", 
                       "var.diag", "OR.SE", "df", "trait")



results$model <- "population"
results_within_exc$model <- "sib_comparison"

results_all_exc <- rbind(results, results_within_exc)
write.csv(results_all_exc, "glm_sib_EA_results_excl11_20220715.csv")


#graph
results_plot_exc <- results_within_exc[which(results_within_exc$variable == "EA_within" | results_within_exc$variable == "EA_between"),]
dotCOLS= c("#a6d8f0", "#f9b282")
barCOLS= c("#008fd5", "#de6b35")
p_exc <- ggplot(results_plot_exc, aes(x=trait, y=OR, ymin=OR-OR.SE*1.96, ymax=OR+OR.SE*1.96, col= variable, fill=variable)) +
  geom_linerange(size= 5, position=position_dodge(width=0.5))+ 
  geom_hline( yintercept=1, lty=2)+
  geom_point(size=3, shape=21, colour="white", stroke=0.5, position=position_dodge(width=0.5)) +
  scale_fill_manual(values=barCOLS)+
  scale_color_manual(values=dotCOLS) + 
  scale_x_discrete(name= "Disorders") + 
  scale_y_continuous(name= "Odds ratio - yrs 11 excluded") +
  coord_flip()+ 
  theme_minimal()

p_exc
p



## 4.2. Remove yrs= 2 ##########################
# it is in theory not possible to not finish primary school in the NL 

sib_EA_alldata_exc2 <- sib_EA_alldata[!which(sib_EA_alldata$yrs == 2), ] #1723367 #sib_all_data 1743032
# exclude families that are now incomplete
count <- as.data.frame(table(sib_EA_alldata_exc2$FID)) #get number of occurences of FID
sib_EA_alldata_exc2$Fsize_excl2 <- count[[2]][match(sib_EA_alldata_exc2$FID, count[[1]])]  # match data fid with value in count, if true takes value of Freq 
sib_EA_alldata_exc2 <- sib_EA_alldata_exc2[sib_EA_alldata_exc2$Fsize_excl2  > 1]  #1712334

#get between and within
mean_of_sib_EA_exc2 <- 
  group_by(sib_EA_alldata_exc2, FID) %>%
  summarize(EA_between = mean(yrs))

final_exc2 <- merge(sib_EA_alldata_exc2, mean_of_sib_EA_exc2, by = "FID")

final_exc2$EA_within <- final_exc2$yrs - final_exc2$EA_between

write.csv(final_exc2, "final_within_between_data_excl2_20220715.csv", row.names=F)


results_within_exc2 <- NULL
for(x in 1:(nrow(codefile)+1)){ #for each trait + Any
  trait <- trait_list[x]
  print(trait)
  model1 <- glm(paste(trait, "~ EA_within + EA_between + GBAGESLACHT + GBAGEBOORTEJAAR + birth_order"),
                family = binomial,
                data = final_exc2) #need paste to solve "variable lengths differ" error
  results_trait <- as.data.frame(summary(model1)$coefficients)
  results_trait <- rownames_to_column(results_trait, "Variables")
  results_trait <- results_trait %>%
    mutate(OR = exp(Estimate), 
           var.diag= diag(vcov(model1)),
           OR.SE = sqrt(OR^2 * var.diag))
  results_trait$df <- model1$df.residual
  results_trait$trait <- paste(trait)
  results_within_exc2 <- rbind(results_within_exc2, results_trait)
} 

results_within_exc2
colnames(results_within_exc2) <- c("variable", "estimate", "SE", "z_value", "p_value","OR", "var.diag", "OR.SE","df", "trait")



results2 <- NULL
for(x in 1:(nrow(codefile)+1)){ #for each trait + Any
  trait <- trait_list[x]
  print(trait)
  model <- glm(paste(trait, "~ yrs + GBAGESLACHT + GBAGEBOORTEJAAR + birth_order"),
               family = binomial,
               data = final_exc2) #need paste to solve "variable lengths differ" error
  results_trait <- as.data.frame(summary(model)$coefficients)
  results_trait <- rownames_to_column(results_trait, "Variables")
  results_trait <- results_trait %>%
    mutate(OR = exp(Estimate), 
           var.diag= diag(vcov(model)),
           OR.SE = sqrt(OR^2 * var.diag))
  results_trait$df <- model$df.residual
  results_trait$trait <- paste(trait)
  results2 <- rbind(results2, results_trait)
} 
colnames(results2) <- c("variable", "estimate", "SE", "z_value", "p_value","OR", 
                       "var.diag", "OR.SE", "df", "trait")



results2$model <- "population"
results_within_exc2$model <- "sib_comparison"

results_all_exc2 <- rbind(results2, results_within_exc2)
write.csv(results_all_exc2, "glm_sib_EA_results_excl2_20220715.csv", row.names=F)


#graph
results_plot_exc <- results_within_exc2[which(results_within_exc2$variable == "EA_within" | results_within_exc2$variable == "EA_between"),]
dotCOLS= c("#a6d8f0", "#f9b282")
barCOLS= c("#008fd5", "#de6b35")
p_exc <- ggplot(results_plot_exc, aes(x=trait, y=OR, ymin=OR-OR.SE*1.96, ymax=OR+OR.SE*1.96, col= variable, fill=variable)) +
  geom_linerange(size= 5, position=position_dodge(width=0.5))+ 
  geom_hline( yintercept=1, lty=2)+
  geom_point(size=3, shape=21, colour="white", stroke=0.5, position=position_dodge(width=0.5)) +
  scale_fill_manual(values=barCOLS)+
  scale_color_manual(values=dotCOLS) + 
  scale_x_discrete(name= "Disorders") + 
  scale_y_continuous(name= "Odds ratio - yrs 2 excluded") +
  coord_flip()+ 
  theme_minimal()

p_exc
p


## 4.3 Quick check Remove if EA within = 0 ###########
# this removes non discordant families 
# but could also remove individual with more than 1 sib, which EA averages to the one sib EA 
final_discordant <- final[!which(final$EA_within == 0),] # 1330628
summary(final_discordant$EA_within)

summary(glm(MDD ~ EA_within + EA_between + GBAGESLACHT + GBAGEBOORTEJAAR, family= binomial, data=final_discordant))
summary(glm(MDD ~ EA_within + EA_between + GBAGESLACHT + GBAGEBOORTEJAAR, family= binomial, data=final))

