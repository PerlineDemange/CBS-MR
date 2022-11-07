# Costs analyses
# Population Estimates and Within sibling Estimates 
# Author: Perline Demange
# Project: CBS-MR EA-MH 
# 

load("data_cost_sib.Rda")
summary(data_cost_sib$cost_mean)

#1. Create within and between variables ######
mean_of_sib_EA <- 
  group_by(data_cost_sib, FID) %>%
  summarize(EA_between = mean(yrs))

data_cost_sib <- merge(data_cost_sib, mean_of_sib_EA, by = "FID")

data_cost_sib$EA_within <- data_cost_sib$yrs - data_cost_sib$EA_between

#2. Analyses linear model ################

modelwithin  <- lm(cost_mean_log ~ EA_within + EA_between + 
       GBAGESLACHT + GBAGEBOORTEJAAR + birth_order,
     data = data_cost_sib)

results_within <- as.data.frame(summary(modelwithin)$coef)
results_within <- rownames_to_column(results_within, "Variables")
results_within$df <- modelwithin$df.residual
results_within$model <- "sib_comparison"

model <- lm(cost_mean_log ~ yrs + 
     GBAGESLACHT + GBAGEBOORTEJAAR + birth_order, 
   data = data_cost_sib)
results <- as.data.frame(summary(model)$coef)
results <- rownames_to_column(results, "Variables")
results$df <- model$df.residual
results$model <- "population"


results_all <- rbind(results, results_within)
write.csv(results_all, "lm_costs_results_20220722.csv", row.names=F, quote=F)
write.xlsx(results_all, "lm_costs_results_20220722.xlsx", 
           row.names=F)


# 3. Analyses Mixed models ####################

summary(lmer(cost_mean_log ~ EA_within + EA_between + 
               GBAGESLACHT + GBAGEBOORTEJAAR + birth_order +(1|FID),
             data = data_cost_sib))

modelwithin  <- lmer(cost_mean_log ~ EA_within + EA_between + 
                     GBAGESLACHT + GBAGEBOORTEJAAR + birth_order +(1|FID),
                   data = data_cost_sib)

results_within <- as.data.frame(summary(modelwithin)$coef)
results_within <- rownames_to_column(results_within, "Variables")
results_within$Z <- results_within$Estimate/results_within$`Std. Error`
results_within$Pvalue <- 2*pnorm(-abs(results_within$Z))
results_within$model <- "sib_comparison"

model <- lmer(cost_mean_log ~ yrs + 
              GBAGESLACHT + GBAGEBOORTEJAAR + birth_order +(1|FID), 
            data = data_cost_sib)
results <- as.data.frame(summary(model)$coef)
results <- rownames_to_column(results, "Variables")
results$Z <- results$Estimate/results$`Std. Error`
results$Pvalue <- 2*pnorm(-abs(results$Z))
results$model <- "population"


results_all <- rbind(results, results_within)
write.csv(results_all, "lmer_costs_results_20220722.csv", row.names=F, quote=F)
write.xlsx(results_all, "lmer_costs_results_20220722.xlsx", 
           row.names=F)
