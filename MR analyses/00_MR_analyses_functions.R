## Project: EA_MH_CBS_MR 
## Script purpose: Function to run all MR analyses with the TwoSampleMR package at once
## Author: Perline Demange 


run_MR_analyses <- function(data){ 

# Run pre-registered MR: #################

res  <- mr(data, method_list=c("mr_ivw", "mr_egger_regression", 
                               "mr_weighted_mode", "mr_weighted_median"))
res.or <- generate_odds_ratios(res)

# Run pleiotropy test and sensitivity variables : #################

pleio <- mr_pleiotropy_test(data)
pleio$isq <- Isq(data$beta.exposure, data$se.exposure)
F.stat <- data$beta.exposure^2/data$se.exposure^2
pleio$F.stat.mean <- mean(F.stat)
pleio$F.stat.min <- min(F.stat)

# We will report the Cochran Q’s-statistic SNP effect heterogeneity #####
heterog <- mr_heterogeneity(data, method_list=c( "mr_ivw","mr_egger_regression"))

result <- list(res.or, pleio, heterog)
return(result)
}

