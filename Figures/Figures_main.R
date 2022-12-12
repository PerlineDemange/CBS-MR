############################################
## Project: EA_MH_CBS_MR 
## Script purpose: Main figure
## Author: Perline Demange 
##########################################
library(data.table)
library(ggplot2)



# Figure 2 prevalence #############





# Figure 3 ########
# Get CBS data 

results_all <- fread("C:/Users/user/Dropbox/CBS - MR/CBS_Ouput_210225/glm_sib_EA_results_20210224.csv")
results_all$OR <- as.numeric(results_all$OR)
results_CBS <- results_all[which(results_all$variable == "EA_within"| results_all$variable == "yrs"),]
results_CBS$OR_lci <- (results_CBS$OR - (results_CBS$OR.SE*1.96))
results_CBS$OR_uci <- (results_CBS$OR + (results_CBS$OR.SE*1.96))

results_CBS <- results_CBS[, c( "trait", "OR", "OR_lci", "OR_uci","variable", "p_value")]
results_CBS$result <- "CBS"

# All MR with sensitivity analysis
data_mr <- fread("../MR analyses/Results_EA_to_MH/all_MR_results_221020.csv")

data_mr$outcome[data_mr$outcome == "ANX"] <- "GAD"
data_mr$outcome[data_mr$outcome == "BIP"] <- "Bipolar"
data_mr$outcome[data_mr$outcome == "bip1"] <- "Bipolar_1"
data_mr$outcome[data_mr$outcome == "bip2"] <- "Bipolar_2"
data_mr$outcome[data_mr$outcome == "ANNO"] <- "Anorexia"
data_mr$outcome[data_mr$outcome == "ptsd"] <- "PTSD"
data_mr$outcome[data_mr$outcome == "scz"] <- "Schizophrenia"
data_mr$outcome[data_mr$outcome == "ocd"] <- "OCD"
data_mr$outcome[data_mr$outcome == "ASD"] <- "Autism"
data_mr$outcome[data_mr$outcome == "alc"] <- "Alcohol"


results.MR.EAtoMH<- data_mr[, c("outcome", "or", "or_lci95", "or_uci95", "method", "pval")]
names(results.MR.EAtoMH) <- c("trait", "OR", "OR_lci","OR_uci", "variable", "p_value")

results.MR.EAtoMH$result <- "EAtoMH"

# MR within fam 
data_mr_wf <- fread("../MR analyses/Results_EA_to_MH/all_MR_wf_results_221020.csv")
data_mr_wf$outcome[data_mr_wf$outcome == "ANX"] <- "GAD"
data_mr_wf$outcome[data_mr_wf$outcome == "BIP"] <- "Bipolar"
data_mr_wf$outcome[data_mr_wf$outcome == "bip1"] <- "Bipolar_1"
data_mr_wf$outcome[data_mr_wf$outcome == "bip2"] <- "Bipolar_2"
data_mr_wf$outcome[data_mr_wf$outcome == "ANNO"] <- "Anorexia"
data_mr_wf$outcome[data_mr_wf$outcome == "ptsd"] <- "PTSD"
data_mr_wf$outcome[data_mr_wf$outcome == "scz"] <- "Schizophrenia"
data_mr_wf$outcome[data_mr_wf$outcome == "ocd"] <- "OCD"
data_mr_wf$outcome[data_mr_wf$outcome == "ASD"] <- "Autism"
data_mr_wf$outcome[data_mr_wf$outcome == "alc"] <- "Alcohol"
data_mr_wf$method <- paste( data_mr_wf$method,"wf", sep="_")
results.MR.EAtoMH_wf<- data_mr_wf[, c("outcome", "or", "or_lci95", "or_uci95", "method", "pval")]
names(results.MR.EAtoMH_wf) <- c("trait", "OR", "OR_lci","OR_uci", "variable", "p_value")

results.MR.EAtoMH_wf$result <- "EAtoMH"
#keep only ivw
results.MR.EAtoMH_wf <- results.MR.EAtoMH_wf[results.MR.EAtoMH_wf$variable == "Inverse variance weighted_wf",]


# Combine 

results_MR <- rbind(results.MR.EAtoMH, results.MR.EAtoMH_wf)
results_plot <- rbind(results_MR, results_CBS)

results <- as.data.frame(results_plot)

#add empty lines 
results[nrow(results) + 1,] <- c("MDD", 0.95, 0.95, 0.95, "empty", 0, "EAtoMH")
results[nrow(results) + 1,] <- c("PTSD", 0.95, 0.95, 0.95, "empty", 0, "EAtoMH")
results[nrow(results) + 1,] <- c("Alcohol", 0.95, 0.95, 0.95, "empty", 0, "EAtoMH")
results[nrow(results) + 1,] <- c("ADHD", 0.95, 0.95, 0.95, "empty", 0, "EAtoMH")
results[nrow(results) + 1,] <- c("Anorexia", 0.95, 0.95, 0.95, "empty", 0, "EAtoMH")
results[nrow(results) + 1,] <- c("GAD", 0.95, 0.95, 0.95, "empty", 0, "EAtoMH")
results[nrow(results) + 1,] <- c("Schizophrenia", 0.95, 0.95, 0.95, "empty", 0, "EAtoMH")
results[nrow(results) + 1,] <- c("Bipolar", 0.95, 0.95, 0.95, "empty", 0, "EAtoMH")
results[nrow(results) + 1,] <- c("OCD", 0.95, 0.95, 0.95, "empty", 0, "EAtoMH")
results[nrow(results) + 1,] <- c("MDD", 0.95, 0.95, 0.95, "empty2", 0, "EAtoMH")
results[nrow(results) + 1,] <- c("PTSD", 0.95, 0.95, 0.95, "empty2", 0, "EAtoMH")
results[nrow(results) + 1,] <- c("Alcohol", 0.95, 0.95, 0.95, "empty2", 0, "EAtoMH")
results[nrow(results) + 1,] <- c("ADHD", 0.95, 0.95, 0.95, "empty2", 0, "EAtoMH")
results[nrow(results) + 1,] <- c("Anorexia", 0.95, 0.95, 0.95, "empty2", 0, "EAtoMH")
results[nrow(results) + 1,] <- c("GAD", 0.95, 0.95, 0.95, "empty2", 0, "EAtoMH")
results[nrow(results) + 1,] <- c("Schizophrenia", 0.95, 0.95, 0.95, "empty2", 0, "EAtoMH")
results[nrow(results) + 1,] <- c("Bipolar", 0.95, 0.95, 0.95, "empty2", 0, "EAtoMH")
results[nrow(results) + 1,] <- c("OCD", 0.95, 0.95, 0.95, "empty2", 0, "EAtoMH")
results$OR <- as.numeric(results$OR)
results$OR_lci <- as.numeric(results$OR_lci)
results$OR_uci <- as.numeric(results$OR_uci)
results$p_value <- as.numeric(results$p_value)

# Plot

results$variable <- factor(results$variable, levels = c("Inverse variance weighted_wf",
                                                        "Weighted median", "Weighted mode", "MR Egger","Inverse variance weighted",
                                                        "empty","empty2",
                                                        "EA_within","yrs" ), 
                           labels = c("Within-sibling IVW",
                                      "Weighted median", "Weighted mode", "MR Egger","Inverse variance weighted", 
                                      "empty","empty2",
                                      "Within-sibling effect","Population effect"))

order_traits <- c("MDD", "PTSD","Alcohol","ADHD","GAD",
                  "Schizophrenia", 
                  "Bipolar", "Anorexia","OCD")

results$trait <- factor(results$trait, levels = order_traits)
p<-ggplot(results, aes(x=trait, y=OR, ymin=ifelse(OR_lci>0.65, OR_lci , 0.65), ymax= ifelse(OR_uci < 1.3,OR_uci , 1.3), col=variable, fill=variable)) + 
  geom_linerange(position=position_dodge2(width=.7, preserve = "single"))+
  geom_hline(yintercept=1, lty=2)+
  geom_point(size=2,shape=21,stroke=0.5, position=position_dodge2(width=0.7, preserve = "single"))+
  scale_x_discrete(limits = rev(order_traits), name="Disorders")+ #
  scale_y_continuous(name= "Odds ratio", minor_breaks = seq(0.7,1.2, 0.1))+
  scale_fill_manual(values=c("black", "#be2596", 
                             "white","white",
                             "#96be25", "#43ae54", 
                             "#009773", "#007e7e",
                             "blue") , 
                    breaks = c("Population effect", "Within-sibling effect",
                               "empty","empty2",
                               "Inverse variance weighted", "MR Egger",
                               "Weighted mode","Weighted median", 
                               "Within-sibling IVW")) +
  scale_color_manual(values=c("black", "#be2596", 
                              "white","white",
                              "#96be25", "#43ae54",
                              "#009773", "#007e7e",
                              "blue"), 
                     breaks = c("Population effect", "Within-sibling effect",
                                "empty","empty2",
                                "Inverse variance weighted", "MR Egger",
                                "Weighted mode","Weighted median", 
                                "Within-sibling IVW"))+
  #scale_size_manual(values=c(2,5))+
  coord_flip()+ 
  ggtitle("CBS + MR MH to EA")+
  theme_minimal(base_size = 15)+
  theme(panel.grid.major=element_blank(), legend.title=element_blank(), axis.title.y=element_blank())
p

ggsave(p, file="Figure_main_sensitivity_wf_EAtoMH.tiff", width=11, height=9)




# Fig 4 MH to EA ###########


data_mr <- fread("../MR analyses/Results_MH_to_EA/all_MR_results_rev_221020.csv")


data_mr$exposure[data_mr$exposure == "ANX"] <- "GAD"
data_mr$exposure[data_mr$exposure == "BIP"] <- "Bipolar"
data_mr$exposure[data_mr$exposure == "bip1"] <- "Bipolar_1"
data_mr$exposure[data_mr$exposure == "bip2"] <- "Bipolar_2"
data_mr$exposure[data_mr$exposure == "ANNO"] <- "Anorexia"
data_mr$exposure[data_mr$exposure == "ptsd"] <- "PTSD"
data_mr$exposure[data_mr$exposure == "SCZ"] <- "Schizophrenia"
data_mr$exposure[data_mr$exposure == "ocd"] <- "OCD"
data_mr$exposure[data_mr$exposure == "ASD"] <- "Autism"
data_mr$exposure[data_mr$exposure == "alc"] <- "Alcohol"


results.MR.MHtoEA<- data_mr[, c("exposure", "or", "or_lci95", "or_uci95", "method", "pval")]
names(results.MR.MHtoEA) <- c("trait", "OR", "OR_lci","OR_uci", "variable", "p_value")

results.MR.MHtoEA$result <- "MHtoEA"


# MR within fam 
data_mr_wf <- fread("../MR analyses/Results_MH_to_EA/all_MR_wf_results_rev_221020.csv")
data_mr_wf$exposure[data_mr_wf$exposure == "ANX"] <- "GAD"
data_mr_wf$exposure[data_mr_wf$exposure == "BIP"] <- "Bipolar"
data_mr_wf$exposure[data_mr_wf$exposure == "bip1"] <- "Bipolar_1"
data_mr_wf$exposure[data_mr_wf$exposure == "bip2"] <- "Bipolar_2"
data_mr_wf$exposure[data_mr_wf$exposure == "ANNO"] <- "Anorexia"
data_mr_wf$exposure[data_mr_wf$exposure == "ptsd"] <- "PTSD"
data_mr_wf$exposure[data_mr_wf$exposure == "SCZ"] <- "Schizophrenia"
data_mr_wf$exposure[data_mr_wf$exposure == "ocd"] <- "OCD"
data_mr_wf$exposure[data_mr_wf$exposure == "ASD"] <- "Autism"
data_mr_wf$exposure[data_mr_wf$exposure == "alc"] <- "Alcohol"
data_mr_wf$method <- paste( data_mr_wf$method,"wf", sep="_")
results.MR.MHtoEA_wf<- data_mr_wf[, c("exposure", "or", "or_lci95", "or_uci95", "method", "pval")]
names(results.MR.MHtoEA_wf) <- c("trait", "OR", "OR_lci","OR_uci", "variable", "p_value")

results.MR.MHtoEA_wf$result <- "MHtoEA"
#keep only ivw
results.MR.MHtoEA_wf <- results.MR.MHtoEA_wf[results.MR.MHtoEA_wf$variable == "Inverse variance weighted_wf",]


# Combine 

results_MR <- rbind(results.MR.MHtoEA, results.MR.MHtoEA_wf)

results <- as.data.frame(results_MR)


# plot

results$variable <- factor(results$variable, levels = c("Inverse variance weighted_wf",
                                                        "Weighted median", "Weighted mode", "MR Egger","Inverse variance weighted" ), 
                           labels = c("Within-sibling IVW",
                                      "Weighted median", "Weighted mode", "MR Egger","Inverse variance weighted"))

order_traits <- c("MDD", "PTSD","Alcohol","ADHD","GAD",
                  "Schizophrenia", 
                  "Bipolar", "Anorexia","OCD")

results$trait <- factor(results$trait, levels = order_traits)
p<-ggplot(results, aes(x=trait, y=OR, ymin=ifelse(OR_lci>0.5, OR_lci , 0.5), 
                       ymax= ifelse(OR_uci < 2.3,OR_uci , 2.3), 
                       col=variable, fill=variable)) + 
  geom_linerange(position=position_dodge2(width=.7, preserve = "single"))+
  geom_hline(yintercept=1, lty=2)+
  geom_point(size=2,shape=21,stroke=0.5, position=position_dodge2(width=0.7, preserve = "single"))+
  scale_x_discrete(limits = rev(order_traits), name="Disorders")+ #
  scale_y_continuous(name= "Odds ratio", minor_breaks = seq(0.5,1.2, 0.1))+
  scale_fill_manual(values=c("#96be25", "#43ae54", 
                             "#009773", "#007e7e",
                             "blue") , 
                    breaks = c("Inverse variance weighted", "MR Egger",
                               "Weighted mode","Weighted median", 
                               "Within-sibling IVW")) +
  scale_color_manual(values=c("#96be25", "#43ae54",
                              "#009773", "#007e7e",
                              "blue"), 
                     breaks = c("Inverse variance weighted", "MR Egger",
                                "Weighted mode","Weighted median", 
                                "Within-sibling IVW"))+
  #scale_size_manual(values=c(2,5))+
  coord_flip()+ 
  ggtitle("MR MH to EA")+
  theme_minimal(base_size = 15)+
  theme(panel.grid.major=element_blank(), legend.title=element_blank(), axis.title.y=element_blank())
p
ggsave(p, file="Figure_main_sensitivity_wf_MHtoEA.tiff", width=11, height=9)


#
results <- results[results$variable == "Inverse variance weighted",]

order_traits <- c("MDD", "PTSD","Alcohol","ADHD","Anorexia","GAD",
                  "Schizophrenia", 
                  "Bipolar", "OCD")

results$trait <- factor(results$trait, levels = order_traits)
p<-ggplot(results, aes(x=trait, y=OR, ymin=ifelse(OR_lci>0.43, OR_lci , 0.43), 
                       ymax= ifelse(OR_uci < 1.3,OR_uci , 1.3), 
                       col=variable, fill=variable)) + 
  geom_linerange(position=position_dodge2(width=.7, preserve = "single"))+
  geom_hline(yintercept=1, lty=2)+
  geom_point(size=4,shape=21,stroke=0.5, position=position_dodge2(width=0.7, preserve = "single"))+
  scale_x_discrete(limits = rev(order_traits), name="Disorders")+ #
  scale_y_continuous(name= "Odds ratio", minor_breaks = seq(0.5,1.2, 0.1))+
  scale_fill_manual(values=c("#96be25") , 
                    breaks = c("Inverse variance weighted")) +
  scale_color_manual(values=c("#96be25"), 
                     breaks = c("Inverse variance weighted"))+
  #scale_size_manual(values=c(2,5))+
  coord_flip()+ 
  ggtitle("MR MH to EA")+
  theme_minimal(base_size = 15)+
  theme(panel.grid.major=element_blank(), legend.title=element_blank(), axis.title.y=element_blank())
p

# Figure 5 average education ######


compare <- fread("../CBS_Output_220908/comparison_family_members_results.csv")

patients <- compare[,c("Diagnoses", "Mean_patients", 
                       "SD_patients", "size_patients")]
names(patients) <- c("Diagnoses", "Mean", "SD", "sample_size")
patients$sample <- "patients"

sib_patients <- compare[,c("Diagnoses", "Mean_sib_patients", 
                           "SD_sib_patients", "size_sib_patients")]
names(sib_patients) <- c("Diagnoses", "Mean", "SD", "sample_size")
sib_patients$sample <- "sib_patients"

no_diagnoses <- compare[,c("Diagnoses", "Mean_no_diagnoses", 
                           "SD_no_diagnoses", "size_no_diagnoses")]
# The average is the same for all disorders as these are individuals without any disorders
# The average is 15.53849 and the SD is 2.65. The CI is so small it doesnt appear, so I just represent it as a line

results <- rbind(patients, sib_patients)
results$lci <-  results$Mean - 1.96*(results$SD/sqrt(results$sample_size))
results$uci <-  results$Mean + 1.96*(results$SD/sqrt(results$sample_size))
results <- results[!which(Diagnoses == "No_disorder"), ]

order_traits <- c("MDD", "PTSD","Alcohol","ADHD","GAD",
                  "Schizophrenia", 
                  "Bipolar","Anorexia", "OCD")

p<-ggplot(data=results, aes(x=Diagnoses, y=Mean,
                            ymin=lci,
                            ymax=uci, 
                            col=factor(sample), 
                            fill=factor(sample))) +
  geom_point(size=3,shape=21,colour="white", stroke=0.5,
             position=position_dodge2(width=0.5, preserve = "single"))+
  geom_linerange(position=position_dodge2(width=.5, preserve = "single"))+ 
  geom_hline(yintercept=15.53849, col="red")+ 
  theme_minimal(base_size = 15)+
  scale_x_discrete(limits = rev(order_traits), name="Disorders")+ 
  coord_flip()+ 
  scale_fill_manual(values=c("#006400", "#D2691E") , 
                    labels = c("Patients", "Healthy siblings of patients"))  +
  scale_color_manual(values=c("#006400", "#D2691E"), 
                     labels = c("Patients", "Healthy siblings of patients")) +
  
  theme(panel.grid.major=element_blank(),
        legend.title=element_blank(),
        axis.title.y=element_blank(), 
        axis.title.x=element_blank())+
  ggtitle("Comparison of the average number of years of education")
p

ggsave(p, file="figure_MAIN_comparison_patients_siblings.tiff", width=11, height=9)


# Move to other scripts ######
# Figure 1 A without MR sensitivity analyses for talks 

results.MR.EAtoMH <- results.MR.EAtoMH[results.MR.EAtoMH$variable == "Inverse variance weighted",]

results_plot <- rbind(results.MR.EAtoMH, results_CBS)

results <- as.data.frame(results_plot)

results$variable <- factor(results$variable, levels = c("Inverse variance weighted",
                                                        "empty","empty2",
                                                        "EA_within","yrs" ), 
                           labels = c("Inverse variance weighted", 
                                      "empty","empty2",
                                      "Within-sibling effect","Population effect"))

order_traits <- c("MDD", "PTSD","Alcohol","ADHD","Anorexia","GAD",
                  "Schizophrenia", 
                  "Bipolar", "OCD")

results$trait <- factor(results$trait, levels = order_traits)
p<-ggplot(results, aes(x=trait, y=OR, ymin=ifelse(OR_lci>0.65, OR_lci , 0.65), ymax= ifelse(OR_uci < 1.3,OR_uci , 1.3), col=variable, fill=variable)) + 
  geom_linerange(position=position_dodge2(width=.7, preserve = "single"))+
  geom_hline(yintercept=1, lty=2)+
  geom_point(size=2,shape=21,stroke=0.5, position=position_dodge2(width=0.7, preserve = "single"))+
  scale_x_discrete(limits = rev(order_traits), name="Disorders")+ #
  scale_y_continuous(name= "Odds ratio", minor_breaks = seq(0.7,1.2, 0.1))+
  scale_fill_manual(values=c("black", "#be2596", 
                             "white","white",
                             "#96be25") , 
                    breaks = c("Population effect", "Within-sibling effect",
                               "empty","empty2",
                               "Inverse variance weighted")) +
  scale_color_manual(values=c("black", "#be2596", 
                              "white","white",
                              "#96be25"), 
                     breaks = c("Population effect", "Within-sibling effect",
                                "empty","empty2",
                                "Inverse variance weighted"))+
  #scale_size_manual(values=c(2,5))+
  coord_flip()+ 
  ggtitle("CBS + MR MH to EA")+
  theme_minimal(base_size = 15)+
  theme(panel.grid.major=element_blank(), legend.title=element_blank(), axis.title.y=element_blank())
p

ggsave(p, file="Figure_main_EAtoMH_ECSR.tiff", width=11, height=9)


