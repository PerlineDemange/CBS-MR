## Project: EA_MH_CBS_MR 
## Script purpose: Main figures
## Author: Perline Demange 

library(data.table)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(corrplot)
source("functions_plot.R")

# when minor grid coincide with major grid that has been disabled it is not appearing 
# workaround, to do everytime 
trace(ggplot2:::guide_grid, edit=TRUE)
# I manually remove 
# x.minor <- setdiff(x.minor, x.major)
# y.minor <- setdiff(y.minor, y.major)
#can re-edit with untrace(ggplot2:::guide_grid)

# Figure 1 #############
order_traits <- c("MDD", "PTSD","Alcohol","ADHD","GAD",
                  "Schizophrenia", "OCD", "ASD",
                  "Bipolar", "Eating")


#polychoric correlations
poly<- fread("../CBS_Output_231113/matrix_polychoriccorrelation_sibEA_20230711.csv")
poly <- as.matrix(poly, rownames=1)
poly_short <- poly[order_traits, order_traits]

#I forgot to remove the value of Schizophrenia x Eating which has a cell count lower than 10
poly_short[6,10] <- NA
poly_short[10,6] <- NA

#genetic correlations 
load("../MR Analyses/LDSCoutput_matrix_psychdisorders_EA.RData")
LDSCoutput$S
rG_traits <- c("MDD_Howard", "PTSD_Nievergelt", "Alc_Walters",
               "ADHD_Demontis_2023", "Anx_Purves_meta", "SCZ_PGC3", 
               "OCD_Arnold", "ASD_Grove", "BIP_all_Mullins", "Ano_Watson")
gen <- as.matrix(cov2cor(LDSCoutput$S))
rownames(gen) <- colnames(gen) 
gen <- gen[rG_traits, rG_traits]
gen["Ano_Watson",] <- -gen["Ano_Watson",] # as we saw Anorexia alleles were flipped 
gen[,"Ano_Watson"] <- -gen[,"Ano_Watson"]


#get pvalue
k<-nrow(LDSCoutput$S)
SE<-matrix(0, k, k)
SE[lower.tri(SE,diag=TRUE)] <-sqrt(diag(LDSCoutput$V))
pvalue <- as.matrix(2*pnorm((abs(LDSCoutput$S)/SE), lower.tail=F))
rownames(pvalue) <- colnames(pvalue) 
#pvalue <- pvalue[rG_traits, rG_traits] # this line is wrong! because the matrix is not symmetric, this shuffles correct values and 0.... 
#need to render the matrix symmetric to make it subsettable
lowertri <- t(pvalue) # need this step because r populate matrix col by col  https://davetang.org/muse/2014/01/22/making-symmetric-matrices-in-r/
pvalue[upper.tri(pvalue)] <- lowertri[upper.tri(lowertri)]
pvalue <- pvalue[rG_traits, rG_traits] #now this should be correct 

#combined cor and rg in one figure 

combined <- poly_short
combined[lower.tri(combined)] <- gen[lower.tri(gen)]

corrplot(combined, method="square",is.corr=F, 
         na.label = "NA",na.label.col = "grey",
         diag=F, col= COL1("YlOrBr"))

pdf(file = "poly_rG_withoutEA.pdf")
corrplot(poly_short, method="color",is.corr=F, 
         na.label = "NA",na.label.col = "grey",
         diag=F, col= COL2("RdBu"),
         type="upper", 
         outline=T, tl.pos = "n", tl.col = "indianred4", 
         addgrid.col = "darkgray", mar = c(8,10,0,0), 
         addCoef.col = 'grey36')
corrplot(gen, method="color",is.corr=F, 
         na.label = "NA",na.label.col = "grey",
         diag=F, col= COL2("PRGn"),
         type="lower", add=T, 
         outline=T, tl.pos = "l", tl.col = "indianred4", 
         addgrid.col = "darkgray", mar = c(0,0,0,0), 
         addCoef.col = 'grey36', 
         p.mat = pvalue, sig.level = 0.001, 
         insig = 'pch', pch.col = "grey", pch.cex=5)
dev.off()

# Figure 2 prevalence #############

diagnoses_by_yrs_sex <- fread("../CBS_Output_231113/diagnoses_by_yrs_sex_20230707_output.csv") # of sibEA
trait_list <- colnames(diagnoses_by_yrs_sex[,5:24]) 

# High prevalence group  

high_prev <- prevalence_plot_third(c("MDD", "PTSD", "Alcohol"), 10.5)

tiff("allplots_highprev.tiff", width=15000, height=4500, res=900)
do.call(grid.arrange, c(high_prev, nrow=1))
dev.off()

# Medium prevalence group  

mid_prev <- prevalence_plot_third(c( "ADHD", "Schizophrenia","ASD"), 5.20)

tiff("allplots_midprev.tiff", width=15000, height=4500, res=900)
do.call(grid.arrange, c(mid_prev, nrow=1))
dev.off()


# Low prevalence 

low_prev <- prevalence_plot_third(c("GAD", "OCD","Bipolar"), 1.5)
tiff("allplots_lowprev.tiff", width=15000, height=4500, res=900)
do.call(grid.arrange, c(low_prev, nrow=1))
dev.off()

# eating
eating <- prevalence_plot_third(c( "Eating","OCD","Bipolar"), 1.5)
tiff("allplots_eating.tiff", width=15000, height=4500, res=900)
do.call(grid.arrange, c(eating, nrow=1))
dev.off()





# Figure 3 ########
# Get CBS data 

results_all <- fread("../CBS_Output_231113/glm_sib_EA_results_20230922.csv")
results_all$OR_robust <- as.numeric(results_all$OR_robust)
results_CBS <- results_all[which(results_all$variable == "EA_within"|
                                   results_all$variable == "yrs"),]
results_CBS$OR_lci_robust <- (results_CBS$OR_robust - (results_CBS$OR.SE_robust*1.96))
results_CBS$OR_uci_robust <- (results_CBS$OR_robust + (results_CBS$OR.SE_robust*1.96))

results_CBS <- results_CBS[, c( "trait", "OR_robust", "OR_lci_robust",
                                "OR_uci_robust","variable", "p_value_robust")]
names(results_CBS ) <- c("trait", "OR", "OR_lci","OR_uci",
                                 "variable", "p_value")
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
data_mr$outcome[data_mr$outcome == "ASD"] <- "ASD"
data_mr$outcome[data_mr$outcome == "alc"] <- "Alcohol"
# to plot adhd 2023 only
data_mr <- data_mr[!(outcome == "ADHD"),] 
data_mr$outcome[data_mr$outcome == "ADHD2023"] <- "ADHD"

results.MR.EAtoMH<- data_mr[, c("outcome", "or", "or_lci95",
                                "or_uci95", "method", "pval")]
names(results.MR.EAtoMH) <- c("trait", "OR", "OR_lci","OR_uci",
                              "variable", "p_value")

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
data_mr_wf$outcome[data_mr_wf$outcome == "ASD"] <- "ASD"
data_mr_wf$outcome[data_mr_wf$outcome == "alc"] <- "Alcohol"
# to plot adhd 2023 only
data_mr_wf <- data_mr_wf[!(outcome == "ADHD"),] 
data_mr_wf$outcome[data_mr_wf$outcome == "ADHD2023"] <- "ADHD"


data_mr_wf$method <- paste( data_mr_wf$method,"wf", sep="_")
results.MR.EAtoMH_wf<- data_mr_wf[, c("outcome", "or", "or_lci95",
                                      "or_uci95", "method", "pval")]
names(results.MR.EAtoMH_wf) <- c("trait", "OR", "OR_lci","OR_uci",
                                 "variable", "p_value")

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
results[nrow(results) + 1,] <- c("ASD", 0.95, 0.95, 0.95, "empty2", 0, "EAtoMH")
results$OR <- as.numeric(results$OR)
results$OR_lci <- as.numeric(results$OR_lci)
results$OR_uci <- as.numeric(results$OR_uci)
results$p_value <- as.numeric(results$p_value)

# Plot

results$variable <- factor(results$variable,
                           levels = c("Inverse variance weighted_wf",
                                      "Weighted median", "Weighted mode",
                                      "MR Egger","Inverse variance weighted",
                                      "empty","empty2",
                                      "EA_within","yrs" ), 
                           labels = c("Within-sibling IVW",
                                      "Weighted median", "Weighted mode",
                                      "MR Egger","Inverse variance weighted", 
                                      "empty","empty2",
                                      "Within-sibling association","Observational association"))

order_traits <- c("MDD", "PTSD","Alcohol","ADHD","GAD",
                  "Schizophrenia", "OCD", "ASD",
                  "Bipolar", "Anorexia")

results$trait <- factor(results$trait, levels = order_traits)
p<-ggplot() + 
  geom_linerange(data=results, aes(x=trait, y=OR,
                                   ymin=ifelse(OR_lci>0.65, OR_lci , 0.65),
                                   ymax= ifelse(OR_uci < 1.3,OR_uci , 1.3),
                                   col=variable),
                 position=position_dodge2(width=.7, preserve = "single"))+
  geom_hline(yintercept=1, lty=2)+
  geom_point(data=results,aes(x=trait, y=OR,
                              col=variable, fill=variable),
             size=3,shape=21,stroke=0.5,
             position=position_dodge2(width=0.7, preserve = "single"))+
  coord_flip()+
  scale_x_discrete(limits = rev(order_traits), name="Disorders")+ #
  scale_y_continuous(name= "Odds ratio", minor_breaks = seq(0.7,1.2, 0.1))+
  scale_fill_manual(values=c("black", "#be2596", 
                             "white","white",
                             "#96be25", "#43ae54", 
                             "#009773", "#007e7e",
                             "blue") , 
                    breaks = c("Observational association", "Within-sibling association",
                               "empty","empty2",
                               "Inverse variance weighted", "MR Egger",
                               "Weighted mode","Weighted median", 
                               "Within-sibling IVW")) +
  scale_color_manual(values=c("black", "#be2596", 
                              "white","white",
                              "#96be25", "#43ae54",
                              "#009773", "#007e7e",
                              "blue"), 
                     breaks = c("Observational association", "Within-sibling association",
                                "empty","empty2",
                                "Inverse variance weighted", "MR Egger",
                                "Weighted mode","Weighted median", 
                                "Within-sibling IVW"))+
  #scale_size_manual(values=c(2,5))+
  coord_flip()+ 
  theme_minimal(base_size = 22)+
  theme(panel.grid.major=element_blank(),
        legend.title=element_blank(),
        axis.title.y=element_blank(), 
        rect = element_rect(fill = "transparent"),
        legend.text = element_text(
          margin = margin(l=5, t=5, b=5, unit = "pt")))


p <- p + geom_rect(aes(ymin = 0.66, ymax = 1.3, xmin = c(1.5, 3.5, 5.5, 7.5, 9.5),
                  xmax = c(2.5, 4.5, 6.5, 8.5, 10.5)), fill = "gray", alpha = 0.2)

p <- p + guides(color = FALSE, fill = FALSE) # to remove both color and fill legends

p
ggsave(p,
       file="Figure_main_sensitivity_wf_EAtoMH.tiff",
       width=12, height=12, 
       bg = "transparent")
ggsave(p,
       file="Figure_main_sensitivity_wf_EAtoMH.png",
       width=12, height=12, 
       bg = "transparent")
ggsave(p,
       file="Figure_main_sensitivity_wf_EAtoMH.pdf",
       width=12, height=12, 
       bg = "transparent")





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
data_mr$exposure[data_mr$exposure == "ASD"] <- "ASD"
data_mr$exposure[data_mr$exposure == "alc"] <- "Alcohol"
# to plot adhd 2023 only
data_mr <- data_mr[!(exposure == "ADHD"),] 
data_mr$exposure[data_mr$exposure == "ADHD2023"] <- "ADHD"

results.MR.MHtoEA<- data_mr[, c("exposure", "b", "lo_ci",
                                "up_ci", "method", "pval")]
names(results.MR.MHtoEA) <- c("trait", "Beta", "Beta_lci","Beta_uci",
                              "variable", "p_value")

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
data_mr_wf$exposure[data_mr_wf$exposure == "ASD"] <- "ASD"
data_mr_wf$exposure[data_mr_wf$exposure == "alc"] <- "Alcohol"
# to plot adhd 2023 only
data_mr_wf <- data_mr_wf[!(exposure == "ADHD"),] 
data_mr_wf$exposure[data_mr_wf$exposure == "ADHD2023"] <- "ADHD"

data_mr_wf$method <- paste( data_mr_wf$method,"wf", sep="_")
results.MR.MHtoEA_wf<- data_mr_wf[,
                                  c("exposure","b", "lo_ci",
                                    "up_ci", "method", "pval")]
names(results.MR.MHtoEA_wf) <- c("trait", "Beta", "Beta_lci","Beta_uci",
                                 "variable", "p_value")

results.MR.MHtoEA_wf$result <- "MHtoEA"
#keep only ivw
results.MR.MHtoEA_wf <- results.MR.MHtoEA_wf[results.MR.MHtoEA_wf$variable == "Inverse variance weighted_wf",]


# Combine 

results_MR <- rbind(results.MR.MHtoEA, results.MR.MHtoEA_wf)

results <- as.data.frame(results_MR)


# plot

results$variable <- factor(results$variable,
                           levels = c("Inverse variance weighted_wf",
                                      "Weighted median", "Weighted mode",
                                      "MR Egger","Inverse variance weighted" ), 
                           labels = c("Within-sibling IVW",
                                      "Weighted median", "Weighted mode",
                                      "MR Egger","Inverse variance weighted"))

order_traits <- c("MDD", "PTSD","Alcohol","ADHD","GAD",
                  "Schizophrenia", 
                  "OCD", "ASD", "Bipolar", "Anorexia")

results$trait <- factor(results$trait, levels = order_traits)
p<-ggplot() + 
  geom_linerange(data=results, aes(x=trait, y=Beta, ymin= Beta_lci, 
                              ymax= Beta_uci , 
                              col=variable, fill=variable), 
                 position=position_dodge2(width=.7, preserve = "single"))+
  geom_hline(yintercept=0, lty=2)+
  geom_point(data=results, aes(x=trait, y=Beta, ymin= Beta_lci, 
                          ymax= Beta_uci , 
                          col=variable, fill=variable),
             size=3,shape=21,stroke=0.5,
             position=position_dodge2(width=0.7, preserve = "single"))+
  scale_x_discrete(limits = rev(order_traits), name="Disorders")+ #
  scale_y_continuous(name= "Beta", minor_breaks = seq(-0.5,0.3, 0.1))+
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
  coord_flip(ylim=c(-0.55,0.3))+ 
  #ggtitle("MR MH to EA")+
  theme_minimal(base_size = 22)+
  theme(panel.grid.major=element_blank(),
        legend.title=element_blank(),
        axis.title.y=element_blank(), 
        rect = element_rect(fill = "transparent"),
        legend.text = element_text(
          margin = margin(l=5, t=5, b=5, unit = "pt")))
p

p <- p + geom_rect(aes(ymin = -0.6, ymax = 0.345, xmin = c(1.5, 3.5, 5.5, 7.5, 9.5),
                       xmax = c(2.5, 4.5, 6.5, 8.5, 10.5)), fill = "gray", alpha = 0.2)

p <- p + guides(color = FALSE, fill = FALSE) # to remove both color and fill legends
p <- p +
  annotate("segment", x = 1.15, xend = 1.15, y = 0.3, yend = 0.34, color = "#43ae54", size = 0.2,
           arrow = arrow(length = unit(0.08, "inches"), type = "closed"))



ggsave(p, file="Figure_main_sensitivity_wf_MHtoEA.tiff",
       width=12, height=12, bg = "transparent")
ggsave(p, file="Figure_main_sensitivity_wf_MHtoEA.png",
       width=12, height=12, 
       bg = "transparent")

ggsave(p, file="Figure_4.pdf",
       width=12, height=12, 
       bg = "transparent")


# Figure 5 average education ######


compare <- fread("../CBS_Output_231113/comparison_family_members_results_20231013.csv")

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
# The average is 15.5402 and the SD is 2.65. The CI is so small it doesnt appear, so I just represent it as a line

results <- rbind(patients, sib_patients)
results$lci <-  results$Mean - 1.96*(results$SD/sqrt(results$sample_size))
results$uci <-  results$Mean + 1.96*(results$SD/sqrt(results$sample_size))
results <- results[!which(Diagnoses == "No_disorder"), ]

order_traits <- c("MDD", "PTSD","Alcohol","ADHD","GAD",
                  "Schizophrenia", "OCD", "ASD",
                  "Bipolar","Anorexia")

p<-ggplot(data=results, aes(x=Diagnoses, y=Mean,
                            ymin=lci,
                            ymax=uci, 
                            col=factor(sample), 
                            fill=factor(sample))) +
  geom_hline(aes(yintercept=15.5402, linetype="Healthy sibships"),
             size=1.1,col="#3c8900")+ 
  geom_point(size=5,shape=21,colour="white", stroke=0.5,
             position=position_dodge2(width=0.5, preserve = "single"))+
  geom_linerange(size=1.1,
                 position=position_dodge2(width=.5, preserve = "single"))+ 
  
  theme_minimal(base_size = 25)+
  scale_x_discrete(limits = rev(order_traits), name="Disorders")+ 
  coord_flip()+ 
  scale_fill_manual(values=c("#ee1517", "#f2a212") , 
                    labels = c("Patients", "Healthy siblings\nof patients"))  +
  scale_color_manual(values=c("#ee1517", "#f2a212"), 
                     labels = c("Patients", "Healthy siblings\nof patients")) +
  labs(y="Average number of years of education")+
  theme(panel.grid.major=element_blank(),
        legend.title=element_blank(),
        axis.title.y=element_blank(), 
        axis.title.x=element_text(),
        legend.text = element_text(
        margin = margin(l=10, t=5, b=5, unit = "pt")))+  #legend.key.size = unit(3, "cm")
  guides(fill=guide_legend(
        keyheight=10,
        default.unit="pt"))
p


ggsave(p, file="figure_MAIN_comparison_patients_siblings.tiff",
       width=11, height=9)
ggsave(p, file="Figure_5.pdf",
       width=11, height=9)




