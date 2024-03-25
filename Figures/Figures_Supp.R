## Project: EA_MH_CBS_MR 
## Script purpose: All figures in SI 
## Author: Perline Demange 

# Set up #####
library(data.table)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(corrplot)

# when minor grid coincide with major grid that has been disabled it is not appearing 
# workaround, to do everytime 
trace(ggplot2:::guide_grid, edit=TRUE)
# I manually remove 
# x.minor <- setdiff(x.minor, x.major)
# y.minor <- setdiff(y.minor, y.major)
#can re-edit with untrace(ggplot2:::guide_grid)


# Supp. Fig.2 Cooccurrence, polychoric, and genetic correlations ######

order_traits <- c("MDD", "PTSD","Alcohol","ADHD","GAD",
                  "Schizophrenia", "OCD",
                  "Bipolar", "Eating")


cooc <- fread("../CBS_Output_231113/matrix_cooccurrence_freq_sibEA_20230711.txt")
cooc <- t(as.matrix(cooc, rownames=1))
cooc_short <- cooc[order_traits, order_traits]

#supp fig2 A 
corrplot(cooc, method="square",is.corr=F, 
         na.label = "NA", diag=T,col= COL1("OrRd"))

corrplot(cooc_short, method="square",is.corr=F, 
         na.label = "NA", diag=F, col= COL1("OrRd"))

poly<- fread("../CBS_Output_231113/matrix_polychoriccorrelation_sibEA_20230711.csv")
poly <- as.matrix(poly, rownames=1)
poly_short <- poly[order_traits, order_traits]

#supp fig 2 B
corrplot(poly, method="square",is.corr=F, 
         na.label = "NA",diag=F, col= COL2("RdBu"), 
         na.label.col = "grey", type = "lower")
corrplot(poly_short, method="square",is.corr=F, 
         na.label = "NA", diag=F, col= COL1("OrRd"))

# genetic correlation matrix 
load("../MR Analyses/LDSCoutput_matrix_psychdisorders_EA.RData")
LDSCoutput$S
rG_traits <- c("ASD_Grove", "ADHD_Demontis", "ADHD_Demontis_2023", "Alc_Walters","SCZ_PGC3",
               "MDD_Howard","BIP_all_Mullins","BIP_I_Mullins", "BIP_II_Mullins",
               "Anx_Purves_meta",  "OCD_Arnold","PTSD_Nievergelt",
              "Ano_Watson", "EA", "EAwf2209") #removed "ADDH_Demontis" for the plot
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

#save rg data 
write.table(gen, "Supp.Table.geneticcorrelations.csv", row.names=F, quote=F)
write.table(pvalue, "Supp.Table.geneticcorrelations.pvalue.csv", row.names=F, quote=F)

#Supp. Fig 2. C. 
corrplot(gen, method="square",is.corr=F, 
         na.label = "NA", diag=F, col= COL1("OrRd"))


# Supp. Fig.3 Prevalence ######

#gba
#diagnoses_by_yrs_sex <- fread("diagnoses_by_yrs_sex_gba_20210225_output.csv")
#sibEA
diagnoses_by_yrs_sex <- fread("../CBS_Output_231113/diagnoses_by_yrs_sex_20230707_output.csv")
trait_list <- colnames(diagnoses_by_yrs_sex[,5:24]) #24

# All in one figure 
# Prevalence per years of education and per sex 
plot_list <- list()
plot_list_grid <- list()
for (trait in trait_list) {
  trait_sub <- subset(diagnoses_by_yrs_sex, select=c("yrs", "sex", "measure", "sample_size", trait))
  trait_sub_prev <- trait_sub[trait_sub$measure == "prevalence",]
  trait_sub_prev_CI <- trait_sub[trait_sub$measure == "prevalence_CI",]
  trait_sub_prev$CI <- trait_sub_prev_CI[[trait]]
  trait_sub_prev$sex <- as.factor(trait_sub_prev$sex)
  trait_sub_prev$sex <- factor(trait_sub_prev$sex, c(2,1))
  trait_sub_prev$CI_min <- trait_sub_prev[[trait]] - trait_sub_prev$CI 
  trait_sub_prev$CI_max <- trait_sub_prev[[trait]] + trait_sub_prev$CI 
  
  # Plot for grid of all plots
  plot_list_grid[[trait]] <- ggplot(data=trait_sub_prev, aes_string(x="yrs", y=trait, color="sex"))+
    geom_point()+ #size=6, position=position_dodge(width=0.5)
    geom_errorbar(aes(x=yrs, ymin=CI_min, ymax=CI_max))+ #, size=1, position=position_dodge(width=0.5)
    #ylim(0,0.1)+
    #labs(x= "Years of Education", color= "Sex")+
    scale_color_discrete(labels=c("Women","Men"))+
    theme_minimal()+ #base_size = 22
    theme(panel.grid.minor=element_blank(), 
          legend.position="none",
          axis.title.x= element_blank())
  
} 


jpeg("allplots_sibEA.jpg", width=7000, height=8000, res=900)
do.call(grid.arrange, c(plot_list_grid, ncol=3))
dev.off()


# Supp. Fig.4. All MR results (bipolar 1 and 2  included) #########

# EA to MH 

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
data_mr_wf$outcome[data_mr_wf$outcome == "ASD"] <- "Autism"
data_mr_wf$outcome[data_mr_wf$outcome == "alc"] <- "Alcohol"
# to plot adhd 2023 only
data_mr_wf <- data_mr_wf[!(outcome == "ADHD"),] 
data_mr_wf$outcome[data_mr_wf$outcome == "ADHD2023"] <- "ADHD"
data_mr_wf$method <- paste( data_mr_wf$method,"wf", sep="_")
results.MR.EAtoMH_wf<- data_mr_wf[,
                                  c("outcome", "or", "or_lci95",
                                    "or_uci95", "method", "pval")]
names(results.MR.EAtoMH_wf) <- c("trait", "OR", "OR_lci","OR_uci",
                                 "variable", "p_value")

results.MR.EAtoMH_wf$result <- "EAtoMH"
#keep only ivw
results.MR.EAtoMH_wf <- results.MR.EAtoMH_wf[results.MR.EAtoMH_wf$variable == "Inverse variance weighted_wf",]


# Combine 

results_MR <- rbind(results.MR.EAtoMH, results.MR.EAtoMH_wf)

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
                  "Bipolar","Bipolar_1", "Bipolar_2",
                  "Anorexia","OCD", "ASD")

results$trait <- factor(results$trait, levels = order_traits)
p<-ggplot(results, aes(x=trait, y=OR, ymin=ifelse(OR_lci>0.5, OR_lci , 0.5), 
                       ymax= ifelse(OR_uci < 2.3,OR_uci , 2.3), 
                       col=variable, fill=variable)) + 
  geom_linerange(position=position_dodge2(width=.7, preserve = "single"))+
  geom_hline(yintercept=1, lty=2)+
  geom_point(size=2,shape=21,stroke=0.5,
             position=position_dodge2(width=0.7, preserve = "single"))+
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
  coord_flip(ylim=c(0.65,1.3))+ 
  ggtitle("MR EA to MH")+
  theme_minimal(base_size = 15)+
  theme(panel.grid.major=element_blank(),
        legend.title=element_blank(),
        axis.title.y=element_blank())
p
ggsave(p, file="Figure_supp_sensitivity_wf_EAtoMH.tiff", width=11, height=9)




# MH to EA 
data_mr <- fread("../MR analyses/Results_MH_to_EA/all_MR_results_rev_221020.csv")


data_mr$exposure[data_mr$exposure == "ANX"] <- "GAD"
data_mr$exposure[data_mr$exposure == "BIP"] <- "Bipolar"
data_mr$exposure[data_mr$exposure == "BIP1"] <- "Bipolar_1"
data_mr$exposure[data_mr$exposure == "BIP2"] <- "Bipolar_2"
data_mr$exposure[data_mr$exposure == "ANNO"] <- "Anorexia"
data_mr$exposure[data_mr$exposure == "ptsd"] <- "PTSD"
data_mr$exposure[data_mr$exposure == "SCZ"] <- "Schizophrenia"
data_mr$exposure[data_mr$exposure == "ocd"] <- "OCD"
data_mr$exposure[data_mr$exposure == "ASD"] <- "ASD"
data_mr$exposure[data_mr$exposure == "alc"] <- "Alcohol"
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
data_mr_wf$exposure[data_mr_wf$exposure == "BIP1"] <- "Bipolar_1"
data_mr_wf$exposure[data_mr_wf$exposure == "BIP2"] <- "Bipolar_2"
data_mr_wf$exposure[data_mr_wf$exposure == "ANNO"] <- "Anorexia"
data_mr_wf$exposure[data_mr_wf$exposure == "ptsd"] <- "PTSD"
data_mr_wf$exposure[data_mr_wf$exposure == "SCZ"] <- "Schizophrenia"
data_mr_wf$exposure[data_mr_wf$exposure == "ocd"] <- "OCD"
data_mr_wf$exposure[data_mr_wf$exposure == "ASD"] <- "ASD"
data_mr_wf$exposure[data_mr_wf$exposure == "alc"] <- "Alcohol"
data_mr_wf <- data_mr_wf[!(exposure == "ADHD"),] 
data_mr_wf$exposure[data_mr_wf$exposure == "ADHD2023"] <- "ADHD"
data_mr_wf$method <- paste( data_mr_wf$method,"wf", sep="_")
results.MR.MHtoEA_wf<- data_mr_wf[,
                                  c("exposure", "b", "lo_ci",
                                    "up_ci", "method", "pval")]
names(results.MR.MHtoEA_wf) <- c("trait", "Beta", "Beta_lci","Beta_uci",
                                 "variable", "p_value")

results.MR.MHtoEA_wf$result <- "MHtoEA"
#keep only ivw
results.MR.MHtoEA_wf <- results.MR.MHtoEA_wf[
  results.MR.MHtoEA_wf$variable == "Inverse variance weighted_wf",]


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
                  "Bipolar","Bipolar_1", "Bipolar_2",
                  "Anorexia","OCD", "ASD")

results$trait <- factor(results$trait, levels = order_traits)
p<-ggplot(results, aes(x=trait, y=Beta, ymin= Beta_lci, 
                       ymax= Beta_uci,
                       col=variable, fill=variable)) + 
  geom_linerange(position=position_dodge2(width=.7, preserve = "single"))+
  geom_hline(yintercept=0, lty=2)+
  geom_point(size=2,shape=21,stroke=0.5,
             position=position_dodge2(width=0.7, preserve = "single"))+
  scale_x_discrete(limits = rev(order_traits), name="Disorders")+ #
  scale_y_continuous(name= "Beta", minor_breaks = seq(-0.60,0.3, 0.1))+
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
  coord_flip(ylim=c(-0.65,0.3))+ 
  ggtitle("MR MH to EA")+
  theme_minimal(base_size = 15)+
  theme(panel.grid.major=element_blank(),
        legend.title=element_blank(),
        axis.title.y=element_blank())
p
ggsave(p, file="Figure_supp_sensitivity_wf_MHtoEA.tiff", width=11, height=9)

# Supp. Fig.5. Descriptives of expenditures ##########
cost_sex_gba <- fread("../CBS_Output_220908/costs_range_by_sex_gbaEA_2209005.csv")
head(cost_sex_gba)
cost_sex_sib <- fread("../CBS_Output_220908/costs_range_by_sex_sibEA_2209005.csv")
head(cost_sex_sib)

levels_range <- unique(cost_sex_gba$range) # get the order of the ranges as it should 
cost_sex_sib$range <- as.factor(cost_sex_sib$range) 
cost_sex_gba$range <- as.factor(cost_sex_gba$range)

## For figure with total numbers 

cost_tot_gba <- cost_sex_gba %>% 
  group_by(range) %>%
  summarise(Freq = sum(Freq))

cost_tot_sib <- cost_sex_sib %>% 
  group_by(range) %>%
  summarise(Freq = sum(Freq))


cost_tot_sib$frequency <- cost_tot_sib$Freq/sum(cost_tot_sib$Freq)
cost_tot_sib$range <- factor(cost_tot_sib$range, levels=levels_range)
cost_tot_sib$sample <- "sib"

cost_tot_gba$frequency <- cost_tot_gba$Freq/sum(cost_tot_gba$Freq)
cost_tot_gba$range <- factor(cost_tot_gba$range, levels=levels_range)
cost_tot_gba$sample <- "gba"

cost_tot <- rbind(cost_tot_gba, cost_tot_sib)

## Figure with both genders and both samples 
cost_sex_gba %>% 
  group_by(GBAGESLACHT) %>%
  summarise(Freq = sum(Freq))
#sample size of 1 is 1553683 and of 2 is 1627824

cost_sex_sib%>% 
  group_by(GBAGESLACHT) %>%
  summarise(Freq = sum(Freq))
#sample size of 1 is 842105 and of 2 is 846248

cost_sex_sib$frequency <- ifelse(cost_sex_sib$GBAGESLACHT == 1, 
                                 cost_sex_sib$Freq/842105,
                                 cost_sex_sib$Freq/846248)
cost_sex_sib$range <- factor(cost_sex_sib$range, levels=levels_range)
cost_sex_sib$sample <- "sib"

cost_sex_gba$frequency <- ifelse(cost_sex_gba$GBAGESLACHT == 1, 
                                 cost_sex_gba$Freq/1553683,
                                 cost_sex_gba$Freq/1627824)
cost_sex_gba$range <- factor(cost_sex_gba$range, levels=levels_range)
cost_sex_gba$sample <- "gba"

cost_sex <- rbind(cost_sex_gba, cost_sex_sib)

p<-ggplot(data=cost_sex, aes(x=range, y=frequency, fill=factor(GBAGESLACHT))) +
  geom_bar(stat= "identity", position="dodge")+ 
  facet_wrap(~sample)+
  theme_minimal(base_size = 15)+
  scale_x_discrete(labels = c(
    "(-1,0]"= "0", 
    "(100,200]" = "(100,200]", 
    "(900,1e+03]"= "(900,1000]", 
    "(3e+03,3e+05]"= "(3000,+]"))+ #limits = c("(-1,0]","(100,200]","(900,1e+03]", "(3e+03,3e+05]") only show, while I want to remove ticks 
  theme(panel.grid.major=element_blank(),
        legend.title=element_blank(),
        axis.title.y=element_blank(), 
        axis.title.x=element_blank(), 
        axis.ticks.x=element_line(color="black"), 
        axis.text.x= element_text(angle=90))+
  ggtitle("Frequency of mental health expenditures")
tiff("figure_descriptive_costs_sex_gbasib.tiff", width=700, height=480)
p
dev.off()
# clean manually on powerpoint


# Supp. Fig.6. Within-sib in same sex sibships ##########################
order_traits_all_CBS<- c("Any", "ASD", "ADHD", 
                         "Alcohol","Schizos", "Schizophrenia",
                         "MDD","Bipolar", "Bipolar_1", "Bipolar_2", 
                         "Anxiety", "GAD", "Panic", "Phobia", "OCD", "PTSD",
                         "Anorexia", "Bulimia", 
                         "Personality", "ClusterA", "ClusterB", "ClusterC")


results_CBSwomen <- fread("../CBS_Output_231113/glm_sib_EA_women_results_20231108.csv")
results_CBSwomen$OR_robust <- as.numeric(results_CBSwomen$OR_robust)
results_CBSwomen <- results_CBSwomen[which(results_CBSwomen$variable == "EA_within"| results_CBSwomen$variable == "yrs"),]
results_CBSwomen$OR_lci_robust <- (results_CBSwomen$OR_robust - (results_CBSwomen$OR.SE_robust*1.96))
results_CBSwomen$OR_uci_robust <- (results_CBSwomen$OR_robust + (results_CBSwomen$OR.SE_robust*1.96))
results_CBSwomen <- results_CBSwomen[, c( "trait", "OR_robust", "OR_lci_robust",
                                "OR_uci_robust","variable", "p_value_robust")]
names(results_CBSwomen ) <- c("trait", "OR", "OR_lci","OR_uci",
                         "variable", "p_value")

results_CBSwomen$result <- "Women"

results_CBSmen <- fread("../CBS_Output_231113/glm_sib_EA_men_results_20231108.csv")
results_CBSmen$OR_robust <- as.numeric(results_CBSmen$OR_robust)
results_CBSmen <- results_CBSmen[which(results_CBSmen$variable == "EA_within"| results_CBSmen$variable == "yrs"),]
results_CBSmen$OR_lci_robust <- (results_CBSmen$OR_robust - (results_CBSmen$OR.SE_robust*1.96))
results_CBSmen$OR_uci_robust <- (results_CBSmen$OR_robust + (results_CBSmen$OR.SE_robust*1.96))
results_CBSmen <- results_CBSmen[, c( "trait", "OR_robust", "OR_lci_robust",
                                          "OR_uci_robust","variable", "p_value_robust")]
names(results_CBSmen ) <- c("trait", "OR", "OR_lci","OR_uci",
                              "variable", "p_value")
results_CBSmen$result <- "Men"

results_CBSsamesex <- rbind(results_CBSwomen,results_CBSmen)

results <- results_CBSsamesex
order_traits <- order_traits_all_CBS 
results$variable <- factor(results$variable, levels = c("EA_within",
                                                        "yrs" ), 
                           labels = c("Within-sibling effect",
                                      "Population effect"))
results$trait <- factor(results$trait, levels = order_traits)

p<-ggplot(results, aes(x=trait, y=OR,
                       ymin=ifelse(OR_lci>0.69, OR_lci , 0.69),
                       ymax= ifelse(OR_uci < 1.3,OR_uci , 1.3),
                       col=variable,
                       fill=variable)) + 
  geom_linerange(position=position_dodge2(width=.5, preserve = "single"))+
  geom_hline(yintercept=1, lty=1)+
  geom_point(size=3,shape=21,colour="white", stroke=0.5,
             position=position_dodge2(width=0.5, preserve = "single"))+
  scale_x_discrete(limits = rev(order_traits), name="Disorders")+ #
  scale_y_continuous(name= "Odds ratio", minor_breaks = seq(0.6,1.2, 0.1))+
  coord_flip()+ 
  facet_wrap(~result)+
  theme_minimal(base_size = 15)+
  scale_fill_manual(values=c("black", "#be2596") , 
                    breaks = c("Population effect",
                               "Within-sibling effect"))  +
  scale_color_manual(values=c("black", "#be2596"), 
                     breaks = c("Population effect",
                                "Within-sibling effect")) +
  theme(panel.grid.major=element_blank(),
        legend.title=element_blank(),
        axis.title.y=element_blank())
tiff("figure_allCBS_women_men.tiff", width=700, height=480)
p
dev.off()

# Supp. Fig.7. Within-sib, all, excluding 11 and 2 years EA ############
order_traits_all_CBS<- c("Any","ASD","ADHD", 
                         "Alcohol","Schizos", "Schizophrenia",
                         "MDD","Bipolar", "Bipolar_1", "Bipolar_2", 
                         "Anxiety", "GAD", "Panic", "Phobia", "OCD", "PTSD",
                         "Anorexia", "Bulimia", 
                         "Personality", "ClusterA", "ClusterB", "ClusterC")


results_CBSa <- fread("../CBS_Output_231113/glm_sib_EA_results_20230922.csv")
results_CBSa$OR_robust <- as.numeric(results_CBSa$OR_robust)
results_CBSa <- results_CBSa[which(results_CBSa$variable == "EA_within"|
                                     results_CBSa$variable == "yrs"),]
results_CBSa$OR_lci_robust <- (results_CBSa$OR_robust - (results_CBSa$OR.SE_robust*1.96))
results_CBSa$OR_uci_robust <- (results_CBSa$OR_robust + (results_CBSa$OR.SE_robust*1.96))

results_CBSa <- results_CBSa[, c("trait", 
                                  "OR_robust", "OR_lci_robust", "OR_uci_robust",
                                  "variable", "p_value_robust")]
names(results_CBSa ) <- c("trait", "OR", "OR_lci","OR_uci",
                            "variable", "p_value")
results_CBSa$result <- "All years of education"

results_CBS11 <- fread("../CBS_Output_231113/glm_sib_EA_results_excl11_20230922.csv")
results_CBS11$OR_robust <- as.numeric(results_CBS11$OR_robust)
results_CBS11 <- results_CBS11[which(results_CBS11$variable == "EA_within"|
                                     results_CBS11$variable == "yrs"),]
results_CBS11$OR_lci_robust <- (results_CBS11$OR_robust - (results_CBS11$OR.SE_robust*1.96))
results_CBS11$OR_uci_robust <- (results_CBS11$OR_robust + (results_CBS11$OR.SE_robust*1.96))

results_CBS11 <- results_CBS11[, c("trait", 
                                 "OR_robust", "OR_lci_robust", "OR_uci_robust",
                                 "variable", "p_value_robust")]
names(results_CBS11 ) <- c("trait", "OR", "OR_lci","OR_uci",
                          "variable", "p_value")
results_CBS11$result <- "Excluding 11 years of education"


results_CBS2 <- fread("../CBS_Output_231113/glm_sib_EA_results_excl2_20230922.csv")
results_CBS2$OR_robust <- as.numeric(results_CBS2$OR_robust)
results_CBS2 <- results_CBS2[which(results_CBS2$variable == "EA_within"|
                                     results_CBS2$variable == "yrs"),]
results_CBS2$OR_lci_robust <- (results_CBS2$OR_robust - (results_CBS2$OR.SE_robust*1.96))
results_CBS2$OR_uci_robust <- (results_CBS2$OR_robust + (results_CBS2$OR.SE_robust*1.96))

results_CBS2 <- results_CBS2[, c("trait", 
                                 "OR_robust", "OR_lci_robust", "OR_uci_robust",
                                 "variable", "p_value_robust")]
names(results_CBS2 ) <- c("trait", "OR", "OR_lci","OR_uci",
                          "variable", "p_value")
results_CBS2$result <- "Excluding 2 years of education"

results_CBSEA <- rbind(results_CBSa,results_CBS11)
results_CBSEA <- rbind(results_CBSEA,results_CBS2)

results <- results_CBSEA
order_traits <- order_traits_all_CBS 
results$variable <- factor(results$variable, levels = c("EA_within",
                                                        "yrs" ), 
                           labels = c("Within-sibling effect",
                                      "Population effect"))
results$trait <- factor(results$trait, levels = order_traits)

p<-ggplot(results, aes(x=trait, y=OR,
                       ymin=ifelse(OR_lci>0.69, OR_lci , 0.69),
                       ymax= ifelse(OR_uci < 1.3,OR_uci , 1.3),
                       col=variable,
                       fill=variable)) + 
  geom_linerange(position=position_dodge2(width=.5, preserve = "single"))+
  geom_hline(yintercept=1, lty=1)+
  geom_point(size=3,shape=21,colour="white", stroke=0.5,
             position=position_dodge2(width=0.5, preserve = "single"))+
  scale_x_discrete(limits = rev(order_traits), name="Disorders")+ #
  scale_y_continuous(name= "Odds ratio", minor_breaks = seq(0.6,1.2, 0.1))+
  coord_flip()+ 
  facet_wrap(~result)+
  theme_minimal(base_size = 15)+
  scale_fill_manual(values=c("black", "#be2596") , 
                    breaks = c("Population effect",
                               "Within-sibling effect"))  +
  scale_color_manual(values=c("black", "#be2596"), 
                     breaks = c("Population effect",
                                "Within-sibling effect")) +
  theme(panel.grid.major=element_blank(),
        legend.title=element_blank(),
        axis.title.y=element_blank())

tiff("figure_allCBS_excl11_excl2.tiff", width=1000, height=480)
p
dev.off()


# Supp. Fig.8. compare EA sib and patients #####################

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

p<-ggplot(data=results, aes(x=reorder(Diagnoses, Mean), y=Mean,
                            ymin=lci,
                            ymax=uci, 
                            col=factor(sample), 
                            fill=factor(sample))) +
  geom_point(size=3,shape=21,colour="white", stroke=0.5,
             position=position_dodge2(width=0.5, preserve = "single"))+
  geom_linerange(position=position_dodge2(width=.5, preserve = "single"))+ 
  geom_hline(yintercept=15.5402, col="red")+ 
  theme_minimal(base_size = 15)+
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
tiff("figure_comparison_patients_siblings.tiff", width=700, height=700)
p
dev.off()




# Additional non presented figures ######
## descriptives of education data   

diploma <- fread("C:/Users/user/Dropbox/CBS - MR/CBS_Ouput_210225/descriptive_table_count_diploma_20210223.csv",
                 header=F)
diploma
yrs <- fread("C:/Users/user/Dropbox/CBS - MR/CBS_Ouput_210225/descriptive_table_count_yrs_20210223.csv", 
             header=F)
yrs_sib <- as.data.frame(t(yrs[,1:10]))
colnames(yrs_sib) <- c( "yrs", "gba", "sib")
barplot(yrs_sib$sib, 
        names.arg = yrs_sib$yrs, 
        ylim=c(0,9e+05), 
        main="Sample size per years of education - in siblings")
barplot(yrs_sib$gba, 
        names.arg = yrs_sib$yrs,
        ylim=c(0,9e+05), 
        main="Sample size per years of education - in global population")

# do frequency instead
yrs_sib$sib_freq <- yrs_sib$sib/sum(yrs_sib$sib)
barplot(yrs_sib$sib_freq,  names.arg = yrs_sib$yrs, ylim=c(0,0.5), 
        main="Frequency per years of education - in siblings", 
        cex.main=1.5,
        cex.axis=1.5, 
        cex.names=1.5)


yrs_sib$gba_freq <- yrs_sib$gba/sum(yrs_sib$gba)
barplot(yrs_sib$gba_freq,
        names.arg = yrs_sib$yrs,
        ylim=c(0,0.5),
        main="Frequency per years of education - in global population")

summary(yrs_sib)
#mean
sum(yrs_sib$yrs * yrs_sib$gba)/sum(yrs_sib$gba)
sum(yrs_sib$yrs * yrs_sib$sib)/sum(yrs_sib$sib)
