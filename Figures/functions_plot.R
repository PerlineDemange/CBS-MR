#Function to make it easier to plot \

prevalence_plot_third <- function(trait_list, ylim_max){
  plot_list <- list()
  plot_list_grid <- list()
  for (trait in trait_list) {
    trait_sub <- subset(diagnoses_by_yrs_sex, select=c("yrs", "sex", "measure", "sample_size", trait))
    trait_sub_prev <- trait_sub[trait_sub$measure == "prevalence",]
    trait_sub_prev_CI <- trait_sub[trait_sub$measure == "prevalence_CI",]
    trait_sub_prev$CI <- trait_sub_prev_CI[[trait]]
    trait_sub_prev$sex <- as.factor(trait_sub_prev$sex)
    trait_sub_prev$sex <- factor(trait_sub_prev$sex, c(2,1))
    trait_sub_prev[[trait]] <- trait_sub_prev[[trait]]*100
    trait_sub_prev$CI_min <- trait_sub_prev[[trait]] - trait_sub_prev$CI*100 
    trait_sub_prev$CI_max <- trait_sub_prev[[trait]] + trait_sub_prev$CI*100 
    
    plot_list[[trait]] <- ggplot(data=trait_sub_prev, 
                                 aes_string(x="yrs", y=trait, color="sex"))+
      geom_point(size=4, position=position_dodge(width=0.5))+ #
      geom_errorbar(aes(x=yrs, ymin=CI_min, ymax=CI_max), 
                    size=1, position=position_dodge(width=0.5))+ # 
      ylim(0, ylim_max)+
      labs(x= "Years of Education", color= "Sex")+
      scale_color_discrete(labels=c("Women","Men" ))+
      theme_minimal(base_size = 22)+ 
      ggtitle(trait)+
      theme(panel.grid.minor=element_blank(), 
            legend.position="none",
            axis.title.x= element_blank(), 
            axis.title.y= element_blank(), 
            plot.title = element_text(hjust = 0.5))
  }
  return(plot_list) 
} 
