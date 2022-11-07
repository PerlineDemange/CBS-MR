# Health data
# Costs 
# Author: Perline Demange
# Project: CBS-MR EA-MH 
# 

# Costs #############
# Load gba6585 to get ID of individuals 
gba <- fread("gba_6585_20201015.csv")
gba <- gba[,2 ]

#1. Get all costs for all years #####
data_list <- as.data.frame(fread("H:/Data/Costs/Costs_file_list.csv", header=T))
gba_allyears <- NULL 
year <- data_list$year
file_list_cost <- data_list$costfile
for (y in 1:nrow(data_list)) { #for each year
  gba_year <- gba
  gba_year$year <- year[y] #get year
  print(year[y])
  costdata <- fread(file_list_cost[y], dec=",") #load data cost for the correct year
  costdata <- costdata[,c("RINPERSOON", "ZVWKEERSTELIJNSPSYCHO","ZVWKGGZ",
                          "ZVWKGENBASGGZ", "ZVWKSPECGGZ") ] #first two are NA fr60 2014 and last tow are NA before 
  costdata$cost_psy_year <- rowSums(costdata[, c("ZVWKEERSTELIJNSPSYCHO",
                                                 "ZVWKGGZ","ZVWKGENBASGGZ",
                                                  "ZVWKSPECGGZ")], na.rm=T)#sum of 1st line and 2nd line 
  cost_year <- left_join(gba_year, costdata, by = "RINPERSOON") # left join instead of merge as some individual miss only some years
  gba_allyears <- rbind(gba_allyears, cost_year)
  } 
    
#write.csv(gba_allyears, "gba_costs_allyears_20210307.csv", row.names=F)
gba_allyears <- fread("gba_costs_allyears_20210307.csv")

head(gba_allyears)

#2. Check if some people have no data and exclude them #####
missing_data <- gba_allyears %>% 
  group_by(RINPERSOON) %>%
  summarize_if(is.numeric, mean, na.rm = TRUE)  #ignore Na so that only people with only NA get NA (if one NA do the average of the rest)

head(missing_data)
summary(missing_data$cost_psy_year) # NAs = 1486364

gba_allyears_cl <- gba_allyears[!(gba_allyears$RINPERSOON %in% 
                                    missing_data[is.na(missing_data$cost_psy_year)==T,]$RINPERSOON), ]
#65397670 - 50524030 =14873640 and there is ten years so it fits 
head(gba_allyears_cl)


rm(gba)
rm(gba_year)
rm(missing_data)
rm(costdata)
rm(gba_allyears)
rm(cost_year)

#3. Exclude people with negative costs ######
# Negative costs is not supposed to happen. 
# It is unclear what negative costs means (paid too much the year before? 
#but sometimes doesnt fit, etc.)
# CBS gets data from Vektis, and Vektis checks the data but does not make 
#any correction, which is why sometimes there is negative amount
# Because of this, and the low sample of people with negative costs, 
#we exclude everyone who has negative costs once
summary(gba_allyears) # there are people with negative costs in the total gba6585 data

nrow(gba_allyears_cl[gba_allyears_cl$ZVWKEERSTELIJNSPSYCHO <0, ])#60
nrow(gba_allyears_cl[gba_allyears_cl$ZVWKGGZ <0, ])#116
nrow(gba_allyears_cl[gba_allyears_cl$ZVWKGENBASGGZ <0, ]) #91
nrow(gba_allyears_cl[gba_allyears_cl$ZVWKSPECGGZ <0, ]) #63
nrow(gba_allyears_cl[gba_allyears_cl$ZVWKEERSTELIJNSPSYCHO <0|
                       gba_allyears_cl$ZVWKGGZ <0 | 
                       gba_allyears_cl$ZVWKGENBASGGZ <0 |
                       gba_allyears_cl$ZVWKGENBASGGZ <0|
                       gba_allyears_cl$ZVWKSPECGGZ <0 ,]) #330

gba_allyears_cl <- gba_allyears_cl[!which(gba_allyears_cl$ZVWKEERSTELIJNSPSYCHO <0|
                                                         gba_allyears_cl$ZVWKGGZ <0 | 
                                                        gba_allyears_cl$ZVWKGENBASGGZ <0 |
                                                         gba_allyears_cl$ZVWKGENBASGGZ <0|
                                                         gba_allyears_cl$ZVWKSPECGGZ <0), ]
#exclude 330 people out of 5 millions

#4. Exclude people with more than 50% missing years in 2009-2013 and 2014-2018 #### 
# Check people with missing some years # to check after merging with siblings 
# possible reasons for missing years are living abroad, death(?),
# eventually request to be out of the system (not likely)
summary(gba_allyears_cl$cost_psy_year) #3941690  missing years 

gba_allyears_cl_0913 <- gba_allyears_cl[gba_allyears_cl$year < 2014, 1:4]
gba_allyears_cl_1418 <- gba_allyears_cl[gba_allyears_cl$year >= 2014, -c(3,4,7)]

#Identify people with more than half of the years missing to exclude 
# 2009-2013: 5 years, 2014-2018: 5 years 
missing_data_0913 <- gba_allyears_cl_0913 %>% 
  group_by(RINPERSOON) %>%
  summarize_at(2:3, .funs = funs('NA' = sum(is.na(.)))) 
summary(missing_data_0913)
nrow(missing_data_0913[missing_data_0913$ZVWKGGZ_NA > 2.5,]) #423960
nrow(missing_data_0913[missing_data_0913$ZVWKEERSTELIJNSPSYCHO_NA > 2.5,])#423960

missing_data_0913_id <- missing_data_0913[missing_data_0913$ZVWKEERSTELIJNSPSYCHO_NA > 2.5,]$RINPERSOON
  
missing_data_1418 <- gba_allyears_cl_1418 %>% 
  group_by(RINPERSOON) %>%
  summarize_at(2:3, .funs = funs('NA' = sum(is.na(.)))) 

summary(missing_data_1418)
nrow(missing_data_1418[missing_data_1418$ZVWKGENBASGGZ_NA > 2.5,]) #383902
nrow(missing_data_1418[missing_data_1418$ZVWKSPECGGZ_NA  > 2.5,]) #383902

missing_data_1418_id <- missing_data_1418[missing_data_1418$ZVWKSPECGGZ_NA  > 2.5,]$RINPERSOON

gba_allyears_cl <- gba_allyears_cl[!(gba_allyears_cl$RINPERSOON %in% missing_data_0913_id |gba_allyears_cl$RINPERSOON %in% missing_data_1418_id), ]
#  50533700  - 44836628 = 5697072
length(unique(c(missing_data_0913_id, missing_data_1418_id))) #569708
length(unique(gba_allyears_cl$RINPERSOON)) #4483695 #lost 569708 individuals 
# data is missing 8 rows (from one person or several )
# check what is going on 

missing_rows <- gba_allyears_cl %>% 
  group_by(RINPERSOON) %>%
  summarize(year, sum(year)) 
head(missing_rows)
missing_rows[missing_rows$`sum(year)` < 20135,]


gba_allyears_cl[gba_allyears_cl$RINPERSOON == 155,]

#write.csv(gba_allyears_cl, "gba_costs_allyears_cl_20220307.csv", row.names=F)
gba_allyears_cl <- fread("gba_costs_allyears_cl_20220307.csv")
rm(missing_data_0913)
rm(missing_data_1418)
rm(gba_allyears_cl_0913)
rm(gba_allyears_cl_1418)


#5. Histogram of the population ######

summary(gba_allyears_cl$cost_psy_year)
hist(gba_allyears_cl$cost_psy_year)

#6. Data transformation #######
## 6.1 Log transformation #####
gba_allyears_cl$cost_psy_year_log <- log(gba_allyears_cl$cost_psy_year + 1)
hist(gba_allyears_cl$cost_psy_year_log)
#write.csv(gba_allyears_cl, "gba_costs_allyears_cl_20220307.csv", row.names=F)

summary(gba_allyears_cl$cost_psy_year_log)
summary(gba_allyears_cl$cost_psy_year)

## 6.2 Sum and mean over years ####

gba_allyears_sum <- gba_allyears_cl %>% 
  group_by(RINPERSOON) %>%
  summarise_at(vars(cost_psy_year, cost_psy_year_log), sum, na.rm = TRUE)
summary(gba_allyears_sum)
gba_allyears_mean <- gba_allyears_cl %>% 
  group_by(RINPERSOON) %>%
  summarise_at(vars(cost_psy_year, cost_psy_year_log), mean, na.rm = TRUE)
summary(gba_allyears_mean)

gba_allyears <- merge(gba_allyears_mean, gba_allyears_sum, by = "RINPERSOON") # 4483695
colnames(gba_allyears) <- c("RINPERSOON", "cost_mean","cost_log_mean", "cost_sum",  "cost_log_sum")
hist(gba_allyears$cost_mean)
hist(gba_allyears$cost_sum_log)

## 6.3 Log transformation after sum and mean #####
gba_allyears$cost_mean_log <- log(gba_allyears$cost_mean +1) #this is the variable we use in the analyses
gba_allyears$cost_sum_log <- log(gba_allyears$cost_sum +1)

#write.csv(gba_allyears, "gba_costs_allyears_summean_20220307.csv", row.names=F)
gba_allyears<- fread("gba_costs_allyears_summean_20220307.csv")

# 7. Descriptives in full population ######
gba_mean <- c(mean(gba_allyears$cost_mean),median(gba_allyears$cost_mean),
                 sd(gba_allyears$cost_mean), nrow(gba_allyears))
gba_sum <- c(mean(gba_allyears$cost_sum),median(gba_allyears$cost_sum),
                sd(gba_allyears$cost_sum), nrow(gba_allyears))

## 7.1 Create categories following price range of De Zeeuw. ######
# because the higher bound is lower than De zeeuw, I pick a smaller range, and 
# I reduce the highest range to 3000 because the mean cost is lower than in de zeeuw (only MH costs here)
gba_allyears$range <- cut(gba_allyears$cost_mean, 
                          breaks=c(-1,0, 100, 200, 300, 400, 500,
                                   600, 700, 800, 900,1000,
                                   1100, 1200, 1300, 1400, 1500,
                                   1600, 1700, 1800, 1900, 2000,
                                   2100, 2200, 2300, 2400, 2500,
                                   2600, 2700, 2800, 2900, 3000,300000))

#boudaries are (,]

gba_table <- table(gba_allyears$range)
barplot(table(gba_allyears$range)) # shoudl remove 0 to have a visible plot 


# 8. Descriptives in full population gba with EA ##########
## 8.1 Merge with gba EA data #######
gba_EA <- fread("gba6585_edu_subset_20210219.csv")
gba_allyears_EA <- inner_join(gba_allyears, gba_EA, by= "RINPERSOON") #3181507

## 8.2 Create categories following price range of De Zeeuw. ######
gba_mean_EA <- c(mean(gba_allyears_EA$cost_mean),median(gba_allyears_EA$cost_mean),
                 sd(gba_allyears_EA$cost_mean), nrow(gba_allyears_EA))
gba_sum_EA <- c(mean(gba_allyears_EA$cost_sum), median(gba_allyears_EA$cost_sum),
                sd(gba_allyears_EA$cost_sum), nrow(gba_allyears_EA))

gba_allyears_EA$range <- cut(gba_allyears_EA$cost_mean, 
                          breaks=c(-1,0,100, 200, 300, 400, 500,
                                   600, 700, 800, 900,1000,
                                   1100, 1200, 1300, 1400, 1500,
                                   1600, 1700, 1800, 1900, 2000,
                                   2100, 2200, 2300, 2400, 2500,
                                   2600, 2700, 2800, 2900, 3000,300000))


gba_table_EA <- table(gba_allyears_EA$range)
barplot(table(gba_allyears_EA$range))

## 8.3 Histogram of average costs, per sex #################

gba_table_EA_sex <- gba_allyears_EA %>% 
  group_by(range, GBAGESLACHT) %>%
  summarize(Freq=n())

ggplot(gba_table_EA_sex, aes(x=factor(range), y=Freq, fill=factor(GBAGESLACHT))) +
  geom_col(position="dodge")

write.table(gba_table_EA_sex, "costs_range_by_sex_gbaEA_2209005.csv",
            row.names=F, quote=F, sep=" ")

## 8.4 Histogram of average costs, per EA in years ####
gba_table_EA_yrs <- gba_allyears_EA %>% 
  group_by(range, yrs) %>%
  summarize(Freq=n())

#many yrs 22 are <10 so I dont output this now. #not that useful 
ggplot(gba_table_EA_yrs, aes(x=factor(range), y=Freq, fill=factor(yrs))) +
  geom_col(position="dodge")


# 9. Get subset with Sib Edu ######
sib_EA <- fread("sib_edu_subset_20210219.csv") #1743032
sib_allyears <- gba_allyears[which(gba_allyears$RINPERSOON %in% sib_EA$RINPERSOON),]
# 1709851

save(sib_allyears, file= "sib_allyears_costs.Rda")

## 9.1 Merge with siblings data######
load("sib_allyears_costs.Rda")
sib_EA <- fread("sib_edu_subset_20210219.csv")
head(sib_EA)


data_cost_sib <- merge(sib_EA, sib_allyears, by = "RINPERSOON")


## 9.2 Get new number of family size and remove families < 2 #####
count <- as.data.frame(table(data_cost_sib$FID)) #get number of occurences of FID
data_cost_sib$Fsize_cost <- count[[2]][match(data_cost_sib$FID, count[[1]])] # match data fid with value in count, if true takes value of Freq 
data_cost_sib <- data_cost_sib[data_cost_sib$Fsize_cost > 1,] #1688353

length(unique(data_cost_sib$FID) ) #744 428 families 
nrow(data_cost_sib)/ length(unique(data_cost_sib$FID)) #Average of 2.27 sibling per sibships

## 9.3. Save #####
save(data_cost_sib, file= "data_cost_sib.Rda")

load("data_cost_sib.Rda")

#10. Descriptives in siblings EA pop ##################
hist(data_cost_sib$cost_mean)
hist(data_cost_sib$cost_mean_log)
hist(data_cost_sib$cost_log_mean)

## 10.1 Create categories following price range of De Zeeuw. ######
sib_mean_EA<- c(mean(data_cost_sib$cost_mean),median(data_cost_sib$cost_mean),
                 sd(data_cost_sib$cost_mean), nrow(data_cost_sib))
sib_sum_EA <- c(mean(data_cost_sib$cost_sum), median(data_cost_sib$cost_sum),
                sd(data_cost_sib$cost_sum), nrow(data_cost_sib))

data_cost_sib$range <- cut(data_cost_sib$cost_mean, 
                             breaks=c(-1,0,100, 200, 300, 400, 500,
                                      600, 700, 800, 900,1000,
                                      1100, 1200, 1300, 1400, 1500,
                                      1600, 1700, 1800, 1900, 2000,
                                      2100, 2200, 2300, 2400, 2500,
                                      2600, 2700, 2800, 2900, 3000,300000))


sib_table_EA <- table(data_cost_sib$range)
barplot(table(data_cost_sib$range))

## 10.2 Histogram of average costs, per sex #################

sib_table_EA_sex <- data_cost_sib %>% 
  group_by(range, GBAGESLACHT) %>%
  summarize(Freq=n())

ggplot(sib_table_EA_sex, aes(x=factor(range), y=Freq, fill=factor(GBAGESLACHT))) +
  geom_col(position="dodge")

write.table(sib_table_EA_sex, "costs_range_by_sex_sibEA_2209005.csv",
            row.names=F, quote=F, sep=" ")

## 10.3 Histogram of average costs, per birthorder ####
data_cost_sib$birth_order_cat <- data_cost_sib$birth_order
data_cost_sib$birth_order_cat[data_cost_sib$birth_order >= 3] <- "3+"
table(data_cost_sib$birth_order_cat)

sib_table_EA_birthorder <- data_cost_sib %>% 
  group_by(range, birth_order_cat) %>%
  summarize(Freq=n())

ggplot(sib_table_EA_birthorder, aes(x=factor(range), y=Freq, fill=factor(birth_order_cat))) +
  geom_col(position="dodge")

write.table(sib_table_EA_birthorder, "costs_range_by_birthorder_sibEA_2209005.csv",
            row.names=F, quote=F, sep=" ")

## 10.4 Histogram of average costs, per year of birth ####

data_cost_sib$birthyear <- data_cost_sib$GBAGEBOORTEJAAR
data_cost_sib$birthyear[data_cost_sib$GBAGEBOORTEJAAR >= 1981] <- "81-85"
data_cost_sib$birthyear[data_cost_sib$birthyear >= 1976 & 
                               data_cost_sib$birthyear < 1981] <- "76-80"
data_cost_sib$birthyear[data_cost_sib$birthyear >= 1971& 
                               data_cost_sib$birthyear < 1976] <- "71-75"
data_cost_sib$birthyear[data_cost_sib$birthyear >= 1965& 
                               data_cost_sib$birthyear < 1971] <- "65-70"
table(data_cost_sib$birthyear)


sib_table_EA_birthyear <- data_cost_sib %>% 
  group_by(range, birthyear) %>%
  summarize(Freq=n())

ggplot(sib_table_EA_birthyear, aes(x=factor(range), y=Freq, fill=factor(birthyear))) +
  geom_col(position="dodge")

write.table(sib_table_EA_birthyear, "costs_range_by_birthyear_sibEA_2209005.csv",
            row.names=F, quote=F, sep=" ")


# 11. Get basic descriptive in full sib sample ####
sib <- fread("gba_6585_FID_fullsibship_final_20210215.csv")
sib_allyears <- inner_join(gba_allyears, sib, by= "RINPERSOON") #3004885

## 11.1 Get new number of family size and remove families < 2 #####
count <- as.data.frame(table(sib_allyears$FID)) #get number of occurences of FID
sib_allyears$Fsize_cost <- count[[2]][match(sib_allyears$FID, count[[1]])] # match data fid with value in count, if true takes value of Freq 
sib_allyears <- sib_allyears[sib_allyears$Fsize_cost > 1,] #2886902

length(unique(sib_allyears$FID) ) #1 223 144 families 
nrow(sib_allyears)/ length(unique(sib_allyears$FID)) #Average of 2.36 sibling per sibships

## 11.2 Get descriptives ######

sib_mean <- c(mean(sib_allyears$cost_mean),median(sib_allyears$cost_mean),
              sd(sib_allyears$cost_mean), nrow(sib_allyears))
sib_sum <- c(mean(sib_allyears$cost_sum),median(sib_allyears$cost_sum),
             sd(sib_allyears$cost_sum), nrow(sib_allyears))

sib_allyears$range <- cut(sib_allyears$cost_mean, 
                          breaks=c(-1,0, 100, 200, 300, 400, 500,
                                   600, 700, 800, 900,1000,
                                   1100, 1200, 1300, 1400, 1500,
                                   1600, 1700, 1800, 1900, 2000,
                                   2100, 2200, 2300, 2400, 2500,
                                   2600, 2700, 2800, 2900, 3000,300000))

#boudaries are (,]

sib_table <- table(sib_allyears$range)

## 11.2 Save basic Descriptives for all #####

desc_mean <- as.data.frame(rbind(gba_mean, gba_mean_EA, sib_mean, sib_mean_EA))
colnames(desc_mean) <- c("mean", "median", "SD", "sample_size")
desc_mean <- rownames_to_column(desc_mean, "population")

desc_sum <- as.data.frame(rbind(gba_sum, gba_sum_EA, sib_sum, sib_sum_EA))
colnames(desc_sum) <- c("mean_sum", "median_sum", "SD_sum", "sample_size_sum")
desc <- cbind(desc_mean, desc_sum)
desc <- cbind(desc, rbind(gba_table, gba_table_EA, sib_table, sib_table_EA) )

write.table(desc, "descriptive_table_costs_2209005.csv",
            row.names=F, quote=F, sep=" ")



