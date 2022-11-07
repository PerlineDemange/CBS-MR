# Script to identify siblings in our subset datasets 
# Author: Perline Demange
# Project: CBS-MR EA-MH 
# 

# 1. Clean data #############

##* 1.1 Open the datasets gbapersoon and kindouder ------------

# Open Gbapersoontab
# GBAPERSOON2019TABV1.csv was saved as cvs in spss from the file GBAPERSOON2019TABV1.sav
gba <- fread("../Data/GBA/GBAPERSOON2019TABV1.csv")
head(gba) #N=25489305

# Open Kindoudertab
# KINDOUDER2019TABV2.csv was saved as csv from spss from the file KINDOUDER2019TABV2.sav
kindouder <- fread("../Data/GBA/KINDOUDER2019TABV2.csv")
head(kindouder) # N=16979273


##* 1.2 Select only people born between 1965 and 1985 from gbapersoontab ----------

gba6585 <- gba[which(gba$GBAGEBOORTEJAAR >= 1965 & gba$GBAGEBOORTEJAAR <= 1985), ]
summary(gba6585) #N=6539767
#hist(gba6585$GBAGEBOORTEJAAR) # give incorrect breaks rather do: 
barplot(table(gba6585$GBAGEBOORTEJAAR))
rm(gba)

##** 1.2.1 Save GBA 1965-1985 for later use-----

write.csv(gba6585, "gba_6585_20201015.csv", row.names = F)

##* 1.3 Remove individual with unknown parents (parents not in GBA ) from kindoudertab -----

# R is in BGA, G is not in GBA 
kindpa <- kindouder[which(kindouder$RINPERSOONSpa == 'R' & kindouder$RINPERSOONSMa == "R"),] # N=14073108 (-3millions)

##* 1.4 Remove individuals not in GBA (R is in GBA, D is dead at birth and N is unknown) from kindoudertab -----------

kindclean <- kindpa[which(kindpa$RINPERSOONS == 'R'),] # N=13855944
head(kindclean)
rm(kindpa)
rm(kindouder)

##* 1.5 Merge selected gbapersoon and kindouder -----

data <- merge(gba6585, kindclean, by=c('RINPERSOONS', "RINPERSOON"))
head(data) #N = 4063765
barplot(table(data$GBAGEBOORTEJAAR))


# 2. Identify family with several children #########

##* 2.1  Identify individuals with same parents based on RINPERSOON -----

data$FID <- paste(data$RINPERSOONpa, data$RINPERSOONMa, sep="_")
data <- data[order(data$FID),]
head(data)

##* 2.2 Get number of full siblings in family: Fsize ------

uniquefid <- unique(data$FID) #all different FID #2154525 different full family 
count <- as.data.frame(table(data$FID)) #get number of occurences of FID
head(count)
data$Fsize <- count[[2]][match(data$FID, count[[1]])] # match data fid with value in count, if true takes value of Freq 
table(data$Fsize)
tail(data)

##* 2.3 Remove families with only one child -----

data <- data[data$Fsize > 1] #3282626

##* 2.4 Get number of individual with same father SamePa or same mother SameMa ------

nbfather <- as.data.frame(table(data$RINPERSOONpa))
data$SamePa <- nbfather[[2]][match(data$RINPERSOONpa, nbfather[[1]])] 
head(data)
hist(data$SamePa)
hist(data$Fsize)
nbmother <- as.data.frame(table(data$RINPERSOONMa))
data$SameMa <- nbmother[[2]][match(data$RINPERSOONMa, nbmother[[1]])]
hist(data$SameMa)

##* 2.5 Identify individuals who have half siblings ------

data$halfPa <- data$Fsize < data$SamePa # TRUE if the father had children with other mother, halfsibling exist
summary(data$halfPa) # true= 21755
head(data[data$halfPa == T])

data$halfMa <- data$Fsize < data$SameMa # TRUE if the mother had children with other father, halfsibling exist
summary(data$halfMa) #true= 6978

data$halfsib <- data$halfPa == T | data$halfMa == T
summary(data$halfsib) # T= 27725
head(data[data$halfsib == T])

##* 2.6 Remove families whose parents had more than one child with other (identical) partner ---

data <- data[data$halfsib == F]
summary(data) #N=3254901

##* 2.7 Exclude homoparental families ----

data$GBAGESLACHTMOEDER  <- as.factor(data$GBAGESLACHTMOEDER)
data$GBAGESLACHTVADER  <- as.factor(data$GBAGESLACHTVADER)
summary(data$GBAGESLACHTMOEDER) #sex of the mother: 1 is man and 2 is woman and - missing data
summary(data$GBAGESLACHTVADER) #sex of the father: 1 is man and 2 is woman and - missing data

data <- data[data$GBAGESLACHTMOEDER == 2] #N=3254848 #keep individual whose mother is a woman 
data <- data[data$GBAGESLACHTVADER == 1] #N=3254608 #keep individual whose father is a man 

rm(kindclean)
rm(nbfather)
rm(nbmother)

##* 2.8 Exclude families with mismatching parental information #####
# Some parents have the same RINDPERSON but appear to have different birthday and country of birth 

###** 2.8.1 Exclude parents within families with different birthday -----
####*** 2.8.1.1 Fathers ----
table(data$GBAGEBOORTEJAARVADER) #1814 missing data, range 1873 to 1970
table(data$GBAGEBOORTEMAANDVADER) # 1814 missing data 
data$GBAGEBOORTEJAARVADER <- as.numeric(data$GBAGEBOORTEJAARVADER)
table(data$GBAGEBOORTEJAARVADER) # doesn't print NA if numeric 
data <- data[data$GBAGEBOORTEJAARVADER > 1870,] #3252794
table(data$GBAGEBOORTEMAANDVADER) # we removed all missing data

data$birthdayPa <- paste(data$GBAGEBOORTEJAARVADER, data$GBAGEBOORTEMAANDVADER, data$GBAGEBOORTEDAGVADER, sep ="-")

data_birthdayPa <- data %>% 
  group_by(FID) %>%
  mutate(SameBirthdayPa = ifelse(duplicated(birthdayPa)| duplicated(birthdayPa, fromLast = TRUE), T, F))
data_birthdayPa <- as.data.frame(data_birthdayPa)
summary(data_birthdayPa) #5368 fathers with a different birthday indicated
tail(data_birthdayPa[data_birthdayPa$SameBirthdayPa == F,]) # Xkoppelnumber C6 is coming often 

# several are mismatched by only one month, maybe should only sort on year? 
# either i remove all family with a mismatch father or just the mismatch father (keep fam of >2 in the data)
# I decide to do: birthday has to match entirely, only mismatched father is excluded 

data_birthdayPa <- data_birthdayPa[data_birthdayPa$SameBirthdayPa == T, ] #3247426

# Get new number of family size and remove families < 2 (should not happen)
count <- as.data.frame(table(data_birthdayPa$FID)) #get number of occurences of FID
data_birthdayPa$Fsize_final <- count[[2]][match(data_birthdayPa$FID, count[[1]])] # match data fid with value in count, if true takes value of Freq 
data_birthdayPa <- data_birthdayPa[data_birthdayPa$Fsize_final > 1,] #3247426

data <- data_birthdayPa #3247426

####*** 2.8.1.2 Mothers ----------
table(data$GBAGEBOORTEJAARMOEDER) #495 missing data from 1911 to 1971 
table(data$GBAGEBOORTEMAANDMOEDER) #495 missing data
data$GBAGEBOORTEJAARMOEDER<- as.numeric(data$GBAGEBOORTEJAARMOEDER)
data <- data[which(is.na(data$GBAGEBOORTEJAARMOEDER) == F),] #3246931

data$birthdayMa <- paste(data$GBAGEBOORTEJAARMOEDER, data$GBAGEBOORTEMAANDMOEDER, data$GBAGEBOORTEDAGMOEDER, sep ="-")

data_birthdayMa <- data %>% 
  group_by(FID) %>%
  mutate(SameBirthdayMa = ifelse(duplicated(birthdayMa)| duplicated(birthdayMa, fromLast = TRUE), T, F))
data_birthdayMa <- as.data.frame(data_birthdayMa)
summary(data_birthdayMa) # 4066 mothers with a different birthday indicated

data_birthdayMa <- data_birthdayMa[data_birthdayMa$SameBirthdayMa == T, ] #3242865

# Get new number of family size and remove families < 2 (should not happen)
count <- as.data.frame(table(data_birthdayMa$FID)) #get number of occurences of FID
data_birthdayMa$Fsize_final <- count[[2]][match(data_birthdayMa$FID, count[[1]])] # match data fid with value in count, if true takes value of Freq 
data_birthdayMa <- data_birthdayMa[data_birthdayMa$Fsize_final > 1,] #3242865

data <- data_birthdayMa #3242865

##** 2.8.2 Exclude parents with mismatching land of origin within same family####
###*** 2.8.2.1 Fathers ---- 
table(data$GBAGEBOORTELANDVADER)
table(is.na(data$GBAGEBOORTELANDVADER)) # no missing data 

data_landPa <- data %>% 
  group_by(FID) %>%
  mutate(SameLandPa = ifelse(duplicated(GBAGEBOORTELANDVADER)| duplicated(GBAGEBOORTELANDVADER, fromLast = TRUE), T, F))
data_landPa <- as.data.frame(data_landPa) 
summary(data_landPa) #4328 fathers with different lands of origin 

data_landPa <- data_landPa[data_landPa$SameLandPa == T, ] #3238537

# Get new number of family size and remove families < 2 (should not happen)
count <- as.data.frame(table(data_landPa$FID)) #get number of occurences of FID
data_landPa$Fsize_final <- count[[2]][match(data_landPa$FID, count[[1]])] # match data fid with value in count, if true takes value of Freq 
data_landPa <- data_landPa[data_landPa$Fsize_final > 1,] #3238537

data <- data_landPa #3238537

###*** 2.8.2.2 Mothers -----
table(data$GBAGEBOORTELANDMOEDER)
table(is.na(data$GBAGEBOORTELANDMOEDER)) # no missing data 

data_landMa <- data %>% 
  group_by(FID) %>%
  mutate(SameLandMa = ifelse(duplicated(GBAGEBOORTELANDMOEDER)| duplicated(GBAGEBOORTELANDMOEDER, fromLast = TRUE), T, F))
data_landMa <- as.data.frame(data_landMa) 
summary(data_landMa) # 1222 mothers with different lands of origin 

data_landMa <- data_landMa[data_landMa$SameLandMa == T, ] #3237315

# Get new number of family size and remove families < 2 (should not happen)
count <- as.data.frame(table(data_landMa$FID)) #get number of occurences of FID
data_landMa$Fsize_final <- count[[2]][match(data_landMa$FID, count[[1]])] # match data fid with value in count, if true takes value of Freq 
data_landMa <- data_landMa[data_landMa$Fsize_final > 1,] #3237315

data <- data_landMa #3237315

rm(data_birthdayPa)
rm(data_birthdayMa)
rm(data_landPa)
rm(data_landMa)

##* 2.9 Exclude siblings born less 6months apart, but keep twins----
###** 2.9.1 Identify twins (born on the same day) ----
data$birthday <- paste(data$GBAGEBOORTEJAAR, data$GBAGEBOORTEMAAND, data$GBAGEBOORTEDAG, sep ="-")
barplot(table(data$GBAGEBOORTEMAAND))
barplot(table(gba6585$GBAGEBOORTEDAG)) # Day of birth is set at 1 by security. 

data_twin <- data %>% 
  group_by(FID) %>%
  mutate(multiple = ifelse(duplicated(birthday)| duplicated(birthday, fromLast = TRUE), T, F))
data_twin <- as.data.frame(data_twin)
summary(data_twin$multiple) #76672 multiples 


# Get number of multiples within family to investigate problems
summary <- data_twin %>%
  group_by(FID) %>%
  tally(multiple == T)
table(summary$n)

#I checked families with 7-8 multiples, nothing looks weird. 

# Save data
write.csv(data_twin, "gba_6585_FID_fullsibship_twin_20201016.csv", row.names = F)
#data <- fread("gba_6585_FID_fullsibship_twin_20201016.csv")
data <- data_twin

###** 2.9.2 Remove siblings born less 6months apart----
library(lubridate)
data$birthday <- ymd(data$birthday)


# Function returning T if a pair of sib is born less than 6 months apart wihtin the family family
is_lower_6_months <- function(x,y){ 
  person1 <- data[FID == family, ][x,]
  person2 <- data[FID == family, ][y,]
  ifelse(abs(person1$birthday - person2$birthday) < months(7), T, F)
}
# use months(7) and not months(6) or person born on exactly at 6 months return F, while we want T 

# Run a loop on all families looking at every siblings combinations 
# This loop takes more than a day to run! 
families  <- unique(data$FID)
closer_than_6_months <- NULL 
for (family in families) { 
  number_family_members <- nrow(data[FID == family, ])
  birthday_occurrence <- matrix(0, nrow=number_family_members, ncol=number_family_members)
  siblings <- data[FID == family, ]$RINPERSOON
  possible_combi <- combn(length(siblings), 2)
  for (combi in 1:ncol(possible_combi)) { 
    i <- possible_combi[1, combi]
    j <- possible_combi[2, combi]
    birthday_occurrence[i,j] <- is_lower_6_months(i,j)
    birthday_occurrence[j,i] <- birthday_occurrence[i,j]
  }
  for (individual in 1:ncol(birthday_occurrence)) {
    closer_than_6_months <- c(closer_than_6_months, ifelse(sum(birthday_occurrence[individual,]) > 0, T, F))
  }
}
data$closer_than_6months <- closer_than_6_months

write.csv(data, "gba_6585_FID_fullsibship_twin_6months_20201030.csv", row.names = F)

summary(data$closer_than_6months) #T 78755
# 76672 multiples so 2083 to exclude

###** 2.9.3 Exclude siblings born less than 6 months apart and keep twins ----

data$closer_than_6months_not_multiple <- ifelse(data$closer_than_6months == T & data$multiple == F, T, F)
summary(data$closer_than_6months_not_multiple) #2083
head(data[data$closer_than_6months_not_multiple == T,])

data_final <- data[data$closer_than_6months_not_multiple == F, ] #3235232

#Exclude families with only one sib left 

count <- as.data.frame(table(data_final$FID)) #get number of occurences of FID
data_final$Fsize_final <- count[[2]][match(data_final$FID, count[[1]])] # match data fid with value in count, if true takes value of Freq 
data_final <- data_final[data_final$Fsize_final > 1,] #3234923

write.csv(data_final, "gba_6585_FID_fullsibship_final_20201101.csv", row.names = F)

## *3 Add birth order variable #####
data_final <- fread("gba_6585_FID_fullsibship_final_20201101.csv")
head(data_final)
data_final$birthday <- ymd(data_final$birthday)

data_birth_order <- data_final %>% 
  group_by(FID) %>%
  mutate(birth_order = dense_rank(birthday))

data_birth_order <- as.data.frame(data_birth_order)


# Save data  All siblings available #############

summary(data_final)

uniquefid <- unique(data_final$FID) #1 355 165

write.csv(data_birth_order, "gba_6585_FID_fullsibship_final_20210215.csv", row.names = F)





