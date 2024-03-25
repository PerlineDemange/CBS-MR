# Set up libraries
# Author: Perline Demange
# Project: CBS-MR EA-MH 
# 
  
library(data.table)
library(tidyverse)
set.seed(42)
library(lmerTest) #instead of lme4, to be able to get df 
library(lubridate)
library(dplyr)
library("xlsx")
library(nlme)
library(sandwich)
library(lmtest)

#read function 03 
source("03_0_functions.R")

