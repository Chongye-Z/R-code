# -----------------------------------------------------------------------------
# SCRIPT DESCRIPTION
# -----------------------------------------------------------------------------
# 
# This script performs how to generate high-CO2 events in the manuscript
# A high-CO2 event is defined relative to each participant's own CO2 distribution 
# e.g., above the 97.5th percentile for above 2 consecutive days
#
# -----------------------------------------------------------------------------
# DATA STRUCTURES
# -----------------------------------------------------------------------------
# This script relies on these datasets:
#
# 1. 'com_data_new' (loaded from 'com_data_new_example.csv')
#    The population data used in the study.
#    - yr (integer): year of visit
#    - longitude (numeric): longitude of addresses provided by participants
#    - latitude (Numeric): longitude of addresses provided by participants
#    - where (integer): 1 for Taiwan and 2 for Hong Kong
#
# 2. 'taiwan_CO2' (loaded from 'taiwan_CO2.xlsx')
#    Taiwan CO2 daily concentration data from 2003 to 2020 (wide format).
#    column V3335:V3679 are 13 grids of CO2 data covering the Taiwan island.
#    'long_lat_TW' (loaded from 'long_lat_TW.csv') linking 'taiwan_CO2' to 'com_data_new' dataset.
#    column longitude and latitude are unique values of 'com_data_new' in Taiwan
#    column grid_id are id of column V3335:V3679 in 'taiwan_CO2'
#
# 3. 'hk_CO2' (loaded from 'hk_CO2.xlsx')
#    HK CO2 daily concentration data from 2003 to 2020 (wide format).
#    column V3669:V3670 are 2 grids of CO2 data covering Hong Kong.
#    'long_lat_HK' (loaded from 'long_lat_HK.csv') linking 'hk_CO2' to 'com_data_new' dataset.
#    column longitude and latitude are unique values of 'com_data_new' in Hong Kong
#    column grid_id are id of column V3669:V3670 in 'hk_CO2'
#
# -----------------------------------------------------------------------------
# DERIVED VARIABLE NAMING CONVENTION AND OUTPUT DATA
# -----------------------------------------------------------------------------
# The final output contains 4 datasets: 
# 'CumE.csv','duration.csv','frequency.csv','intensity.csv'
# Corresponding to the specific metric calculated:
#     - CumE: Cumulative exposure, sum of CO2 concentrations across all high-CO2 event days.
#     - duration: Total number of days within all high-CO2 event days.
#     - frequency: Total number of distinct high-CO2 events.
#     - intensity: Average CO2 concentrations across all high-CO2 event days.
#     
# All output datasets have the same column names in patterns `CO2_over_`
#
# EXAMPLE: `CO2_over_p97.5_2D` in 'CumE.csv' means:
#   The cumulative exposure in high-CO2 events defined by daily CO2 
#   concentration exceeding the 97th percentile (_p97.5) for at least
#   2 consecutive days (_2D).
#   If `CO2_over_p97.5_2D` in 'duration.csv', that means the total number of days in
#   high-CO2 events defined by ... (the same as above)
#
# -----------------------------------------------------------------------------

library(tidyverse)

# population data
com_data_new <- data.table::fread("data/example_data.csv")

# CO2 data
{
  ## TW 
  {
    long_lat_TW <- read.csv("data/long_lat_TW.csv")
    
    taiwan_CO2 <- rio::import("data/taiwan_CO2.xlsx")
    ndays <- as.Date("2020-12-31") - as.Date("2003-01-01") + 1
    dates <- seq(from=as.Date("2003-01-01"), by=1, length.out=ndays)
    taiwan_CO2$date <- dates
    taiwan_CO2$year <- year(taiwan_CO2$date)
    taiwan_CO2$month <- month(taiwan_CO2$date)
  }
  
  ## HK
  {
    long_lat_HK <- read.csv("data/long_lat_HK.csv")
    
    hk_CO2 <- rio::import("data/hk_CO2.xlsx")
    ndays <- as.Date("2020-12-31") - as.Date("2003-01-01") + 1
    dates <- seq(from=as.Date("2003-01-01"), by=1, length.out=ndays)
    hk_CO2$date <- dates
    hk_CO2$year <- year(hk_CO2$date)
    hk_CO2$month <- month(hk_CO2$date)
  }
}


# label names
xname <- c("CO2_over_p90","CO2_over_p92.5","CO2_over_p95","CO2_over_p97.5")

## 4 type exposure generation (CumE, duration, frequency and intensity) 
## exposure consecutive 2/3/4 days
{
  ### gen_consecutive() is the function to identify consecutive days
  ### if CO2 level was larger than the percentile on someday, it would be 1 in 'highdata'
  ### rle() can compute the lengths of equal values in 'highdata', then can find 2/3/4 consecutive patterns
  ### inverse.rle() can add the pattern to 'highdata'
  ### columns H2DCO2/H3DCO2/H4DCO2 are high-CO2 event defined on 2/3/4 consecutive days
  ### columns H2DCO2_f/H3DCO2_f/H4DCO2_f are high-CO2 event frequency
  gen_consecutive <- function(data,threshold){
    highdata <- data |>  mutate(high = ifelse(value>threshold,yes = 1,no = 0))
    tempc <- rle(highdata$high)
    tempc2D <- tempc3D <- tempc4D <- tempc
    tempc2D$values <- ifelse((tempc$lengths>=2 & tempc$values == 1),yes = 1, no = 0)
    tempc3D$values <- ifelse((tempc$lengths>=3 & tempc$values == 1),yes = 1, no = 0)
    tempc4D$values <- ifelse((tempc$lengths>=4 & tempc$values == 1),yes = 1, no = 0)
    
    highdata$H2DCO2 <- inverse.rle(tempc2D)
    highdata$H3DCO2 <- inverse.rle(tempc3D)
    highdata$H4DCO2 <- inverse.rle(tempc4D)
    
    highdata$H2DCO2_f <- sum(tempc$lengths>=2 & tempc$values == 1)
    highdata$H3DCO2_f <- sum(tempc$lengths>=3 & tempc$values == 1)
    highdata$H4DCO2_f <- sum(tempc$lengths>=4 & tempc$values == 1)
    
    return(highdata)
  }
  
  ### CumE
  {
    ## matching CO2 data for Taiwan cohort ##
    {
      tempd <- com_data_new |>  filter(where == 1)
      tempd <- left_join(tempd, long_lat_TW)
      keydata1 <- tempd |>  select(where,yr,grid_id) |>  distinct()
      keydata1 <- as.data.frame(keydata1)
      na_matrix <- matrix(NA, nrow = nrow(keydata1), ncol = 12)
      colnames(na_matrix) <- paste0(xname,rep(c("_2D","_3D","_4D"),each=4))
      keydata1 <- cbind(keydata1, na_matrix)
      s1 <- Sys.time()
      for (i in 1:nrow(keydata1)) {
        tempyr <- keydata1[i,"yr"]
        
        tempdate1 <- paste0(tempyr,"-",12,"-31") |>  as.Date()
        tempdate2 <- paste0(tempyr-1,"-",1,"-01") |>  as.Date()
        
        tempCO2 <- taiwan_CO2[,c(keydata1$grid_id[i],14:16)] |>  
          filter(date >= tempdate2 & date <= tempdate1)
        colnames(tempCO2)[1] <- "value"
        
        temphigh <- gen_consecutive(data = tempCO2,threshold = quantile(tempCO2$value,0.9))
        keydata1$CO2_over_p90_2D[i] <- sum(temphigh$value[which(temphigh$H2DCO2 == 1)])
        keydata1$CO2_over_p90_3D[i] <- sum(temphigh$value[which(temphigh$H3DCO2 == 1)])
        keydata1$CO2_over_p90_4D[i] <- sum(temphigh$value[which(temphigh$H4DCO2 == 1)])
        
        temphigh <- gen_consecutive(data = tempCO2,threshold = quantile(tempCO2$value,0.925))
        keydata1$CO2_over_p92.5_2D[i] <- sum(temphigh$value[which(temphigh$H2DCO2 == 1)])
        keydata1$CO2_over_p92.5_3D[i] <- sum(temphigh$value[which(temphigh$H3DCO2 == 1)])
        keydata1$CO2_over_p92.5_4D[i] <- sum(temphigh$value[which(temphigh$H4DCO2 == 1)])
        
        temphigh <- gen_consecutive(data = tempCO2,threshold = quantile(tempCO2$value,0.95))
        keydata1$CO2_over_p95_2D[i] <- sum(temphigh$value[which(temphigh$H2DCO2 == 1)])
        keydata1$CO2_over_p95_3D[i] <- sum(temphigh$value[which(temphigh$H3DCO2 == 1)])
        keydata1$CO2_over_p95_4D[i] <- sum(temphigh$value[which(temphigh$H4DCO2 == 1)])
        
        temphigh <- gen_consecutive(data = tempCO2,threshold = quantile(tempCO2$value,0.975))
        keydata1$CO2_over_p97.5_2D[i] <- sum(temphigh$value[which(temphigh$H2DCO2 == 1)])
        keydata1$CO2_over_p97.5_3D[i] <- sum(temphigh$value[which(temphigh$H3DCO2 == 1)])
        keydata1$CO2_over_p97.5_4D[i] <- sum(temphigh$value[which(temphigh$H4DCO2 == 1)])
      }
      Sys.time() - s1
    }
    
    ## matching CO2 data for HK cohort ##
    {
      tempd <- com_data_new |>  filter(where == 2)
      tempd <- left_join(tempd, long_lat_HK)
      keydata2 <- tempd |>  select(where,yr,grid_id) |>  distinct()
      keydata2 <- as.data.frame(keydata2)
      na_matrix <- matrix(NA, nrow = nrow(keydata2), ncol = 12)
      colnames(na_matrix) <- paste0(xname,rep(c("_2D","_3D","_4D"),each=4))
      keydata2 <- cbind(keydata2, na_matrix)
      for (i in 1:nrow(keydata2)) {
        tempyr <- keydata2[i,"yr"]
        
        tempdate1 <- paste0(tempyr,"-",12,"-31") |>  as.Date()
        tempdate2 <- paste0(tempyr-1,"-",1,"-01") |>  as.Date()
        
        tempCO2 <- hk_CO2[,c(keydata2$grid_id[i],3:5)] |>  
          filter(date >= tempdate2 & date <= tempdate1)
        colnames(tempCO2)[1] <- "value"

        temphigh <- gen_consecutive(data = tempCO2,threshold = quantile(tempCO2$value,0.9))
        keydata2$CO2_over_p90_2D[i] <- sum(temphigh$value[which(temphigh$H2DCO2 == 1)])
        keydata2$CO2_over_p90_3D[i] <- sum(temphigh$value[which(temphigh$H3DCO2 == 1)])
        keydata2$CO2_over_p90_4D[i] <- sum(temphigh$value[which(temphigh$H4DCO2 == 1)])
        
        temphigh <- gen_consecutive(data = tempCO2,threshold = quantile(tempCO2$value,0.925))
        keydata2$CO2_over_p92.5_2D[i] <- sum(temphigh$value[which(temphigh$H2DCO2 == 1)])
        keydata2$CO2_over_p92.5_3D[i] <- sum(temphigh$value[which(temphigh$H3DCO2 == 1)])
        keydata2$CO2_over_p92.5_4D[i] <- sum(temphigh$value[which(temphigh$H4DCO2 == 1)])
        
        temphigh <- gen_consecutive(data = tempCO2,threshold = quantile(tempCO2$value,0.95))
        keydata2$CO2_over_p95_2D[i] <- sum(temphigh$value[which(temphigh$H2DCO2 == 1)])
        keydata2$CO2_over_p95_3D[i] <- sum(temphigh$value[which(temphigh$H3DCO2 == 1)])
        keydata2$CO2_over_p95_4D[i] <- sum(temphigh$value[which(temphigh$H4DCO2 == 1)])
        
        temphigh <- gen_consecutive(data = tempCO2,threshold = quantile(tempCO2$value,0.975))
        keydata2$CO2_over_p97.5_2D[i] <- sum(temphigh$value[which(temphigh$H2DCO2 == 1)])
        keydata2$CO2_over_p97.5_3D[i] <- sum(temphigh$value[which(temphigh$H3DCO2 == 1)])
        keydata2$CO2_over_p97.5_4D[i] <- sum(temphigh$value[which(temphigh$H4DCO2 == 1)])
      }
    }
    
    keydata <- rbind(keydata1,keydata2)
    write.csv(keydata,"data/CumE.csv")
  }
  ### duration
  {
    ## matching CO2 data for Taiwan cohort ##
    {
      tempd <- com_data_new |>  filter(where == 1)
      tempd <- left_join(tempd, long_lat_TW)
      keydata1 <- tempd |>  select(where,yr,grid_id) |>  distinct()
      keydata1 <- as.data.frame(keydata1)
      na_matrix <- matrix(NA, nrow = nrow(keydata1), ncol = 12)
      colnames(na_matrix) <- paste0(xname,rep(c("_2D","_3D","_4D"),each=4))
      keydata1 <- cbind(keydata1, na_matrix)
      s1 <- Sys.time()
      for (i in 1:nrow(keydata1)) {
        tempyr <- keydata1[i,"yr"]
        
        tempdate1 <- paste0(tempyr,"-",12,"-31") |>  as.Date()
        tempdate2 <- paste0(tempyr-1,"-",1,"-01") |>  as.Date()
        
        tempCO2 <- taiwan_CO2[,c(keydata1$grid_id[i],14:16)] |>  
          filter(date >= tempdate2 & date <= tempdate1)
        colnames(tempCO2)[1] <- "value"

        temphigh <- gen_consecutive(data = tempCO2,threshold = quantile(tempCO2$value,0.9))
        keydata1$CO2_over_p90_2D[i] <- sum(temphigh$H2DCO2 == 1)
        keydata1$CO2_over_p90_3D[i] <- sum(temphigh$H3DCO2 == 1)
        keydata1$CO2_over_p90_4D[i] <- sum(temphigh$H4DCO2 == 1)
        
        temphigh <- gen_consecutive(data = tempCO2,threshold = quantile(tempCO2$value,0.925))
        keydata1$CO2_over_p92.5_2D[i] <- sum(temphigh$H2DCO2 == 1)
        keydata1$CO2_over_p92.5_3D[i] <- sum(temphigh$H3DCO2 == 1)
        keydata1$CO2_over_p92.5_4D[i] <- sum(temphigh$H4DCO2 == 1)
        
        temphigh <- gen_consecutive(data = tempCO2,threshold = quantile(tempCO2$value,0.95))
        keydata1$CO2_over_p95_2D[i] <- sum(temphigh$H2DCO2 == 1)
        keydata1$CO2_over_p95_3D[i] <- sum(temphigh$H3DCO2 == 1)
        keydata1$CO2_over_p95_4D[i] <- sum(temphigh$H4DCO2 == 1)
        
        temphigh <- gen_consecutive(data = tempCO2,threshold = quantile(tempCO2$value,0.975))
        keydata1$CO2_over_p97.5_2D[i] <- sum(temphigh$H2DCO2 == 1)
        keydata1$CO2_over_p97.5_3D[i] <- sum(temphigh$H3DCO2 == 1)
        keydata1$CO2_over_p97.5_4D[i] <- sum(temphigh$H4DCO2 == 1)
      }
      Sys.time() - s1
    }
    
    ## matching CO2 data for HK cohort ##
    {
      tempd <- com_data_new |>  filter(where == 2)
      tempd <- left_join(tempd, long_lat_HK)
      keydata2 <- tempd |>  select(where,yr,grid_id) |>  distinct()
      keydata2 <- as.data.frame(keydata2)
      na_matrix <- matrix(NA, nrow = nrow(keydata2), ncol = 12)
      colnames(na_matrix) <- paste0(xname,rep(c("_2D","_3D","_4D"),each=4))
      keydata2 <- cbind(keydata2, na_matrix)
      for (i in 1:nrow(keydata2)) {
        tempyr <- keydata2[i,"yr"]
        
        tempdate1 <- paste0(tempyr,"-",12,"-31") |>  as.Date()
        tempdate2 <- paste0(tempyr-1,"-",1,"-01") |>  as.Date()
        
        tempCO2 <- hk_CO2[,c(keydata2$grid_id[i],3:5)] |>  
          filter(date >= tempdate2 & date <= tempdate1)
        colnames(tempCO2)[1] <- "value"

        temphigh <- gen_consecutive(data = tempCO2,threshold = quantile(tempCO2$value,0.9))
        keydata2$CO2_over_p90_2D[i] <- sum(temphigh$H2DCO2 == 1)
        keydata2$CO2_over_p90_3D[i] <- sum(temphigh$H3DCO2 == 1)
        keydata2$CO2_over_p90_4D[i] <- sum(temphigh$H4DCO2 == 1)
        
        temphigh <- gen_consecutive(data = tempCO2,threshold = quantile(tempCO2$value,0.925))
        keydata2$CO2_over_p92.5_2D[i] <- sum(temphigh$H2DCO2 == 1)
        keydata2$CO2_over_p92.5_3D[i] <- sum(temphigh$H3DCO2 == 1)
        keydata2$CO2_over_p92.5_4D[i] <- sum(temphigh$H4DCO2 == 1)
        
        temphigh <- gen_consecutive(data = tempCO2,threshold = quantile(tempCO2$value,0.95))
        keydata2$CO2_over_p95_2D[i] <- sum(temphigh$H2DCO2 == 1)
        keydata2$CO2_over_p95_3D[i] <- sum(temphigh$H3DCO2 == 1)
        keydata2$CO2_over_p95_4D[i] <- sum(temphigh$H4DCO2 == 1)
        
        temphigh <- gen_consecutive(data = tempCO2,threshold = quantile(tempCO2$value,0.975))
        keydata2$CO2_over_p97.5_2D[i] <- sum(temphigh$H2DCO2 == 1)
        keydata2$CO2_over_p97.5_3D[i] <- sum(temphigh$H3DCO2 == 1)
        keydata2$CO2_over_p97.5_4D[i] <- sum(temphigh$H4DCO2 == 1)
      }
    }
    
    keydata <- rbind(keydata1,keydata2)
    write.csv(keydata,"data/duration.csv")
  }
  ### frequency
  {
    ## matching CO2 data for Taiwan cohort ##
    {
      tempd <- com_data_new |>  filter(where == 1)
      tempd <- left_join(tempd, long_lat_TW)
      keydata1 <- tempd |>  select(where,yr,grid_id) |>  distinct()
      keydata1 <- as.data.frame(keydata1)
      na_matrix <- matrix(NA, nrow = nrow(keydata1), ncol = 12)
      colnames(na_matrix) <- paste0(xname,rep(c("_2D","_3D","_4D"),each=4))
      keydata1 <- cbind(keydata1, na_matrix)
      s1 <- Sys.time()
      for (i in 1:nrow(keydata1)) {
        tempyr <- keydata1[i,"yr"]
        
        tempdate1 <- paste0(tempyr,"-",12,"-31") |>  as.Date()
        tempdate2 <- paste0(tempyr-1,"-",1,"-01") |>  as.Date()
        
        tempCO2 <- taiwan_CO2[,c(keydata1$grid_id[i],14:16)] |>  
          filter(date >= tempdate2 & date <= tempdate1)
        colnames(tempCO2)[1] <- "value"

        temphigh <- gen_consecutive(data = tempCO2,threshold = quantile(tempCO2$value,0.9))
        keydata1$CO2_over_p90_2D[i] <- unique(temphigh$H2DCO2_f)
        keydata1$CO2_over_p90_3D[i] <- unique(temphigh$H3DCO2_f)
        keydata1$CO2_over_p90_4D[i] <- unique(temphigh$H4DCO2_f)
        
        temphigh <- gen_consecutive(data = tempCO2,threshold = quantile(tempCO2$value,0.925))
        keydata1$CO2_over_p92.5_2D[i] <- unique(temphigh$H2DCO2_f)
        keydata1$CO2_over_p92.5_3D[i] <- unique(temphigh$H3DCO2_f)
        keydata1$CO2_over_p92.5_4D[i] <- unique(temphigh$H4DCO2_f)
        
        temphigh <- gen_consecutive(data = tempCO2,threshold = quantile(tempCO2$value,0.95))
        keydata1$CO2_over_p95_2D[i] <- unique(temphigh$H2DCO2_f)
        keydata1$CO2_over_p95_3D[i] <- unique(temphigh$H3DCO2_f)
        keydata1$CO2_over_p95_4D[i] <- unique(temphigh$H4DCO2_f)
        
        temphigh <- gen_consecutive(data = tempCO2,threshold = quantile(tempCO2$value,0.975))
        keydata1$CO2_over_p97.5_2D[i] <- unique(temphigh$H2DCO2_f)
        keydata1$CO2_over_p97.5_3D[i] <- unique(temphigh$H3DCO2_f)
        keydata1$CO2_over_p97.5_4D[i] <- unique(temphigh$H4DCO2_f)
      }
      Sys.time() - s1
    }
    
    ## matching CO2 data for HK cohort ##
    {
      tempd <- com_data_new |>  filter(where == 2)
      tempd <- left_join(tempd, long_lat_HK)
      keydata2 <- tempd |>  select(where,yr,grid_id) |>  distinct()
      keydata2 <- as.data.frame(keydata2)
      na_matrix <- matrix(NA, nrow = nrow(keydata2), ncol = 12)
      colnames(na_matrix) <- paste0(xname,rep(c("_2D","_3D","_4D"),each=4))
      keydata2 <- cbind(keydata2, na_matrix)
      for (i in 1:nrow(keydata2)) {
        tempyr <- keydata2[i,"yr"]
        
        tempdate1 <- paste0(tempyr,"-",12,"-31") |>  as.Date()
        tempdate2 <- paste0(tempyr-1,"-",1,"-01") |>  as.Date()
        
        tempCO2 <- hk_CO2[,c(keydata2$grid_id[i],3:5)] |>  
          filter(date >= tempdate2 & date <= tempdate1)
        colnames(tempCO2)[1] <- "value"

        temphigh <- gen_consecutive(data = tempCO2,threshold = quantile(tempCO2$value,0.9))
        keydata2$CO2_over_p90_2D[i] <- unique(temphigh$H2DCO2_f)
        keydata2$CO2_over_p90_3D[i] <- unique(temphigh$H3DCO2_f)
        keydata2$CO2_over_p90_4D[i] <- unique(temphigh$H4DCO2_f)
        
        temphigh <- gen_consecutive(data = tempCO2,threshold = quantile(tempCO2$value,0.925))
        keydata2$CO2_over_p92.5_2D[i] <- unique(temphigh$H2DCO2_f)
        keydata2$CO2_over_p92.5_3D[i] <- unique(temphigh$H3DCO2_f)
        keydata2$CO2_over_p92.5_4D[i] <- unique(temphigh$H4DCO2_f)
        
        temphigh <- gen_consecutive(data = tempCO2,threshold = quantile(tempCO2$value,0.95))
        keydata2$CO2_over_p95_2D[i] <- unique(temphigh$H2DCO2_f)
        keydata2$CO2_over_p95_3D[i] <- unique(temphigh$H3DCO2_f)
        keydata2$CO2_over_p95_4D[i] <- unique(temphigh$H4DCO2_f)
        
        temphigh <- gen_consecutive(data = tempCO2,threshold = quantile(tempCO2$value,0.975))
        keydata2$CO2_over_p97.5_2D[i] <- unique(temphigh$H2DCO2_f)
        keydata2$CO2_over_p97.5_3D[i] <- unique(temphigh$H3DCO2_f)
        keydata2$CO2_over_p97.5_4D[i] <- unique(temphigh$H4DCO2_f)
      }
    }
    
    keydata <- rbind(keydata1,keydata2)
    write.csv(keydata,"data/frequency.csv")
  }
  ### intensity
  {
    ## matching CO2 data for Taiwan cohort ##
    {
      tempd <- com_data_new |>  filter(where == 1)
      tempd <- left_join(tempd, long_lat_TW)
      keydata1 <- tempd |>  select(where,yr,grid_id) |>  distinct()
      keydata1 <- as.data.frame(keydata1)
      na_matrix <- matrix(NA, nrow = nrow(keydata1), ncol = 12)
      colnames(na_matrix) <- paste0(xname,rep(c("_2D","_3D","_4D"),each=4))
      keydata1 <- cbind(keydata1, na_matrix)
      s1 <- Sys.time()
      for (i in 1:nrow(keydata1)) {
        tempyr <- keydata1[i,"yr"]
        
        tempdate1 <- paste0(tempyr,"-",12,"-31") |>  as.Date()
        tempdate2 <- paste0(tempyr-1,"-",1,"-01") |>  as.Date()
        
        tempCO2 <- taiwan_CO2[,c(keydata1$grid_id[i],14:16)] |>  
          filter(date >= tempdate2 & date <= tempdate1)
        colnames(tempCO2)[1] <- "value"
        
        temphigh <- gen_consecutive(data = tempCO2,threshold = quantile(tempCO2$value,0.9))
        keydata1$CO2_over_p90_2D[i] <- ifelse(sum(temphigh$H2DCO2)!=0,mean(temphigh$value[which(temphigh$H2DCO2 == 1)]),0)
        keydata1$CO2_over_p90_3D[i] <- ifelse(sum(temphigh$H3DCO2)!=0,mean(temphigh$value[which(temphigh$H3DCO2 == 1)]),0)
        keydata1$CO2_over_p90_4D[i] <- ifelse(sum(temphigh$H4DCO2)!=0,mean(temphigh$value[which(temphigh$H4DCO2 == 1)]),0)
        
        temphigh <- gen_consecutive(data = tempCO2,threshold = quantile(tempCO2$value,0.925))
        keydata1$CO2_over_p92.5_2D[i] <- ifelse(sum(temphigh$H2DCO2)!=0,mean(temphigh$value[which(temphigh$H2DCO2 == 1)]),0)
        keydata1$CO2_over_p92.5_3D[i] <- ifelse(sum(temphigh$H3DCO2)!=0,mean(temphigh$value[which(temphigh$H3DCO2 == 1)]),0)
        keydata1$CO2_over_p92.5_4D[i] <- ifelse(sum(temphigh$H4DCO2)!=0,mean(temphigh$value[which(temphigh$H4DCO2 == 1)]),0)
        
        temphigh <- gen_consecutive(data = tempCO2,threshold = quantile(tempCO2$value,0.95))
        keydata1$CO2_over_p95_2D[i] <- ifelse(sum(temphigh$H2DCO2)!=0,mean(temphigh$value[which(temphigh$H2DCO2 == 1)]),0)
        keydata1$CO2_over_p95_3D[i] <- ifelse(sum(temphigh$H3DCO2)!=0,mean(temphigh$value[which(temphigh$H3DCO2 == 1)]),0)
        keydata1$CO2_over_p95_4D[i] <- ifelse(sum(temphigh$H4DCO2)!=0,mean(temphigh$value[which(temphigh$H4DCO2 == 1)]),0)
        
        temphigh <- gen_consecutive(data = tempCO2,threshold = quantile(tempCO2$value,0.975))
        keydata1$CO2_over_p97.5_2D[i] <- ifelse(sum(temphigh$H2DCO2)!=0,mean(temphigh$value[which(temphigh$H2DCO2 == 1)]),0)
        keydata1$CO2_over_p97.5_3D[i] <- ifelse(sum(temphigh$H3DCO2)!=0,mean(temphigh$value[which(temphigh$H3DCO2 == 1)]),0)
        keydata1$CO2_over_p97.5_4D[i] <- ifelse(sum(temphigh$H4DCO2)!=0,mean(temphigh$value[which(temphigh$H4DCO2 == 1)]),0)
      }
      Sys.time() - s1
    }
    
    ## matching CO2 data for HK cohort ##
    {
      tempd <- com_data_new |>  filter(where == 2)
      tempd <- left_join(tempd, long_lat_HK)
      keydata2 <- tempd |>  select(where,yr,grid_id) |>  distinct()
      keydata2 <- as.data.frame(keydata2)
      na_matrix <- matrix(NA, nrow = nrow(keydata2), ncol = 12)
      colnames(na_matrix) <- paste0(xname,rep(c("_2D","_3D","_4D"),each=4))
      keydata2 <- cbind(keydata2, na_matrix)
      for (i in 1:nrow(keydata2)) {
        tempyr <- keydata2[i,"yr"]
        
        tempdate1 <- paste0(tempyr,"-",12,"-31") |>  as.Date()
        tempdate2 <- paste0(tempyr-1,"-",1,"-01") |>  as.Date()
        
        tempCO2 <- hk_CO2[,c(keydata2$grid_id[i],3:5)] |>  
          filter(date >= tempdate2 & date <= tempdate1)
        colnames(tempCO2)[1] <- "value"
        
        temphigh <- gen_consecutive(data = tempCO2,threshold = quantile(tempCO2$value,0.9))
        keydata2$CO2_over_p90_2D[i] <- ifelse(sum(temphigh$H2DCO2)!=0,mean(temphigh$value[which(temphigh$H2DCO2 == 1)]),0)
        keydata2$CO2_over_p90_3D[i] <- ifelse(sum(temphigh$H3DCO2)!=0,mean(temphigh$value[which(temphigh$H3DCO2 == 1)]),0)
        keydata2$CO2_over_p90_4D[i] <- ifelse(sum(temphigh$H4DCO2)!=0,mean(temphigh$value[which(temphigh$H4DCO2 == 1)]),0)
        
        temphigh <- gen_consecutive(data = tempCO2,threshold = quantile(tempCO2$value,0.925))
        keydata2$CO2_over_p92.5_2D[i] <- ifelse(sum(temphigh$H2DCO2)!=0,mean(temphigh$value[which(temphigh$H2DCO2 == 1)]),0)
        keydata2$CO2_over_p92.5_3D[i] <- ifelse(sum(temphigh$H3DCO2)!=0,mean(temphigh$value[which(temphigh$H3DCO2 == 1)]),0)
        keydata2$CO2_over_p92.5_4D[i] <- ifelse(sum(temphigh$H4DCO2)!=0,mean(temphigh$value[which(temphigh$H4DCO2 == 1)]),0)
        
        temphigh <- gen_consecutive(data = tempCO2,threshold = quantile(tempCO2$value,0.95))
        keydata2$CO2_over_p95_2D[i] <- ifelse(sum(temphigh$H2DCO2)!=0,mean(temphigh$value[which(temphigh$H2DCO2 == 1)]),0)
        keydata2$CO2_over_p95_3D[i] <- ifelse(sum(temphigh$H3DCO2)!=0,mean(temphigh$value[which(temphigh$H3DCO2 == 1)]),0)
        keydata2$CO2_over_p95_4D[i] <- ifelse(sum(temphigh$H4DCO2)!=0,mean(temphigh$value[which(temphigh$H4DCO2 == 1)]),0)
        
        temphigh <- gen_consecutive(data = tempCO2,threshold = quantile(tempCO2$value,0.975))
        keydata2$CO2_over_p97.5_2D[i] <- ifelse(sum(temphigh$H2DCO2)!=0,mean(temphigh$value[which(temphigh$H2DCO2 == 1)]),0)
        keydata2$CO2_over_p97.5_3D[i] <- ifelse(sum(temphigh$H3DCO2)!=0,mean(temphigh$value[which(temphigh$H3DCO2 == 1)]),0)
        keydata2$CO2_over_p97.5_4D[i] <- ifelse(sum(temphigh$H4DCO2)!=0,mean(temphigh$value[which(temphigh$H4DCO2 == 1)]),0)
      }
    }
    
    keydata <- rbind(keydata1,keydata2)
    write.csv(keydata,"data/intensity.csv")
  }
}

