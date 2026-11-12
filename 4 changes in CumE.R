# -----------------------------------------------------------------------------
# SCRIPT DESCRIPTION
# -----------------------------------------------------------------------------
# 
# This script is designed to assess health impacts of 
# changes in cumulative CO2 exposure
#
# -----------------------------------------------------------------------------
# DATA STRUCTURES
# -----------------------------------------------------------------------------
# The script loads a primary dataset 'com_data_new'.
#
# We divided the participants into two groups based on 
# the overall change in CumE between baseline and the final visit: 
# an increased exposure group (> 0 ppm) and 
# a stable or decreased exposure group (≤ 0 ppm)
#
# We then matched the two groups using generalised propensity scores with 
# 1:1 nearest-neighbour matching and adjustment for all covariates in the baseline
#
#   - Covariates: `gender`,`age`,`bmi`,`n_disease`,`edu`,`smoke1`,`drink1`,
#                 `fruit_vege`,`sport1`,`yr`,`season`,`rural_ubran`,`work`
# -----------------------------------------------------------------------------


library(tidyverse)
library(lme4)
library(MatchIt)


## data
{
  long_lat_CO2 <- read.csv("data/long_lat (TW and HK).csv")
  
  com_data_new <- data.table::fread("data/example_data.csv")
  CumE_data <- read.csv("high_exposure/consecutive day for one event/CumE.csv")
  tempdata <- left_join(com_data_new,long_lat_CO2)
  tempdata <- left_join(tempdata, CumE_data)
  tempdata <- as.data.frame(tempdata)
}


## data processing (CumE change)
{
  move_data_cume <- tempdata |>  
    group_by(pid) |>  
    mutate(n = length(pid)) |>  
    ungroup() |>  
    filter(n>1)
  
  move_data_cume <- move_data_cume |>  
    mutate(time = yr * 12 + mon) |>  
    group_by(pid) |>  
    arrange(time, .by_group = T) |>  
    mutate(first_end = case_when(
      time == first(time) ~ "first",
      time == last(time) ~ "last",
      .default = "other"
    )) |>  
    ungroup()
  
  move_data_higher_cume <- move_data_cume |>  
    filter(first_end != "other")
  move_data_higher_cume <- move_data_higher_cume |>  
    group_by(pid) |>  
    mutate(first = CO2_over_p97.5_2D[first_end == "first"], 
           last = CO2_over_p97.5_2D[first_end == "last"])
  ## 116,152 participants 232,304 observation
  
  move_data_higher_new_cume <- move_data_higher_cume |>  
    group_by(pid) |>  
    mutate(Cume_change = last - first) |>  
    ungroup() |>  
    mutate(group = ifelse(Cume_change > 0,1,0)) |>  
    dplyr::select(pid,pf_fvc:pf_mmf,ffratio,Impaired_fev1:Abnormal_1,
                  gender,age,edu:sport1,fruit_vege,n_disease,yr,season,rural_ubran,work,
                  first_end,Cume_change,group) |>  
    distinct() |>  
    pivot_wider(names_from = first_end,
                values_from = c(pf_fvc:pf_mmf,ffratio,Impaired_fev1:Abnormal_1,
                                gender,age,edu:sport1,fruit_vege,n_disease,yr,season,rural_ubran,work))
  
  ## quantile(move_data_higher_new_cume$Cume_change)
  ##      0%        25%        50%        75%       100% 
  ## -4481.7932  -799.5917   190.3581  1129.8858  6929.3789 
}
### rename the covariate for following presentation of plots
{
  move_data_higher_new_cume <- move_data_higher_new_cume |> 
    mutate(Gender = gender_first,
           Age = age_first,
           Education = as.factor(edu_first),
           BMI = bmi_first,
           Smoke = as.factor(smoke1_first),
           Drink = as.factor(drink1_first),
           `Fruit vegetable intake` = fruit_vege_first,
           `Physical activity level` = as.factor(sport1_first),
           `Number of existing diseases` = as.factor(n_disease_first),
           `Year of visit` = yr_first,
           `Season of visit` = as.factor(season_first),
           `Address area` = rural_ubran_first,
           Occupation = work_first)
}

# --------------------------------- 1:1 match ---------------------------------
matching <- matchit(
  formula = group ~
    Gender+Age+Education+BMI+
    Smoke+Drink+`Fruit vegetable intake`+`Physical activity level`+
    `Number of existing diseases`+`Year of visit`+`Season of visit`+
    `Address area`+ Occupation, 
  method = "nearest",ratio=1, caliper = 0.2, replace = T, 
  distance='logit',distance.options=list(size=5),
  data = move_data_higher_new_cume)
dta_m1_cume <- get_matches(matching, data = move_data_higher_new_cume)

## Supplementary Figure 5: cov balance after matching
{
  cobalt::love.plot(
    matching,
    stats = "mean.diffs",
    binary = "std",
    abs = TRUE,
    un = TRUE,
    thresholds = c(m = 0.015),
    var.order = "unadjusted",
    line = TRUE,
    drop.distance = TRUE,
    sample.names = c("Unadjusted", "Adjusted"),
    colors = c("#fc8d58", "#4474b3"),
    shapes = c(16, 16),
    size = 2,
    title = "Covariate Balance") +
    labs(x = "Absolute Standardized Mean Differences",y = NULL) +
    theme_bw() +
    theme(text = element_text(family = "serif"),
          panel.grid = element_blank(),
          plot.title = element_text(size = 12,hjust = 0.5),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 11),
          legend.position = "bottom")

  ggsave("result/cov balance after matching.tiff",
         width = 10,height = 8,dpi = 350,compression = "lzw")
}

# --------------------------- plot for CumE change -----------------------------

plotdata <- move_data_higher_new_cume |>  
  group_by(group) |>  
  summarise(across(.cols = c(pf_fvc_last,pf_fev1_last,pf_mmf_last,
                             ffratio_last,Impaired_fev1_last,
                             Impaired_fvc_last,Impaired_ratio_last,
                             Abnormal_1_last),.fns = mean))
plotdata$group_name <- factor(plotdata$group,levels = c(0,1),
                              labels = c("Not increased","Increased"))
plotdata <- plotdata |>  
  pivot_longer(cols = pf_fvc_last:Abnormal_1_last,names_to = "outcome",
               values_to = "value")
plotdata$ylabel <- rep(c("FVC","FEV1","MMEF","FEV1/FVC",
                         "Impaired FEV1","Impaired FVC",
                         "Impaired FEV1/FVC","Abnormal pattern"),2)
plotdata <- plotdata |>  
  mutate(value = ifelse(!ylabel %in% c("FVC","FEV1","MMEF","FEV1/FVC"),
                        yes = value * 100, no =  value))

### effect estimation
est_res <- data.frame(yname = c("pf_fvc_last","pf_fev1_last",
                                "ffratio_last","pf_mmf_last",
                                "Impaired_fvc_last","Impaired_fev1_last",
                                "Impaired_ratio_last","Abnormal_1_last"),
                      ylabel = c("FVC","FEV1","FEV1/FVC","MMEF",
                                 "Impaired FVC","Impaired FEV1",
                                 "Impaired FEV1/FVC","Abnormal pattern"),
                      est = NA)
s1 <- Sys.time()
for (i in 1:nrow(est_res)) {
  tempformula <- as.formula(paste0(est_res$yname[i],
  "~ group +
  Gender+Age+Education+BMI+
  Smoke+Drink+`Fruit vegetable intake`+`Physical activity level`+
  `Number of existing diseases`+`Year of visit`+`Season of visit`+
  `Address area`+ Occupation + (1|subclass)"))
  if(est_res$ylabel[i] %in% c("FVC","FEV1","MMEF","FEV1/FVC")){
    tempm <- lmer(tempformula,data=dta_m1_cume,weights=weights)
    a <- summary(tempm)$coefficients
    e <- a["group","Estimate"] |>  round(2)
    l <- (e - 1.96*a["group","Std. Error"]) |>  round(2)
    u <- (e + 1.96*a["group","Std. Error"]) |>  round(2)
    est <- paste0(e," (",l,",",u,")")
  }else{
    tempm <- glmer(tempformula,data=dta_m1_cume,weights=weights,
                   family = binomial,nAGQ = 0)
    a <- summary(tempm)$coefficients
    e <- a["group","Estimate"] |>  exp() |>  round(2)
    l <- exp(a["group","Estimate"] - 1.96*a["group","Std. Error"]) |>  round(2)
    u <- exp(a["group","Estimate"] + 1.96*a["group","Std. Error"]) |>  round(2)
    est <- paste0(e," (",l,",",u,")")
  }
  
  est_res$est[i] <- est
}
s2 <- Sys.time()
s2-s1 ### Time difference of 2.318213 mins
rio::export(est_res, "result/Figure 3 estimate results.xlsx")


plotlist <- list()
for (i in 1:nrow(est_res)) {
  if(est_res$ylabel[i] %in% c("FVC","FEV1","MMEF")){
    ymin <- 2600
    ymax <- 4000
    ytitle <- paste0(est_res$ylabel[i]," (ml)")
    est_text <- paste0("Estimate coefficient: \n",est_res$est[i])
  }else if(est_res$ylabel[i] == "FEV1/FVC"){
    ymin <- 89
    ymax <- 90
    ytitle <- paste0(est_res$ylabel[i]," (%)")
    est_text <- paste0("Estimate coefficient: \n",est_res$est[i])
  }else if(est_res$ylabel[i] == "Abnormal pattern"){
    ymin <- 0
    ymax <- 24
    ytitle <- "Abnormal proportion (%)"
    est_text <- paste0("Estimate OR: \n",est_res$est[i])
  }else{
    ymin <- 0
    ymax <- 24
    ytitle <- "Impaired proportion (%)"
    est_text <- paste0("Estimate OR: \n",est_res$est[i])
  }
  
  tempplot <- plotdata |>  
    filter(ylabel == est_res$ylabel[i]) |>  
    ggplot(aes(x = group,y = value))+
    theme_classic()+
    geom_segment(aes(xend = group, y = ymin, yend = value), 
                 color = "#C9CACA",linewidth = 1.5) +
    geom_line(color = "black", linetype = "dotted")+
    geom_point(aes(colour = group_name),size = 5)+
    scale_color_manual(values = c("#42B540FF","#ED0000FF"))+
    scale_x_continuous(breaks = c(0,1),
                       labels = c("Not increased","Increased"),
                       limits = c(-0.25,1.25))+ 
    geom_text(x = 0.5, y = ymax - (ymax-ymin)*0.05,label = est_text,
              family = "serif",fontface = "plain")+
    coord_cartesian(ylim = c(ymin,ymax))+
    labs(y = ytitle,title = est_res$ylabel[i])+
    theme(text = element_text(family = "serif"),
          plot.title = element_text(size = 14,face = "bold",hjust = 0.5),
          axis.title.y = element_text(size = 14,face = "bold"),
          axis.title.x = element_blank(),
          axis.text = element_text(size = 12),
          legend.position = "none")
  
  plotlist[[i]] <- tempplot
}
cowplot::plot_grid(plotlist = plotlist,nrow = 2)
ggsave("result/Figure 3.tiff",
       width = 13,height = 5.5,dpi = 350,compression = "lzw")


# ----------------------------- sensitivity analyse ----------------------------

sen_data <- move_data_higher_cume |>  
  group_by(pid) |>  
  mutate(temptime = yr*12+mon,
         timediff = temptime[first_end == "last"] - temptime[first_end == "first"])

sen_data$timediff |>  quantile() 
##  0%  25%  50%  75% 100% 
##  1   25   50   92  215


move_data_higher_new_cume_sen <- sen_data |>  
  filter(timediff >= 24) |>  ## 2 years (24 months)
  mutate(Cume_change = last - first) |>  
  ungroup() |>  
  mutate(group = ifelse(Cume_change > 0,1,0)) |>  
  dplyr::select(pid,pf_fvc:pf_mmf,ffratio,Impaired_fev1:Abnormal_1,
                gender,age,edu:sport1,fruit_vege,n_disease,yr,season,rural_ubran,work,
                first_end,group) |>  
  distinct() |>  
  pivot_wider(names_from = first_end,
              values_from = c(pf_fvc:pf_mmf,ffratio,Impaired_fev1:Abnormal_1,
                              gender,age,edu:sport1,fruit_vege,n_disease,yr,season,rural_ubran,work))

move_data_higher_new_cume_sen |>  nrow()
## 92,826 participants

### rename the covariate 
{
  move_data_higher_new_cume_sen <- move_data_higher_new_cume_sen |> 
    mutate(Gender = gender_first,
           Age = age_first,
           Education = as.factor(edu_first),
           BMI = bmi_first,
           Smoke = as.factor(smoke1_first),
           Drink = as.factor(drink1_first),
           `Fruit vegetable intake` = fruit_vege_first,
           `Physical activity level` = as.factor(sport1_first),
           `Number of existing diseases` = as.factor(n_disease_first),
           `Year of visit` = yr_first,
           `Season of visit` = as.factor(season_first),
           `Address area` = rural_ubran_first,
           Occupation = work_first)
}


### 1:1 match
matching <- matchit(
  formula = group ~
    Gender+Age+Education+BMI+
    Smoke+Drink+`Fruit vegetable intake`+`Physical activity level`+
    `Number of existing diseases`+`Year of visit`+`Season of visit`+
    `Address area`+ Occupation, 
  method = "nearest",ratio=1, caliper = 0.2, replace = T, 
  distance='logit',distance.options=list(size=5),
  data = move_data_higher_new_cume_sen)
dta_m1_cume <- get_matches(matching, data = move_data_higher_new_cume_sen)


est_res <- data.frame(yname = c("pf_fvc_last","pf_fev1_last",
                                "ffratio_last","pf_mmf_last",
                                "Impaired_fvc_last","Impaired_fev1_last",
                                "Impaired_ratio_last","Abnormal_1_last"),
                      ylabel = c("FVC","FEV1","FEV1/FVC","MMEF",
                                 "Impaired FVC","Impaired FEV1",
                                 "Impaired FEV1/FVC","Abnormal pattern"),
                      est = NA)
s1 <- Sys.time()
for (i in 1:nrow(est_res)) {
  tempformula <- as.formula(paste0(est_res$yname[i],"~ group +
                gender_first+age_first+edu_first+bmi_first+
                smoke1_first+drink1_first+fruit_vege_first+sport1_first+
                n_disease_first+season_first+rural_ubran_first+work_first + 
                (1|subclass)"))
  if(est_res$ylabel[i] %in% c("FVC","FEV1","MMEF","FEV1/FVC")){
    tempm <- lmer(tempformula,data=dta_m1_cume,weights=weights)
    a <- summary(tempm)$coefficients
    e <- a["group","Estimate"] |>  round(2)
    l <- (e - 1.96*a["group","Std. Error"]) |>  round(2)
    u <- (e + 1.96*a["group","Std. Error"]) |>  round(2)
    est <- paste0(e," (",l,",",u,")")
  }else{
    tempm <- glmer(tempformula,data=dta_m1_cume,weights=weights,
                   family = binomial,nAGQ = 0)
    a <- summary(tempm)$coefficients
    e <- a["group","Estimate"] |>  exp() |>  round(2)
    l <- exp(a["group","Estimate"] - 1.96*a["group","Std. Error"]) |>  round(2)
    u <- exp(a["group","Estimate"] + 1.96*a["group","Std. Error"]) |>  round(2)
    est <- paste0(e," (",l,",",u,")")
  }
  
  est_res$est[i] <- est
}
s2 <- Sys.time()
s2-s1  ### Time difference of 1.444795 mins 
rio::export(est_res, "sensitivity analyses/estimate result of sensitivity analyse for Figure 3.xlsx")


plotdata <- move_data_higher_new_cume_sen |>  
  group_by(group) |>  
  summarise(across(.cols = c(pf_fvc_last,pf_fev1_last,pf_mmf_last,
                             ffratio_last,Impaired_fev1_last,
                             Impaired_fvc_last,Impaired_ratio_last,
                             Abnormal_1_last),.fns = mean))
plotdata$group_name <- factor(plotdata$group,levels = c(0,1),
                              labels = c("Not increased","Increased"))
plotdata <- plotdata |>  
  pivot_longer(cols = pf_fvc_last:Abnormal_1_last,names_to = "outcome",
               values_to = "value")
plotdata$ylabel <- rep(c("FVC","FEV1","MMEF","FEV1/FVC",
                         "Impaired FEV1","Impaired FVC",
                         "Impaired FEV1/FVC","Abnormal pattern"),2)
plotdata <- plotdata |>  
  mutate(value = ifelse(!ylabel %in% c("FVC","FEV1","MMEF","FEV1/FVC"),
                        yes = value * 100, no =  value))


plotlist <- list()
for (i in 1:nrow(est_res)) {
  if(est_res$ylabel[i] %in% c("FVC","FEV1","MMEF")){
    ymin <- 2600
    ymax <- 4000
    ytitle <- paste0(est_res$ylabel[i]," (ml)")
    est_text <- paste0("Estimate coefficient: \n",est_res$est[i])
  }else if(est_res$ylabel[i] == "FEV1/FVC"){
    ymin <- 89
    ymax <- 90.05
    ytitle <- paste0(est_res$ylabel[i]," (%)")
    est_text <- paste0("Estimate coefficient: \n",est_res$est[i])
  }else if(est_res$ylabel[i] == "Abnormal pattern"){
    ymin <- 0
    ymax <- 24
    ytitle <- "Abnormal proportion (%)"
    est_text <- paste0("Estimate OR: \n",est_res$est[i])
  }else{
    ymin <- 0
    ymax <- 24
    ytitle <- "Impaired proportion (%)"
    est_text <- paste0("Estimate OR: \n",est_res$est[i])
  }
  
  tempplot <- plotdata |>  
    filter(ylabel == est_res$ylabel[i]) |>  
    ggplot(aes(x = group,y = value))+
    theme_classic()+
    geom_segment(aes(xend = group, y = ymin, yend = value), 
                 color = "#C9CACA",linewidth = 1.5) +
    geom_line(color = "black", linetype = "dotted")+
    geom_point(aes(colour = group_name),size = 5)+
    scale_color_manual(values = c("#42B540FF","#ED0000FF"))+
    scale_x_continuous(breaks = c(0,1),
                       labels = c("Not increased","Increased"),
                       limits = c(-0.25,1.25))+ 
    geom_text(x = 0.5, y = ymax - (ymax-ymin)*0.05,label = est_text,
              family = "serif",fontface = "plain")+
    coord_cartesian(ylim = c(ymin,ymax))+
    labs(y = ytitle,title = est_res$ylabel[i])+
    theme(text = element_text(family = "serif"),
          plot.title = element_text(size = 14,face = "bold",hjust = 0.5),
          axis.title.y = element_text(size = 14,face = "bold"),
          axis.title.x = element_blank(),
          axis.text = element_text(size = 12),
          legend.position = "none")
  
  plotlist[[i]] <- tempplot
}
cowplot::plot_grid(plotlist = plotlist,nrow = 2)
ggsave("sensitivity analyses/SF6 (sensitivity analyse for Figure 3).tiff",
       width = 13,height = 5.5,dpi = 350,compression = "lzw")


