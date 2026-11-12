# -----------------------------------------------------------------------------
# SCRIPT DESCRIPTION
# -----------------------------------------------------------------------------
# 
# This script is designed to check effect modification, aiming to determine 
# if the effect of CumE to high-CO2 event on the outcome differs across
# various subgroups of the population.
#
# The process:
# 1.  Interaction Testing: tests for statistical interaction between 
#     CumE and potential effect modifiers using LMMs or MELRMs.
# 2.  Data Stratification: Based on the potential modifiers,
#     it creates clean, stratified subsets of the data.
# 3.  Stratified Analysis: It runs the main regression models separately within
#     each stratum.
# 4.  Visualization (Forest Plots): The results of the stratified analyses
#     are visualized using forest plots to compare the magnitude and 
#     significance of the effects across different subgroups.

# -----------------------------------------------------------------------------
# DATA STRUCTURES
# -----------------------------------------------------------------------------
# The script loads a primary dataset 'com_data_new'.
#
# OUTCOME VARIABLE:
#   - pf_fvc (Numeric): forced vital capacity (FVC)
#   - pf_fev1 (Numeric): forced expiratory volume in 1 second (FEV1)
#   - pf_mmf (Numeric): maximum mid-expiratory flow (MMEF) 
#   - ffratio (Numeric): the FEV1/FVC ratio
#   - Impaired_fev1 (binary): Impaired FEV1 (1: FEV1 < LLN or 0: not)
#   - Impaired_fvc (binary): Impaired FVC (1: FVC < LLN or 0: not)
#   - Impaired_ratio (binary): Impaired FEV1/FVC (1: FEV1/FVC < LLN or 0: not)
#   - Abnormal_1 (binary): 1: FEV1/FVC ≥ LLN and FVC < LLN or 0: not 
#
# IDENTIFIERS:
#   - pid (Character/Numeric): A unique participant identifier, used for the
#     random intercept `(1|pid)` in the mixed-effects models.
#
# FORMULA:
#   - tempsubf1: remove `gender`
#   - tempsubf2: remove `age`
#   - tempsubf3: remove `rural_ubran`
# 
# -----------------------------------------------------------------------------


library(tidyverse)
library(lme4)
library(lmerTest)


# population data
com_data_new <- data.table::fread("data/example_data.csv")

# label names
yname <- c("pf_fvc","pf_fev1","pf_mmf","ffratio",
           "Impaired_fev1","Impaired_fvc","Impaired_ratio","Abnormal_1")
ylabel <- c("FVC","FEV1","MMEF","FEV1/FVC",
            "Impaired FEV1","Impaired FVC","Impaired FEV1/FVC","Abnormal pattern")

## function to extract model results
coeffun <- function(model, row, mtype) {
  if(mtype == 1){
    coef <- as.data.frame(summary(model)$coefficients)
    coef$CI_low <- coef$Estimate - 1.96 * coef$`Std. Error`
    coef$CI_high <- coef$Estimate + 1.96 * coef$`Std. Error`
    coef <- coef[, c("Estimate", "CI_low", "CI_high", "Pr(>|t|)")]
    colnames(coef) <- c("Estimate", "low", "high", "P")
    coef$effect <- paste0(sprintf("%.3f", coef$Estimate), " (", sprintf("%.3f", coef$low), ",", sprintf("%.3f", coef$high), ")")
    coef$effect <- ifelse(coef$P < 0.05, paste0(coef$effect, "*"), coef$effect)
    coef <- coef[row,]
  }else{
    coef <- as.data.frame(summary(model)$coefficients)
    coef$CI_low <- exp(coef$Estimate - 1.96 * coef$`Std. Error`)
    coef$CI_high <- exp(coef$Estimate + 1.96 * coef$`Std. Error`)
    coef$Estimate <- exp(coef$Estimate)
    coef <- coef[, c("Estimate", "CI_low", "CI_high", "Pr(>|z|)")]
    colnames(coef) <- c("Estimate", "low", "high", "P")
    colnames(coef) <- c("Estimate", "low", "high", "P")
    coef$effect <- paste0(sprintf("%.3f", coef$Estimate), " (", sprintf("%.3f", coef$low), ",", sprintf("%.3f", coef$high), ")")
    coef$effect <- ifelse(coef$P < 0.05, paste0(coef$effect, "*"), coef$effect)
    coef <- coef[row,]
  }
  return(coef)
}

## location id of CO2 data
long_lat_CO2 <- read.csv("data/long_lat (TW and HK).csv")


# ---------------------------- subgroup analysis -------------------------------
subg <- c("gender_g","age_g","resident_g")

tempsubf1 <- "+age+bmi+n_disease+edu+smoke1+drink1+fruit_vege+sport1+yr+season+rural_ubran+work+(1|pid)"
tempsubf2 <- "+gender+bmi+n_disease+edu+smoke1+drink1+fruit_vege+sport1+yr+season+rural_ubran+work+(1|pid)"
tempsubf3 <- "+gender+age+bmi+n_disease+edu+smoke1+drink1+fruit_vege+sport1+yr+season+work+(1|pid)"


## CumE
CumE_data <- read.csv("data/CumE.csv")
tempdata <- left_join(com_data_new,long_lat_CO2)
tempdata <- left_join(tempdata, CumE_data)
tempdata <- as.data.frame(tempdata)

tempdata <- tempdata |>  
  mutate(gender_g = gender,
         age_g = ifelse(age > 60, yes = 1, no = 2),
         resident_g = ifelse(rural_ubran == 0,yes = 1,no = 2))
iq <- IQR(tempdata[,"CO2_over_p97.5_2D"], na.rm = TRUE)
tempdata$xiqr <- tempdata[,"CO2_over_p97.5_2D"] / iq # Scale exposure by IQR


# loop
{
  para_comb <- expand.grid(subg=subg,y=yname)
  para_comb$subg <- as.character(para_comb$subg)
  out <- data.frame(outcome=NULL,subg=NULL,group=NULL,p_inter=NULL,
                    Estimate=NULL,low=NULL,high=NULL,P=NULL,effect=NULL)
  # Loop through each outcome and subgroup
  s1 <- Sys.time()
  for (i in 1:nrow(para_comb)) {
    tempsubf <- case_when(
      para_comb$subg[i] == "gender_g" ~ tempsubf1,
      para_comb$subg[i] == "age_g" ~ tempsubf2,
      para_comb$subg[i] == "resident_g" ~ tempsubf3
    )
    
    tempf <- paste(para_comb$y[i],"~ xiqr", tempsubf)
    
    inetrvar <- paste0(para_comb$subg[i],":xiqr")
    tempf_inter <- paste0(para_comb$y[i],"~",inetrvar,"+xiqr+",
                          para_comb$subg[i], tempsubf)
    
    if(para_comb$y[i] %in% c("pf_fvc","pf_fev1","pf_mmf","ffratio")){
      tempm1 <- lmer(tempf, data = tempdata |>  filter(.data[[para_comb$subg[i]]] == 1))
      tempm2 <- lmer(tempf, data = tempdata |>  filter(.data[[para_comb$subg[i]]] == 2))
      tempm3 <- lmer(tempf_inter, data = tempdata)
      p_inter <- parameters::p_value(tempm3)$p[length(parameters::p_value(tempm3)$p)]
      res1 <- cbind(data.frame(outcome=para_comb$y[i],subg=para_comb$subg[i],group=1,p_inter=p_inter), 
                    coeffun(model = tempm1, row = 2, mtype = 1))
      res2 <- cbind(data.frame(outcome=para_comb$y[i],subg=para_comb$subg[i],group=2,p_inter=p_inter), 
                    coeffun(model = tempm2,row = 2, mtype = 1))
    }else{
      tempm1 <- glmer(tempf,family=binomial,nAGQ=0, data=tempdata |>  filter(.data[[para_comb$subg[i]]] == 1))
      tempm2 <- glmer(tempf,family=binomial,nAGQ=0, data=tempdata |>  filter(.data[[para_comb$subg[i]]] == 2))
      tempm3 <- glmer(tempf_inter,family=binomial,nAGQ=0,data = tempdata)
      p_inter <- parameters::p_value(tempm3)$p[length(parameters::p_value(tempm3)$p)]
      res1 <- cbind(data.frame(outcome=para_comb$y[i],subg=para_comb$subg[i],group=1,p_inter=p_inter), 
                    coeffun(model = tempm1, row = 2, mtype = 2))
      res2 <- cbind(data.frame(outcome=para_comb$y[i],subg=para_comb$subg[i],group=2,p_inter=p_inter), 
                    coeffun(model = tempm2,row = 2, mtype = 2))
    }
    
    out <- rbind(out,res1,res2)
    print(i)
  }
  s2 <- Sys.time()
  s2 - s1
  rio::export(out,"result/subg_CumE.xlsx")
  
  ## Time difference of 1.516835 hours
}


## plot
{
  plotdata <- rio::import("result/subg_CumE.xlsx")
  
  plotdata <- plotdata |>  
    mutate(sub = case_when(
      subg == "gender_g" & group == 1 ~ "Male",
      subg == "gender_g" & group == 2 ~ "Female",
      subg == "age_g" & group == 1 ~ "> 60 years",
      subg == "age_g" & group == 2 ~ "≤ 60 years",
      subg == "resident_g" & group == 1 ~ "Rural",
      subg == "resident_g" & group == 2 ~ "Urban",
    ))
  sub_name <- c("Rural","Urban","> 60 years","≤ 60 years","Male","Female")
  
  
  plotdata$xmin <- rep(rev(seq(from = 0.7, to = 4.7, by = 2)), each = 2)
  plotdata$xmax <- rep(rev(seq(from = 2.3, to = 6.3, by = 2)), each = 2)
  plotdata$sub <- factor(plotdata$sub,levels = sub_name)
  plotdata$outcome <- factor(plotdata$outcome,levels = yname,labels = ylabel)
  plotdata$ref <- ifelse(plotdata$outcome %in% c("FVC","FEV1","MMEF","FEV1/FVC"),
                         yes = 0, no = 1)
  plotdata$sig <- ifelse(plotdata$p_inter < 0.1, yes = "yes", no = "no")
  
  
  plotdata |>  
    ggplot(aes(x = sub, y = Estimate))+
    theme_test()+ 
    geom_rect(aes(xmin = xmin,xmax = xmax,ymin = -Inf,ymax = Inf,fill = sig), 
              alpha = 0.5)+
    scale_fill_manual(values = c("white","grey95"))+
    geom_point()+
    geom_errorbar(aes(ymin = low, ymax = high), width=.2)+
    geom_hline(aes(yintercept = ref), linetype = 2)+
    coord_flip()+
    labs(y = "CumE effect estimate (95% CI)")+
    facet_wrap(~outcome,ncol=2,scales = "free") +
    theme(text = element_text(family = "serif"),
          title = element_text(size = 12,face = "bold"),
          axis.title = element_blank(),
          axis.text = element_text(size = 12),
          legend.position = "none",
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(size = 12,face = "bold"))
  
  ggsave("result/sub effect (3 groups).tiff",
         width = 12,height = 15,dpi = 350,compression = "lzw")
}

