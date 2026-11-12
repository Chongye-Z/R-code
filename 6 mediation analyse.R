# -----------------------------------------------------------------------------
# SCRIPT DESCRIPTION
# -----------------------------------------------------------------------------
# 
# This script is designed to evaluate whether 
# 2-year average ambient temperature mediated the associations between 
# high-CO2 events and respiratory health.
#
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
# -----------------------------------------------------------------------------


library(tidyverse)
library(lme4)
library(lmerTest)
library(mediation)
library(RMediation)

# CO2 grid
long_lat_CO2 <- read.csv("grid_id_CO2 (TW and HK).csv")
# data
CumE_data <- read.csv("data/CumE.csv")
tempdata <- left_join(com_data_new,long_lat_CO2)
tempdata <- left_join(tempdata, CumE_data)
tempdata <- as.data.frame(tempdata)

# iqr increase
iq <- IQR(tempdata[, "CO2_over_p97.5_2D"])
tempdata$xiqr <- tempdata[, "CO2_over_p97.5_2D"] / iq


# label names
yname <- c("pf_fvc","pf_fev1","pf_mmf","ffratio",
           "Impaired_fev1","Impaired_fvc","Impaired_ratio","Abnormal_1")
ylabel <- c("FVC","FEV1","MMEF","FEV1/FVC",
            "Impaired FEV1","Impaired FVC","Impaired FEV1/FVC","Abnormal pattern")

# covariates
cov3 <- "+gender+age+bmi+n_disease+edu+smoke1+drink1+fruit_vege+sport1+yr+season+rural_ubran+work+(1|pid)"

# function to extract model results
coeffun <- function(model, row, mtype) {
  if(mtype == 1){
    coef <- as.data.frame(summary(model)$coefficients)
    coef$CI_low <- coef$Estimate - 1.96 * coef$`Std. Error`
    coef$CI_high <- coef$Estimate + 1.96 * coef$`Std. Error`
    coef <- coef[, c("Estimate","Std. Error","CI_low", "CI_high", "Pr(>|t|)")]
    colnames(coef) <- c("Estimate", "SE", "low", "high", "P")
    coef$effect <- paste0(sprintf("%.3f", coef$Estimate), " (", sprintf("%.3f", coef$low), ",", sprintf("%.3f", coef$high), ")")
    coef$effect <- ifelse(coef$P < 0.05, paste0(coef$effect, "*"), coef$effect)
    coef <- coef[row,]
  }else{
    coef <- as.data.frame(summary(model)$coefficients)
    coef$CI_low <- exp(coef$Estimate - 1.96 * coef$`Std. Error`)
    coef$CI_high <- exp(coef$Estimate + 1.96 * coef$`Std. Error`)
    coef$Estimate <- exp(coef$Estimate)
    coef <- coef[, c("Estimate","Std. Error","CI_low", "CI_high", "Pr(>|z|)")]
    colnames(coef) <- c("Estimate", "SE", "low", "high", "P")
    coef$effect <- paste0(sprintf("%.3f", coef$Estimate), " (", sprintf("%.3f", coef$low), ",", sprintf("%.3f", coef$high), ")")
    coef$effect <- ifelse(coef$P < 0.05, paste0(coef$effect, "*"), coef$effect)
    coef <- coef[row,]
  }
  return(coef)
}


# ------------------------- effect estimation ----------------------------------

library(future)
library(future.apply)

## mediator: 2 year mean temperature
{
  ### path c direct effect CumE to outcome ("outcome_res")
  ### path b effect temperature to outcome ("M_Y_res")
  plan(multisession, workers = 10)
  s1 <- Sys.time()
  out_list <- future_lapply(
    X = 1:length(yname),
    FUN = function(i) {

      outcome_f <- paste0(yname[i], "~ xiqr","+temperature_2y", cov3) %>% as.formula()
      if (yname[i] %in% c("pf_fvc", "pf_fev1", "pf_mmf", "ffratio")) {
        outcome_model <- lmer(outcome_f, data = tempdata)
        mtype <- 1
      } else {
        outcome_model <- glmer(outcome_f,data = tempdata,family = binomial,nAGQ = 0)
        mtype <- 2
      }
      outcome_res <- coeffun(outcome_model, row = 2, mtype = mtype)
      colnames(outcome_res) <- paste0(colnames(outcome_res),"_X_Y")
      
      M_Y_res <- coeffun(outcome_model, row = 3, mtype = mtype)
      colnames(M_Y_res) <- paste0(colnames(M_Y_res),"_M_Y")
      
      tempresult <- cbind(M_Y_res,outcome_res)
      
      tempresult
    },
    future.seed = TRUE
  )
  ## stop parallel after running
  plan(sequential)
  
  out_iqr <- do.call(rbind, out_list)
  
  ### path a effect CumE to temperature ("X_M_res")
  X_M_f <- paste0("temperature_2y", "~ xiqr", cov3) %>% as.formula()
  X_M_model <- lmer(X_M_f, data = tempdata)
  X_M_res <- coeffun(X_M_model, row = 2, mtype = 1)
  colnames(X_M_res) <- paste0(colnames(X_M_res),"_X_M")
  
  res_all <- cbind(out_iqr,X_M_res)
  res_all$outcome <- yname
  
  ### Monte Carlo approach to compute the CI 
  res_all_with_CI <- res_all
  set.seed(888)
  for (i in 1:nrow(res_all)) {
    if(res_all$outcome[i] %in% c("pf_fvc", "pf_fev1", "pf_mmf", "ffratio")){
      tempci <- medci(
        mu.x = res_all$Estimate_X_M[i],
        mu.y = res_all$Estimate_M_Y[i],
        se.x = res_all$SE_X_M[i],
        se.y = res_all$SE_M_Y[i],
        rho = 0,
        type = "MC",
        n.mc = 1e+05)
      res_all_with_CI$ab[i] <- res_all$Estimate_X_M[i] * res_all$Estimate_M_Y[i]
      res_all_with_CI$ab_low[i] <- tempci$`95% CI`[1]
      res_all_with_CI$ab_high[i] <- tempci$`95% CI`[2]
      res_all_with_CI$popM[i] <- res_all_with_CI$ab[i] / 
        (res_all_with_CI$Estimate_X_Y[i] + res_all_with_CI$ab[i])
    }else{
      tempci <- medci(
        mu.x = res_all$Estimate_X_M[i],
        mu.y = log(res_all$Estimate_M_Y[i]),
        se.x = res_all$SE_X_M[i],
        se.y = res_all$SE_M_Y[i],
        rho = 0,
        type = "MC",
        n.mc = 1e+05)
      res_all_with_CI$ab[i] <- res_all$Estimate_X_M[i] * log(res_all$Estimate_M_Y[i])
      res_all_with_CI$ab_low[i] <- tempci$`95% CI`[1]
      res_all_with_CI$ab_high[i] <- tempci$`95% CI`[2]
      res_all_with_CI$popM[i] <- res_all_with_CI$ab[i] / 
        (log(res_all_with_CI$Estimate_X_Y[i]) + res_all_with_CI$ab[i])
      ### this mediated proportions were calculated in the logit scale
    }
  }
  res_all_with_CI <- res_all_with_CI %>% 
    mutate(ab_exp = exp(ab),
           ab_low_exp = exp(ab_low),
           ab_high_exp = exp(ab_high))
  res_all_with_CI[1:4,c("ab_exp","ab_low_exp","ab_high_exp")] <- NA
  
  rio::export(res_all_with_CI,"result/mediator_temperature.xlsx", row.names = F)
}
