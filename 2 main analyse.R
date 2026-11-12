# -----------------------------------------------------------------------------
# SCRIPT DESCRIPTION
# -----------------------------------------------------------------------------
# 
# This script is designed to assess the association between various high-CO2 
# event exposure definitions and lung function using 
# linear mixed-effects model (LMM) and mixed-effects logistic regression model (MELRM).
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
# PREDICTOR VARIABLES:
#   - High-CO2 event exposures: A series of variables like `CO2_over_p97.5_2D`.
#   - Covariates: `gender`,`age`,`bmi`,`n_disease`,`edu`,`smoke1`,`drink1`,
#                 `fruit_vege`,`sport1`,`yr`,`season`,`rural_ubran`,`work`
# -----------------------------------------------------------------------------


library(tidyverse)
library(lme4)
library(lmerTest)


# population data
com_data_new <- data.table::fread("data/example_data.csv")

# label names
type <- c("CumE", "duration", "frequency", "intensity")
yname <- c("pf_fvc","pf_fev1","pf_mmf","ffratio",
           "Impaired_fev1","Impaired_fvc","Impaired_ratio","Abnormal_1")
ylabel <- c("FVC (ml)","FEV1 (ml)","FEV1/FVC (%)","MMEF (ml/s)",
           "Impaired FVC", "Impaired FEV1","Impaired FEV1/FVC","Abnormal pattern")
xname <- c("CO2_over_p90","CO2_over_p92.5","CO2_over_p95","CO2_over_p97.5")

# covariate formulas
cov1 <- "+gender+age+bmi+n_disease+yr+season+(1|pid)"
cov2 <- "+gender+age+bmi+n_disease+yr+season+edu+smoke1+drink1+fruit_vege+sport1+(1|pid)"
cov3 <- "+gender+age+bmi+n_disease+edu+smoke1+drink1+fruit_vege+sport1+yr+season+rural_ubran+work+(1|pid)"


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


# ------------------------------ Main analysis ---------------------------------
#### linear model (one IQR increase)
{
  library(future)
  library(future.apply)
  
  xname_new <- paste0(xname,rep(c("_2D","_3D","_4D"),each=4))
  para_comb <- expand.grid(x=xname_new,y=yname,type=type)
  para_comb <- para_comb |>  mutate(across(.cols = x:type,.fns = as.character))
  
  plan(multisession, workers = 10)
  s1 <- Sys.time()
  out_list <- future_lapply(
    X = seq_len(nrow(para_comb)),
    FUN = function(i) {
      
      tempd <- read.csv(paste0("data/",para_comb$type[i],".csv"))
      
      tempdata <- left_join(com_data_new, long_lat_CO2)
      tempdata <- left_join(tempdata, tempd)
      tempdata <- as.data.frame(tempdata)
      
      iq <- IQR(tempdata[, para_comb$x[i]], na.rm = TRUE)
      
      ## if no one interquartile range is zero then return NULL
      if (is.na(iq) || iq == 0) {return(NULL)}
      
      tempdata$xiqr <- tempdata[, para_comb$x[i]] / iq # Scale exposure by IQR
      tempformula <- paste0(para_comb$y[i], "~ xiqr", cov3) |>  as.formula()
      
      if (para_comb$y[i] %in% c("pf_fvc", "pf_fev1", "pf_mmf", "ffratio")) {
        tempmodel <- lmer(tempformula, data = tempdata)
        mtype <- 1
      } else {
        tempmodel <- glmer(tempformula,data = tempdata,family = binomial,nAGQ = 0)
        mtype <- 2
      }
      
      tempresult <- cbind(
        data.frame(
          outcome = para_comb$y[i],
          exposure = para_comb$x[i],
          type = para_comb$type[i],
          AIC = AIC(tempmodel),
          BIC = BIC(tempmodel)
        ),
        coeffun(tempmodel, row = 2, mtype = mtype)
      )
      
      tempresult
    },
    future.seed = TRUE,
    future.packages = c("dplyr", "lme4", "lmerTest")
  )
  
  out_iqr <- do.call(rbind, out_list)
  
  s2 <- Sys.time()
  s2-s1 ### Time difference of 1.052402 hours
  
  write.csv(out_iqr,"result/IQR (model 3).csv", row.names = F)
  
  ## stop parallel after running
  plan(sequential)
}


## the optional threshold was 97.5 (2 days)
### adjust cov on different model
### effect estimate (cov1 and cov2)
{
  para_comb <- expand.grid(x="CO2_over_p97.5_2D",y=yname, covarate = c(1:2))
  para_comb <- para_comb |>  mutate(across(.cols = x:y,.fns = as.character))
  
  out_iqr <- data.frame(outcome=NULL,exposure=NULL,covarate=NULL,AIC=NULL,BIC=NULL,
                        Estimate=NULL,low=NULL,high=NULL,P=NULL,effect=NULL)
  
  CumE_data <- read.csv("data/CumE.csv")
  tempdata <- left_join(com_data_new,long_lat_CO2)
  tempdata <- left_join(tempdata, CumE_data)
  tempdata <- as.data.frame(tempdata)
  
  s1 <- Sys.time()
  for (i in 1:nrow(para_comb)) {
    iq <- IQR(tempdata[,para_comb$x[i]], na.rm = TRUE) # Calculate IQR 
    if(iq == 0){next}
    tempdata$xiqr <- tempdata[,para_comb$x[i]] / iq # Scale exposure by IQR
    
    if (para_comb$covarate[i] == 1) {
      tempformula <- paste0(para_comb$y[i],"~ xiqr",cov1) |>  as.formula()
    }else{
      tempformula <- paste0(para_comb$y[i],"~ xiqr",cov2) |>  as.formula()
    }
    
    if(para_comb$y[i] %in% c("pf_fvc","pf_fev1","pf_mmf","ffratio")){
      tempmodel <- lmer(tempformula,data = tempdata)
      mtype <- 1
    }else{
      tempmodel <- glmer(tempformula,data = tempdata,family = binomial,nAGQ = 0)
      mtype <- 2
    }
    
    tempresult <- cbind(data.frame(
      outcome = para_comb$y[i], exposure = para_comb$x[i],
      covarate = para_comb$covarate[i],
      AIC = AIC(tempmodel),BIC = BIC(tempmodel)), 
      coeffun(tempmodel, row = 2, mtype = mtype))
    
    out_iqr <- rbind(out_iqr,tempresult)
    
    print(i)
  }
  s2 <- Sys.time()
  s2 - s1
  ## Time difference of 7.902658 mins
  
  rio::export(out_iqr,"result/cov1and2_CumE_iqr.xlsx")
  
}
# Figure 1: effect plot 97.5 two day
{
  tempres <- read.csv("result/IQR (model 3).csv")
  plotdata <- tempres |> filter(exposure == "CO2_over_p97.5_2D" & type == "CumE")
  cov_12 <- rio::import("result/cov1and2_CumE_iqr.xlsx")
  plotdata$covarate <- 3
  plotdata <- rbind(plotdata |>  select(-type), cov_12)
  
  plotdata$yname <- ylabel
  plotdata$yname <- factor(plotdata$yname,levels = ylabel)
  
  plotdata <- plotdata |>  
    mutate(Model = case_when(
      covarate == 1 ~ "Model 1",
      covarate == 2 ~ "Model 2",
      covarate == 3 ~ "Model 3"
    ))
  
  
  plotdata |> 
    ggplot(aes(x = yname, y = Estimate,colour = Model))+
    theme_test()+
    scale_color_manual(values = c("#42B540FF","#00468BFF","#ED0000FF"))+
    geom_errorbar(aes(ymin = low, ymax = high),linewidth = 0.8,
                  position = position_dodge(0.3),width=0.15)+
    geom_point(position = position_dodge(0.3),shape = 18,size = 2)+
    coord_cartesian(clip = "off") +
    facet_wrap(~yname, scales = "free", nrow = 2)+
    labs(y = "Effect (95% CI)") +
    theme(text = element_text(family = "serif"),
          axis.title.y = element_text(size = 12,face = "bold"),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_blank(),
          legend.position = "bottom",
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          strip.background = element_blank(),
          strip.text = element_text(size = 10,face = "bold"))
  
  ggsave("result/effect plot.tiff",
         width = 10,height = 4,dpi = 350,compression = "lzw")
}

