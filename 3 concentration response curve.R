# -----------------------------------------------------------------------------
# SCRIPT DESCRIPTION
# -----------------------------------------------------------------------------
# This script conducts the analysis of the concentration-response
# relationship between CumE and lung function.
#
# Linearity Assessment: 
#     fit both a linear model and a non-linear model (using a natural cubic spline).
#     Then use a Likelihood Ratio Test (LRT) to determine if the non-linear
#     fit is statistically superior.
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
library(dlnm)


# population data
com_data_new <- data.table::fread("data/example_data.csv")

# label names
type <- c("CumE", "duration", "frequency", "intensity")
yname <- c("pf_fvc","pf_fev1","pf_mmf","ffratio",
           "Impaired_fev1","Impaired_fvc","Impaired_ratio","Abnormal_1")
ylabel <- c("FVC (ml)","FEV1 (ml)","MMEF (ml/s)","FEV1/FVC (%)",
            "Impaired FEV1","Impaired FVC","Impaired FEV1/FVC","Abnormal pattern")

# covariate formulas
cov3 <- "+gender+age+bmi+n_disease+edu+smoke1+drink1+fruit_vege+sport1+yr+season+rural_ubran+work+(1|pid)"

## location id of CO2 data
long_lat_CO2 <- read.csv("data/long_lat (TW and HK).csv")


# ---------------------------- nonlinear curve ---------------------------------
# nonlinear check (only for 97.5 two days CumE)
{
  xname_new <- "CO2_over_p97.5_2D"
  para_comb <- expand.grid(x=xname_new,y=yname,type=type[1])
  para_comb <- para_comb |>  mutate(across(.cols = x:type,.fns = as.character))
  
  out <- data.frame(outcome=NULL,exposure=NULL,type=NULL,ns=NULL,AIC=NULL,
                    Chisq=NULL,P=NULL)
  CumE_data <- read.csv("data/CumE.csv")
  tempdata <- left_join(com_data_new,long_lat_CO2)
  tempdata <- left_join(tempdata, CumE_data)
  tempdata <- as.data.frame(tempdata)
  
  s1 <- Sys.time()
  for (i in 1:nrow(para_comb)) {
    for(j in 2:5){
      tempformula_b <- paste0(para_comb$y[i],"~",para_comb$x[i],cov3) |>  as.formula()
      if(para_comb$y[i] %in% c("pf_fvc","pf_fev1","pf_mmf","ffratio")){
        base_model <- lmer(tempformula_b,data = tempdata)
      }else{
        base_model <- glmer(tempformula_b,data = tempdata,family = binomial,nAGQ = 0)
      }
      
      if(median(tempdata[,para_comb$x[i]])!=0){
        tempformula <- as.formula(paste0(para_comb$y[i],"~","ns(",
                                         para_comb$x[i],",df=",j,")",cov3))
      }else{
        tempformula <- as.formula(paste0(para_comb$y[i],"~","ns(",
                                         para_comb$x[i],
                                         ",df=",j,",Boundary.knots=c(-1,max(tempdata[,para_comb$x[i]])))",
                                         cov3))
      }
      
      if(para_comb$y[i] %in% c("pf_fvc","pf_fev1","pf_mmf","ffratio")){
        tempmodel <- lmer(tempformula,data = tempdata)
      }else{
        tempmodel <- glmer(tempformula,data = tempdata,family = binomial,nAGQ = 0)
      }
      
      a <- anova(base_model, tempmodel)
      
      tempresult <- data.frame(
        outcome = para_comb$y[i], exposure = para_comb$x[i],type = para_comb$type[i],
        ns = j,AIC = AIC(tempmodel),Chisq = a$Chisq[2],P = a$`Pr(>Chisq)`[2])
      
      out <- rbind(out,tempresult)
      
      print(c(i,j))
    }
  }
  s2 <- Sys.time()
  s2 - s1
  # Time difference of 1.638018 hours
  
  rio::export(out,"result/nonlinear check.xlsx")
  
  
  check_nonlinear <- rio::import("result/nonlinear check.xlsx")
  check_nonlinear |>  group_by(outcome) |>  arrange(AIC) |>  view()
  
}
# Figure 2: Concentration–response associations plot
{
  ## save the estimate results
  {
    xname_new <- "CO2_over_p97.5_2D"
    para_comb <- expand.grid(x=xname_new,y=yname,type="CumE")
    para_comb <- para_comb |>  mutate(across(.cols = x:type,.fns = as.character))
    
    CumE_data <- read.csv("data/CumE.csv")
    tempdata <- left_join(com_data_new,long_lat_CO2)
    tempdata <- left_join(tempdata, CumE_data)
    tempdata <- as.data.frame(tempdata)
    ob_t <- onebasis(tempdata[,"CO2_over_p97.5_2D"], fun = "ns", df = 2)
    
    nonl_est <- data.frame(outcome=NULL,CumE=NULL,est=NULL,low=NULL,high=NULL)
    s1 <- Sys.time()
    for (i in 1:nrow(para_comb)) {
      
      formula <- as.formula(paste(para_comb$y[i]," ~ ob_t", cov3, sep = ""))
      # Fit the model
      if(para_comb$y[i] %in% c("pf_fvc","pf_fev1","pf_mmf","ffratio")){
        model <- lmer(formula, data = tempdata)
        pr_t <- crosspred(ob_t, model, cen = min(tempdata[,"CO2_over_p97.5_2D"]))
        tempest <- data.frame(outcome = ylabel[i], 
                              CumE = pr_t$predvar, pr_t$matfit,
                              pr_t$matlow, pr_t$mathigh)
        colnames(tempest)[3:5] <- c("est","low","high")
      }else{
        model <- glmer(formula, data = tempdata,family = binomial(),nAGQ = 0)
        pr_t <- crosspred(ob_t, model, cen = min(tempdata[,"CO2_over_p97.5_2D"]))
        tempest <- data.frame(outcome = ylabel[i], 
                              CumE = pr_t$predvar, pr_t$matRRfit,
                              pr_t$matRRlow, pr_t$matRRhigh)
        colnames(tempest)[3:5] <- c("est","low","high")
      }
      
      nonl_est <- rbind(nonl_est, tempest)
    }
    s2 <- Sys.time()
    s2 - s1 ## Time difference of 6.737471 mins
    
    rio::export(nonl_est,"result/nonl_est.xlsx")
  }
  ## draw plot
  {
    CumE_data <- read.csv("data/CumE.csv")
    tempdata <- left_join(com_data_new,long_lat_CO2)
    tempdata <- left_join(tempdata, CumE_data)
    tempdata <- as.data.frame(tempdata)
    
    hist_plot <- ggplot(tempdata, aes(x = CO2_over_p97.5_2D)) +
      geom_histogram(fill = "grey", bins = 26,alpha = 0.3) +
      theme_void()
    
    mycolor <- c("#db346e")
    
    plotdata <- rio::import("result/nonl_est.xlsx")
    plotdata$ref <- ifelse(plotdata$outcome %in% c("FVC (ml)","FEV1 (ml)","MMEF (ml/s)","FEV1/FVC (%)"),
                           yes = 0, no = 1)
    plotdata$outcome <- factor(plotdata$outcome,levels = ylabel)
    
    nonl_plotlist <- list()
    for (i in 1:length(ylabel)) {
      y_label <- ifelse(ylabel[i] %in% c("FVC (ml)","FEV1 (ml)","MMEF (ml/s)","FEV1/FVC (%)"),
                        yes = "Coefficient (95%CI)", no = "Odd ratio (95%CI)")
      tmepplot <- plotdata |>  
        filter(outcome == ylabel[i]) |>  
        ggplot(aes(x = CumE, y = est))+
        theme_classic()+
        geom_line(linewidth=0.8,colour = mycolor)+
        geom_ribbon(aes(ymin = low, ymax = high),
                    colour = mycolor, fill = mycolor, alpha = 0.2)+
        geom_hline(aes(yintercept = ref),linetype = 2)+
        geom_vline(xintercept = quantile(tempdata$CO2_over_p97.5_2D,0.975),linetype = 1)+
        geom_vline(xintercept = quantile(tempdata$CO2_over_p97.5_2D,0.025),linetype = 1)+
        facet_wrap(~outcome)+
        labs(x = "CumE (ppm)", y = y_label) +
        theme(text = element_text(family = "serif"),
              axis.title = element_text(size = 12,face = "bold"),
              axis.text = element_text(size = 12),
              legend.position = "none",
              strip.background = element_blank(),
              strip.text = element_text(size = 12,face = "bold"))
      tempcomplot <- tmepplot + patchwork::inset_element(
        hist_plot, left = 0, bottom = 0, right = 1, top = 1)
      nonl_plotlist[[i]] <- tempcomplot
    }
    
    tiff("result/Figure 2.tiff", 
         width = 4000, height = 1700, res = 350, compression = "lzw")
    cowplot::plot_grid(plotlist = nonl_plotlist,nrow = 2)
    dev.off()
  }
}


# ---------------------------- quantile effect ---------------------------------
### effect estimate
{
  CumE_data <- read.csv("data/CumE.csv")
  tempdata <- left_join(com_data_new,long_lat_CO2)
  tempdata <- left_join(tempdata, CumE_data)
  
  q4 <- quantile(tempdata$CO2_over_p97.5_2D,c(0.25,0.5,0.75))
  tempdata <- tempdata |>  
    mutate(
      CO2_1 = ifelse(CO2_over_p97.5_2D > q4[1] & CO2_over_p97.5_2D <= q4[2], 1, 0),
      CO2_2 = ifelse(CO2_over_p97.5_2D > q4[2] & CO2_over_p97.5_2D <= q4[3], 1, 0),
      CO2_3 = ifelse(CO2_over_p97.5_2D > q4[3], 1, 0)
    )
  tempdata <- as.data.frame(tempdata)
  
  para_comb <- expand.grid(x="CO2_over_p97.5_2D",y=yname)
  para_comb <- para_comb |>  mutate(across(.cols = x:y,.fns = as.character))
  
  out_4q <- data.frame(outcome=NULL,exposure=NULL,AIC=NULL,Q=NULL,varname=NULL,
                       Estimate=NULL,low=NULL,high=NULL,P=NULL,effect=NULL)
  s1 <- Sys.time()
  for (i in 1:nrow(para_comb)) {
    tempformula <- paste0(para_comb$y[i],"~ CO2_1+CO2_2+CO2_3",cov3) |>  as.formula()
    if (para_comb$y[i] %in% c("pf_fvc", "pf_fev1", "pf_mmf", "ffratio")) {
      tempmodel <- lmer(tempformula, data = tempdata)
      mtype <- 1
    } else {
      tempmodel <- glmer(tempformula,data = tempdata,family = binomial,nAGQ = 0)
      mtype <- 2
    }
    tempresult <- cbind(data.frame(
      outcome = para_comb$y[i], exposure = para_comb$x[i],
      AIC = AIC(tempmodel),Q = c("CO2_25","CO2_50","CO2_75"),
      varname=rownames(summary(tempmodel)$coefficients)[c(2,3,4)]), 
      coeffun(tempmodel, row = c(2,3,4), mtype = mtype))
    out_4q <- rbind(out_4q,tempresult)
    
    print(i)
  }
  s2 <- Sys.time()
  s2 - s1
  
  rio::export(out_4q,"result/quantile_effect.xlsx")
  
}
### Supplemental Figure 4:
{
  plotdata <- rio::import("result/quantile_effect.xlsx")
  plotdata <- plotdata |> 
    mutate(Quartile = case_when(
      Q == "CO2_25" ~ "Second",
      Q == "CO2_50" ~ "Third",
      Q == "CO2_75" ~ "Fourth"
    ))
  plotdata$Quartile <- factor(plotdata$Quartile,levels = c("Second","Third","Fourth"))
  
  plotdata$yname <- plotdata$outcome
  plotdata$yname <- factor(plotdata$yname,levels = yname, labels = ylabel)
  
  plotdata <- plotdata |> 
    mutate(g = case_when(
      outcome %in% c("pf_fvc","pf_fev1","pf_mmf") ~ "Lung function parameters",
      outcome == "ffratio" ~ "FEV1/FVC",
      .default = "Adverse respiratory outcomes"
    ))
  
  
  plotdata |> 
    ggplot(aes(x = yname, y = Estimate,colour = Quartile))+
    theme_test()+
    scale_color_manual(values = c("#ff6f91","#d65db1","#845ec2"))+
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
  
  ggsave("result/Supplemental Figure 4.tiff",
         width = 10,height = 4,dpi = 350,compression = "lzw")
}