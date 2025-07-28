# =============================================================================
# Script: 04_plots_tables.R
# Project: Linear and Generalized Linear Models with Poisson Distribution
#          Applied to Count Data in Neuropsychology: A Comparative Analysis
# Author: Diego Rivera
# Date: July 27, 2025
# Purpose: Generate PIP tables/plots, CV comparison plots, HC percentile comparisons,
#          KS tests, ROC curves, and violin plots by CDR severity
# =============================================================================

# ---- 0. Load dependencies ----
library(tidyverse)   # ggplot2, dplyr, tidyr, readr, purrr
library(patchwork)   # for arranging plots
library(cowplot)     # for plot_grid
library(pROC)        # ROC analysis

# =============================================================================
# 1. Posterior Inclusion Probabilities (PIP)
# =============================================================================
# Load fitted PIP models
load(file.path("models", "PIP_models.Rdata"))

# Extract unique covariates
all_covariates <- unique(unlist(
  map(models, ~ names(summary(.x[[1]])[,1]))
))

# Function to extract PIPs
extract_pips <- function(models_list, type = c("lm","glm")) {
  type  <- match.arg(type)
  index <- if(type == "lm") 1 else 2
  mat   <- map_df(models_list, function(m) {
    pip <- summary(m[[index]])[,1]
    full <- set_names(rep(0, length(all_covariates)), all_covariates)
    full[names(pip)] <- pip
    full
  }, .id = "model")
  mat
}

# Create and save PIP tables
pip_lm  <- extract_pips(models, "lm")
pip_glm <- extract_pips(models, "glm")
write_csv(pip_lm,  file.path("results","PIP_linear_models.csv"))
write_csv(pip_glm, file.path("results","PIP_glm_models.csv"))

# Plot settings
covariates <- c(
  "age1","age2","log(edu)","sexFemale",
  "age1:log(edu)","age1:sexFemale",
  "age2:log(edu)","age2:sexFemale",
  "log(edu):sexFemale"
)
labels_x <- c(
  "Age", expression(Age^2), expression(ln(Edu)), "Sex",
  expression(Age %*% ln(Edu)), expression(Age %*% Sex),
  expression(Age^2 %*% ln(Edu)), expression(Age^2 %*% Sex),
  expression(ln(Edu) %*% Sex)
)
titles <- c(
  "Stroop Palabra","Stroop Color","Stroop PC","SDMT",
  "Letter F","Letter A","Letter S","Letter M",
  "Animals","Fruits","Professions"
)

create_pip_plot <- function(name,title) {
  df <- tibble(
    covariate = factor(covariates, levels=covariates),
    LM  = summary(models[[name]]$model.lm)[covariates,1],
    GLM = summary(models[[name]]$model.glm)[covariates,1]
  ) %>% pivot_longer(c(LM,GLM), names_to="model", values_to="PIP")
  ggplot(df,aes(covariate,PIP,fill=model))+
    geom_col(position=position_dodge(width=0.7),width=0.6)+
    geom_hline(yintercept=0.5,linetype="dashed")+
    scale_x_discrete(labels=labels_x)+
    scale_y_continuous(limits=c(0,1.05))+
    labs(title=title,y="PIP",x=NULL)+
    theme_minimal(base_size=14)+
    theme(axis.text.x=element_text(angle=45,hjust=1),legend.title=element_blank())
}
# Save PIP plots
pdir <- file.path("figures","plots_pip")
if(!dir.exists(pdir)) dir.create(pdir,recursive=TRUE)
plots <- map2(names(models),titles, create_pip_plot)
walk2(plots,names(models),~ ggsave(file.path(pdir,paste0(.y,"_pip.jpg")),.x,width=8,height=5,dpi=300))
# Grid
ggsave(file.path(pdir,"pip_grid.pdf"), wrap_plots(plots,ncol=3),width=14,height=18,device=cairo_pdf)

# =============================================================================
# 2. Cross-Validation Comparison
# =============================================================================
cv_lm  <- read_csv(file.path("results","cross_validation_LM_summary_parallel.csv")) %>% filter(type=="LM")
cv_glm <- read_csv(file.path("results","cross_validation_GLM_summary_parallel.csv")) %>% filter(type=="GLM")
cv <- left_join(
  cv_lm  %>% select(model,mean_RMSE_lm=mean_RMSE,mean_MAE_lm=mean_MAE),
  cv_glm %>% select(model,mean_RMSE_glm=mean_RMSE,mean_MAE_glm=mean_MAE),
  by="model"
) %>% mutate(diff_rmse=mean_RMSE_lm-mean_RMSE_glm,diff_mae=mean_MAE_lm-mean_MAE_glm)

# Plot diff
p1 <- cv %>% ggplot(aes(reorder(model,diff_rmse),diff_rmse,fill=diff_rmse>0))+geom_col()+coord_flip()+labs(y="RMSE diff")+theme_minimal()
p2 <- cv %>% ggplot(aes(reorder(model,diff_mae),diff_mae,fill=diff_mae>0))+geom_col()+coord_flip()+labs(y="MAE diff")+theme_minimal()
# Save
cvdir <- file.path("figures","cv_comparison")
if(!dir.exists(cvdir)) dir.create(cvdir,recursive=TRUE)
ggsave(file.path(cvdir,"cv_diff.pdf"),p1/p2,width=10,height=12)

# =============================================================================
# 3. Healthy Controls Percentile Comparison
# =============================================================================
combined_hc <- combined %>% filter(Grupo.y=="HC")
vars <- names(combined_hc) %>% str_subset("_lm$") %>% str_remove("_lm$")
titles <- titles
hcdir <- file.path("figures","hc_percentile"); if(!dir.exists(hcdir)) dir.create(hcdir,recursive=TRUE)
plot_hc <- function(var,title){
  df <- combined_hc %>% transmute(mean=(.data[[paste0(var,"_lm")]]+.data[[paste0(var,"_glm")]])/2,
                                  diff=(.data[[paste0(var,"_lm")]]-.data[[paste0(var,"_glm")]]),id=row_number())
  m<-mean(df$diff);s<-sd(df$diff)
  p_ba <- ggplot(df,aes(mean,diff))+geom_point(alpha=0.5)+geom_hline(yintercept=m)+geom_hline(yintercept=m+1.96*s,linetype="dashed")+geom_hline(yintercept=m-1.96*s,linetype="dashed")+labs(title=paste(title,"BA"))+theme_minimal()
  df2<-df %>% pivot_longer(c(.data[[paste0(var,"_lm")]],.data[[paste0(var,"_glm")]]),names_to="Model",values_to="P")
  p_ln<-ggplot(tibble(Participant=row_number(),lm=.data[[paste0(var,"_lm")]],glm=.data[[paste0(var,"_glm")]]),aes(Participant,P,color=Model))+geom_line(aes(group=Participant),alpha=0.5)+geom_point()+labs(title=paste(title,"Pct"))+theme_minimal()+theme(axis.text.x=element_blank())
  p_ba|p_ln
}
# Save 1-4,5-8,9-11
i1<-wrap_plots(map2(vars[1:4],titles[1:4],plot_hc),ncol=1);ggsave(file.path(hcdir,"hc_1to4.pdf"),i1,width=10,height=16)
i2<-wrap_plots(map2(vars[5:8],titles[5:8],plot_hc),ncol=1);ggsave(file.path(hcdir,"hc_5to8.pdf"),i2,width=10,height=16)
i3<-wrap_plots(map2(vars[9:11],titles[9:11],plot_hc),ncol=1);ggsave(file.path(hcdir,"hc_9to11.pdf"),i3,width=10,height=12)

# =============================================================================
# 4. KS Tests
# =============================================================================
ks <- map_df(vars, ~{
  x <- combined_hc[[paste0(.x,"_lm")]]
  y <- combined_hc[[paste0(.x,"_glm")]]
  t<-ks.test(x,y)
  tibble(var=.x,stat=t$statistic,pvalue=t$p.value)
}) %>% arrange(pvalue)
write_csv(ks,file.path(hcdir,"ks_results.csv"))

# =============================================================================
# 5. ROC Curves
# =============================================================================
rocdir <- file.path("figures","roc_curves"); if(!dir.exists(rocdir)) dir.create(rocdir,recursive=TRUE)
roc_plots<-map2(vars,titles,~{
  r1<-roc(combined$Grupo.y,combined[[paste0(.x,"_lm")]],quiet=TRUE)
  r2<-roc(combined$Grupo.y,combined[[paste0(.x,"_glm")]],quiet=TRUE)
  df1<-ggroc(r1)$data%>%mutate(Model="LM")
  df2<-ggroc(r2)$data%>%mutate(Model="GLM")
  ggplot(bind_rows(df1,df2),aes(1-specificity,sensitivity,color=Model))+geom_line()+geom_abline(linetype="dashed")+labs(title=.y)+theme_minimal()
})
pdf(file.path(rocdir,"ROC_all.pdf"),width=12,height=16);print(plot_grid(plotlist=roc_plots,ncol=3));dev.off()

# =============================================================================
# 6. Violin Plots by CDR (Alz sample)
# =============================================================================
alzdir<-file.path("figures","violin_plots"); if(!dir.exists(alzdir)) dir.create(alzdir,recursive=TRUE)
combined_alz<-combined %>% filter(Grupo.y=="Alz") %>% mutate(CDR=factor(data_dx_al$CDR,levels=c("Very Mild / Mild","Moderate","Severe")))
violins<-map2(vars,titles,~{
  df1<-combined_alz %>% transmute(CDR,percentil=.data[[paste0(.x,"_lm")]],Model="LM")
  df2<-combined_alz %>% transmute(CDR,percentil=.data[[paste0(.x,"_glm")]],Model="GLM")
  ggplot(bind_rows(df1,df2),aes(CDR,percentil,fill=Model))+geom_violin(position=position_dodge(0.8),alpha=0.7)+geom_boxplot(width=0.1,position=position_dodge(0.8),outlier.shape=NA)+labs(title=.y)+theme_minimal()
})
walk(list(1:4,5:8,9:11),~{
  pdf(file.path(alzdir,paste0("violins_",min(.x),"_to_",max(.x),".pdf")),width=12,height=10)
  print(plot_grid(plotlist=violins[.x],ncol=2));dev.off()
})

message("All figures and tables generated in 'results/' and 'figures/' directories.")
