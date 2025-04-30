## -------------------------------------------------------------
##
## Script name: BirdNET Sierra Validation and Thresholds
##
## Script purpose:
##
## Author: Spencer R Keyser & Connor M. Wood
##
## Date Created: 2025-04-23
##
## Email: srk252@cornell.edu
##
## Github: https://github.com/skeyser
##
## -------------------------------------------------------------
##
## Notes:
##
##
## -------------------------------------------------------------

## Defaults
options(scipen = 10, digits = 10)

## -------------------------------------------------------------

## Package Loading
library(dplyr)
library(ggplot2)
library(here)
library(PresenceAbsence)
library(precrec)

## -------------------------------------------------------------

## -------------------------------------------------------------
##
## Begin Section: 2021 Validation and Thresholding
##
## -------------------------------------------------------------

## load species data
species.list_2021=list.files(path=here('./Data/Species_Thresholds/2021 Validation Results/'),
                             recursive=T,pattern='csv$',full.names = T)

## Pre-configured params
output_2021=NULL
score.pred_2021.c=seq(100,1000,1)
score.pred_2021.r=seq(-3,7,.1)

## Expected probability
desired_probability <- 0.9

## Set plot logicals
histogram="no"
logistic.plot="yes"
confidence.plot="yes"
prauc.plot="no"
par(mfrow=c(2,4))

## Loops for species
for(s in 1:length(species.list_2021)){
  
  ## Read species data in
  dt=read.csv(species.list_2021[s])
  
  ## Generate logit-scale scores
  dt$raw=log((dt$score/1000)/(1-(dt$score/1000)))
  
  ## Pull species name from the files
  species=gsub(".*/(.*)\\.csv", "\\1", species.list_2021[s])
  
  ## Pull matching species
  output_2021$species[s]=species
  
  ## Run glm for the validation data
  tmp.mod.c=glm(correct~score,dt,family = 'binomial')
  tmp.mod.r=glm(correct~raw,dt,family = 'binomial')
  
  ## Predict to the desired scoring threshold
  predictions.c=predict(tmp.mod.c,list(score=score.pred_2021.c),type='r')
  predictions.r=predict(tmp.mod.r,list(raw=score.pred_2021.r),type='r')
  
  ## Extract coefficients from the model
  output_2021$intercept.c[s]=tmp.mod.c$coefficients[1]
  output_2021$beta.c[s]=tmp.mod.c$coefficients[2]
  output_2021$intercept.r[s]=tmp.mod.r$coefficients[1]
  output_2021$beta.r[s]=tmp.mod.r$coefficients[2]
  
  ## Generate the cutoffs across different corresponding Pr(TP)
  cutoff85.c=(log(.85/(1-.85))-tmp.mod.c$coefficients[1])/tmp.mod.c$coefficients[2]
  cutoff90.c=(log(.90/(1-.90))-tmp.mod.c$coefficients[1])/tmp.mod.c$coefficients[2]
  cutoff95.c=(log(.95/(1-.95))-tmp.mod.c$coefficients[1])/tmp.mod.c$coefficients[2]
  cutoff975.c=(log(.975/(1-.975))-tmp.mod.c$coefficients[1])/tmp.mod.c$coefficients[2]
  cutoff99.c=(log(.99/(1-.99))-tmp.mod.c$coefficients[1])/tmp.mod.c$coefficients[2]
  
  cutoff85.r=(log(.85/(1-.85))-tmp.mod.r$coefficients[1])/tmp.mod.r$coefficients[2]
  cutoff90.r=(log(.90/(1-.90))-tmp.mod.r$coefficients[1])/tmp.mod.r$coefficients[2]
  cutoff95.r=(log(.95/(1-.95))-tmp.mod.r$coefficients[1])/tmp.mod.r$coefficients[2]
  cutoff975.r=(log(.975/(1-.975))-tmp.mod.r$coefficients[1])/tmp.mod.r$coefficients[2]
  cutoff99.r=(log(.99/(1-.99))-tmp.mod.r$coefficients[1])/tmp.mod.r$coefficients[2]
  
  #output_2021$th.950[s]=round(predict(tmp.mod,newdata=data.frame(score=950),type='r'),3)
  #output_2021$th.1000[s]=round(predict(tmp.mod,newdata=data.frame(score=1000),type='r'),3)
  
  ## Visualize
  if(logistic.plot=="yes"){
    plot(correct~raw,dt,main=paste(species, "Raw"),
         xlim=range(score.pred_2021.r),pch=16,cex=1.5,col=rgb(0,0,0,.2))
    lines(predictions.r~score.pred_2021.r,lwd=4,col=rgb(0,.75,1,.5))
    abline(v=cutoff85.r,col='yellow',lwd=4)
    abline(v=cutoff90.r,col='orange',lwd=4)
    abline(v=cutoff95.r,col='red',lwd=4)
    abline(v=cutoff975.r, col = 'green', lwd=4)
    abline(v=cutoff99.r,col='magenta',lwd=4)
  }
  if(confidence.plot=="yes"){
    plot(correct~score,dt,main=paste(species, "Conf"),
         xlim=range(score.pred_2021.c),pch=16,cex=1.5,col=rgb(0,0,0,.2))
    lines(predictions.c~score.pred_2021.c,lwd=4,col=rgb(0,.75,1,.5))
    abline(v=cutoff85.c,col='yellow',lwd=4)
    abline(v=cutoff90.c,col='orange',lwd=4)
    abline(v=cutoff95.c,col='red',lwd=4)
    abline(v=cutoff975.c, col='green',lwd=4)
    abline(v=cutoff99.c,col='magenta',lwd=4)
  }
    
  
  if(histogram=="yes"){
    hist(dt$score,main=species,
         xlim=c(0,100),breaks=seq(0,100,10),lty='blank',xlab='Score')
    abline(v=cutoff95,col=rgb(1,.5,0,.7),lwd=5)
    abline(v=cutoff99,col=rgb(1,0,0,.7),lwd=5)
  }
  
  {
    output_2021$cutoff_85.r=(log(.85/(1-.85))-output_2021$intercept.r)/output_2021$beta.r
    output_2021$cutoff_90.r=(log(.90/(1-.90))-output_2021$intercept.r)/output_2021$beta.r
    output_2021$cutoff_95.r=(log(.95/(1-.95))-output_2021$intercept.r)/output_2021$beta.r
    output_2021$cutoff_975.r=(log(.975/(1-.975))-output_2021$intercept.r)/output_2021$beta.r
    output_2021$cutoff_99.r=(log(.99/(1-.99))-output_2021$intercept.r)/output_2021$beta.r
    
    output_2021$cutoff_90.r_conf=1/(1+exp(-output_2021$cutoff_90.r))
    output_2021$cutoff_95.r_conf=1/(1+exp(-output_2021$cutoff_95.r))
    output_2021$cutoff_975.r_conf=1/(1+exp(-output_2021$cutoff_975.r))
    output_2021$cutoff_99.r_conf=1/(1+exp(-output_2021$cutoff_99.r))
  }
  
  
  ## Check precision and recall curves
  ## We need to threshold at the values of interest
  dt$conf_pred_90 <- ifelse(dt$score >= cutoff90.c, 1, 0)
  dt$conf_pred_95 <- ifelse(dt$score >= cutoff95.c, 1, 0)
  dt$conf_pred_975 <- ifelse(dt$score >= cutoff975.c, 1, 0)
  dt$conf_pred_99 <- ifelse(dt$score >= cutoff99.c, 1, 0)
  dt$r_pred_90 <- ifelse(dt$raw >= cutoff90.r, 1, 0)
  dt$r_pred_95 <- ifelse(dt$raw >= cutoff95.r, 1, 0)
  dt$r_pred_975 <- ifelse(dt$raw >= cutoff975.r, 1, 0)
  dt$r_pred_99 <- ifelse(dt$raw >= cutoff99.r, 1, 0)
  
  ## Join the scores for precrec object
  scores <- join_scores(dt$conf_pred_90, dt$conf_pred_95, dt$conf_pred_975, dt$conf_pred_99)
  labels <- join_scores(dt$correct, dt$correct, dt$correct, dt$correct)
  dt.c <- mmdata(scores, labels, modnames = c("90", "95", "975", "99"))
  
  scores.r <- join_scores(dt$r_pred_90, dt$r_pred_95, dt$r_pred_975, dt$r_pred_99)
  labels.r <- join_scores(dt$correct, dt$correct, dt$correct, dt$correct)
  dt.r <- mmdata(scores.r, labels.r, modnames = c("90", "95", "975", "99"))
  
  ## Calculate PRAUC curve (focuses on positive class)
  em.conf <- precrec::evalmod(dt.c)
  pr_auc.conf <- precrec::auc(em.conf) |> 
    filter(curvetypes == "PRC") |> 
    pull(aucs)
  
  em.r <- precrec::evalmod(dt.r)
  pr_auc.r <- precrec::auc(em.r) |> 
    filter(curvetypes == "PRC") |> 
    pull(aucs)
  
  if(prauc.plot == "yes"){
    autoplot(em.conf) + ggtitle(species)
  }
  
  ## PRAUC values
  output_2021$PRAUC_90.c[s] <- pr_auc.conf[1]
  output_2021$PRAUC_95.c[s] <- pr_auc.conf[2]
  output_2021$PRAUC_975.c[s] <- pr_auc.conf[3]
  output_2021$PRAUC_99.c[s] <- pr_auc.conf[4]
  
  output_2021$PRAUC_90.r[s] <- pr_auc.r[1]
  output_2021$PRAUC_95.r[s] <- pr_auc.r[2]
  output_2021$PRAUC_975.r[s] <- pr_auc.r[3]
  output_2021$PRAUC_99.r[s] <- pr_auc.r[4]
  
  ## Calculate a variety of prediction metrics
  pa_metrics.c <- dt |> 
    mutate(ID = 1:nrow(dt), score = score/1000) |> 
    select(ID, correct, score) |> 
    PresenceAbsence::presence.absence.accuracy(na.rm = TRUE, 
                                               st.dev = FALSE, 
                                               threshold = c(0.9, 0.95, 0.975, 0.99),
                                               )
  
  pa_metrics.r <- dt |> 
    mutate(ID = 1:nrow(dt), score = raw) |> 
    select(ID, correct, score) |> 
    PresenceAbsence::presence.absence.accuracy(na.rm = TRUE, 
                                               st.dev = FALSE, 
                                               threshold = c(0.9, 0.95, 0.975, 0.99),
    )
  
  ## Generate various traditional test metrics for each threshold from BirdNet Conf Score
  output_2021$PCC_90.c[s] <- pa_metrics.c$PCC[pa_metrics.c$threshold == 0.90]
  output_2021$PCC_95.c[s] <- pa_metrics.c$PCC[pa_metrics.c$threshold == 0.95]
  output_2021$PCC_975.c[s] <- pa_metrics.c$PCC[pa_metrics.c$threshold == 0.975]
  output_2021$PCC_99.c[s] <- pa_metrics.c$PCC[pa_metrics.c$threshold == 0.99]
  output_2021$PCC_90.r[s] <- pa_metrics.r$PCC[pa_metrics.r$threshold == 0.90]
  output_2021$PCC_95.r[s] <- pa_metrics.r$PCC[pa_metrics.r$threshold == 0.95]
  output_2021$PCC_975.r[s] <- pa_metrics.r$PCC[pa_metrics.r$threshold == 0.975]
  output_2021$PCC_99.r[s] <- pa_metrics.r$PCC[pa_metrics.r$threshold == 0.99]
  
  output_2021$Sensitivity_90.c[s] <- pa_metrics.c$sensitivity[pa_metrics.c$threshold == 0.90]
  output_2021$Sensitivity_95.c[s] <- pa_metrics.c$sensitivity[pa_metrics.c$threshold == 0.95]
  output_2021$Sensitivity_99.c[s] <- pa_metrics.c$sensitivity[pa_metrics.c$threshold == 0.99]
  output_2021$Sensitivity_90.r[s] <- pa_metrics.r$sensitivity[pa_metrics.r$threshold == 0.90]
  output_2021$Sensitivity_95.r[s] <- pa_metrics.r$sensitivity[pa_metrics.r$threshold == 0.95]
  output_2021$Sensitivity_99.r[s] <- pa_metrics.r$sensitivity[pa_metrics.r$threshold == 0.99]
  
  output_2021$Specificity_90.c[s] <- pa_metrics.c$specificity[pa_metrics.c$threshold == 0.90]
  output_2021$Specificity_95.c[s] <- pa_metrics.c$specificity[pa_metrics.c$threshold == 0.95]
  output_2021$Specificity_99.c[s] <- pa_metrics.c$specificity[pa_metrics.c$threshold == 0.99]
  output_2021$Specificity_90.r[s] <- pa_metrics.r$specificity[pa_metrics.r$threshold == 0.90]
  output_2021$Specificity_95.r[s] <- pa_metrics.r$specificity[pa_metrics.r$threshold == 0.95]
  output_2021$Specificity_99.r[s] <- pa_metrics.r$specificity[pa_metrics.r$threshold == 0.99]
  
  output_2021$Kappa_90.c[s] <- pa_metrics.c$Kappa[pa_metrics.c$threshold == 0.90]
  output_2021$Kappa_95.c[s] <- pa_metrics.c$Kappa[pa_metrics.c$threshold == 0.95]
  output_2021$Kappa_99.c[s] <- pa_metrics.c$Kappa[pa_metrics.c$threshold == 0.99]
  output_2021$Kappa_90.r[s] <- pa_metrics.r$Kappa[pa_metrics.r$threshold == 0.90]
  output_2021$Kappa_95.r[s] <- pa_metrics.r$Kappa[pa_metrics.r$threshold == 0.95]
  output_2021$Kappa_99.r[s] <- pa_metrics.r$Kappa[pa_metrics.r$threshold == 0.99]
  
  output_2021$AUC.c[s] <- pa_metrics.c$AUC[pa_metrics.c$threshold == 0.90]
  output_2021$AUC.r[s] <- pa_metrics.r$AUC[pa_metrics.r$threshold == 0.90]
  
  # Create complete 2x2 confusion matrices with all possible combinations
  ConfMat90.c <- table(factor(dt$correct, levels=c(0,1)), 
                     factor(dt$conf_pred_90, levels=c(0,1)))
  ConfMat95.c <- table(factor(dt$correct, levels=c(0,1)), 
                     factor(dt$conf_pred_95, levels=c(0,1)))
  ConfMat975.c <- table(factor(dt$correct, levels=c(0,1)), 
                       factor(dt$conf_pred_975, levels=c(0,1)))
  ConfMat99.c <- table(factor(dt$correct, levels=c(0,1)), 
                     factor(dt$conf_pred_99, levels=c(0,1)))
  
  ConfMat90.r <- table(factor(dt$correct, levels=c(0,1)), 
                       factor(dt$r_pred_90, levels=c(0,1)))
  ConfMat95.r <- table(factor(dt$correct, levels=c(0,1)), 
                       factor(dt$r_pred_95, levels=c(0,1)))
  ConfMat975.r <- table(factor(dt$correct, levels=c(0,1)), 
                       factor(dt$r_pred_975, levels=c(0,1)))
  ConfMat99.r <- table(factor(dt$correct, levels=c(0,1)), 
                       factor(dt$r_pred_99, levels=c(0,1)))
  
  ## Calculate precision with safe division (if denominator is 0, return 0)
  output_2021$Precision_90.c[s] <- ifelse(
    (ConfMat90.c[2,2] + ConfMat90.c[1,2]) == 0,
    NA,
    ConfMat90.c[2,2] / (ConfMat90.c[2,2] + ConfMat90.c[1,2])
  )
  output_2021$Precision_95.c[s] <- ifelse(
    (ConfMat95.c[2,2] + ConfMat95.c[1,2]) == 0,
    NA,
    ConfMat95.c[2,2] / (ConfMat95.c[2,2] + ConfMat95.c[1,2])
  )
  output_2021$Precision_975.c[s] <- ifelse(
    (ConfMat975.c[2,2] + ConfMat975.c[1,2]) == 0,
    NA,
    ConfMat975.c[2,2] / (ConfMat975.c[2,2] + ConfMat975.c[1,2])
  )
  output_2021$Precision_99.c[s] <- ifelse(
    (ConfMat99.c[2,2] + ConfMat99.c[1,2]) == 0,
    NA,
    ConfMat99.c[2,2] / (ConfMat99.c[2,2] + ConfMat99.c[1,2])
  )
  
  output_2021$Precision_90.r[s] <- ifelse(
    (ConfMat90.r[2,2] + ConfMat90.r[1,2]) == 0,
    NA,
    ConfMat90.r[2,2] / (ConfMat90.r[2,2] + ConfMat90.r[1,2])
  )
  output_2021$Precision_95.r[s] <- ifelse(
    (ConfMat95.r[2,2] + ConfMat95.r[1,2]) == 0,
    NA,
    ConfMat95.r[2,2] / (ConfMat95.r[2,2] + ConfMat95.r[1,2])
  )
  output_2021$Precision_975.r[s] <- ifelse(
    (ConfMat975.r[2,2] + ConfMat975.r[1,2]) == 0,
    NA,
    ConfMat975.r[2,2] / (ConfMat975.r[2,2] + ConfMat975.r[1,2])
  )
  output_2021$Precision_99.r[s] <- ifelse(
    (ConfMat99.r[2,2] + ConfMat99.r[1,2]) == 0,
    NA,
    ConfMat99.r[2,2] / (ConfMat99.r[2,2] + ConfMat99.r[1,2])
  )
  
  ## Calculate recall with safe division
  output_2021$Recall_90.c[s] <- ifelse(
    (ConfMat90.c[2,2] + ConfMat90.c[2,1]) == 0,
    NA,
    ConfMat90.c[2,2] / (ConfMat90.c[2,2] + ConfMat90.c[2,1])
  )
  output_2021$Recall_95.c[s] <- ifelse(
    (ConfMat95.c[2,2] + ConfMat95.c[2,1]) == 0,
    NA,
    ConfMat95.c[2,2] / (ConfMat95.c[2,2] + ConfMat95.c[2,1])
  )
  output_2021$Recall_975.c[s] <- ifelse(
    (ConfMat975.c[2,2] + ConfMat975.c[2,1]) == 0,
    NA,
    ConfMat975.c[2,2] / (ConfMat975.c[2,2] + ConfMat975.c[2,1])
  )
  output_2021$Recall_99.c[s] <- ifelse(
    (ConfMat99.c[2,2] + ConfMat99.c[2,1]) == 0,
    NA,
    ConfMat99.c[2,2] / (ConfMat99.c[2,2] + ConfMat99.c[2,1])
  )
  
  output_2021$Recall_90.r[s] <- ifelse(
    (ConfMat90.r[2,2] + ConfMat90.r[2,1]) == 0,
    NA,
    ConfMat90.r[2,2] / (ConfMat90.r[2,2] + ConfMat90.r[2,1])
  )
  output_2021$Recall_95.r[s] <- ifelse(
    (ConfMat95.r[2,2] + ConfMat95.r[2,1]) == 0,
    NA,
    ConfMat95.r[2,2] / (ConfMat95.r[2,2] + ConfMat95.r[2,1])
  )
  output_2021$Recall_975.r[s] <- ifelse(
    (ConfMat975.r[2,2] + ConfMat975.r[2,1]) == 0,
    NA,
    ConfMat975.r[2,2] / (ConfMat975.r[2,2] + ConfMat975.r[2,1])
  )
  output_2021$Recall_99.r[s] <- ifelse(
    (ConfMat99.r[2,2] + ConfMat99.r[2,1]) == 0,
    NA,
    ConfMat99.r[2,2] / (ConfMat99.r[2,2] + ConfMat99.r[2,1])
  )
  
  ## MCC Score and threshold
  mcc_f1.c <- mccf1::mccf1(
    response = dt$correct,
    predictor = dt$score/1000
  )
  
  mcc_f1.r <- mccf1::mccf1(
    response = dt$correct,
    predictor = dt$raw
  )
  
  ## Scores for MCCF1 and identified "best threshold"
  output_2021$MCCF1.c[s] <- summary(mcc_f1.c)$mccf1_metric
  output_2021$MCCF1_ScoreThresh.c[s] <- summary(mcc_f1.c)$best_threshold
  
  output_2021$MCCF1.r[s] <- summary(mcc_f1.r)$mccf1_metric
  output_2021$MCCF1_ScoreThresh.r[s] <- summary(mcc_f1.r)$best_threshold
  
  
  if(s==length(species.list_2021)){
    output_2021=data.frame(output_2021)
  }
  
}

glimpse(output_2021)

output_2021_long <- output_2021 %>%
  # First, select columns that contain threshold information
  select(species, matches("cutoff_\\d+\\.r_|Precision_|Recall_|PRAUC_|PCC_|Sensitivity_|Specificity_|Kappa_|AUC_")) %>%
  pivot_longer(
    cols = -species,
    names_to = "temp",
    values_to = "value"
  ) %>%
  # Split the temp column into metric and threshold
  separate(temp, 
           into = c("metric", "threshold"), 
           sep = "_",
           extra = "merge") %>%
  separate(threshold,
           into = c("threshold", "type"),
           sep = "\\.",
           extra = "merge") |> 
  # Clean up the threshold values (remove any remaining text)
  mutate(threshold = gsub("\\D", "", threshold))

## Pull the rows where precision is greatest but recall is still above 0.5
output_2021_MaxPrec.r <- output_2021_long  |> 
  # Create temporary dataframe with precision and recall paired
  filter(type == "r") |> 
  group_by(species, threshold) |> 
  mutate(
    recall_value = value[metric == "Recall"],
    precision_value = value[metric == "Precision"]
  ) |> 
  # Filter for recall >= 0.5
  filter(recall_value > 0) |> 
  filter(threshold %in% c("975", "99")) |> 
  # For each species, keep only the threshold with highest precision
  group_by(species) |> 
  filter(precision_value == max(precision_value)) |>
  filter(recall_value == max(recall_value)) |> 
  # Remove the temporary columns and maintain original format
  select(-recall_value, -precision_value) %>%
  # Sort if desired
  arrange(species, metric, threshold)

output_2021_MaxPrec <- output_2021_long  |> 
  # Create temporary dataframe with precision and recall paired
  filter(type == "c") |> 
  group_by(species, threshold) |> 
  mutate(
    recall_value = value[metric == "Recall"],
    precision_value = value[metric == "Precision"]
  ) |> 
  # Filter for recall >= 0.5
  filter(recall_value > 0)  |>
  filter(threshold %in% c("975", "99")) |> 
  # For each species, keep only the threshold with highest precision
  group_by(species) |> 
  filter(precision_value == max(precision_value)) |>  
  filter(recall_value == max(recall_value)) |> 
  # Remove the temporary columns and maintain original format
  #select(-recall_value, -precision_value) %>%
  # Sort if desired
  arrange(species, metric, threshold) |>
  bind_rows(output_2021_MaxPrec.r) |> 
  arrange(species) |> 
  select(-value, -metric) |> 
  distinct(species, type, threshold, recall_value, precision_value) |> 
  select(-recall_value, -precision_value) |> 
  group_by(species, type) |> 
  filter(threshold == min(as.numeric(threshold))) |> 
  mutate(type = ifelse(type == "c", "r_conf", "r")) |> 
  pivot_wider(names_from = type, values_from = threshold, values_fill = NA) |> 
  ungroup()

## Add in the best threshold for Raw and Conf
output_2021_best_thresh <- output_2021 |> 
  left_join(output_2021_MaxPrec) |> 
  select(species:cutoff_99.r_conf, BestThresh.r_conf = r_conf, BestThresh.r = r)

## Take the 99th scores and see what we are missing
out99 <- output_2021_long  |> 
  # Create temporary dataframe with precision and recall paired
  filter(threshold %in% c("95", "99") & metric %in% c("Precision", "Recall") & type == "r")


par(mfrow=c(1,2))
hist(output_2021$cutoff90.r_conf, xlab="Confidence Score", xlim=c(0,1),
     main="Confidence score yielding pr(correct)=0.90")
hist(output_2021$cutoff90.r[-74], xlab="Logit Score", #xlim=c(0,1),
     main="Logit score yielding pr(correct)=0.90")

write.csv(output_2021, here('./Data/Thresholds_2021_20230309_AllMetrics.csv'),row.names = F)
write.csv(output_2021_MaxPrec, here('./Data/Thresholds_2021_20230309_BestPrecRec.csv'),row.names = F)
write.csv(output_2021_best_thresh, here('./Data/Thresholds_2021_20230309_BestThreshold_975min.csv'),row.names = F)
output_2021[output_2021$species=="Olive-sided Flycatcher",]

output_2021_Sequential.r <- output_2021_long  |> 
  # Create temporary dataframe with precision and recall paired
  filter(type == "r") |> 
  group_by(species, threshold) |> 
  mutate(
    recall_value = value[metric == "Recall"],
    precision_value = value[metric == "Precision"]
  ) |> 
  # First try threshold 99
  group_by(species) |> 
  mutate(
    selected_threshold = case_when(
      any(threshold == "99" & recall_value >= 0) ~ "99",
      any(threshold == "975" & recall_value >= 0) ~ "975",
      any(threshold == "95" & recall_value >= 0) ~ "95",
      any(threshold == "90" & recall_value >= 0) ~ "90",
      TRUE ~ NA# default if none of the above meet criteria
    )
  ) |> 
  filter(threshold == selected_threshold) |> 
  # Rest of your original code
  arrange(species, metric, threshold)

output_2021_Sequential <- output_2021_long  |> 
  # Create temporary dataframe with precision and recall paired
  filter(type == "c") |> 
  group_by(species, threshold) |> 
  mutate(
    recall_value = value[metric == "Recall"],
    precision_value = value[metric == "Precision"]
  ) |> 
  # First try threshold 99
  group_by(species) |> 
  mutate(
    selected_threshold = case_when(
      any(threshold == "99" & recall_value >= 0.1) ~ "99",
      any(threshold == "975" & recall_value >= 0.1) ~ "975",
      any(threshold == "95" & recall_value >= 0.1) ~ "95",
      any(threshold == "90" & recall_value >= 0.1) ~ "90",
      TRUE ~ NA# default if none of the above meet criteria
    )
  ) |> 
  filter(threshold == selected_threshold) |> 
  # Rest of your original code
  arrange(species, metric, threshold) |>
  bind_rows(output_2021_Sequential.r) |> 
  arrange(species) |> 
  select(-value, -metric) |> 
  distinct(species, type, threshold, recall_value, precision_value) |> 
  select(-recall_value, -precision_value) |> 
  group_by(species, type) |> 
  filter(threshold == min(as.numeric(threshold))) |> 
  mutate(type = ifelse(type == "c", "r_conf", "r")) |> 
  pivot_wider(names_from = type, values_from = threshold, values_fill = NA) |> 
  ungroup()




#################################################################################################




save.image("BirdNET_thresholds.RData")


library(pROC)

# Create ROC object
roc_obj <- roc(dt$correct, dt$score)

# Find optimal threshold balancing sensitivity/specificity
coords <- coords(roc_obj, "best", ret="threshold")

# Plot ROC curve
plot(roc_obj)
points(coords[1], coords[2], pch=19, col="red")

# Get various thresholds based on different criteria
thresholds <- coords(roc_obj, x="all", ret=c("threshold", "specificity", "sensitivity"))

# Find threshold that gives desired specificity
desired_spec <- thresholds[which.min(abs(thresholds$specificity - 0.99)),]



