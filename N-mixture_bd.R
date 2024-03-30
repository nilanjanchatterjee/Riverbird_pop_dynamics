library(unmarked)
library(dplyr)
library(ggplot2)
library(MuMIn)
library(AICcmodavg)


######################################################################################################
####################### Plumbeous Water Redstart
#setwd("E:/work folder/manuscripts/river bird co-occurrence and encounter fate/Bird Study")
bd<-read.csv("./Data/BD_multi.csv",header=T)
head(bd)
dim(bd)
site_cov <- data.frame(alt=bd$Alt,flow=(bd$Flow),veg=as.factor(bd$Bank_veg),wor=bd$WoR )

###naive estimates for year 2014 modelling with habitat factors
bd_nm_14<-unmarkedFramePCount(y=bd[,2:4], siteCovs = site_cov)
head(bd_nm_14)
summary(bd_nm_14)
bd_14_x1<- pcount(~1~ flow,bd_nm_14, K=20)
bd_14_x1
bd_14_x2<- pcount(~1~ alt+ wor +flow ,bd_nm_14, K=20)
bd_14_x2
bd_14_x3<- pcount(~1~ alt+ flow,bd_nm_14, K=20) ## Best fit model
bd_14_x3
bd_14_x4<- pcount(~1~ alt,bd_nm_14, K=20)
bd_14_x4

#### dreadge to get the global estimate and compare the models
Nmix_dredge <- dredge(global.model = bd_14_x2,rank = "AICc")# the rank argument can use AICc, QAICc, or others (see help page)
                      
Nmix_dredge ##the top model has a weightage of 0.923 and the differance of delta >5

##### predict estimate 
nmix_modavg_lam_predict_bd14 <- modavgPred(list(nmix_full = bd_14_x1, nmix_sub = bd_14_x3),
                                            # c.hat = 1, # to change variance inflation factor, default = 1) 
                                            parm.type = "lambda", # psi = occupancy, can also be "detect" for detection probability
                                            newdata = bd_nm_14@siteCovs)[c("mod.avg.pred",  "lower.CL",  "upper.CL")]


## Put predictions, CI, and all site covariates into one data frame
nmix_modavg_lam_predict_df_bd14 <- data.frame(Predicted = nmix_modavg_lam_predict_bd14$mod.avg.pred,
                                               lower = nmix_modavg_lam_predict_bd14$lower.CL,
                                               upper = nmix_modavg_lam_predict_bd14$upper.CL,
                                               site_cov)

# Look at first values
head(nmix_modavg_lam_predict_df_bd14)

##### plot the nmixture model estimate
# First, set-up a new dataframe to predict along a sequence of the covariate.
# Predicting requires all covariates, so let's hold the other covariates constant at their mean value
nmix_alt_newdata <- data.frame(alt = seq(300, 3300, by = 50),
                               wor = mean(site_cov$wor), # hold other variables constant
                               flow = mean(site_cov$flow)) # hold other variables constant

# Model-averaged prediction of occupancy and confidence interval
nmix_alt_pred <- modavgPred(list(nmix_full = bd_14_x2, nmix_sub = bd_14_x3),
                            # c.hat =    # to change variance inflation factor, default = 1) 
                            parm.type = "lambda", # psi = occupancy
                            newdata = nmix_alt_newdata)[c("mod.avg.pred",
                                                          "lower.CL",
                                                          "upper.CL")]

# Put prediction, confidence interval, and covariate values together in a data frame
nmix_alt_pred_df_bd14 <- data.frame(Predicted = nmix_alt_pred$mod.avg.pred,
                                     lower = nmix_alt_pred$lower.CL,
                                     upper = nmix_alt_pred$upper.CL,
                                     nmix_alt_newdata)

# Plot the relationship
nmix_alt_pred_plot <- ggplot(nmix_alt_pred_df_bd14, aes(x = alt, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5, linetype = "dashed") +
  geom_path(size = 1) +
  labs(x = "Altitude", y = "Abundance Estimate") +
  theme_bw() +
  #coord_cartesian(ylim = c(0,1)) +
  theme(text = element_text(family = "HelveticaNeue", colour = "black", size=18),
        axis.text = element_text(colour = "black", size=14))
nmix_alt_pred_plot

#######################################################################################################
#####estimates for year 2015 modelling with habitat factors
bd_nm_15<-unmarkedFramePCount(y=bd[,5:7], siteCovs = site_cov)
head(bd_nm_15)
summary(bd_nm_15)
bd_15_x1<- pcount(~1~ flow,bd_nm_15, K=20)
bd_15_x1
bd_15_x2<- pcount(~1~ alt+ wor +flow ,bd_nm_15, K=20)
bd_15_x2
bd_15_x3<- pcount(~1~ flow+ wor,bd_nm_15, K=20) ## Best fit model
bd_15_x3
bd_15_x4<- pcount(~1~ alt+flow,bd_nm_15, K=20)
bd_15_x4

#### dreadge to get the global estimate and compare the models
Nmix_dredge <- dredge(global.model = bd_15_x2, rank = "AICc") # the rank argument can use AICc, QAICc, or others (see help page)
                      
Nmix_dredge ##the top model has a weightage of 0.923 and the differance of delta >5

##### predict estimate 
nmix_modavg_lam_predict_bd15 <- modavgPred(list(nmix_full = bd_15_x1, nmix_sub = bd_15_x3, nmix_sub1 = bd_15_x4),
                                            # c.hat = 1, # to change variance inflation factor, default = 1) 
                                            parm.type = "lambda", # psi = occupancy, can also be "detect" for detection probability
                                            newdata = bd_nm_15@siteCovs)[c("mod.avg.pred",  "lower.CL",  "upper.CL")]


## Put predictions, CI, and all site covariates into one data frame
nmix_modavg_lam_predict_df_bd15 <- data.frame(Predicted = nmix_modavg_lam_predict_bd15$mod.avg.pred,
                                               lower = nmix_modavg_lam_predict_bd15$lower.CL,
                                               upper = nmix_modavg_lam_predict_bd15$upper.CL,
                                               site_cov)

# Look at first values
head(nmix_modavg_lam_predict_df_bd15)

##### plot the nmixture model estimate
# First, set-up a new dataframe to predict along a sequence of the covariate.
# Predicting requires all covariates, so let's hold the other covariates constant at their mean value
nmix_alt_newdata <- data.frame(alt = mean(site_cov$alt),
                               wor = mean(site_cov$wor), # hold other variables constant
                               flow = seq(0,3, by=0.1)) # hold other variables constant

# Model-averaged prediction of occupancy and confidence interval
nmix_alt_pred_bd15 <- modavgPred(list(nmix_full = bd_15_x1, nmix_sub = bd_15_x3, nmix_sub1 = bd_15_x4),
                                  # c.hat =    # to change variance inflation factor, default = 1) 
                                  parm.type = "lambda", # psi = occupancy
                                  newdata = nmix_alt_newdata)[c("mod.avg.pred",  "lower.CL", "upper.CL")]

# Put prediction, confidence interval, and covariate values together in a data frame
nmix_alt_pred_df_bd15 <- data.frame(Predicted = nmix_alt_pred_bd15$mod.avg.pred,
                                     lower = nmix_alt_pred_bd15$lower.CL,
                                     upper = nmix_alt_pred_bd15$upper.CL,
                                     nmix_alt_newdata)

# Plot the relationship
nmix_alt_pred_bd15_plot <- ggplot(nmix_alt_pred_df_bd15, aes(x = flow, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5, linetype = "dashed") +
  geom_path(size = 1) +
  labs(x = "Flow-score", y = "Abundance Estimate") +
  theme_bw() +
  #coord_cartesian(ylim = c(0,1)) +
  theme(text = element_text(family = "HelveticaNeue", colour = "black", size=18),
        axis.text = element_text(colour = "black", size=14))
nmix_alt_pred_bd15_plot

##### naive estimates for year 2016
bd_nm_16<-unmarkedFramePCount(y=bd[,8:10], siteCovs = site_cov)
head(bd_nm_16)
summary(bd_nm_16)
bd_16_x1<- pcount(~1~ flow,bd_nm_16, K=20)
bd_16_x1
bd_16_x2<- pcount(~1~ alt+ wor +flow ,bd_nm_16, K=20)
bd_16_x2
bd_16_x3<- pcount(~1~ flow+ wor,bd_nm_16, K=20) ## Best fit model
bd_16_x3
bd_16_x4<- pcount(~1~ alt+flow,bd_nm_16, K=20)
bd_16_x4

#### dreadge to get the global estimate and compare the models
Nmix_dredge_16 <- dredge(global.model = bd_16_x2,rank = "AICc") # the rank argument can use AICc, QAICc, or others (see help page)

Nmix_dredge_16 ##the top model has a weightage of 0.923 and the differance of delta >5

##### predict estimate 
nmix_modavg_lam_predict_bd16 <- modavgPred(list(nmix_full = bd_16_x2, nmix_sub = bd_16_x3, nmix_sub1 = bd_15_x4),
                                            # c.hat = 1, # to change variance inflation factor, default = 1) 
                                            parm.type = "lambda", # psi = occupancy, can also be "detect" for detection probability
                                            newdata = bd_nm_16@siteCovs)[c("mod.avg.pred",  "lower.CL",  "upper.CL")]


## Put predictions, CI, and all site covariates into one data frame
nmix_modavg_lam_predict_df_bd16 <- data.frame(Predicted = nmix_modavg_lam_predict_bd16$mod.avg.pred,
                                               lower = nmix_modavg_lam_predict_bd16$lower.CL,
                                               upper = nmix_modavg_lam_predict_bd16$upper.CL,
                                               site_cov)

# Look at first values
head(nmix_modavg_lam_predict_df_bd16)

##### plot the nmixture model estimate
# First, set-up a new dataframe to predict along a sequence of the covariate.
# Predicting requires all covariates, so let's hold the other covariates constant at their mean value
nmix_alt_newdata <- data.frame(alt = mean(site_cov$alt),
                               wor = mean(site_cov$wor), # hold other variables constant
                               flow = seq(0,3, by=0.1)) # hold other variables constant

# Model-averaged prediction of occupancy and confidence interval
nmix_alt_pred_bd16 <- modavgPred(list(nmix_full = bd_16_x2, nmix_sub = bd_16_x3),
                                  # c.hat =    # to change variance inflation factor, default = 1) 
                                  parm.type = "lambda", # psi = occupancy
                                  newdata = nmix_alt_newdata)[c("mod.avg.pred",  "lower.CL", "upper.CL")]

# Put prediction, confidence interval, and covariate values together in a data frame
nmix_alt_pred_df_bd16 <- data.frame(Predicted = nmix_alt_pred_bd16$mod.avg.pred,
                                     lower = nmix_alt_pred_bd16$lower.CL,
                                     upper = nmix_alt_pred_bd16$upper.CL,
                                     nmix_alt_newdata)

# Plot the relationship
nmix_alt_pred_bd16_plot <- ggplot(nmix_alt_pred_df_bd16, aes(x = flow, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5, linetype = "dashed") +
  geom_path(size = 1) +
  labs(x = "Flow-score", y = "Abundance Estimate") +
  theme_bw() +
  #coord_cartesian(ylim = c(0,1)) +
  theme(text = element_text(family = "HelveticaNeue", colour = "black", size=18),
        axis.text = element_text(colour = "black", size=14))
nmix_alt_pred_bd16_plot
