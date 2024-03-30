library(unmarked)
library(dplyr)
library(ggplot2)
library(MuMIn)
library(AICcmodavg)
library(cowplot)

######################################################################################################
####################### Plumbeous Water Redstart
#setwd("E:/work folder/manuscripts/river bird co-occurrence and encounter fate/Bird Study")
pwr<-read.csv("./Data/PWR_multi.csv",header=T)
head(pwr)
dim(pwr)
site_cov <- data.frame(alt=pwr$Alt,flow=(pwr$Flow),veg=as.factor(pwr$Bank_veg),wor=pwr$WoR )

###naive estimates for year 2014 modelling with habitat factors
pwr_nm_14<-unmarkedFramePCount(y=pwr[,2:4], siteCovs = site_cov)
head(pwr_nm_14)
summary(pwr_nm_14)
pwr_14_x1<- pcount(~1~ wor,pwr_nm_14, K=20)
pwr_14_x1
pwr_14_x2<- pcount(~1~ alt+ wor +flow ,pwr_nm_14, K=20)
pwr_14_x2
pwr_14_x3<- pcount(~1~ alt+ wor,pwr_nm_14, K=20) ## Best fit model
pwr_14_x3
pwr_14_x4<- pcount(~1~ alt,pwr_nm_14, K=20)
pwr_14_x4
pwr_14_n<- pcount(~1~ 1,pwr_nm_14, K=20)
pwr_14_n

# fl <- fitList(Null=pwr_14_n, rw=pwr_14_x1, All=pwr_14_x2, 
#               rwelev = pwr_14_x3, elev=pwr_14_x4)
# fl
# ms <- modSel(fl, nullmod="Null")
# ms
#### dreadge to get the global estimate and compare the models
Nmix_dredge <- dredge(global.model = pwr_14_x2,
                      # the rank argument can use AICc, QAICc, or others (see help page)
                      rank = "AICc")
Nmix_dredge ##the top model has a weightage of 0.923 and the differance of delta >5

# ####site-wise estimates
# ranef(pwr_14_x3)
# 
# backTransform(pwr_14_x1,type="det")
# backTransform(pwr_14_x1,type="state") #does not run
# 
# lc <- linearComb(pwr_14_x3, c(1, x1= 2000, x4= 10), type="state") # Estimate abundance on the log scale when forest=0
# backTransform(lc) 
# plogis(coef(pwr_14_x1, type="det"))

##### predict estimate 
nmix_modavg_lam_predict_pwr14 <- modavgPred(list(nmix_full = pwr_14_x2, nmix_sub = pwr_14_x3),
                                      # c.hat = 1, # to change variance inflation factor, default = 1) 
                                      parm.type = "lambda", # psi = occupancy, can also be "detect" for detection probability
                                      newdata = pwr_nm_14@siteCovs)[c("mod.avg.pred",  "lower.CL",  "upper.CL")]


## Put predictions, CI, and all site covariates into one data frame
nmix_modavg_lam_predict_df_pwr14 <- data.frame(Predicted = nmix_modavg_lam_predict_pwr14$mod.avg.pred,
                                         lower = nmix_modavg_lam_predict_pwr14$lower.CL,
                                         upper = nmix_modavg_lam_predict_pwr14$upper.CL,
                                         site_cov)

# Look at first values
head(nmix_modavg_lam_predict_df_pwr14)

##### plot the nmixture model estimate
# First, set-up a new dataframe to predict along a sequence of the covariate.
# Predicting requires all covariates, so let's hold the other covariates constant at their mean value
nmix_alt_newdata <- data.frame(alt = seq(300, 3300, by = 50),
                                  wor = mean(site_cov$wor), # hold other variables constant
                                  flow = mean(site_cov$flow)) # hold other variables constant

# Model-averaged prediction of occupancy and confidence interval
nmix_alt_pred <- modavgPred(list(nmix_full = pwr_14_x2, nmix_sub = pwr_14_x3),
                               # c.hat =    # to change variance inflation factor, default = 1) 
                               parm.type = "lambda", # psi = occupancy
                               newdata = nmix_alt_newdata)[c("mod.avg.pred",
                                                                "lower.CL",
                                                                "upper.CL")]

# Put prediction, confidence interval, and covariate values together in a data frame
nmix_alt_pred_df_pwr14 <- data.frame(Predicted = nmix_alt_pred$mod.avg.pred,
                                  lower = nmix_alt_pred$lower.CL,
                                  upper = nmix_alt_pred$upper.CL,
                                  nmix_alt_newdata)

# Plot the relationship
nmix_alt_pred_pwr14_plot <- ggplot(nmix_alt_pred_df_pwr14, aes(x = alt, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5, linetype = "dashed") +
  geom_path(size = 1) +
  labs(x = "Altitude", y = "Abundance Estimate", title = "PWR_summer_2014") +
  theme_bw() +
  #coord_cartesian(ylim = c(0,1)) +
  theme(text = element_text(family = "HelveticaNeue", colour = "black", size=18),
        axis.text = element_text(colour = "black", size=14))
nmix_alt_pred_pwr14_plot

#######################################################################################################
#####estimates for year 2015 modelling with habitat factors
pwr_nm_15<-unmarkedFramePCount(y=pwr[,5:7], siteCovs = site_cov)
head(pwr_nm_15)
summary(pwr_nm_15)
pwr_15_x1<- pcount(~1~ wor,pwr_nm_15, K=20)
pwr_15_x1
pwr_15_x2<- pcount(~1~ alt+ wor +flow ,pwr_nm_15, K=20)
pwr_15_x2
pwr_15_x3<- pcount(~1~ alt+ wor,pwr_nm_15, K=20) ## Best fit model
pwr_15_x3
pwr_15_x4<- pcount(~1~ alt+flow,pwr_nm_15, K=20)
pwr_15_x4
pwr_15_n<- pcount(~1~ 1,pwr_nm_15, K=20)
pwr_15_n

#### dreadge to get the global estimate and compare the models
Nmix_dredge <- dredge(global.model = pwr_15_x2,
                      # the rank argument can use AICc, QAICc, or others (see help page)
                      rank = "AICc")
Nmix_dredge ##the top model has a weightage of 0.923 and the differance of delta >5

##### predict estimate 
nmix_modavg_lam_predict_pwr15 <- modavgPred(list(nmix_full = pwr_15_x4, nmix_sub = pwr_15_x3),
                                            # c.hat = 1, # to change variance inflation factor, default = 1) 
                                            parm.type = "lambda", # psi = occupancy, can also be "detect" for detection probability
                                            newdata = pwr_nm_15@siteCovs)[c("mod.avg.pred",  "lower.CL",  "upper.CL")]


## Put predictions, CI, and all site covariates into one data frame
nmix_modavg_lam_predict_df_pwr15 <- data.frame(Predicted = nmix_modavg_lam_predict_pwr15$mod.avg.pred,
                                               lower = nmix_modavg_lam_predict_pwr15$lower.CL,
                                               upper = nmix_modavg_lam_predict_pwr15$upper.CL,
                                               site_cov)

# Look at first values
head(nmix_modavg_lam_predict_df_pwr15)

##### plot the nmixture model estimate
# First, set-up a new dataframe to predict along a sequence of the covariate.
# Predicting requires all covariates, so let's hold the other covariates constant at their mean value
nmix_alt_newdata <- data.frame(alt = seq(300, 3300, by = 50),
                               wor = mean(site_cov$wor), # hold other variables constant
                               flow = mean(site_cov$flow)) # hold other variables constant

# Model-averaged prediction of occupancy and confidence interval
nmix_alt_pred_pwr15 <- modavgPred(list(nmix_full = pwr_15_x4, nmix_sub = pwr_15_x3),
                            # c.hat =    # to change variance inflation factor, default = 1) 
                            parm.type = "lambda", # psi = occupancy
                            newdata = nmix_alt_newdata)[c("mod.avg.pred",  "lower.CL", "upper.CL")]

# Put prediction, confidence interval, and covariate values together in a data frame
nmix_alt_pred_df_pwr15 <- data.frame(Predicted = nmix_alt_pred_pwr15$mod.avg.pred,
                                     lower = nmix_alt_pred_pwr15$lower.CL,
                                     upper = nmix_alt_pred_pwr15$upper.CL,
                                     nmix_alt_newdata)

# Plot the relationship
nmix_alt_pred_pwr15_plot <- ggplot(nmix_alt_pred_df_pwr15, aes(x = alt, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5, linetype = "dashed") +
  geom_path(size = 1) +
  labs(x = "Altitude", y = "Abundance Estimate", title = "PWR_summer_2015") +
  theme_bw() +
  coord_cartesian(ylim = c(0,5)) +
  theme(text = element_text(family = "HelveticaNeue", colour = "black", size=18),
        axis.text = element_text(colour = "black", size=14))
nmix_alt_pred_pwr15_plot

##### naive estimates for year 2016
pwr_nm_16<-unmarkedFramePCount(y=pwr[,8:10], siteCovs = site_cov)
head(pwr_nm_16)
summary(pwr_nm_16)
pwr_16_x1<- pcount(~1~ wor,pwr_nm_16, K=20)
pwr_16_x1
pwr_16_x2<- pcount(~1~ alt+ wor +flow ,pwr_nm_16, K=20)
pwr_16_x2
pwr_16_x3<- pcount(~1~ alt+ wor,pwr_nm_16, K=20) ## Best fit model
pwr_16_x3
pwr_16_x4<- pcount(~1~ alt+flow,pwr_nm_16, K=20)
pwr_16_x4

#### dreadge to get the global estimate and compare the models
Nmix_dredge_16 <- dredge(global.model = pwr_16_x2,rank = "AICc") # the rank argument can use AICc, QAICc, or others (see help page)
                      
Nmix_dredge_16 ##the top model has a weightage of 0.923 and the differance of delta >5

##### predict estimate 
nmix_modavg_lam_predict_pwr16 <- modavgPred(list(nmix_full = pwr_16_x2, nmix_sub = pwr_16_x3),
                                            # c.hat = 1, # to change variance inflation factor, default = 1) 
                                            parm.type = "lambda", # psi = occupancy, can also be "detect" for detection probability
                                            newdata = pwr_nm_16@siteCovs)[c("mod.avg.pred",  "lower.CL",  "upper.CL")]


## Put predictions, CI, and all site covariates into one data frame
nmix_modavg_lam_predict_df_pwr16 <- data.frame(Predicted = nmix_modavg_lam_predict_pwr16$mod.avg.pred,
                                               lower = nmix_modavg_lam_predict_pwr16$lower.CL,
                                               upper = nmix_modavg_lam_predict_pwr16$upper.CL,
                                               site_cov)

# Look at first values
head(nmix_modavg_lam_predict_df_pwr16)

##### plot the nmixture model estimate
# First, set-up a new dataframe to predict along a sequence of the covariate.
# Predicting requires all covariates, so let's hold the other covariates constant at their mean value
nmix_alt_newdata <- data.frame(alt = seq(300, 3300, by = 50),
                               wor = mean(site_cov$wor), # hold other variables constant
                               flow = mean(site_cov$flow)) # hold other variables constant

# Model-averaged prediction of occupancy and confidence interval
nmix_alt_pred_pwr16 <- modavgPred(list(nmix_full = pwr_16_x2, nmix_sub = pwr_16_x3),
                                  # c.hat =    # to change variance inflation factor, default = 1) 
                                  parm.type = "lambda", # psi = occupancy
                                  newdata = nmix_alt_newdata)[c("mod.avg.pred",  "lower.CL", "upper.CL")]

# Put prediction, confidence interval, and covariate values together in a data frame
nmix_alt_pred_df_pwr16 <- data.frame(Predicted = nmix_alt_pred_pwr16$mod.avg.pred,
                                     lower = nmix_alt_pred_pwr16$lower.CL,
                                     upper = nmix_alt_pred_pwr16$upper.CL,
                                     nmix_alt_newdata)

# Plot the relationship
nmix_alt_pred_pwr16_plot <- ggplot(nmix_alt_pred_df_pwr16, aes(x = alt, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5, linetype = "dashed") +
  geom_path(size = 1) +
  labs(x = "Altitude", y = "Abundance Estimate", title = "PWR_summer_2016") +
  theme_bw() +
  coord_cartesian(ylim = c(0,5)) +
  theme(text = element_text(family = "HelveticaNeue", colour = "black", size=18),
        axis.text = element_text(colour = "black", size=14))
nmix_alt_pred_pwr16_plot

##### naive estimates for year 2017

pwr_nm_17<-unmarkedFramePCount(y=pwr[,11:13], siteCovs = site_cov)
head(pwr_nm_17)
summary(pwr_nm_17)
pwr_17_x1<- pcount(~1~ wor,pwr_nm_17, K=20)
pwr_17_x1
pwr_17_x2<- pcount(~1~ alt+ wor +flow ,pwr_nm_17, K=20)
pwr_17_x2
pwr_17_x3<- pcount(~1~ alt+ wor,pwr_nm_17, K=20) ## Best fit model
pwr_17_x3
pwr_17_x4<- pcount(~1~ alt+flow,pwr_nm_17, K=20)
pwr_17_x4

#### dreadge to get the global estimate and compare the models
Nmix_dredge_17 <- dredge(global.model = pwr_17_x2,rank = "AICc") # the rank argument can use AICc, QAICc, or others (see help page)

Nmix_dredge_17 ##the top model has a weightage of 0.923 and the differance of delta >5

##### predict estimate 
nmix_modavg_lam_predict_pwr17 <- modavgPred(list(nmix_full = pwr_17_x4, nmix_sub = pwr_17_x3),
                                            # c.hat = 1, # to change variance inflation factor, default = 1) 
                                            parm.type = "lambda", # psi = occupancy, can also be "detect" for detection probability
                                            newdata = pwr_nm_17@siteCovs)[c("mod.avg.pred",  "lower.CL",  "upper.CL")]


## Put predictions, CI, and all site covariates into one data frame
nmix_modavg_lam_predict_df_pwr17 <- data.frame(Predicted = nmix_modavg_lam_predict_pwr17$mod.avg.pred,
                                               lower = nmix_modavg_lam_predict_pwr17$lower.CL,
                                               upper = nmix_modavg_lam_predict_pwr17$upper.CL,
                                               site_cov)

# Look at first values
head(nmix_modavg_lam_predict_df_pwr17)

##### plot the nmixture model estimate
# First, set-up a new dataframe to predict along a sequence of the covariate.
# Predicting requires all covariates, so let's hold the other covariates constant at their mean value
nmix_alt_newdata_pwr <- data.frame(alt = seq(300, 3300, by = 50),
                               wor = mean(site_cov$wor), # hold other variables constant
                               flow = mean(site_cov$flow)) # hold other variables constant

# Model-averaged prediction of occupancy and confidence interval
nmix_alt_pred_pwr17 <- modavgPred(list(nmix_full = pwr_17_x4, nmix_sub = pwr_17_x3),
                                  # c.hat =    # to change variance inflation factor, default = 1) 
                                  parm.type = "lambda", # psi = occupancy
                                  newdata = nmix_alt_newdata_pwr)[c("mod.avg.pred",  "lower.CL", "upper.CL")]

# Put prediction, confidence interval, and covariate values together in a data frame
nmix_alt_pred_df_pwr17 <- data.frame(Predicted = nmix_alt_pred_pwr17$mod.avg.pred,
                                     lower = nmix_alt_pred_pwr17$lower.CL,
                                     upper = nmix_alt_pred_pwr17$upper.CL,
                                     nmix_alt_newdata_pwr)

# Plot the relationship
nmix_alt_pred_pwr17_plot <- ggplot(nmix_alt_pred_df_pwr17, aes(x = alt, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5, linetype = "dashed") +
  geom_path(size = 1) +
  labs(x = "Altitude", y = "Abundance Estimate", title = "PWR_summer_2017") +
  theme_bw() +
  coord_cartesian(ylim = c(0,5)) +
  theme(text = element_text(family = "HelveticaNeue", colour = "black", size=18),
        axis.text = element_text(colour = "black", size=14))
nmix_alt_pred_pwr17_plot

plot_grid(nmix_alt_pred_pwr14_plot,nmix_alt_pred_pwr15_plot,nmix_alt_pred_pwr16_plot,nmix_alt_pred_pwr17_plot)

ggsave("PWR_years.jpeg", width = 12, height = 8, units = "in", dpi=300)
##### naive estimates for year 2018

head(pwr[,14:16])
pwr_nm_18<-unmarkedFramePCount(y=pwr[,14:16])
summary(pwr_nm_18)
pwr_18<- pcount(~1 ~1, pwr_nm_18, K=20)
pwr_18

backTransform(pwr_18, type= "det")
backTransform(pwr_18, type= "state")



####Grey Wagtail
gw<-read.csv("GW_multi.csv",header = T)
head(gw)
gw_h<-left_join(gw,hab,by="Site.ID")
head(gw_h)
#####estimates for year 2014
head(gw[,2:4])
gw_nm_14<-unmarkedFramePCount(y=gw[,2:4])
summary(gw_nm_14)
gw_14 <- pcount(~1 ~1, gw_nm_14, K=20)
gw_14
backTransform(gw_14,type="det")
backTransform(gw_14,type = "state")

#####estimates for year 2015
head(gw[,5:7])
gw_nm_15<-unmarkedFramePCount(y=gw[,5:7])
summary(gw_nm_15)
gw_15 <- pcount(~1 ~1, gw_nm_15, K=20)
gw_15

backTransform(gw_15,type="det")
backTransform(gw_15,type = "state")

####estimates for year 2016
head(gw[,7:9])
gw_nm_16<-unmarkedFramePCount(y=gw[,7:10])
summary(gw_nm_16)
gw_16 <- pcount(~1 ~1, gw_nm_16, K=20)
gw_16
backTransform(gw_16,type="det")
backTransform(gw_16,type = "state")


#####estimates for year 2017
head(gw[,11:13])
gw_nm_17<-unmarkedFramePCount(y=gw[,11:13])
summary(gw_nm_17)
gw_17 <- pcount(~1 ~1, gw_nm_17, K=20)
gw_17

backTransform(gw_17,type="det")
backTransform(gw_17,type = "state")

####estimates for year 2018
head(gw[,14:16])
gw_nm_18<-unmarkedFramePCount(y=gw[,14:16])
summary(gw_nm_18)
gw_18 <- pcount(~1 ~1, gw_nm_17, K=20)
gw_18
backTransform(gw_18,type="det")
backTransform(gw_18,type = "state")

#####White-capped Redstart
wcr<-read.csv("WCR_multi.csv",header = T)

#### estimates for year 2014
head(wcr[,2:4])
wcr_nm_14<-unmarkedFramePCount(y=wcr[,2:4])
summary(wcr_nm_14)
wcr_14 <- pcount(~1 ~1, wcr_nm_14, K=20)
wcr_14
backTransform(wcr_14,type="det")
backTransform(wcr_14,type = "state")

#### estimates for year 2015
head(wcr[,5:7])
wcr_nm_15<-unmarkedFramePCount(y=wcr[,5:7])
summary(wcr_nm_15)
wcr_15 <- pcount(~1 ~1, wcr_nm_15, K=20)
wcr_15
backTransform(wcr_15,type="det")
backTransform(wcr_15,type = "state")

#### estimates for year 2016
head(wcr[,8:10])
wcr_nm_16<-unmarkedFramePCount(y=wcr[,8:10])
summary(wcr_nm_16)
wcr_16 <- pcount(~1 ~1, wcr_nm_16, K=20)
wcr_16
backTransform(wcr_16,type="det")
backTransform(wcr_16,type = "state")

#### estimates for year 2017
head(wcr[,11:13])
wcr_nm_17<-unmarkedFramePCount(y=wcr[,11:13])
summary(wcr_nm_16)
wcr_17 <- pcount(~1 ~1, wcr_nm_16, K=20)
wcr_17
backTransform(wcr_17,type="det")
backTransform(wcr_17,type = "state")

#### estimates for year 2018
head(wcr[,14:16])
wcr_nm_18<-unmarkedFramePCount(y=wcr[,14:16])
summary(wcr_nm_18)
wcr_18 <- pcount(~1 ~1, wcr_nm_18, K=20)
wcr_18
backTransform(wcr_18,type="det")
backTransform(wcr_18,type = "state")

#####Brown Dipper
bd<-read.csv("BD_multi.csv",header = T)

#### estimates for year 2014
head(bd[,2:4])
bd_nm_14<-unmarkedFramePCount(y=bd[,2:4])
summary(bd_nm_14)
bd_14 <- pcount(~1 ~1, bd_nm_14, K=20)
bd_14
backTransform(bd_14,type="det")
backTransform(bd_14,type = "state")

#### estimates for year 2015
head(bd[,5:7])
bd_nm_15<-unmarkedFramePCount(y=bd[,5:7])
summary(bd_nm_15)
bd_15 <- pcount(~1 ~1, bd_nm_15, K=20)
bd_15
backTransform(bd_15,type="det")
backTransform(bd_15,type = "state")

#### estimates for year 2016
head(bd[,8:10])
bd_nm_16<-unmarkedFramePCount(y=bd[,8:10])
summary(bd_nm_16)
bd_16 <- pcount(~1 ~1, bd_nm_16, K=20)
bd_16
backTransform(bd_16,type="det")
backTransform(bd_16,type = "state")

#### estimates for year 2017
head(bd[,11:13])
bd_nm_17<-unmarkedFramePCount(y=bd[,11:13])
summary(bd_nm_17)
bd_17 <- pcount(~1 ~1, bd_nm_17, K=20)
bd_17
backTransform(bd_17,type="det")
backTransform(bd_17,type = "state")

#### estimates for year 2018
head(bd[,14:16])
bd_nm_18<-unmarkedFramePCount(y=bd[,14:16])
summary(bd_nm_18)
bd_18 <- pcount(~1 ~1, bd_nm_18, K=20)
bd_18
backTransform(bd_18,type="det")
backTransform(bd_18,type = "state")


######Crested Kingfisher
ck<-read.csv("CK_multi.csv",header = T)

#### estimates for year 2014
head(ck[,2:4])
ck_nm_14<-unmarkedFramePCount(y=ck[,2:4])
summary(ck_nm_14)
ck_14 <- pcount(~1 ~1, ck_nm_14, K=20)
ck_14
backTransform(ck_14,type="det")
backTransform(ck_14,type = "state")

ck_14_abd<-(ranef(ck_14))
ck_14_site_mn<-bup(ck_14_abd,stat = "mean")
str(ck_14_abd)
head(ck_14_site_mn)
length(ck_14_site_mn)

#### estimates for year 2015
head(ck[,5:7])
ck_nm_15<-unmarkedFramePCount(y=ck[,5:7])
summary(ck_nm_15)
ck_15 <- pcount(~1 ~1, ck_nm_15, K=20)
ck_15
backTransform(ck_15,type="det")
backTransform(ck_15,type = "state")

#### estimates for year 2016
head(ck[,8:10])
ck_nm_16<-unmarkedFramePCount(y=ck[,8:10])
summary(ck_nm_16)
ck_16 <- pcount(~1 ~1, ck_nm_16, K=20)
ck_16
backTransform(ck_16,type="det")
backTransform(ck_16,type = "state")

#### estimates for year 2017
head(ck[,11:13])
ck_nm_17<-unmarkedFramePCount(y=ck[,11:13])
summary(ck_nm_17)
ck_17<- pcount(~1 ~1, ck_nm_17, K=20)
ck_17
backTransform(ck_17,type="det")
backTransform(ck_17,type = "state")

#### estimates for year 2018
head(ck[,14:16])
ck_nm_18<-unmarkedFramePCount(y=ck[,14:16])
summary(ck_nm_18)
ck_18 <- pcount(~1 ~1, ck_nm_18, K=20)
ck_18
backTransform(ck_18,type="det")
backTransform(ck_18,type = "state")
