library(dplyr) # for data organization
library(magrittr) # for piping
library(ggplot2) # for plotting
# install.packages("unmarked") # first time only
library(unmarked) # for occupancy models
library(AICcmodavg) # For GOF tests

# Load detection history (43 sites with 3 visits each for 5 years)
pwr_det_hist <- read.csv("./Data/PWR_multi.csv",header=T) 
bd_det_hist <- read.csv("./Data/BD_multi.csv",header=T) 
wcr_det_hist <- read.csv("./Data/WCR_multi.csv",header=T) 
bwt_det_hist <- read.csv("./Data/BWT_multi.csv",header=T) 
ck_det_hist <- read.csv("./Data/CK_multi.csv",header=T) 
sf_det_hist <- read.csv("./Data/SF_multi.csv",header=T) 
cmk_det_hist <- read.csv("./Data/Common_Kingfisher_multi.csv",header=T) 
wbw_det_hist <- read.csv("./Data/White browed Wagtail_multi.csv",header=T) 
ww_det_hist <- read.csv("./Data/White wagtail_multi.csv",header=T) 
wtk_det_hist <- read.csv("./Data/WTK_multi.csv",header=T) 
gw_det_hist <- read.csv("./Data/GW_multi.csv",header=T) 
lf_det_hist <- read.csv("./Data/LF_multi.csv",header=T) 

# Examine data
head(pwr_det_hist)

# Load covariate data
site_cov <- data.frame(alt=pwr$Alt,flow=(pwr$Flow),veg=as.factor(pwr$Bank_veg),wor=pwr$WoR )

# Build an unmarkedmultFramOccu
pwr_unmarkedMultFrame <- unmarkedFramePCO( # y is a matrix with observed detection history 
  # 0's and 1's, one row per site, one column per survey
  y = as.matrix(pwr_det_hist[,2:16]),
  # numPrimary is the number of primary surveys; 
  # function derives # of secondary surveys based on
  # ncol(y) / numPrimary
  numPrimary = 5,
  # obsCovs = observation covariates in a list,  # each variable has site rows x survey columns
  #obsCovs = list(effort = effort),
  # siteCovs = dataframe with site rows x column variables..the first column of site_cov is the site number
  siteCovs = site_cov,
  ) 

# S4 class for colext occupancy model data
summary(pwr_unmarkedMultFrame)

## ----buildcolext---------------------------------------------------------
dynamic_nmix_m1 <- pcountOpen(
  # Psi depends on initial habitat, climate
  lambdaformula= ~alt+flow+wor,
  # colonization depends on wetland change
  gammaformula = ~1, 
  # extinction depends on wetland change
  omegaformula = ~1,
  # detection depends on survey effort
  pformula = ~1,
  # data must be a unmarkedMultFrame
  data = pwr_unmarkedMultFrame, dynamics = "trend",
  # method is optim method, leave as "BFGS"
  method = "BFGS", K=50)

#### dreadge to get the global estimate and compare the models
Nmixopen_dredge <- dredge(global.model = dynamic_nmix_m1,
                      # the rank argument can use AICc, QAICc, or others (see help page)
                      rank = "AIC")
Nmixopen_dredge ##the top model has a weightage of 0.923 and the differance of delta >5

## ----modelsummary--------------------------------------------------------
summary(dynamic_nmix_m1)
lam <- coef(backTransform(dynamic_nmix_m1, "lambda")) # or
lam <- exp(coef(dynamic_nmix_m1, type="lambda"))
gam <- exp(coef(dynamic_nmix_m1, type="gamma"))
om <- plogis(coef(dynamic_nmix_m1, type="omega"))
p <- plogis(coef(dynamic_nmix_m1, type="det"))

pwr_abd <- data_frame(year = 2014:2018)

abd_estimates <- ranef(dynamic_nmix_m1)
pwr_abd$mean<- colSums(bup(abd_estimates))
pwr_abd[,3:4] <- t(colSums(confint(abd_estimates)))

## ----plotderivedoccupancy----
sf_mult <- ggplot(pwr_abd, aes(x = year, y = mean)) +
  geom_errorbar(aes(ymin = V1,
                    ymax = V2),
                width = 0) +
  geom_point(size = 3) +
  geom_line() +
  labs(x = "Year", y = "Estimated Abundance", title = "Spotted Forktail") +
  theme_bw()

plot_grid(pwr_mult, bd_mult, wcr_mult, bwt_mult, ck_mult, sf_mult)

ggsave("Pop_trends_riverbirds.jpeg", width = 12, height = 8, units = "in", dpi=600)


################################################################################
data_files <-list.files("./Data/")

sp_sum <- list()
sp_abd <- list()
sp_sum_zip <- list()
all_sp_poptrend <- list()

setwd("./Data/")
# Load covariate data
site_cov <- data.frame(alt=pwr_det_hist$Alt,flow=(pwr_det_hist$Flow),wor=pwr_det_hist$WoR )
dynamic_nmix_mod_p <-list()

for(i in 1:length(data_files))
{
  spdat <- read.csv(data_files[i])
  # Build an unmarkedmultFramOccu
 unmarked_Frame <- unmarkedFramePCO(y = as.matrix(spdat[1:43,2:16]), ##two files common_kingfisher and lf has one row short
                                             numPrimary = 5,
                                             siteCovs = site_cov
  ) 
  
  ## ----build open nmixture---------------------------------------------------------
  dynamic_nmix_mod_p[[i]] <- pcountOpen(lambdaformula= ~scale(alt)+flow+wor,
                                 gammaformula = ~1,
                                 omegaformula = ~1,
                                 pformula = ~1,
                                 data = unmarked_Frame,    method = "BFGS", K=50, mixture = c("P"))
 
 #goftest <-Nmix.gof.test(dynamic_nmix_mod_p, nsim = 1000, parallel = TRUE, ncores = 4) ##gof test
 #print(goftest)
 # dynamic_nmix_mod_zip <- pcountOpen(lambdaformula= ~scale(alt)+flow+wor, 
 #                                    gammaformula = ~1,
 #                                    omegaformula = ~1,
 #                                    pformula = ~1,
 #                                    data = unmarked_Frame,    method = "BFGS", K=50, mixture = c("ZIP"))
 # 
 # print(AIC(dynamic_nmix_mod_p, dynamic_nmix_mod_zip))
 sp_unmk_sum <- summary(dynamic_nmix_mod_p[[i]])
 print(backTransform(dynamic_nmix_mod_p[[i]],type="det"))
 sp_sum[[i]] <- cbind(paste(data_files[i]), do.call(rbind, sp_unmk_sum))
 
 # sp_unmk_sum_zip <- summary(dynamic_nmix_mod_zip)
 # sp_sum_zip[[i]] <- cbind(paste(data_files[i]), do.call(rbind, sp_unmk_sum_zip))
 # 
 #### Abundance estimate
 abd_est <- data_frame(year = 2014:2018, species = paste(data_files[i]))
 abd_est$mean<- colSums(bup(ranef(dynamic_nmix_mod_p[[i]])))
 abd_est[,4:5] <- t(colSums(confint(ranef(dynamic_nmix_mod_p[[i]]))))
 sp_abd[[i]] <- abd_est
 
 tempdat <- bup(ranef(dynamic_nmix_mod_p[[i]]))
 trend_est<- data.frame(trend = (tempdat[,5]/tempdat[,1])^0.2,
                       species = paste(data_files[i]))
  all_sp_poptrend[[i]]<- trend_est
  }

setwd("~/Documents/GitHub/Riverbird_Nmixture")

abd_cmbnd1 <-do.call(rbind, sp_abd)
abd_cmbnd1$species<- gsub("_multi.csv", "",abd_cmbnd1$species)
head(abd_cmbnd)
write.csv(abd_cmbnd, "./Results/All_species_abundance_combined_1906.csv", row.names = F)

abd_cmbnd <-read.csv( "./Results/All_species_abundance_combined.csv")

coef_cmbnd <-do.call(rbind, sp_sum)
colnames(coef_cmbnd)[1] <- c("species")
coef_cmbnd$species<- gsub("_multi.csv", "",coef_cmbnd$species)
coef_cmbnd$variable <- rep(rownames(coef_cmbnd)[1:8], 12)
coef_cmbnd$variable <- gsub("det", "det_prob",coef_cmbnd$variable)
coef_cmbnd$variable <- gsub("lambda.(Intercept)", "lambda_intercept",coef_cmbnd$variable, fixed = TRUE)
coef_cmbnd$variable <- gsub("lambda_elevation", "lambda_elev",coef_cmbnd$variable, fixed = TRUE)
coef_cmbnd$species[coef_cmbnd$species=="BD"] <- "Brown_Dipper"
coef_cmbnd$species[coef_cmbnd$species=="WTK"] <- "White-Throated\nKingfisher"
coef_cmbnd$species[coef_cmbnd$species=="BWT"] <- "Blue-Whistling\nThrush"
coef_cmbnd$species[coef_cmbnd$species=="Common_Kingfisher"] <- "Crested\nKingfisher"
coef_cmbnd$species[coef_cmbnd$species=="CK"] <- "Common\nKingfisher"
coef_cmbnd$species[coef_cmbnd$species=="GW"] <- "Grey_Wagtail"
coef_cmbnd$species[coef_cmbnd$species=="LF"] <- "Little_Forktail"
coef_cmbnd$species[coef_cmbnd$species=="PWR"] <- "Plumbeous\nWater Redstart"
coef_cmbnd$species[coef_cmbnd$species=="SF"] <- "Spotted\nForktail"
coef_cmbnd$species[coef_cmbnd$species=="WCR"] <- "White-capped\nRedstart"
coef_cmbnd$species[coef_cmbnd$species=="White browed Wagtail"] <- "White-browed\nWagtail"
write.csv(coef_cmbnd, "./Results/All_species_coef_combined.csv", row.names = F)

coef_cmbnd_zip <-do.call(rbind, sp_sum_zip)
colnames(coef_cmbnd_zip)[1] <- c("species")
coef_cmbnd_zip$species<- gsub("_multi.csv", "",coef_cmbnd_zip$species)
write.csv(coef_cmbnd_zip, "./Results/All_species_coef_zip_combined.csv", row.names = F)


## ----plots for publication ---
  ggplot(abd_cmbnd, aes(x = year, y = mean1/43)) +
  geom_errorbar(aes(ymin = V1/43,
                    ymax = V2/43),
                width = 0) +
  geom_point(size = 2) +
  geom_line() +
    facet_wrap(~factor(species, levels =c("Brown_Dipper"   ,   "Crested_Kingfisher" , "Little_Forktail" ,         
                                          "Plumbeous_Water_Redstart" , "Spotted_Forktail" ,  "White_capped_Redstart" ,     
                                          "Blue_Whistling_Thrush" , "Common_Kingfisher" , "Grey_Wagtail" ,                
                                          "White browed Wagtail","White Wagtail" , "White_Thraoted_Kingfisher"
    )), scales = "free", ncol=3)+
  labs(x = "Year", y = "Estimated Density per 500m") +
  theme_bw() + theme(axis.title = element_text(size=16), axis.text = element_text(size=10),
                   strip.text.x = element_text(size = 10))

ggsave("All_species_abundance_new_2804.jpeg", width=9, height=9, units = "in", dpi=300)  
  
  ggplot(coef_cmbnd, aes(x= variable, y = Estimate))+
    geom_errorbar(aes(ymin = Estimate -1.96*SE,
                      ymax = Estimate +1.96*SE),
                  width = 0) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed")+
    facet_wrap(~factor(species, levels =c("Brown_Dipper"   ,   "Crested_Kingfisher" , "Little_Forktail" ,         
                                          "Plumbeous_Water_Redstart" , "Spotted_Forktail" ,  "White_capped_Redstart" ,     
                                          "Blue_Whistling_Thrush" , "Common_Kingfisher" , "Grey_Wagtail" ,                
                                          "White browed Wagtail","White Wagtail" , "White_Thraoted_Kingfisher"
                                          )), scales = "free", ncol=3)+
    labs(x = "Variable", y = "Coefficient Estimate") +
    coord_cartesian(ylim = c(-4,4))+
    theme_bw()+
    theme(axis.title = element_text(size=16), axis.text = element_text(size=10),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0),
          strip.text.x = element_text(size = 10))
  
  ggsave("All_species_coef_new_2804.jpeg", width=9, height=12, units = "in", dpi=300)  
  
  coef_cmbnd |> filter(grepl("lambda", variable)) |> mutate(dummy = "lambda") |>
    ggplot( aes(y= dummy, x = Estimate))+
    geom_errorbar(aes(xmin = Estimate -1.96*SE,
                      xmax = Estimate +1.96*SE),
                  width = 0) +
    geom_point(size = 2) +
    geom_vline(xintercept = 0, linetype = "dashed")+
    facet_grid(factor(species, levels =c("Brown_Dipper"   ,   "Crested\nKingfisher" , "Little_Forktail" ,         
                                          "Plumbeous\nWater Redstart" , "Spotted\nForktail" ,  "White-capped\nRedstart" ,     
                                          "Blue-Whistling\nThrush" , "Common\nKingfisher" , "Grey_Wagtail" ,                
                                          "White-browed\nWagtail","White Wagtail" , "White-Throated\nKingfisher"
    )) ~variable, space = "fixed", scales = "free", switch = "y",
    labeller = labeller(variable = c("lambda_intrcpt" = "Intercept",
                             "lambda.scale(alt)" = "Elevation",
                             "lambda.flow" = "Flow-score",
                             "lambda.wor" = "River width")))+
    labs(x = "Coefficient Estimate", y = "Species") +
    coord_cartesian(ylim = c(-4,4))+
    theme_bw()+
    theme(axis.title = element_text(size=16), axis.text.y = element_blank(),
          axis.text.x = element_text(size=8),
          strip.text.x = element_text(size = 10))
  ggsave("All_species_coef_new_2606.jpeg", width=8, height=12, units = "in", dpi=300)  
  
  #####################################################
  ####### ---- Prediction plot ----
  
  # Predicting requires all covariates, so let's hold the other covariates constant at their mean value
  nmix_alt_newdata <- data.frame(alt = seq(300, 3300, by = 50),
                                 wor = mean(site_cov$wor), # hold other variables constant
                                 flow = mean(site_cov$flow)) # hold other variables constant
  
  nmix_alt_pred_df <- list()
  for(i in 1:12){
  # Model-averaged prediction of occupancy and confidence interval
  nmix_alt_pred_ck <- modavgPred(list(nmix_full = dynamic_nmix_mod_p[[i]]),
                                    # c.hat =    # to change variance inflation factor, default = 1) 
                                    parm.type = "lambda", # psi = occupancy
                                    newdata = nmix_alt_newdata)[c("mod.avg.pred",  "lower.CL", "upper.CL")]
  
  # Put prediction, confidence interval, and covariate values together in a data frame
  nmix_alt_pred_df[[i]] <- data.frame(Predicted = nmix_alt_pred_ck$mod.avg.pred,
                                       lower = nmix_alt_pred_ck$lower.CL,
                                       upper = nmix_alt_pred_ck$upper.CL,
                                       nmix_alt_newdata, 
                                      species = data_files[i])
  }
  
  nmix_alt_pred_all <- do.call(rbind, nmix_alt_pred_df)
  nmix_alt_pred_all$species<- gsub("_multi.csv", "",nmix_alt_pred_all$species)
  nmix_alt_pred_all$species[nmix_alt_pred_all$species=="BD"] <- "Brown_Dipper"
  nmix_alt_pred_all$species[nmix_alt_pred_all$species=="WTK"] <- "White_Throated_Kingfisher"
  nmix_alt_pred_all$species[nmix_alt_pred_all$species=="BWT"] <- "Blue_Whistling_Thrush"
  nmix_alt_pred_all$species[nmix_alt_pred_all$species=="CK"] <- "Crested_Kingfisher"
  nmix_alt_pred_all$species[nmix_alt_pred_all$species=="GW"] <- "Grey_Wagtail"
  nmix_alt_pred_all$species[nmix_alt_pred_all$species=="LF"] <- "Little_Forktail"
  nmix_alt_pred_all$species[nmix_alt_pred_all$species=="PWR"] <- "Plumbeous_Water_Redstart"
  nmix_alt_pred_all$species[nmix_alt_pred_all$species=="SF"] <- "Spotted_Forktail"
  nmix_alt_pred_all$species[nmix_alt_pred_all$species=="WCR"] <- "White_capped_Redstart"
  # Plot the relationship
  nmix_alt_pred_plot <- ggplot(nmix_alt_pred_all, aes(x = alt, y = Predicted)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5, linetype = "dashed") +
    geom_path(size = 1) +
    facet_wrap(~factor(species, levels =c("Brown_Dipper"  ,   "Crested_Kingfisher" , "Little_Forktail" ,         
                                          "Plumbeous_Water_Redstart" , "Spotted_Forktail" ,  "White_capped_Redstart" ,     
                                          "Blue_Whistling_Thrush" , "Common_Kingfisher" , "Grey_Wagtail" ,                
                                          "White browed Wagtail","White Wagtail" , "White_Throated_Kingfisher"
    )), scales = "free", ncol=3)+
    labs(x = "Altitude (in meters)", y = "Abundance Estimate per site") +
    theme_bw() +
    coord_cartesian(ylim = c(0,8)) +
    theme(text = element_text(family = "HelveticaNeue", colour = "black", size=18),
          axis.text = element_text(colour = "black", size=14))
  nmix_alt_pred_plot

  #cowplot::plot_grid(nmix_alt_pred_pwr_plot, nmix_alt_pred_sf_plot, nmix_alt_pred_ck_plot)  
ggsave(nmix_alt_pred_plot, filename ="All_Sp_elev_plot.jpeg", width = 9, height = 12, units = "in", dpi=300)  

##### ---- New abundance trend plot ----
abd_cmbnd <-read.csv( "./Results/All_species_abundance_combined.csv")

unique(abd_cmbnd$species)

dat <- list()

for(i in 1:12)
{
  tempabd <- bup(ranef(dynamic_nmix_mod_p[[12]]))
  (tempabd[,5]/tempabd[,1])^0.2
}
for(i in unique(abd_cmbnd1$species))
{
  tempdat <- abd_cmbnd1 |> filter(species ==i)
  dat[[i]] <- data.frame(species =i, 
                        pop_ratio = (tempdat$mean[5]/tempdat$mean[1]), 
                        pop_trend = (tempdat$mean[5]/tempdat$mean[1])^0.25,
                        ucl_trend = (tempdat$V2[5]/tempdat$V2[1])^0.25,
                        lcl_trend = (tempdat$V1[5]/tempdat$V1[1])^0.25
  ) #
}

mergd_dat1 <- do.call(rbind, dat)
mergd_dat$species <-gsub("_multi.csv", "",mergd_dat$species)
mergd_dat

ggplot(mergd_dat, aes(y = species)) +
  geom_point(aes(x = pop_trend))+
  geom_linerange(aes(xmin = lcl_trend, xmax = ucl_trend)) +
  theme_bw()
  
all_pop_trend <- do.call(rbind, all_sp_poptrend)

all_pop_trend %>%
  #left_join(sample_size) %>%
  #mutate(myaxis = paste0(name, "\n", "n=", num)) %>%
  ggplot( aes(y=species, x=trend), fill +"grey") +
  geom_violin(width=1.4) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  viridis::scale_fill_viridis(discrete = TRUE) +
  theme_bw() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("A Violin wrapping a boxplot") +
  xlab("")
