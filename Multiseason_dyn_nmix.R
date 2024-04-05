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

setwd("./Data/")
# Load covariate data
site_cov <- data.frame(alt=pwr_det_hist$Alt,flow=(pwr_det_hist$Flow),wor=pwr_det_hist$WoR )

for(i in 1:length(data_files))
{
  spdat <- read.csv(data_files[i])
  # Build an unmarkedmultFramOccu
 unmarked_Frame <- unmarkedFramePCO(y = as.matrix(spdat[1:43,2:16]), ##two files common_kingfisher and lf has one row short
                                             numPrimary = 5,
                                             siteCovs = site_cov
  ) 
  
  ## ----build open nmixture---------------------------------------------------------
  dynamic_nmix_mod <- pcountOpen(lambdaformula= ~scale(alt)+flow+wor,
                                 gammaformula = ~1,
                                 omegaformula = ~1,
                                 pformula = ~1,
                                 data = unmarked_Frame,    method = "BFGS", K=50)
 
 sp_unmk_sum <- summary(dynamic_nmix_mod)
 sp_sum[[i]] <- cbind(paste(data_files[i]), do.call(rbind, sp_unmk_sum))
 
 #### Abundance estimate
 abd_est <- data_frame(year = 2014:2018, species = paste(data_files[i]))
 abd_est$mean<- colSums(bup(ranef(dynamic_nmix_mod)))
 abd_est[,4:5] <- t(colSums(confint(ranef(dynamic_nmix_mod))))
 sp_abd[[i]] <- abd_est
}

setwd("~/Documents/GitHub/Riverbird_Nmixture")

abd_cmbnd <-do.call(rbind, sp_abd)
abd_cmbnd$species<- gsub("_multi.csv", "",abd_cmbnd$species)
head(abd_cmbnd)
write.csv(abd_cmbnd, "./Results/All_species_abundance_combined.csv", row.names = F)

coef_cmbnd <-do.call(rbind, sp_sum)
colnames(coef_cmbnd)[1] <- c("species")
coef_cmbnd$species<- gsub("_multi.csv", "",coef_cmbnd$species)
coef_cmbnd$variable <- rep(rownames(coef_cmbnd)[1:7], 12)
coef_cmbnd$variable <- gsub("det", "det_prob",coef_cmbnd$variable)
coef_cmbnd$variable <- gsub("lambda.(Intercept)", "lambda_intrcpt",coef_cmbnd$variable, fixed = TRUE)
coef_cmbnd$variable <- gsub("lambda_elevation", "lambda_elev",coef_cmbnd$variable, fixed = TRUE)
coef_cmbnd$species[coef_cmbnd$species=="BD"] <- "Brown_Dipper"
coef_cmbnd$species[coef_cmbnd$species=="WTK"] <- "White_Thraoted_Kingfisher"
coef_cmbnd$species[coef_cmbnd$species=="BWT"] <- "Blue_Whistling_Thrush"
coef_cmbnd$species[coef_cmbnd$species=="CK"] <- "Crested_Kingfisher"
coef_cmbnd$species[coef_cmbnd$species=="GW"] <- "Grey_Wagtail"
coef_cmbnd$species[coef_cmbnd$species=="LF"] <- "Little_Forktail"
coef_cmbnd$species[coef_cmbnd$species=="PWR"] <- "Plumbeous_Water_Redstart"
coef_cmbnd$species[coef_cmbnd$species=="SF"] <- "Spotted_Forktail"
coef_cmbnd$species[coef_cmbnd$species=="WCR"] <- "White_capped_Redstart"
write.csv(coef_cmbnd, "./Results/All_species_coef_combined.csv", row.names = F)

##plot
  ggplot(abd_cmbnd, aes(x = year, y = mean1/43)) +
  geom_errorbar(aes(ymin = V1/43,
                    ymax = V2/43),
                width = 0) +
  geom_point(size = 2) +
  geom_line() +
    facet_wrap(~species, scales = "free")+
  labs(x = "Year", y = "Estimated Density per 500m") +
  theme_bw()

ggsave("./Results/All_species_abundance_new.jpeg", width=12, height=9, units = "in", dpi=300)  
  
  ggplot(coef_cmbnd, aes(x= variable, y = Estimate))+
    geom_errorbar(aes(ymin = Estimate -1.96*SE,
                      ymax = Estimate +1.96*SE),
                  width = 0) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed")+
    facet_wrap(~species, scales = "free")+
    labs(x = "Variable", y = "Coefficient Estimate") +
    coord_cartesian(ylim = c(-4,4))+
    theme_bw()+
    theme(axis.title = element_text(size=16), axis.text = element_text(size=10),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0))
  
  ggsave("./Results/All_species_coef_new.jpeg", width=12, height=9, units = "in", dpi=300)  
  