# Created by: AF
# Updated to include year
# All missings are zeros, and when a year is not present it is all zeros
# Updated 2019 12 06
# exploring demographics 
rm(list = ls())


# --- READ IN DATA ---
class_bisph <- read.csv(file = 'Z:\\Pharma_Influence\\Data\\Updated_data_10_17_2019\\class_Bisphosphonates.csv', header = T)
names(class_bisph)
summary(class_bisph)



# --- CLEAN DATA ---
# put year period into ref dataset
ref <- class_bisph[class_bisph$NPI == 0,]
ref <- ref[,2, drop = F]
# remove first empty rows
class_bisph <- class_bisph[class_bisph$NPI != 0,]
# pull all ids
id <- as.data.frame(unique(class_bisph$NPI))
colnames(id) <- 'NPI'
# merge with year to get complete reference
cref <- merge(ref, id)
# merge cref onto data to have instances for every year for every id
dat <- merge(x=class_bisph, y=cref, by = c('NPI', 'Year'), all = T)
# turn all missing to zeros
dat[is.na(dat)] <- 0




# --- HISTOGRAMS ---
# histogram excluding the zeros for all the variables
myhists <- function(var){
  hist(subset(dat[,var], dat[,var] !=0), breaks = 100, main = var)
}
# myhists(var = "Pay_atelvia")
# myhists(var = "Pay_boniva") #all zeros or missings
# myhists(var = "Pay_divigel")
# myhists(var = "Pay_elestrin") #interesting two peaks
# myhists(var = "Pay_fosamax") 
# myhists(var = "Pay_prolia")
# myhists(var = "Pay_risedronate")
# myhists(var = "Pre_atelvia")
# myhists(var = "Pre_boniva") 
# myhists(var = "Pre_divigel")
# myhists(var = "Pre_elestrin")
# myhists(var = "Pre_fosamax")
# myhists(var = "Pre_prolia") 
# myhists(var = "Pre_risedronate")





# --- DEMOS ---
dem <- read.csv('Z:\\Pharma_Influence\\Data\\Physician_demographics\\GOBA_unique.csv')
# dem2 <- read.csv('Z:\\Pharma_Influence\\Data\\Physician_demographics\\GOBA another version.csv')
# dem3 <- read.csv('Z:\\Pharma_Influence\\Data\\Physician_demographics\\GOBA version 3.csv')

# remove NA NPI or NA Year
dem <- dem[!is.na(dem$NPI_Match),]

demid <- as.data.frame(unique(dem$NPI_Match))
# there are 53533 in dem, but 52605 unique NPI, so 928 duplicates?

# just merge a subset of the demo data to see how they merge
demsub <- dem[, c('NPI_Match', 'GOBA_State')]
colnames(demsub) <- c('NPI', 'State')
#there are no missing states
statetest <- demsub[is.na(demsub$State),]
rm(statetest)



# --- PREP DATA FOR MODEL ---

# wide to long 
datl <- reshape(dat, direction = "long", idvar = c("NPI", "Year"), sep = '_', timevar = 'drug',
                varying = c("Pay_atelvia", "Pay_boniva", "Pay_divigel", "Pay_elestrin", "Pay_fosamax", "Pay_prolia", "Pay_risedronate",
                            "Pre_atelvia", "Pre_boniva", "Pre_divigel", "Pre_elestrin", "Pre_fosamax", "Pre_prolia", "Pre_risedronate"))
datl$NPI2 <- as.factor(datl$NPI)


# SORT and add CUMULATIVE PAYMENT column
datl <- datl[order(datl$NPI2, datl$drug, datl$Year),]
datl$Cpay <- ave(datl$Pay, datl$NPI2, datl$drug, FUN=cumsum)


# MERGE ON DEMOS
datdem <- merge(datl, demsub, all = T)# removing all = T will just merge where IDs are in both datasets
test <- datdem[is.na(datdem$State),]
length(unique(test$NPI))
# 1165 are not in the dem data (some in dem are not in prescribing data, but that doesnt matter)
length(unique(datl$NPI))
# there are 9743 total, so that is ~12% of the sample
datdem$NPI2 <- as.factor(datdem$NPI)


# look at percent that is 0 for Pre and Cpay
subcols <- datl[,c(5,7)]
colSums(subcols==0)/nrow(subcols)*100




# --- MODEL --- 

library(MCMCglmm)


# for reference, it is about 10 minutes per every 1000 iterations
# 10,000 takes 1.7 hours
# 20,000 takes 3.3 hours
# 100,000 would take 16.7 hours



# Mutivariate error structures (i.e. using the term 'trait') are required for multinomial data with more than 2 categories,
# or zero-infalted/altered/hurdle models.

testrun5000.100000.100 <- MCMCglmm(data = datl, Pre~ 1 + drug + Cpay + Year, random = ~NPI2, rcov = ~us(trait):units, family = "zipoisson",
                 burnin = 15000, nitt = 100000, thin = 100)
summary(testrun5000.100000.100)
plot(testrun5000.100000.100$Sol) # fixed effects trace plots
plot(testrun5000.100000.100$VCV) # random effects trace plots





# Year is not looking significant, but is it treating year as numeric,
# or do we need year to be an ordered/ranked factor? or could just be a factor
# I think the intercept is so low (doesn't make sense to be negative)
# because it's treating year as numeric, so for 2013 you need to add
# 2.762*2013 and add that to the -5.5 intercept (=5554.406)


# SAVE MODEL
# saveRDS(testrun4000.20000.100, "Z:\\Pharma_Influence\\Code\\AF\\prelim results\\testrun4000.20000.100.rds")
# saveRDS(test2000.10000.100, "Z:\\Pharma_Influence\\Code\\AF\\prelim results\\test2000.10000.100.rds")
# saveRDS(testrun4000.20000.10, "Z:\\Pharma_Influence\\Code\\AF\\prelim results\\testrun4000.20000.10.rds")
# saveRDS(testrun5000.100000.100, "Z:\\Pharma_Influence\\Code\\AF\\prelim results\\testrun4000.20000.100000.rds")


# LOAD IN PREVIOUS MODELS - make sure the MCMCglmm package is loaded when using this
testrun4000.20000.100 <- readRDS("Z:\\Pharma_Influence\\Code\\AF\\prelim results\\testrun4000.20000.100.rds")
test2000.10000.100 <- readRDS("Z:\\Pharma_Influence\\Code\\AF\\prelim results\\test2000.10000.100.rds")
testrun5000.100000.100 <- readRDS("Z:\\Pharma_Influence\\Code\\AF\\prelim results\\testrun4000.20000.100000.rds")
# this last one I believe was was 100,000 iterations, 5000 burn in and 100 thinning, labeling got messed up


# QUESTIONS
# could some of the convergence issues be due to overdispersion?
# could any of the drugs be correlated?

# are the demos even worth it? it sounded like this is all the demo data he has and there's very little



# NEXT STEPS------------------------------
# remove year
# make covariance structure simplier (compound symmetry (cs) or autoregressive (ar1)?)
# reduce iterations to 10,000
# 



model.idh.1000.10000.10 <- MCMCglmm(data = datl, Pre~ 1 + drug + Cpay, random = ~NPI2, rcov = ~idh(trait):units, family = "zipoisson",
                                   burnin = 1000, nitt = 10000, thin = 10)


summary(model.idh.1000.10000.10)
plot(model.idh.1000.10000.10$Sol) # fixed effects trace plots
plot(model.idh.1000.10000.10$VCV) # random effects trace plots


model.idh.1000.10000.10.trait <- MCMCglmm(data = datl, Pre~  trait - 1 + at.level(trait,1):drug + at.level(trait,1):Cpay, random = ~NPI2,rcov=~idh(trait):units,  family = "zipoisson",
                                   burnin = 1000, nitt = 10000, thin = 10)
summary(model.idh.1000.10000.10.trait)
plot(model.idh.1000.10000.10.trait$Sol) # fixed effects trace plots
plot(model.idh.1000.10000.10.trait$VCV) # random effects trace plots

model.idv.1000.10000.10 <- MCMCglmm(data = datl, Pre~ 1 + drug + Cpay, random = ~NPI2, rcov = ~idv(trait):units, family = "zipoisson",
                                    burnin = 1000, nitt = 10000, thin = 10)
summary(model.idv.1000.10000.10)
plot(model.idv.1000.10000.10$Sol) # fixed effects trace plots
plot(model.idv.1000.10000.10$VCV) # random effects trace plots

# idh - diagonal heterogeneous covariance


# SAVE MODELS
saveRDS(model.idh.1000.10000.10, "Z:\\Pharma_Influence\\Code\\AF\\prelim results\\model.idh.1000.10000.10.rds")
saveRDS(model.idh.1000.10000.10.trait, "Z:\\Pharma_Influence\\Code\\AF\\prelim results\\model.idh.1000.10000.10.trait.rds")

saveRDS(model.idv.1000.10000.10, "Z:\\Pharma_Influence\\Code\\AF\\prelim results\\model.idv.1000.10000.10.rds")


model.idv.1000.10000.10 <- readRDS("Z:\\Pharma_Influence\\Code\\AF\\prelim results\\model.idv.1000.10000.10.rds")



# ----------------------------------
# compare to trait model

# john ran (and then with 100,000 iterations and 100 burn in even though it says 10)
model.idh.1000.10000.10.trait<-MCMCglmm(data=datl,Pre~trait-1+at.level(trait,1):drug + at.level(trait,1):Cpay, random=~NPI2,rcov=~idh(trait):units,
                                        family="zipoisson",burnin=1000,nitt=10000,thin=10) 

model.idh.1000.10000.10 <- readRDS("Z:\\Pharma_Influence\\Code\\AF\\prelim results\\model.idh.1000.10000.10.rds")
model.idh.1000.10000.10.trait <- readRDS("Z:\\Pharma_Influence\\Code\\AF\\prelim results\\model.idh.1000.10000.10.trait.rds")
model.idh.1000.100000.10.trait <- readRDS("Z:\\Pharma_Influence\\Code\\AF\\prelim results\\model.idh.1000.100000.10.trait.rds")


summary(model.idh.1000.10000.10)
summary(model.idh.1000.10000.10.trait)
summary(model.idh.1000.100000.10.trait)
plot(model.idh.1000.10000.10$Sol) # fixed effects trace plots
plot(model.idh.1000.10000.10$VCV) # random effects trace plots
plot(model.idh.1000.10000.10.trait$Sol) # fixed effects trace plots
plot(model.idh.1000.10000.10.trait$VCV) # random effects trace plots
plot(model.idh.1000.100000.10.trait$Sol) # fixed effects trace plots
plot(model.idh.1000.100000.10.trait$VCV) # random effects trace plots




# ----------------------------------
# 01.17.2020 try with just poisson family


model.idh.1000.10000.10.poisson<-MCMCglmm(data=datl, Pre~ drug + Cpay, random = ~NPI2, 
                                        family="poisson",burnin=1000,nitt=10000,thin=10) 
summary(model.idh.1000.10000.10.poisson)
plot(model.idh.1000.10000.10.poisson$Sol) # fixed effects trace plots
plot(model.idh.1000.10000.10.poisson$VCV) # random effects trace plots
saveRDS(model.idh.1000.10000.10.poisson, "Z:\\Pharma_Influence\\Code\\AF\\prelim results\\model.idh.1000.10000.10.poisson.rds")


model.idh.5000.50000.10.poisson<-MCMCglmm(data=datl, Pre~ drug + Cpay, random = ~NPI2, 
                                          family="poisson",burnin=5000,nitt=50000,thin=10) 
summary(model.idh.5000.50000.10.poisson)
plot(model.idh.5000.50000.10.poisson$Sol) # fixed effects trace plots
plot(model.idh.5000.50000.10.poisson$VCV) # random effects trace plots
saveRDS(model.idh.5000.50000.10.poisson, "Z:\\Pharma_Influence\\Code\\AF\\prelim results\\model.idh.5000.50000.10.poisson.rds")


model.idh.1000.10000.10.poisson <- readRDS("Z:\\Pharma_Influence\\Code\\AF\\prelim results\\model.idh.1000.10000.10.poisson.rds")
model.idh.5000.50000.10.poisson <- readRDS("Z:\\Pharma_Influence\\Code\\AF\\prelim results\\model.idh.5000.50000.10.poisson.rds")



# Notes from call with John 
# --------------------------
# could use glimmix to drop down to just poison and not zero inflated
# lmer?
  
# step down in difficulty

# use same mcmcplmm and change family to poisson and remove trait stuff, run like normal
# start with 5-10 thousand





