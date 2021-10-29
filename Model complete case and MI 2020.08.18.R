# AF
# 08.18.2020
# Model run with complete case vs imputed
rm(list =ls())


library(mice)
library(rqdatatable)
library(MCMCglmm)
library(randomForest)

# NOTE
# Oral and anti-viral don't have any payments so they can't be modeled, won't be included in this code
# Hormone and vaginal are too large to impute using the method I was trying
# Not sure if the imputing method is the best or if I should be trying a method that works for longitudinal
# (right now using the sum of the prescriptions as the outcome for the imputation)




# READ IN DEMO DATA------------------------

dem <- read.csv("Z:/Pharma_Influence/Data/Updated data 7_7_2020/T1.csv")
dem_old <- read.csv('Z:\\Pharma_Influence\\Data\\Physician_demographics\\GOBA_unique.csv')

# CLEAN OLDER DEMO DATA---------------------

# remove NA NPI or NA Year
dem_old <- dem_old[!is.na(dem_old$NPI_Match),]
dem_oldid <- as.data.frame(unique(dem_old$NPI_Match))
# there are 53533 in dem, but 52605 unique NPI, so 928 duplicates?

# only need to merge a subset of the demos
dem_oldsub <- dem_old[, c('NPI_Match', 'GOBA_Region', 'HG_Gender', 'HG_Age', 'HG_Medical.School.1...Year')]
colnames(dem_oldsub) <- c('NPI', 'Region', 'Gender', 'Age', 'SchoolYear')

dem_oldsubna <- dem_oldsub[(!is.na(dem_oldsub$Region)) & (dem_oldsub$Gender!='') & (!is.na(dem_oldsub$Age)) & (!is.na(dem_oldsub$SchoolYear)),]
# only 27058 when remove all missing


# CLEAN NEW DEMO DATA----------------------------

# remove NA NPI or NA Year
dem <- dem[!is.na(dem$NPI),]  # no missing npi number
demid <- as.data.frame(unique(dem$NPI))
# in this demo there are 12055 but 12047 unique NPI, so 8 duplicates? SHOULD LOOK INTO WHICH ONES TO BE REMOVING OF THE DUPS

# only need to merge a subset of the demos
demsub <- dem[, c('NPI', 'GOBA_Region', 'PCND_gender', 'Age', 'PCND_GradYr')]
colnames(demsub) <- c('NPI', 'Region', 'Gender', 'Age', 'SchoolYear')



# MERGE TWO DEMO DATASETS------------------------
# KEEPING WHATS NON MISSING PRIORITIZING THE NEW DEMOS

demfinal <- natural_join(a=demsub, b=dem_oldsub, 
                         by = "NPI",
                         jointype = "FULL")
# clean gender to be coded the same way 
demfinal$Gender <- ifelse(demfinal$Gender=='3', 'M', demfinal$Gender)
demfinal$Gender <- ifelse(demfinal$Gender=='2', 'F', demfinal$Gender)
demfinal$Gender <- ifelse(demfinal$Gender=='1', NA, demfinal$Gender)
demfinal <- demfinal[demfinal$NPI!=1306,] #not a real NPI number, removed



# LOAD PAYMENT/PRESCRIBING DATA-----------------------------------
load_anti_inf  <- read.csv(file = 'Z:\\Pharma_Influence\\Data\\Updated_data_10_17_2019\\class_Anti_infective.csv', header = T)
load_antichol  <- read.csv(file = 'Z:\\Pharma_Influence\\Data\\Updated_data_10_17_2019\\class_Anticholinergics_for_overactive_bladder.csv', header = T)
load_bisp      <- read.csv(file = 'Z:\\Pharma_Influence\\Data\\Updated_data_10_17_2019\\class_Bisphosphonates.csv', header = T)
load_hormone   <- read.csv(file = 'Z:\\Pharma_Influence\\Data\\Updated_data_10_17_2019\\class_Hormone_therapy_single_ingredient_therapy.csv', header = T)
load_transderm <- read.csv(file = 'Z:\\Pharma_Influence\\Data\\Updated_data_10_17_2019\\class_Transdermal_estrogen.csv', header = T)
load_vaginal   <- read.csv(file = 'Z:\\Pharma_Influence\\Data\\Updated_data_10_17_2019\\class_Vaginal_Estrogen_Hormone_Therapy.csv', header = T)


# CLEAN AND MERGE TO DEMOS---------------------------------
cleandat <- function(loaddata){
  
  # --- CLEAN DATA ---
  # put year period into ref dataset
  ref <- loaddata[loaddata$NPI == 0,]
  ref <- ref[,2, drop = F]
  # remove first empty rows
  loaddata <- loaddata[loaddata$NPI != 0,]
  # pull all ids
  id <- as.data.frame(unique(loaddata$NPI))
  colnames(id) <- 'NPI'
  # merge with year to get complete reference
  cref <- merge(ref, id)
  # merge cref onto data to have instances for every year for every id
  dat <- merge(x=loaddata, y=cref, by = c('NPI', 'Year'), all = T)
  # turn all missing to zeros
  dat[is.na(dat)] <- 0
  
  # wide to long 
  columns <- grep(c("^P"), names(loaddata))
  datl <- reshape(dat, direction = "long", idvar = c("NPI", "Year"), sep = '_', timevar = 'drug',
                  varying = columns)
  datl$NPI2 <- as.factor(datl$NPI)
  
  # sort and add cumulative payment column
  datl <- datl[order(datl$NPI2, datl$drug, datl$Year),]
  datl$Cpay <- ave(datl$Pay, datl$NPI2, datl$drug, FUN=cumsum)
  
  # MERGE ON DEMOS
  datdem <- merge(x=datl, y=demfinal, by = 'NPI', all.x = T)
  datdem$Region2 <- as.factor(datdem$Region)
  datdem$Gender2 <- as.factor(datdem$Gender)
  
  return(datdem)
}


dat_anti_inf_long  <- cleandat(loaddata=load_anti_inf)
dat_antichol_long  <- cleandat(loaddata=load_antichol)
dat_bisp_long      <- cleandat(loaddata=load_bisp)
dat_hormone_long   <- cleandat(loaddata=load_hormone)
dat_transderm_long <- cleandat(loaddata=load_transderm)
dat_vaginal_long   <- cleandat(loaddata=load_vaginal)



# TAKE THE SUM OF PRESCRIPTIONS OVER ALL THE DIFFERENT DRUGS/YEARS
cleantotals <- function(longdata){
  
  one <- aggregate(cbind(TotPay=longdata$Pay, TotPre=longdata$Pre), by=list(NPI2=longdata$NPI2), FUN=sum)
  two <- merge(x=one, y=subset(longdata, select = -c(Pay,Pre,Cpay,Year,drug)), by = 'NPI2', all.x = T, all.y=F)
  three <- two[!duplicated(two$NPI2),]
  return(three)
  
}
dat_anti_inf <- cleantotals(dat_anti_inf_long)    #5354 physicians
dat_antichol <- cleantotals(dat_antichol_long)    #12002 physicians
dat_bisp <- cleantotals(dat_bisp_long)            #9743 physicians
dat_hormone <- cleantotals(dat_hormone_long)      #14941 physicians
dat_transderm <- cleantotals(dat_transderm_long)  #3953 physicians
dat_vaginal <- cleantotals(dat_vaginal_long)      #21098 physicians




# IMPUTING--------------------------------------------------------------------------------------------------


# Prescriptions is the outcome but because that is longitudinal over several years and several drugs
# I will just be using the total prescriptions but im not sure if that's the best method and I 
# should try an imputing method that is made for longitudinal data


pred <- dat_anti_inf[, c("NPI", "Gender2", "Region2", "Age")]
y <- dat_anti_inf[,c("TotPre")]
dat_anti_inf_impute <- rfImpute(pred, y)

pred <- dat_antichol[, c("NPI", "Gender2", "Region2", "Age")]
y <- dat_antichol[,c("TotPre")]
dat_antichol_impute <- rfImpute(pred, y)

pred <- dat_bisp[, c("NPI", "Gender2", "Region2", "Age")]
y <- dat_bisp[,c("TotPre")]
dat_bisp_impute <- rfImpute(pred, y)

pred <- dat_hormone[, c("NPI", "Gender2", "Region2", "Age")]
y <- dat_hormone[,c("TotPre")]
dat_hormone_impute <- rfImpute(pred, y)

pred <- dat_transderm[, c("NPI", "Gender2", "Region2", "Age")]
y <- dat_transderm[,c("TotPre")]
dat_transderm_impute <- rfImpute(pred, y)

pred <- dat_vaginal[, c("NPI", "Gender2", "Region2", "Age")]
y <- dat_vaginal[,c("TotPre")]
dat_vaginal_impute <- rfImpute(pred, y)

# NOTE
# hormone and vaginal did not work "error: cannot allocate vector of size..."
# running a second time, also antichol had error




# MERGE IMPUTED DEMOS ONTO FULL DATA---------------------------------
cleandat <- function(imputedata, longdata){
  
  # remove columns not needed or that will be joined from imputed demos
  sub <- subset(longdata, select=-c(Age, Gender, Region, Region2, Gender2, SchoolYear))
  # MERGE ON DEMOS
  datdem <- merge(x=sub, y=imputedata, by = 'NPI', all.x = T)
  datdem$NPI2 <- as.factor(datdem$NPI)
  
  return(datdem)
}

dat_anti_inf_final  <- cleandat(imputedata=dat_anti_inf_impute,  longdata= dat_anti_inf_long)
dat_antichol_final  <- cleandat(imputedata=dat_antichol_impute,  longdata= dat_antichol_long)
dat_bisp_final      <- cleandat(imputedata=dat_bisp_impute,      longdata= dat_bisp_long)
dat_hormone_final   <- cleandat(imputedata=dat_hormone_impute,   longdata= dat_hormone_long)
dat_transderm_final <- cleandat(imputedata=dat_transderm_impute, longdata= dat_transderm_long)
dat_vaginal_final   <- cleandat(imputedata=dat_vaginal_impute,   longdata= dat_vaginal_long)




# CREATE COMPLETE CASE DATASETS for age, region and gender-----------
dat_anti_inf_finalcc <- dat_anti_inf_long[(!is.na(dat_anti_inf_long$Region)) & (dat_anti_inf_long$Gender!='') & (!is.na(dat_anti_inf_long$Age)),]
dat_antichol_finalcc <- dat_antichol_long[(!is.na(dat_antichol_long$Region)) & (dat_antichol_long$Gender!='') & (!is.na(dat_antichol_long$Age)),]
dat_bisp_finalcc     <- dat_bisp_long[(!is.na(dat_bisp_long$Region)) & (dat_bisp_long$Gender!='') & (!is.na(dat_bisp_long$Age)),]
dat_hormone_finalcc  <- dat_hormone_long[(!is.na(dat_hormone_long$Region)) & (dat_hormone_long$Gender!='') & (!is.na(dat_hormone_long$Age)),]
dat_transderm_finalcc<- dat_transderm_long[(!is.na(dat_transderm_long$Region)) & (dat_transderm_long$Gender!='') & (!is.na(dat_transderm_long$Age)),]
dat_vaginal_finalcc  <- dat_vaginal_long[(!is.na(dat_vaginal_long$Region)) & (dat_vaginal_long$Gender!='') & (!is.na(dat_vaginal_long$Age)),]






# RUN MODELS WITH IMPUTED DEMOS AND WITH COMPLETE CASE-------------------------------------------
# NOTE SCHOOL YEAR IS REMOVED AND GENDER AND REGION ARE CATEGORICAL
modeldat <- function(loaddata){
  
  model <<-MCMCglmm(data=loaddata, Pre~ drug + Cpay + Gender2 + Age + Region2, random = ~NPI2, 
                    family="poisson",burnin=5000,nitt=50000,thin=10) 
  return(model)
}
# RUN FOR IMPUTED DATA
model_anti_inf_impute  <- modeldat(loaddata=dat_anti_inf_final)
model_antichol_impute  <- modeldat(loaddata=dat_antichol_final)
model_bisp_impute      <- modeldat(loaddata=dat_bisp_final)
model_hormone_impute   <- modeldat(loaddata=dat_hormone_final)
model_transderm_impute <- modeldat(loaddata=dat_transderm_final)
model_vaginal_impute   <- modeldat(loaddata=dat_vaginal_final)

# RUN FOR COMPLETE CASE DATA
model_anti_inf_cc    <- modeldat(loaddata=dat_anti_inf_finalcc)
model_antichol_cc    <- modeldat(loaddata=dat_antichol_finalcc)
model_bisp_cc        <- modeldat(loaddata=dat_bisp_finalcc)
model_hormone_cc     <- modeldat(loaddata=dat_hormone_finalcc)
model_transderm_cc   <- modeldat(loaddata=dat_transderm_finalcc)
model_vaginal_cc     <- modeldat(loaddata=dat_vaginal_finalcc)



# SAVE MODEL RESULTS-----------------------------
saveRDS(model_anti_inf_impute, "Z:\\Pharma_Influence\\Results\\Model impute results 08.18.2020\\anti_inf_impute.5000.50000.10.poisson.rds")
saveRDS(model_antichol_impute, "Z:\\Pharma_Influence\\Results\\Model impute results 08.18.2020\\antichol_impute.5000.50000.10.poisson.rds")
saveRDS(model_bisp_impute, "Z:\\Pharma_Influence\\Results\\Model impute results 08.18.2020\\bisp_impute.5000.50000.10.poisson.rds")
saveRDS(model_hormone_impute, "Z:\\Pharma_Influence\\Results\\Model impute results 08.18.2020\\hormone_impute.5000.50000.10.poisson.rds")
saveRDS(model_transderm_impute, "Z:\\Pharma_Influence\\Results\\Model impute results 08.18.2020\\transderm_impute.5000.50000.10.poisson.rds")
saveRDS(model_vaginal_impute, "Z:\\Pharma_Influence\\Results\\Model impute results 08.18.2020\\vaginal_impute.5000.50000.10.poisson.rds")

saveRDS(model_anti_inf_cc, "Z:\\Pharma_Influence\\Results\\Model impute results 08.18.2020\\anti_inf_cc.5000.50000.10.poisson.rds")
saveRDS(model_antichol_cc, "Z:\\Pharma_Influence\\Results\\Model impute results 08.18.2020\\antichol_cc.5000.50000.10.poisson.rds")
saveRDS(model_antiviral_cc, "Z:\\Pharma_Influence\\Results\\Model impute results 08.18.2020\\antiviral_cc.5000.50000.10.poisson.rds")
saveRDS(model_bisp_cc, "Z:\\Pharma_Influence\\Results\\Model impute results 08.18.2020\\bisp_cc.5000.50000.10.poisson.rds")
saveRDS(model_hormone_cc, "Z:\\Pharma_Influence\\Results\\Model impute results 08.18.2020\\hormone_cc.5000.50000.10.poisson.rds")
saveRDS(model_transderm_cc, "Z:\\Pharma_Influence\\Results\\Model impute results 08.18.2020\\transderm_cc.5000.50000.10.poisson.rds")
saveRDS(model_vaginal_cc, "Z:\\Pharma_Influence\\Results\\Model impute results 08.18.2020\\vaginal_cc.5000.50000.10.poisson.rds")


# WRITE OUT COEFFICIENTS FROM EACH MODEL--------------------------------------
save_anti_inf_impute <- data.frame(summary(model_anti_inf)$solutions) 
write.csv(save_anti_inf_impute, "Z:\\Pharma_Influence\\Results\\Model impute results 08.18.2020\\anti_inf_impute_coeff.csv")
save_anti_inf_cc <- data.frame(summary(model_anti_infna)$solutions)
write.csv(save_anti_inf_cc, "Z:\\Pharma_Influence\\Results\\Model impute results 08.18.2020\\anti_inf_cc_coeff.csv")

save_anti_inf_impute <- data.frame(summary(model_anti_inf)$solutions) 
write.csv(save_anti_inf_impute, "Z:\\Pharma_Influence\\Results\\Model impute results 08.18.2020\\anti_inf_impute_coeff.csv")
save_anti_inf_cc <- data.frame(summary(model_anti_infna)$solutions)
write.csv(save_anti_inf_cc, "Z:\\Pharma_Influence\\Results\\Model impute results 08.18.2020\\anti_inf_cc_coeff.csv")

















