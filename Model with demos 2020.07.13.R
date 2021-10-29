# AF
# 07.13.2020
# New Demos
rm(list =ls())



library(MCMCglmm)
library(rqdatatable)



# LOOK AT NEW DEMOS




StudyGroup <- readRDS("Z:/Pharma_Influence/Data/Updated data 7_7_2020/StudyGroup.rds")
PaySum <- readRDS("Z:/Pharma_Influence/Data/Updated data 7_7_2020/PaySum.rds")
Prescriber <- readRDS("Z:/Pharma_Influence/Data/Updated data 7_7_2020/Prescriber.rds")


dem <- read.csv("Z:/Pharma_Influence/Data/Updated data 7_7_2020/T1.csv")
dem_old <- read.csv('Z:\\Pharma_Influence\\Data\\Physician_demographics\\GOBA_unique.csv')


# CLEAN OLD DATA----------------------------

# remove NA NPI or NA Year
dem_old <- dem_old[!is.na(dem_old$NPI_Match),]
dem_oldid <- as.data.frame(unique(dem_old$NPI_Match))
# there are 53533 in dem, but 52605 unique NPI, so 928 duplicates?

# only need to merge a subset of the demos
dem_oldsub <- dem_old[, c('NPI_Match', 'GOBA_Region', 'HG_Gender', 'HG_Age', 'HG_Medical.School.1...Year')]
colnames(dem_oldsub) <- c('NPI', 'Region', 'Gender', 'Age', 'SchoolYear')

dem_oldsubna <- dem_oldsub[(!is.na(dem_oldsub$Region)) & (dem_oldsub$Gender!='') & (!is.na(dem_oldsub$Age)) & (!is.na(dem_oldsub$SchoolYear)),]
# only 27058 when remove all missing


# CLEAN NEW DATA----------------------------

# remove NA NPI or NA Year
dem <- dem[!is.na(dem$NPI),]  # no missing npi number
demid <- as.data.frame(unique(dem$NPI))
# in this demo there are 12055 but 12047 unique NPI, so 8 duplicates? SHOULD LOOK INTO WHICH ONES TO BE REMOVING OF THE DUPS

# only need to merge a subset of the demos
demsub <- dem[, c('NPI', 'GOBA_Region', 'PCND_gender', 'Age', 'PCND_GradYr')]
colnames(demsub) <- c('NPI', 'Region', 'Gender', 'Age', 'SchoolYear')

# MISSINGS: 2575 missing regions, 2944 missing gender, 2968 missing age, 2968 missing year
test <- demsub[is.na(demsub$Region),]
test <- demsub[is.na(demsub$Gender),]
test <- demsub[is.na(demsub$Age),]
test <- demsub[is.na(demsub$SchoolYear),]
rm(test)

demsubna <- demsub[(!is.na(demsub$Region)) & (demsub$Gender!='') & (!is.na(demsub$Age)) & (!is.na(demsub$SchoolYear)),]
# only 6523 when remove all missing






# look how much missing when we combine ids but dont remove incompletes


# merge the two instead
combinem <- merge(x=demsub, y=dem_oldsub, by = 'NPI', all = T)
# that leaves 55192

# clean gender to match other dataframe


# coalesce type merge
demfinal <- natural_join(a=demsub, b=dem_oldsub, 
                      by = "NPI",
                      jointype = "FULL")
demfinal$Gender <- ifelse(demfinal$Gender=='3', 'M', demfinal$Gender)
demfinal$Gender <- ifelse(demfinal$Gender=='2', 'F', demfinal$Gender)
demfinal$Gender <- ifelse(demfinal$Gender=='1', NA, demfinal$Gender)


# MISSINGS: 1216 missing regions, 15908 missing gender, 20378 missing age, 20069 missing year
test <- demfinal[is.na(demfinal$Region),]
test <- demfinal[is.na(demfinal$Gender),]
test <- demfinal[is.na(demfinal$Age),]
test <- demfinal[is.na(demfinal$SchoolYear),]







# LOAD PAYMENT/PRESCRIBING DATA-----------------------------------
load_anti_inf  <- read.csv(file = 'Z:\\Pharma_Influence\\Data\\Updated_data_10_17_2019\\class_Anti_infective.csv', header = T)
load_antichol  <- read.csv(file = 'Z:\\Pharma_Influence\\Data\\Updated_data_10_17_2019\\class_Anticholinergics_for_overactive_bladder.csv', header = T)
load_antiviral <- read.csv(file = 'Z:\\Pharma_Influence\\Data\\Updated_data_10_17_2019\\class_Antiviral.csv', header = T)
load_bisph     <- read.csv(file = 'Z:\\Pharma_Influence\\Data\\Updated_data_10_17_2019\\class_Bisphosphonates.csv', header = T)
load_hormone   <- read.csv(file = 'Z:\\Pharma_Influence\\Data\\Updated_data_10_17_2019\\class_Hormone_therapy_single_ingredient_therapy.csv', header = T)
load_oral      <- read.csv(file = 'Z:\\Pharma_Influence\\Data\\Updated_data_10_17_2019\\class_Oral_Combined_Estrogen_and_Progestin_Products_for_Hormone_Therapy.csv', header = T)
load_transderm <- read.csv(file = 'Z:\\Pharma_Influence\\Data\\Updated_data_10_17_2019\\class_Transdermal_estrogen.csv', header = T)
load_vaginal   <- read.csv(file = 'Z:\\Pharma_Influence\\Data\\Updated_data_10_17_2019\\class_Vaginal_Estrogen_Hormone_Therapy.csv', header = T)



# CLEAN DATA---------------------------------
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
  
  return(datdem)
}


dat_anti_inf  <- cleandat(loaddata=load_anti_inf)
dat_antichol  <- cleandat(loaddata=load_antichol)
dat_antiviral <- cleandat(loaddata=load_antiviral)
dat_bisp      <- cleandat(loaddata=load_bisph)
dat_hormone   <- cleandat(loaddata=load_hormone)
dat_oral      <- cleandat(loaddata=load_oral)
dat_transderm <- cleandat(loaddata=load_transderm)
dat_vaginal   <- cleandat(loaddata=load_vaginal)


# remove duplicates to just see whats left of the IDs and their demos

dat_anti_inf2 <- dat_anti_inf[!duplicated(dat_anti_inf$NPI),] #5354 physicians
nrow(dat_anti_inf2[is.na(dat_anti_inf2$Region),]) #641 (11%)
nrow(dat_anti_inf2[is.na(dat_anti_inf2$Gender),]) #1468 (27%)
nrow(dat_anti_inf2[is.na(dat_anti_inf2$Age),])    #1678 (31%)
nrow(dat_anti_inf2[is.na(dat_anti_inf2$SchoolYear),]) #1808 (34%)


dat_antichol2 <- dat_antichol[!duplicated(dat_antichol$NPI),] #12002 total physicians
nrow(dat_antichol2[is.na(dat_antichol2$Region),]) #1442 (12%)
nrow(dat_antichol2[is.na(dat_antichol2$Gender),]) #3087 (26%)
nrow(dat_antichol2[is.na(dat_antichol2$Age),])    #3572 (30%)
nrow(dat_antichol2[is.na(dat_antichol2$SchoolYear),]) #3751 (31%)


dat_antiviral2 <- dat_antiviral[!duplicated(dat_antiviral$NPI),] #72 total physicians
nrow(dat_antiviral2[is.na(dat_antiviral2$Region),]) #10 (%)
nrow(dat_antiviral2[is.na(dat_antiviral2$Gender),]) #17 (%)
nrow(dat_antiviral2[is.na(dat_antiviral2$Age),])    #20 (%)
nrow(dat_antiviral2[is.na(dat_antiviral2$SchoolYear),]) #22 (%)


dat_bisp2 <- dat_bisp[!duplicated(dat_bisp$NPI),] #9743 total physicians
nrow(dat_bisp2[is.na(dat_bisp2$Region),]) #1082 (11%)
nrow(dat_bisp2[is.na(dat_bisp2$Gender),]) #2629 (27%)
nrow(dat_bisp2[is.na(dat_bisp2$Age),])    #3027 (31%)
nrow(dat_bisp2[is.na(dat_bisp2$SchoolYear),]) #3184 (33%)


dat_hormone2 <- dat_hormone[!duplicated(dat_hormone$NPI),] #14941 total physicians
nrow(dat_hormone2[is.na(dat_hormone2$Region),]) #1416 (1%)
nrow(dat_hormone2[is.na(dat_hormone2$Gender),]) #3665 (25%)
nrow(dat_hormone2[is.na(dat_hormone2$Age),])    #4277 (29%)
nrow(dat_hormone2[is.na(dat_hormone2$SchoolYear),]) #4525 (30%)


dat_oral2 <- dat_oral[!duplicated(dat_oral$NPI),] #5937 total physicians
nrow(dat_oral2[is.na(dat_oral2$Region),]) #535 (1%)
nrow(dat_oral2[is.na(dat_oral2$Gender),]) #1447 (24%)
nrow(dat_oral2[is.na(dat_oral2$Age),])    #1572 (26%)
nrow(dat_oral2[is.na(dat_oral2$SchoolYear),]) #1770 (30%)



# region is around 10%
# gender is around 25%
# age is around 30%
# school year is around 30%



# look at pattern of missing and if its associated with the variables
# if zero vs non zero instead of correlation

# mice package
# elise


# COMPLETE CASES WITHOUT NAS FOR MODEL
dat_anti_inf_nona <- dat_anti_inf[(!is.na(dat_anti_inf$Region)) & (dat_anti_inf$Gender!='') & (!is.na(dat_anti_inf$Age)) & (!is.na(dat_anti_inf$SchoolYear)),]
dat_antichol_nona <- dat_antichol[(!is.na(dat_antichol$Region)) & (dat_antichol$Gender!='') & (!is.na(dat_antichol$Age)) & (!is.na(dat_antichol$SchoolYear)),]
dat_antiviral_nona <- dat_antiviral[(!is.na(dat_antiviral$Region)) & (dat_antiviral$Gender!='') & (!is.na(dat_antiviral$Age)) & (!is.na(dat_antiviral$SchoolYear)),]
dat_bisp_nona <- dat_bisp[(!is.na(dat_bisp$Region)) & (dat_bisp$Gender!='') & (!is.na(dat_bisp$Age)) & (!is.na(dat_bisp$SchoolYear)),]
dat_hormone_nona <- dat_hormone[(!is.na(dat_hormone$Region)) & (dat_hormone$Gender!='') & (!is.na(dat_hormone$Age)) & (!is.na(dat_hormone$SchoolYear)),]
dat_oral_nona <- dat_oral[(!is.na(dat_oral$Region)) & (dat_oral$Gender!='') & (!is.na(dat_oral$Age)) & (!is.na(dat_oral$SchoolYear)),]
dat_transderm_nona <- dat_transderm[(!is.na(dat_transderm$Region)) & (dat_transderm$Gender!='') & (!is.na(dat_transderm$Age)) & (!is.na(dat_transderm$SchoolYear)),]
dat_vaginal_nona <- dat_vaginal[(!is.na(dat_vaginal$Region)) & (dat_vaginal$Gender!='') & (!is.na(dat_vaginal$Age)) & (!is.na(dat_vaginal$SchoolYear)),]

nrow(dat_anti_inf_nona)/nrow(dat_anti_inf)    #60% of data left after complete case
nrow(dat_antichol_nona)/nrow(dat_antichol)    #63% left
nrow(dat_antiviral_nona)/nrow(dat_antiviral)  #63% left
nrow(dat_bisp_nona)/nrow(dat_bisp)            #62% 
nrow(dat_hormone_nona)/nrow(dat_hormone)      #65% 
nrow(dat_oral_nona)/nrow(dat_oral)            #67%
nrow(dat_transderm_nona)/nrow(dat_transderm)  #70%
nrow(dat_vaginal_nona)/nrow(dat_vaginal)      #63%


# RUN MODELS WITH ALL DEMOS--------------------------------------------------------------------
modeldat <- function(loaddata){
  
  model <<-MCMCglmm(data=loaddata, Pre~ drug + Cpay + Gender + Age + Region + SchoolYear, random = ~NPI2, 
                    family="poisson",burnin=5000,nitt=50000,thin=10) 
  
  return(model)
}


model_anti_inf  <- modeldat(loaddata=dat_anti_inf_nona)
model_antichol  <- modeldat(loaddata=dat_antichol_nona)
model_antiviral <- modeldat(loaddata=dat_antiviral_nona)
model_bisp      <- modeldat(loaddata=dat_bisp_nona)     
model_hormone   <- modeldat(loaddata=dat_hormone_nona)
model_oral      <- modeldat(loaddata=dat_oral_nona)
model_transderm <- modeldat(loaddata=dat_transderm_nona)
model_vaginal   <- modeldat(loaddata=dat_vaginal_nona)
# these all took 24 hours to run



# RESULTS SUMMARY-----------------------------
# all of the models had a warning about some effects not being estimable
summary(model_anti_inf)  # gender sig
summary(model_antichol)  # nothing sig
summary(model_antiviral) # gender close but nothing sig
summary(model_bisp)      # gender sig (last time no results, needed rerun)
summary(model_hormone)   # region sig
summary(model_oral)      # nothing sig (last time region sig, gender close)
summary(model_transderm) # nothing sig (las time gender and age sig)
summary(model_vaginal)   # nothing sig





# SAVE MODEL RESULTS-----------------------------
saveRDS(model_anti_inf, "Z:\\Pharma_Influence\\Results\\Model demo results 07.27.2020\\anti_inf_demo.5000.50000.10.poisson.rds")
saveRDS(model_antichol, "Z:\\Pharma_Influence\\Results\\Model demo results 07.27.2020\\antichol_demo.5000.50000.10.poisson.rds")
saveRDS(model_antiviral, "Z:\\Pharma_Influence\\Results\\Model demo results 07.27.2020\\antiviral_demo.5000.50000.10.poisson.rds")
saveRDS(model_bisp, "Z:\\Pharma_Influence\\Results\\Model demo results 07.27.2020\\bisp_demo.5000.50000.10.poisson.rds")
saveRDS(model_hormone, "Z:\\Pharma_Influence\\Results\\Model demo results 07.27.2020\\hormone_demo.5000.50000.10.poisson.rds")
saveRDS(model_oral, "Z:\\Pharma_Influence\\Results\\Model demo results 07.27.2020\\oral_demo.5000.50000.10.poisson.rds")
saveRDS(model_transderm, "Z:\\Pharma_Influence\\Results\\Model demo results 07.27.2020\\transderm_demo.5000.50000.10.poisson.rds")
saveRDS(model_vaginal, "Z:\\Pharma_Influence\\Results\\Model demo results 07.27.2020\\vaginal_demo.5000.50000.10.poisson.rds")



# WRITE OUT COEFFICIENTS FROM EACH MODEL--------------------------------------
save_bisp <- data.frame(summary(model_bisp)$solutions)  #save coefficents
write.csv(save_bisp, "Z:\\Pharma_Influence\\Results\\Model demo results 07.27.2020\\bisph_coeff.csv")

save_anti_inf <- data.frame(summary(model_anti_inf)$solutions)
write.csv(save_anti_inf, "Z:\\Pharma_Influence\\Results\\Model demo results 07.27.2020\\anti_inf_coeff.csv")

save_antichol <- data.frame(summary(model_antichol)$solutions)
write.csv(save_antichol, "Z:\\Pharma_Influence\\Results\\Model demo results 07.27.2020\\antichol_coeff.csv")

save_antiviral <- data.frame(summary(model_antiviral)$solutions)
write.csv(save_antiviral, "Z:\\Pharma_Influence\\Results\\Model demo results 07.27.2020\\antiviral_coeff.csv")

save_hormone <- data.frame(summary(model_hormone)$solutions)
write.csv(save_hormone, "Z:\\Pharma_Influence\\Results\\Model demo results 07.27.2020\\hormone_coeff.csv")

save_oral <- data.frame(summary(model_oral)$solutions)
write.csv(save_oral, "Z:\\Pharma_Influence\\Results\\Model demo results 07.27.2020\\oral_coeff.csv")

save_transderm <- data.frame(summary(model_transderm)$solutions)
write.csv(save_transderm, "Z:\\Pharma_Influence\\Results\\Model demo results 07.27.2020\\transderm_coeff.csv")

save_vaginal <- data.frame(summary(model_vaginal)$solutions)
write.csv(save_vaginal, "Z:\\Pharma_Influence\\Results\\Model demo results 07.27.2020\\vaginal_coeff.csv")


# READ IN MODEL RESULTS----------------------------

model_anti_inf <- readRDS("Z:/Pharma_Influence/Results/Model demo results 07.27.2020/anti_inf_demo.5000.50000.10.poisson.rds")
model_antichol <- readRDS("Z:/Pharma_Influence/Results/Model demo results 07.27.2020/antichol_demo.5000.50000.10.poisson.rds")
model_antiviral <- readRDS("Z:/Pharma_Influence/Results/Model demo results 07.27.2020/antiviral_demo.5000.50000.10.poisson.rds")
model_bisp <- readRDS("Z:/Pharma_Influence/Results/Model demo results 07.27.2020/bisp_demo.5000.50000.10.poisson.rds")
model_hormone <- readRDS("Z:/Pharma_Influence/Results/Model demo results 07.27.2020/hormone_demo.5000.50000.10.poisson.rds")
model_oral <- readRDS("Z:/Pharma_Influence/Results/Model demo results 07.27.2020/oral_demo.5000.50000.10.poisson.rds")
model_transderm <- readRDS("Z:/Pharma_Influence/Results/Model demo results 07.27.2020/transderm_demo.5000.50000.10.poisson.rds")
model_vaginal <- readRDS("Z:/Pharma_Influence/Results/Model demo results 07.27.2020/vaginal_demo.5000.50000.10.poisson.rds")


# PLOTS-----------------------------------------
# cant do antiviral and transderm because no cumulative payments


cleandrug <- function(data, payN, dataframe){
  
  save <- data.frame(summary(data)$solutions)                      # save coefficents
  num <- nrow(save)-4                                              # only want to keep the drugs and cpay, not the four demos at the end
  numbers <<- setNames(data.frame(t((save[2:num,1]))), row.names(save)[2:num])  #transpose to have column for each drug
  meanage <<- mean(dataframe$Age, na.rm=T)                          # get average age to use for average estimate
  meanyear <<- mean(dataframe$SchoolYear, na.rm=T)                  # get average school year to use for average estimate
  numbersrep <<- numbers[rep(seq_len(nrow(numbers)), each = payN),] + # repeat to get number of instances desired
                                    save[1,1]  +                             # add intercept
                                    save[nrow(save)-2,1]*meanage +           # add mean age
                                    save[nrow(save),1]*meanyear              # add mean school year
  numberspay <<- save[nrow(save)-4,1]*c(1:payN) #Cpay times payments
  numbersadd <<- sweep(numbersrep, 1, numberspay, "+") #add to every drug column
  numbersexp <<- exp(numbersadd)   # exponentiate
  numbersexp$pay <<- c(1:payN)
  return(numbersexp)
}




# bisph: Bisphosphonates---------------------------------
summary(model_bisp)
plot_bisp <- cleandrug(model_bisp, 90000, dat_bisp_nona)
plot(plot_bisp$pay, plot_bisp[,6], type = 'l', ylim = c(0, max(plot_bisp[,6])), 
     ylab = paste(colnames(plot_bisp)[3], "Prescriptions"), xlab = "Cumulative Payments",
     main = paste("Bisphosphonates class:\n",
                  colnames(plot_bisp)[3], "Prescriptions vs Cumulative Payments"))

lines(plot_bisp$pay, plot_bisp[,2], col = "red")
lines(plot_bisp$pay, plot_bisp[,3], col = "blue")

plot(plot_bisp$pay, plot_bisp[,2], type = 'l', ylim = c(0, max(plot_bisp[,2])))
lines(plot_bisp$pay, plot_bisp[,5], col = "green")
lines(plot_bisp$pay, plot_bisp[,3], col = "red")

# so small, and they the different drugs are on way different scales,
# but too small to be meaningful



# anti_inf: Anti infective--------------------------------
summary(model_anti_inf)
plot_anti_inf <- cleandrug(model_anti_inf, 100, dat_anti_inf_nona)
plot(plot_anti_inf$pay, plot_anti_inf[,1], type = 'l', ylim = c(0, max(plot_anti_inf[,1])), 
     ylab = paste(colnames(plot_anti_inf)[3], "Prescriptions"), xlab = "Cumulative Payments",
     main = paste("Anti-infective class:\n",
                  colnames(plot_anti_inf)[3], "Prescriptions vs Cumulative Payments"))
lines(plot_anti_inf$pay, plot_anti_inf[,2], col = "green")

plot(plot_anti_inf$pay, plot_anti_inf[,2], type = 'l', ylim = c(0, max(plot_anti_inf[,2])),
     ylab = paste(colnames(plot_anti_inf)[2], "Prescriptions"), xlab = "Cumulative Payments",
     main = paste("Anti-infective class:\n",
                  colnames(plot_anti_inf)[2], "Prescriptions vs Cumulative Payments"))

# tinidazole does start shooting up to meaningful numbers ($1500, 30 rx)
# but metronidazole is so much larger ($1500 and 150,000,000 Rx)



# antichol:  Anticholinergics for overactive bladder-------
summary(model_antichol)
plot_antichol <- cleandrug(model_antichol, 6000, dat_antichol_nona)
plot(plot_antichol$pay, plot_antichol[,3], type = 'l', ylim = c(0, max(plot_antichol[,3])), 
     ylab = paste(colnames(plot_antichol)[3], "Prescriptions"), xlab = "Cumulative Payments",
     main = paste("Anticholinergics class:\n",
                  colnames(plot_antichol)[3], "Prescriptions vs Cumulative Payments"))

lines(plot_antichol$pay, plot_antichol[,9], col = "blue")
lines(plot_antichol$pay, plot_antichol[,1], col = "green")

plot(plot_antichol$pay, plot_antichol[,1], type = 'l', ylim = c(0, max(plot_antichol[,1])), col = "green")
lines(plot_antichol$pay, plot_antichol[,9], col = "blue")

# seeing most reasonable numbers out of this class, $6000 and 400 rx for 3rd column



# hormone:   Hormone therapy single ingredient ------------
summary(model_hormone)
plot_hormone  <- cleandrug(model_hormone, 1000, dat_hormone_nona)
plot(plot_hormone$pay, plot_hormone[,5], type = 'l', ylim = c(0, max(plot_hormone[,5])))
lines(plot_hormone$pay, plot_hormone[,7], col = "blue")

# Cpay is negative so numbers going down, doesn't look meaningful



# oral: Oral combined estrogen and progestin products for hormone therapy---
summary(model_oral)
plot_oral  <- cleandrug(model_oral, 1000, dat_oral_nona)
plot(plot_oral$pay, plot_oral[,9], type = 'l', ylim = c(0, max(plot_oral[,9])))

# again cpay is negative, and bigger negative effect than the last



# vaginal: Vaginal estrogen hormone therapy-----------------------
summary(model_vaginal)
plot_vaginal  <- cleandrug(model_vaginal, 50000, dat_vaginal_nona)
plot(plot_vaginal$pay, plot_vaginal[,6], type = 'l', ylim = c(0, max(plot_vaginal[,6])), 
     ylab = paste(colnames(plot_vaginal)[6], "Prescriptions"), xlab = "Cumulative Payments",
     main = paste("Vaginal estrogen class:\n",
                  colnames(plot_vaginal)[6], "Prescriptions vs Cumulative Payments"))
lines(plot_vaginal$pay, plot_vaginal[,1], col = "blue")
lines(plot_vaginal$pay, plot_vaginal[,7], col = "green")
lines(plot_vaginal$pay, plot_vaginal[,8], col = "red")

# even biggest effects arent meaningful until getting up really high payments, $50,000 for 300 rx

paste("hello", colnames(plot_vaginal)[6])


