# AF
# 08.04.2020
# Missingness analysis and MI
rm(list =ls())


library(mice)
library(rqdatatable)
library(MCMCglmm)


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


# MERGE TWO DEMO DATASETS------------------------
# KEEPING WHATS NON MISSING, PRIORITIZING THE NEW DEMOS

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
load_antiviral <- read.csv(file = 'Z:\\Pharma_Influence\\Data\\Updated_data_10_17_2019\\class_Antiviral.csv', header = T)
load_bisp     <- read.csv(file = 'Z:\\Pharma_Influence\\Data\\Updated_data_10_17_2019\\class_Bisphosphonates.csv', header = T)
load_hormone   <- read.csv(file = 'Z:\\Pharma_Influence\\Data\\Updated_data_10_17_2019\\class_Hormone_therapy_single_ingredient_therapy.csv', header = T)
load_oral      <- read.csv(file = 'Z:\\Pharma_Influence\\Data\\Updated_data_10_17_2019\\class_Oral_Combined_Estrogen_and_Progestin_Products_for_Hormone_Therapy.csv', header = T)
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
dat_antiviral_long <- cleandat(loaddata=load_antiviral)
dat_bisp_long      <- cleandat(loaddata=load_bisp)
dat_hormone_long   <- cleandat(loaddata=load_hormone)
dat_oral_long      <- cleandat(loaddata=load_oral)
dat_transderm_long <- cleandat(loaddata=load_transderm)
dat_vaginal_long   <- cleandat(loaddata=load_vaginal)



# MISSINGNESS ANALYSIS ON DEMFINAL-------------------------------------------------------------------------------------

# look to see if there is an associate with missing demos (age, region and gender) 
# and the payment or prescribing data (can look at values and just zero vs non zero)

# for significance testing, there is so much skew in the data that testing with 
# the actual payment or prescribing values isn't as robust, but it is more important
# to test if there are more zero vs nonzero payments for missing or non missing demos

# But also do I need to lower my significance because I'm doing multiple testing




# TAKE THE SUM OF PRESCRIPTIONS OVER ALL THE DIFFERENT DRUGS/YEARS
cleantotals <- function(longdata){
  
  one <- aggregate(cbind(TotPay=longdata$Pay, TotPre=longdata$Pre), by=list(NPI2=longdata$NPI2), FUN=sum)
  two <- merge(x=one, y=subset(longdata, select = -c(Pay,Pre,Cpay,Year,drug)), by = 'NPI2', all.x = T, all.y=F)
  three <- two[!duplicated(two$NPI2),]
  return(three)
  
}
dat_anti_inf <- cleantotals(dat_anti_inf_long)    #5354 physicians
dat_antichol <- cleantotals(dat_antichol_long)    #12002 total physicians
dat_antiviral <- cleantotals(dat_antiviral_long)  #72 total physicians
dat_bisp <- cleantotals(dat_bisp_long)            #9743 total physicians
dat_hormone <- cleantotals(dat_hormone_long)      #14941 total physicians
dat_oral <- cleantotals(dat_oral_long)            #5937 total physicians
dat_transderm <- cleantotals(dat_transderm_long)  #3953 total physicians
dat_vaginal <- cleantotals(dat_vaginal_long)      #21098 total physicians




# add non-zero indicators for payments and prescriptions
nonzind <- function(dataz){
  dataz$nonzpay <- ifelse(dataz$TotPay > 0, 1, 0)
  dataz$nonzpre <- ifelse(dataz$TotPre > 0, 1, 0)

  dataz$missGender <- ifelse(is.na(dataz$Gender), 1, 0)
  dataz$missGender <- as.factor(dataz$missGender)
  dataz$missAge    <- ifelse(is.na(dataz$Age), 1, 0)
  dataz$missAge    <- as.factor(dataz$missAge)
  dataz$missRegion <- ifelse(is.na(dataz$Region), 1, 0)
  dataz$missRegion <- as.factor(dataz$missRegion)
  return(dataz)
}

dat_anti_inf  <- nonzind(dataz=dat_anti_inf)
dat_antichol  <- nonzind(dataz=dat_antichol)
dat_antiviral <- nonzind(dataz=dat_antiviral)
dat_bisp      <- nonzind(dataz=dat_bisp)
dat_hormone   <- nonzind(dataz=dat_hormone)
dat_oral      <- nonzind(dataz=dat_oral)
dat_transderm <- nonzind(dataz=dat_transderm)
dat_vaginal   <- nonzind(dataz=dat_vaginal)







# ANTI INF-----------------------------------------------------------------------
plot(dat_anti_inf$Age, dat_anti_inf$TotPay)
plot(dat_anti_inf$Age, dat_anti_inf$TotPre, ylim = c(0,1500))

plot(dat_anti_inf$missAge, dat_anti_inf$nonzpay)
plot(dat_anti_inf$missAge, dat_anti_inf$TotPay)
plot(dat_anti_inf$missAge, dat_anti_inf$nonzpre) 
plot(dat_anti_inf$missAge, dat_anti_inf$TotPre, ylim = c(0,1500))

plot(dat_anti_inf$missRegion, dat_anti_inf$TotPay)
plot(dat_anti_inf$missRegion, dat_anti_inf$TotPre, ylim = c(0,1500))
plot(dat_anti_inf$missRegion, dat_anti_inf$nonzpay) #this was only sig one
plot(dat_anti_inf$missRegion, dat_anti_inf$nonzpre)

plot(dat_anti_inf$missGender, dat_anti_inf$nonzpay)
plot(dat_anti_inf$missGender, dat_anti_inf$TotPay)
plot(dat_anti_inf$missGender, dat_anti_inf$nonzpre)
plot(dat_anti_inf$missGender, dat_anti_inf$TotPre, ylim = c(0,1500))

# 
# 
# # wilcox test on pairs that could be significant---------
# 
part1 <- dat_anti_inf$TotPay[dat_anti_inf$missRegion==0]
part2 <- dat_anti_inf$TotPay[dat_anti_inf$missRegion==1]
t.test(part1, part2)
# 
pay1 <- dat_anti_inf$TotPay[dat_anti_inf$missGender==0]
pay2 <- dat_anti_inf$TotPay[dat_anti_inf$missGender==1]
t.test(pay1, pay2)
pay1 <- dat_anti_inf$TotPre[dat_anti_inf$missGender==0]
pay2 <- dat_anti_inf$TotPre[dat_anti_inf$missGender==1]
t.test(pay1, pay2)


# # tests to really look at 
pay1 <- dat_anti_inf$nonzpay[dat_anti_inf$missGender==0]
pay2 <- dat_anti_inf$nonzpay[dat_anti_inf$missGender==1]
t.test(pay1, pay2)
pay1 <- dat_anti_inf$nonzpre[dat_anti_inf$missGender==0]
pay2 <- dat_anti_inf$nonzpre[dat_anti_inf$missGender==1]
t.test(pay1, pay2) 
pay1 <- dat_anti_inf$nonzpay[dat_anti_inf$missAge==0]
pay2 <- dat_anti_inf$nonzpay[dat_anti_inf$missAge==1]
t.test(pay1, pay2)
pay1 <- dat_anti_inf$nonzpre[dat_anti_inf$missAge==0]
pay2 <- dat_anti_inf$nonzpre[dat_anti_inf$missAge==1]
t.test(pay1, pay2) 
pay1 <- dat_anti_inf$nonzpay[dat_anti_inf$missRegion==0]
pay2 <- dat_anti_inf$nonzpay[dat_anti_inf$missRegion==1]
t.test(pay1, pay2)                          #0.005669 non zero Pay and missing Region
pay1 <- dat_anti_inf$nonzpre[dat_anti_inf$missRegion==0]
pay2 <- dat_anti_inf$nonzpre[dat_anti_inf$missRegion==1]
t.test(pay1, pay2) 




# # ANTI CHOL---------------------------------------------------------------------------
plot(dat_antichol$Age, dat_antichol$TotPay)
plot(dat_antichol$Age, dat_antichol$TotPre, ylim = c(0,1500))

plot(dat_antichol$missAge, dat_antichol$nonzpay)
plot(dat_antichol$missAge, dat_antichol$TotPay)
plot(dat_antichol$missAge, dat_antichol$nonzpre)
plot(dat_antichol$missAge, dat_antichol$TotPre, ylim = c(0,3000))

plot(dat_antichol$missRegion, dat_antichol$TotPay)
plot(dat_antichol$missRegion, dat_antichol$TotPre, ylim = c(0,1500))
plot(dat_antichol$missRegion, dat_antichol$nonzpay)
plot(dat_antichol$missRegion, dat_antichol$nonzpre)

plot(dat_antichol$missGender, dat_antichol$nonzpay)
plot(dat_antichol$missGender, dat_antichol$TotPay)
plot(dat_antichol$missGender, dat_antichol$nonzpre)
plot(dat_antichol$missGender, dat_antichol$TotPre, ylim=c(0,3000))



# test on pairs that could be significant---------

pay1 <- dat_antichol$TotPay[dat_antichol$missRegion==0]
pay2 <- dat_antichol$TotPay[dat_antichol$missRegion==1]
t.test(pay1, pay2) #1.047e-06
pay1 <- dat_antichol$TotPre[dat_antichol$missRegion==0]
pay2 <- dat_antichol$TotPre[dat_antichol$missRegion==1]
t.test(pay1, pay2) 
pay1 <- dat_antichol$TotPay[dat_antichol$missAge==0]
pay2 <- dat_antichol$TotPay[dat_antichol$missAge==1]
t.test(pay1, pay2) #9.615e-08 for pay and missing age
pay1 <- dat_antichol$TotPre[dat_antichol$missAge==0]
pay2 <- dat_antichol$TotPre[dat_antichol$missAge==1]
t.test(pay1, pay2) #1.182e-08

pay1 <- dat_antichol$TotPay[dat_antichol$missGender==0]
pay2 <- dat_antichol$TotPay[dat_antichol$missGender==1]
t.test(pay1, pay2) #1.148e-07 Pay and missing Gender
pay1 <- dat_antichol$TotPre[dat_antichol$missGender==0]
pay2 <- dat_antichol$TotPre[dat_antichol$missGender==1]
t.test(pay1, pay2) #6.688e-08

# tests to really look at
pay1 <- dat_antichol$nonzpay[dat_antichol$missGender==0]
pay2 <- dat_antichol$nonzpay[dat_antichol$missGender==1]
t.test(pay1, pay2) #4.28e-13 non zero Pay and missing Gender
pay1 <- dat_antichol$nonzpre[dat_antichol$missGender==0]
pay2 <- dat_antichol$nonzpre[dat_antichol$missGender==1]
t.test(pay1, pay2)
pay1 <- dat_antichol$nonzpay[dat_antichol$missAge==0]
pay2 <- dat_antichol$nonzpay[dat_antichol$missAge==1]
t.test(pay1, pay2) #4.432e-09 non zero Pay and missing Age
pay1 <- dat_antichol$nonzpre[dat_antichol$missAge==0]
pay2 <- dat_antichol$nonzpre[dat_antichol$missAge==1]
t.test(pay1, pay2)
pay1 <- dat_antichol$nonzpay[dat_antichol$missRegion==0]
pay2 <- dat_antichol$nonzpay[dat_antichol$missRegion==1]
t.test(pay1, pay2) #7.131e-11
pay1 <- dat_antichol$nonzpre[dat_antichol$missRegion==0]
pay2 <- dat_antichol$nonzpre[dat_antichol$missRegion==1]
t.test(pay1, pay2) #3.094e-09




# # ANTI VIRAL---------------------------------------------------------------------------
# # has no payments so not going to bother with this class


# # BISP---------------------------------------------------------------------------
plot(dat_bisp$Age, dat_bisp$TotPay, ylim = c(0,10000))
plot(dat_bisp$Age, dat_bisp$TotPre)

plot(dat_bisp$missAge, dat_bisp$nonzpay)
plot(dat_bisp$missAge, dat_bisp$TotPay, ylim=c(0,200))
plot(dat_bisp$missAge, dat_bisp$nonzpre)
plot(dat_bisp$missAge, dat_bisp$TotPre)

plot(dat_bisp$missRegion, dat_bisp$TotPay)
plot(dat_bisp$missRegion, dat_bisp$TotPre, ylim = c(0,1500))
plot(dat_bisp$missRegion, dat_bisp$nonzpay)
plot(dat_bisp$missRegion, dat_bisp$nonzpre)

plot(dat_bisp$missGender, dat_bisp$nonzpay)
plot(dat_bisp$missGender, dat_bisp$TotPay)
plot(dat_bisp$missGender, dat_bisp$nonzpre)
plot(dat_bisp$missGender, dat_bisp$TotPre)


# test on pairs that could be significant---------

pay1 <- dat_bisp$TotPay[dat_bisp$missAge==0]
pay2 <- dat_bisp$TotPay[dat_bisp$missAge==1]
t.test(pay1,pay2)
pay1 <- dat_bisp$TotPre[dat_bisp$missAge==0]
pay2 <- dat_bisp$TotPre[dat_bisp$missAge==1]
t.test(pay1, pay2)
pay1 <- dat_bisp$TotPay[dat_bisp$missRegion==0]
pay2 <- dat_bisp$TotPay[dat_bisp$missRegion==1]
t.test(pay1, pay2) #0.02735
pay1 <- dat_bisp$TotPre[dat_bisp$missRegion==0]
pay2 <- dat_bisp$TotPre[dat_bisp$missRegion==1]
t.test(pay1, pay2) 

pay1 <- dat_bisp$TotPay[dat_bisp$missGender==0]
pay2 <- dat_bisp$TotPay[dat_bisp$missGender==1]
t.test(pay1, pay2) 
pay1 <- dat_bisp$TotPre[dat_bisp$missGender==0]
pay2 <- dat_bisp$TotPre[dat_bisp$missGender==1]
t.test(pay1, pay2) 


# tests to really look at
pay1 <- dat_bisp$nonzpay[dat_bisp$missGender==0]
pay2 <- dat_bisp$nonzpay[dat_bisp$missGender==1]
t.test(pay1, pay2) #0.00634 but not actually meaningful
pay1 <- dat_bisp$nonzpre[dat_bisp$missGender==0]
pay2 <- dat_bisp$nonzpre[dat_bisp$missGender==1]
t.test(pay1, pay2) 
pay1 <- dat_bisp$nonzpay[dat_bisp$missAge==0]
pay2 <- dat_bisp$nonzpay[dat_bisp$missAge==1]
t.test(pay1, pay2) 
pay1 <- dat_bisp$nonzpre[dat_bisp$missAge==0]
pay2 <- dat_bisp$nonzpre[dat_bisp$missAge==1]
t.test(pay1, pay2) #8.502e-05 but not meaningful
pay1 <- dat_bisp$nonzpay[dat_bisp$missRegion==0]
pay2 <- dat_bisp$nonzpay[dat_bisp$missRegion==1]
t.test(pay1, pay2)
pay1 <- dat_bisp$nonzpre[dat_bisp$missRegion==0]
pay2 <- dat_bisp$nonzpre[dat_bisp$missRegion==1]
t.test(pay1, pay2) #0.0001994 non zero Pre and missing Region

# 
# 
# # HORMONE---------------------------------------------------------------------------
plot(dat_hormone$Age, dat_hormone$TotPay)
plot(dat_hormone$Age, dat_hormone$TotPre)

plot(dat_hormone$missAge, dat_hormone$nonzpay)
plot(dat_hormone$missAge, dat_hormone$TotPay)
plot(dat_hormone$missAge, dat_hormone$nonzpre)
plot(dat_hormone$missAge, dat_hormone$TotPre)

plot(dat_hormone$missRegion, dat_hormone$TotPay)
plot(dat_hormone$missRegion, dat_hormone$TotPre)
plot(dat_hormone$missRegion, dat_hormone$nonzpay)
plot(dat_hormone$missRegion, dat_hormone$nonzpre)

plot(dat_hormone$missGender, dat_hormone$nonzpay)
plot(dat_hormone$missGender, dat_hormone$TotPay)
plot(dat_hormone$missGender, dat_hormone$nonzpre)
plot(dat_hormone$missGender, dat_hormone$TotPre)


# tests to really look at
pay1 <- dat_hormone$nonzpay[dat_hormone$missGender==0]
pay2 <- dat_hormone$nonzpay[dat_hormone$missGender==1]
t.test(pay1, pay2)
pay1 <- dat_hormone$nonzpre[dat_hormone$missGender==0]
pay2 <- dat_hormone$nonzpre[dat_hormone$missGender==1]
t.test(pay1, pay2) #0.021 non zero Pre and missing Gender not meaningful
pay1 <- dat_hormone$nonzpay[dat_hormone$missAge==0]
pay2 <- dat_hormone$nonzpay[dat_hormone$missAge==1]
t.test(pay1, pay2)
pay1 <- dat_hormone$nonzpre[dat_hormone$missAge==0]
pay2 <- dat_hormone$nonzpre[dat_hormone$missAge==1]
t.test(pay1, pay2) #1.284e-08 non zero Pre and missing Age
pay1 <- dat_hormone$nonzpay[dat_hormone$missRegion==0]
pay2 <- dat_hormone$nonzpay[dat_hormone$missRegion==1]
t.test(pay1, pay2) 
pay1 <- dat_hormone$nonzpre[dat_hormone$missRegion==0]
pay2 <- dat_hormone$nonzpre[dat_hormone$missRegion==1]
t.test(pay1, pay2)
# 
# 
# 
# 
# # ORAL---------------------------------------------------------------------------
plot(dat_oral$Age, dat_oral$TotPay)
plot(dat_oral$Age, dat_oral$TotPre)

plot(dat_oral$missAge, dat_oral$nonzpay)
plot(dat_oral$missAge, dat_oral$TotPay)
plot(dat_oral$missAge, dat_oral$nonzpre)
plot(dat_oral$missAge, dat_oral$TotPre)

plot(dat_oral$missRegion, dat_oral$TotPay)
plot(dat_oral$missRegion, dat_oral$TotPre)
plot(dat_oral$missRegion, dat_oral$nonzpay)
plot(dat_oral$missRegion, dat_oral$nonzpre)

plot(dat_oral$missGender, dat_oral$nonzpay)
plot(dat_oral$missGender, dat_oral$TotPay)
plot(dat_oral$missGender, dat_oral$nonzpre)
plot(dat_oral$missGender, dat_oral$TotPre)


# tests to really look at --NONE SIG!
pay1 <- dat_oral$nonzpay[dat_oral$missGender==0]
pay2 <- dat_oral$nonzpay[dat_oral$missGender==1]
t.test(pay1, pay2) #0.01496
pay1 <- dat_oral$nonzpre[dat_oral$missGender==0]
pay2 <- dat_oral$nonzpre[dat_oral$missGender==1]
t.test(pay1, pay2) #0.03207 not meaningful
pay1 <- dat_oral$nonzpay[dat_oral$missAge==0]
pay2 <- dat_oral$nonzpay[dat_oral$missAge==1]
t.test(pay1, pay2) #0.01371
pay1 <- dat_oral$nonzpre[dat_oral$missAge==0]
pay2 <- dat_oral$nonzpre[dat_oral$missAge==1]
t.test(pay1, pay2) #0.02255 not meaningful
pay1 <- dat_oral$nonzpay[dat_oral$missRegion==0]
pay2 <- dat_oral$nonzpay[dat_oral$missRegion==1]
t.test(pay1, pay2)
pay1 <- dat_oral$nonzpre[dat_oral$missRegion==0]
pay2 <- dat_oral$nonzpre[dat_oral$missRegion==1]
t.test(pay1, pay2)
# 
# 
# 
# # TRANSDERM---------------------------------------------------------------------------
# # no payments in this drug group too so not going to bother
# 
# 
# 
# # VAGINAL---------------------------------------------------------------------------
plot(dat_vaginal$Age, dat_vaginal$TotPay)
plot(dat_vaginal$Age, dat_vaginal$TotPre)

plot(dat_vaginal$missAge, dat_vaginal$nonzpay)
plot(dat_vaginal$missAge, dat_vaginal$TotPay)
plot(dat_vaginal$missAge, dat_vaginal$nonzpre)
plot(dat_vaginal$missAge, dat_vaginal$TotPre)

plot(dat_vaginal$missRegion, dat_vaginal$TotPay)
plot(dat_vaginal$missRegion, dat_vaginal$TotPre)
plot(dat_vaginal$missRegion, dat_vaginal$nonzpay)
plot(dat_vaginal$missRegion, dat_vaginal$nonzpre)

plot(dat_vaginal$missGender, dat_vaginal$nonzpay)
plot(dat_vaginal$missGender, dat_vaginal$TotPay)
plot(dat_vaginal$missGender, dat_vaginal$nonzpre)
plot(dat_vaginal$missGender, dat_vaginal$TotPre)


# tests to really look at --ALL SIGNIFICANT BUT THE DIFFERENCES ARE REALLY SMALL
pay1 <- dat_vaginal$nonzpay[dat_vaginal$missGender==0]
pay2 <- dat_vaginal$nonzpay[dat_vaginal$missGender==1]
t.test(pay1, pay2) #0.002058 non zero Pay and missing Gender
pay1 <- dat_vaginal$nonzpre[dat_vaginal$missGender==0]
pay2 <- dat_vaginal$nonzpre[dat_vaginal$missGender==1]
t.test(pay1, pay2) #31.219e-06 non zero Pre and missing Gender
pay1 <- dat_vaginal$nonzpay[dat_vaginal$missAge==0]
pay2 <- dat_vaginal$nonzpay[dat_vaginal$missAge==1]
t.test(pay1, pay2) #0.01602 non zero Pay and missing Age
pay1 <- dat_vaginal$nonzpre[dat_vaginal$missAge==0]
pay2 <- dat_vaginal$nonzpre[dat_vaginal$missAge==1]
t.test(pay1, pay2) #< 2.2e-16 non zero Pre and missing Age
pay1 <- dat_vaginal$nonzpay[dat_vaginal$missRegion==0]
pay2 <- dat_vaginal$nonzpay[dat_vaginal$missRegion==1]
t.test(pay1, pay2) #2.027e-05 non zeroo Pay and missing Region
pay1 <- dat_vaginal$nonzpre[dat_vaginal$missRegion==0]
pay2 <- dat_vaginal$nonzpre[dat_vaginal$missRegion==1]
t.test(pay1, pay2) #0.001162 non zero Pre and missing Region




# IMPUTING--------------------------------------------------------------------------------------------------


library(randomForest)

# Prescriptions are the outcomes but because that is longitudinal over several years and several drugs
# I will just be using the total prescriptions but im not sure if that's okay to do



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

# hormone and vaginal did not work "error: cannot allocate vector of size..."


# MERGE DATA---------------------------------
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
dat_bisp_final      <- cleandat(imputedata=dat_bisp_impute,     longdata= dat_bisp_long)
# dat_hormone_final   <- cleandat(imputedata=dat_hormone_impute,   longdata= dat_hormone_long)
dat_transderm_final <- cleandat(imputedata=dat_transderm_impute, longdata= dat_transderm_long)
# dat_vaginal_final   <- cleandat(imputedata=dat_vaginal_impute,   longdata= dat_vaginal_long)






# RUN MODELS WITH ALL DEMOS--------------------------------------------------------------------
# NOTE SCHOOL YEAR IS REMOVED AND GENDER AND REGION ARE CATEGORICAL
modeldat <- function(loaddata){
  
  model <<-MCMCglmm(data=loaddata, Pre~ drug + Cpay + Gender2 + Age + Region2, random = ~NPI2, 
                    family="poisson",burnin=5000,nitt=50000,thin=10) 
  
  return(model)
}

# model_anti_inf  <- modeldat(loaddata=dat_anti_inf_final)
# model_antichol  <- modeldat(loaddata=dat_antichol_final)
# model_bisp      <- modeldat(loaddata=dat_bisp_final)     
# model_hormone   <- modeldat(loaddata=dat_hormone_final)
# model_transderm <- modeldat(loaddata=dat_transderm_final)
# model_vaginal   <- modeldat(loaddata=dat_vaginal_final)
# # these all took 24 hours to run total last time, but there is more data now with imputation
# 
# 
# # RESULTS SUMMARY-----------------------------
# summary(model_anti_inf)  # gender sig
# summary(model_antichol)  # nothing sig
# summary(model_bisp)      # gender sig (last time no results, needed rerun)
# summary(model_hormone)   # region sig
# summary(model_transderm) # nothing sig (las time gender and age sig)
# summary(model_vaginal)   # nothing sig


# RUN MODEL WITH THE MISSINGS REMOVE FOR SENSITIVITY ANALYSIS


# create complete case
dat_anti_inf_finalcc <- dat_anti_inf_long[(!is.na(dat_anti_inf_long$Region)) & (dat_anti_inf_long$Gender!='') & 
                                             (!is.na(dat_anti_inf_long$Age)) & (!is.na(dat_anti_inf_long$SchoolYear)),]


model_anti_inf  <- modeldat(loaddata=dat_anti_inf_final)
model_anti_infcc  <- modeldat(loaddata=dat_anti_inf_finalcc)


summary(model_anti_inf) 
summary(model_anti_infcc) 

# WRITE OUT COEFFICIENTS FROM EACH MODEL--------------------------------------
save_anti_inf_impute <- data.frame(summary(model_anti_inf)$solutions)  #save coefficents
write.csv(save_anti_inf_impute, "Z:\\Pharma_Influence\\Results\\Model impute results 08.18.2020\\anti_inf_impute_coeff.csv")

save_anti_inf_cc <- data.frame(summary(model_anti_infna)$solutions)  #save coefficents
write.csv(save_anti_inf_cc, "Z:\\Pharma_Influence\\Results\\Model impute results 08.18.2020\\anti_inf_cc_coeff.csv")




