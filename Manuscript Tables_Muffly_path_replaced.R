# JCS
# 08.19.2020
############################################################################## 
#From AF New Demos
rm(list =ls())

library(MCMCglmm)
library(rqdatatable)

# LOOK AT NEW DEMOS
StudyGroup <- readRDS("//rvffls01/STATSOBGYN/Pharma_Influence/Data/Updated data 7_7_2020/StudyGroup.rds")
PaySum <- readRDS("//rvffls01/STATSOBGYN/Pharma_Influence/Data/Updated data 7_7_2020/PaySum.rds")
Prescriber <- readRDS("//rvffls01/STATSOBGYN/Pharma_Influence/Data/Updated data 7_7_2020/Prescriber.rds")


dem <- read.csv("//rvffls01/STATSOBGYN/Pharma_Influence/Data/Updated data 7_7_2020/T1.csv")
dem_old <- read.csv('//rvffls01/STATSOBGYN/Pharma_Influence\\Data\\Physician_demographics\\GOBA_unique.csv')


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
library(reshape)
demsub$Male<-ifelse(demsub$Gender=='M',1,0)
demsub_m<-melt(demsub,id.vars=c('NPI','Region'),measure.vars=c('Male','SchoolYear','Age'))
#demsub_c<-cast(unique(demsub_m[!is.na(demsub_m$value),]),NPI+Region~variable,value='value')

##Two unique School Years, 1987 & 1988 used most recent
#dem_oldsub[dem_oldsub$NPI==1255401733,'SchoolYear']<-1988

dem_oldsub$Male<-ifelse(dem_oldsub$Gender=='M',1,0)
dem_oldsub_m<-melt(dem_oldsub,id.vars=c('NPI','Region'),measure.vars=c('Male','SchoolYear','Age'))
#dem_oldsub_c<-cast(unique(dem_oldsub_m[!is.na(dem_oldsub_m$value),]),NPI+Region~variable,value='value')

demfinal_bind<-unique(rbind(demsub_m,dem_oldsub_m))
demfinal<-cast(unique(demfinal_bind[!is.na(demfinal_bind$value),]),NPI+Region~variable,value='value',fun.aggregate=max,na.rm=T)
# coalesce type merge

# demfinal <- natural_join(a=demsub, b=dem_oldsub, 
#                       by = "NPI",
#                       jointype = "FULL")
# demfinal$Gender <- ifelse(demfinal$Gender=='3', 'M', demfinal$Gender)
# demfinal$Gender <- ifelse(demfinal$Gender=='2', 'F', demfinal$Gender)
# demfinal$Gender <- ifelse(demfinal$Gender=='1', NA, demfinal$Gender)


length(unique(demfinal$NPI))
#53838
length(unique(Prescriber$NPI))
#19158

# MISSINGS: 1216 missing regions, 15908 missing gender, 20378 missing age, 20069 missing year
test <- demfinal[is.na(demfinal$Region),]
test <- demfinal[is.na(demfinal$Gender),]
test <- demfinal[is.na(demfinal$Age),]
test <- demfinal[is.na(demfinal$SchoolYear),]


# LOAD PAYMENT/PRESCRIBING DATA-----------------------------------
load_anti_inf  <- read.csv(file = '//rvffls01/STATSOBGYN/Pharma_Influence\\Data\\Updated_data_10_17_2019\\class_Anti_infective.csv', header = T)
load_antichol  <- read.csv(file = '//rvffls01/STATSOBGYN/Pharma_Influence\\Data\\Updated_data_10_17_2019\\class_Anticholinergics_for_overactive_bladder.csv', header = T)
load_antiviral <- read.csv(file = '//rvffls01/STATSOBGYN/Pharma_Influence\\Data\\Updated_data_10_17_2019\\class_Antiviral.csv', header = T)
load_bisph     <- read.csv(file = '//rvffls01/STATSOBGYN/Pharma_Influence\\Data\\Updated_data_10_17_2019\\class_Bisphosphonates.csv', header = T)
load_hormone   <- read.csv(file = '//rvffls01/STATSOBGYN/Pharma_Influence\\Data\\Updated_data_10_17_2019\\class_Hormone_therapy_single_ingredient_therapy.csv', header = T)
load_oral      <- read.csv(file = '//rvffls01/STATSOBGYN/Pharma_Influence\\Data\\Updated_data_10_17_2019\\class_Oral_Combined_Estrogen_and_Progestin_Products_for_Hormone_Therapy.csv', header = T)
load_transderm <- read.csv(file = '//rvffls01/STATSOBGYN/Pharma_Influence\\Data\\Updated_data_10_17_2019\\class_Transdermal_estrogen.csv', header = T)
load_vaginal   <- read.csv(file = '//rvffls01/STATSOBGYN/Pharma_Influence\\Data\\Updated_data_10_17_2019\\class_Vaginal_Estrogen_Hormone_Therapy.csv', header = T)


length(unique(c(load_anti_inf$NPI,load_antichol$NPI,load_antiviral$NPI,load_bisph$NPI,load_hormone$NPI,load_oral$NPI,load_transderm$NPI,load_vaginal$NPI)))
#23160


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

######################################## End AF

library(reshape)

c_anti_inf  <- cast(dat_anti_inf,NPI~.,value='Pay',fun.aggregate=sum,na.rm=T)
c_antichol  <- cast(dat_antichol,NPI~.,value='Pay',fun.aggregate=sum,na.rm=T)
c_antiviral <- cast(dat_antiviral,NPI~.,value='Pay',fun.aggregate=sum,na.rm=T)
c_bisp      <- cast(dat_bisp,NPI~.,value='Pay',fun.aggregate=sum,na.rm=T)
c_hormone   <- cast(dat_hormone,NPI~.,value='Pay',fun.aggregate=sum,na.rm=T)
c_oral      <- cast(dat_oral,NPI~.,value='Pay',fun.aggregate=sum,na.rm=T)
c_transderm <- cast(dat_transderm,NPI~.,value='Pay',fun.aggregate=sum,na.rm=T)
c_vaginal   <- cast(dat_vaginal,NPI~.,value='Pay',fun.aggregate=sum,na.rm=T)

total_pay<-rbind(c_anti_inf,c_antichol,c_antiviral,c_bisp,c_hormone,c_oral,c_transderm,c_vaginal)
names(total_pay)[2]<-'Pay'
ctotal_pay<-cast(total_pay,NPI~.,value='Pay',fun.aggregate=sum,na.rm=T)
names(ctotal_pay)[2]<-'Pay'
pay_ind<-merge(demfinal,ctotal_pay,by='NPI')
