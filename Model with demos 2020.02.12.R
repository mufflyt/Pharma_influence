# AF
# 2.12.2020
# MCMC multi response Poisson model looking at adding demos to the model
rm(list =ls())


library(MCMCglmm)


# LOAD DATA-----------------------------------
load_anti_inf  <- read.csv(file = 'Z:\\Pharma_Influence\\Data\\Updated_data_10_17_2019\\class_Anti_infective.csv', header = T)
load_antichol  <- read.csv(file = 'Z:\\Pharma_Influence\\Data\\Updated_data_10_17_2019\\class_Anticholinergics_for_overactive_bladder.csv', header = T)
load_antiviral <- read.csv(file = 'Z:\\Pharma_Influence\\Data\\Updated_data_10_17_2019\\class_Antiviral.csv', header = T)
load_bisph     <- read.csv(file = 'Z:\\Pharma_Influence\\Data\\Updated_data_10_17_2019\\class_Bisphosphonates.csv', header = T)
load_hormone   <- read.csv(file = 'Z:\\Pharma_Influence\\Data\\Updated_data_10_17_2019\\class_Hormone_therapy_single_ingredient_therapy.csv', header = T)
load_oral      <- read.csv(file = 'Z:\\Pharma_Influence\\Data\\Updated_data_10_17_2019\\class_Oral_Combined_Estrogen_and_Progestin_Products_for_Hormone_Therapy.csv', header = T)
load_transderm <- read.csv(file = 'Z:\\Pharma_Influence\\Data\\Updated_data_10_17_2019\\class_Transdermal_estrogen.csv', header = T)
load_vaginal   <- read.csv(file = 'Z:\\Pharma_Influence\\Data\\Updated_data_10_17_2019\\class_Vaginal_Estrogen_Hormone_Therapy.csv', header = T)



# LOAD AND CLEAN DEMOS------------------------
dem <- read.csv('Z:\\Pharma_Influence\\Data\\Physician_demographics\\GOBA_unique.csv')

# remove NA NPI or NA Year
dem <- dem[!is.na(dem$NPI_Match),]
demid <- as.data.frame(unique(dem$NPI_Match))
# there are 53533 in dem, but 52605 unique NPI, so 928 duplicates?

# only need to merge a subset of the demos
demsub <- dem[, c('NPI_Match', 'GOBA_Region', 'HG_Gender', 'HG_Age', 'HG_Medical.School.1...Year')]
colnames(demsub) <- c('NPI', 'Region', 'Gender', 'Age', 'SchoolYear')

# MISSINGS: 433 missing regions, 0 missing gender, 23250 missing age, 22858 missing year
test <- demsub[is.na(demsub$Region),]
test <- demsub[is.na(demsub$Gender),]
test <- demsub[is.na(demsub$Age),]
test <- demsub[is.na(demsub$SchoolYear),]
rm(test)

demsubna <- demsub[(!is.na(demsub$Region)) & (demsub$Gender!='') & (!is.na(demsub$Age)) & (!is.na(demsub$SchoolYear)),]
# only 27058 when remove all missing
# demsubna2 <- demsub[(!is.na(demsub$Region)) & ((demsub$Gender!='') ),]
# 53100 when remove missing from just region and gender






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
  datdem <- merge(x=datl, y=demsubna, by = 'NPI')
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

colnames(dat_anti_inf_npi_dem)

# (length(unique(load_anti_inf$NPI))-length(unique(dat_anti_inf2$NPI)))/length(unique(load_anti_inf$NPI))*100

# WHEN DEMOS IS SUBSET TO NO MISSING
# anti_inf lost 49%
# anitchol lost 47%
# antiviral lost 49%
# bisph lost 49%
# hormone lost 46%
# oral lost 41%
# transderm lost 39%
# vaginal lost 48%

# WHEN DEMOS IS SUBSET TO ONLY REGION/GENDER MISSING
# anti_inf lost 37%
# anitchol lost 35%
# antiviral lost 37%
# bisph lost 37%
# hormone lost 34%
# oral lost 33%
# transderm lost 30%
# vaginal lost 35%




# OUTPUT THE NPIS THAT HAVE AND DO NOT HAVE DEMO DATA FOR TYLER-------------------------------------

# subset to unique from each class that have complete demos
dat_anti_inf_npi_dem <- data.frame(unique(dat_anti_inf$NPI))
colnames(dat_anti_inf_npi_dem) <- 'NPI'
dat_antichol_npi_dem <- data.frame(unique(dat_antichol$NPI))
colnames(dat_antichol_npi_dem) <- 'NPI'
dat_antiviral_npi_dem <- data.frame(unique(dat_antiviral$NPI))
colnames(dat_antiviral_npi_dem) <- 'NPI'
dat_bisp_npi_dem <- data.frame(unique(dat_bisp$NPI))
colnames(dat_bisp_npi_dem) <- 'NPI'
dat_hormone_npi_dem <- data.frame(unique(dat_hormone$NPI))
colnames(dat_hormone_npi_dem) <- 'NPI'
dat_oral_npi_dem <- data.frame(unique(dat_oral$NPI))
colnames(dat_oral_npi_dem) <- 'NPI'
dat_transderm_npi_dem <- data.frame(unique(dat_transderm$NPI))
colnames(dat_transderm_npi_dem) <- 'NPI'
dat_vaginal_npi_dem <- data.frame(unique(dat_vaginal$NPI))
colnames(dat_vaginal_npi_dem) <- 'NPI'

# stack them all
dat_all_npi_dem <- rbind(dat_anti_inf_npi_dem, dat_antichol_npi_dem, dat_antiviral_npi_dem, dat_bisp_npi_dem,
                         dat_hormone_npi_dem, dat_oral_npi_dem, dat_transderm_npi_dem, dat_vaginal_npi_dem)
# subset to unique across all classes that have demos
npi_with_dem <- data.frame(unique(dat_all_npi_dem$NPI))
colnames(npi_with_dem) <- "NPI"
npi_with_dem$ind <- 1

### make list off all NPIs and to then make list of those not in dem list

# subset to unique from each class from all data
load_anti_inf_npi_dem <- data.frame(unique(load_anti_inf$NPI))
colnames(load_anti_inf_npi_dem) <- 'NPI'
load_antichol_npi_dem <- data.frame(unique(load_antichol$NPI))
colnames(load_antichol_npi_dem) <- 'NPI'
load_antiviral_npi_dem <- data.frame(unique(load_antiviral$NPI))
colnames(load_antiviral_npi_dem) <- 'NPI'
load_bisp_npi_dem <- data.frame(unique(load_bisph$NPI))
colnames(load_bisp_npi_dem) <- 'NPI'
load_hormone_npi_dem <- data.frame(unique(load_hormone$NPI))
colnames(load_hormone_npi_dem) <- 'NPI'
load_oral_npi_dem <- data.frame(unique(load_oral$NPI))
colnames(load_oral_npi_dem) <- 'NPI'
load_transderm_npi_dem <- data.frame(unique(load_transderm$NPI))
colnames(load_transderm_npi_dem) <- 'NPI'
load_vaginal_npi_dem <- data.frame(unique(load_vaginal$NPI))
colnames(load_vaginal_npi_dem) <- 'NPI'

# stack them all
load_all_npi_dem <- rbind(load_anti_inf_npi_dem, load_antichol_npi_dem, load_antiviral_npi_dem, load_bisp_npi_dem,
                         load_hormone_npi_dem, load_oral_npi_dem, load_transderm_npi_dem, load_vaginal_npi_dem)
# subset to unique across all classes
load_all_npi_dem_nodup <- data.frame(unique(load_all_npi_dem$NPI))
colnames(load_all_npi_dem_nodup) <- 'NPI'


# merge total list with demo list and then remove those that are int dem for list of all NPIs missing demos
npi_miss_dem <- merge(x=load_all_npi_dem_nodup, y=npi_with_dem, all=T)
npi_miss_dem$ind <- ifelse(is.na(npi_miss_dem$ind), 0, npi_miss_dem$ind)
npi_miss_dem <- subset(npi_miss_dem, npi_miss_dem$ind ==0 & npi_miss_dem$NPI != 0)
# 11709 with dem
# 23160 in all
# 11451 without dem (there was one NPI0 that got removed)


# output
write.csv(npi_miss_dem, 'Z:\\Pharma_Influence\\Data\\Demo NPIs\\npi_miss_demos.csv')
write.csv(npi_with_dem, 'Z:\\Pharma_Influence\\Data\\Demo NPIs\\npi_with_demos.csv')




# RUN MODELS WITH ALL DEMOS--------------------------------------------------------------------
modeldat <- function(loaddata){
  
  model <<-MCMCglmm(data=loaddata, Pre~ drug + Cpay + Gender + Age + Region + SchoolYear, random = ~NPI2, 
                    family="poisson",burnin=5000,nitt=50000,thin=10) 
  
  return(model)
}


model_anti_inf  <- modeldat(loaddata=dat_anti_inf)
model_antichol  <- modeldat(loaddata=dat_antichol)
model_antiviral <- modeldat(loaddata=dat_antiviral)
model_bisp      <- modeldat(loaddata=dat_bisp)     #need to rerun
model_hormone   <- modeldat(loaddata=dat_hormone)
model_oral      <- modeldat(loaddata=dat_oral)
model_transderm <- modeldat(loaddata=dat_transderm)
model_vaginal   <- modeldat(loaddata=dat_vaginal)




# RESULTS SUMMARY-----------------------------
# all of the models had a warning about some effects not being estimable
summary(model_anti_inf)  # gender sig
summary(model_antichol)  # nothing sig
summary(model_antiviral) # gender close but nothing sig
# summary(model_bisp)
summary(model_hormone)   # region sig
summary(model_oral)      # region sig, gender close
summary(model_transderm) # gender and age sig
summary(model_vaginal)   # nothing sig




# SAVE MODEL RESULTS-----------------------------
saveRDS(model_anti_inf, "Z:\\Pharma_Influence\\Results\\Model demo results 02.13.2020\\anti_inf_demo.5000.50000.10.poisson.rds")
saveRDS(model_antichol, "Z:\\Pharma_Influence\\Results\\Model demo results 02.13.2020\\antichol_demo.5000.50000.10.poisson.rds")
saveRDS(model_antiviral, "Z:\\Pharma_Influence\\Results\\Model demo results 02.13.2020\\antiviral_demo.5000.50000.10.poisson.rds")
saveRDS(model_bisp, "Z:\\Pharma_Influence\\Results\\Model demo results 02.13.2020\\bisp_demo.5000.50000.10.poisson.rds")
saveRDS(model_hormone, "Z:\\Pharma_Influence\\Results\\Model demo results 02.13.2020\\hormone_demo.5000.50000.10.poisson.rds")
saveRDS(model_oral, "Z:\\Pharma_Influence\\Results\\Model demo results 02.13.2020\\oral_demo.5000.50000.10.poisson.rds")
saveRDS(model_transderm, "Z:\\Pharma_Influence\\Results\\Model demo results 02.13.2020\\transderm_demo.5000.50000.10.poisson.rds")
saveRDS(model_vaginal, "Z:\\Pharma_Influence\\Results\\Model demo results 02.13.2020\\vaginal_demo.5000.50000.10.poisson.rds")







