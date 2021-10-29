# JCS
# 08.19.2020
############################################################################## 
#From AF New Demos
rm(list =ls())

library(MCMCglmm)
library(rqdatatable)

# LOOK AT NEW DEMOS
StudyGroup <- readRDS("Z:/Pharma_Influence/Data/Updated data 7_7_2020/StudyGroup.rds")
PaySum <- readRDS("Z:/Pharma_Influence/Data/Updated data 7_7_2020/PaySum.rds")
Prescriber <- readRDS("Z:/Pharma_Influence/Data/Updated data 7_7_2020/Prescriber.rds")


dem <- read.csv("Z:/Pharma_Influence/Data/Updated data 7_7_2020/T1.csv")
dem_old <- read.csv('Z:\\Pharma_Influence\\Data\\Physician_demographics\\GOBA_unique.csv')

dem_old$MD<-ifelse(dem_old$GOBA_Degree.1 %in% c("D. O.","D..","D.O","D.O.","Do","DO"),0,
                  ifelse(!dem_old$GOBA_Degree.1 %in% '',1,NA))
dem$MD<-ifelse(dem$GOBA_title %in% c("DO","DO MMS","DO MS"),0,
                  ifelse(!is.na(dem$GOBA_title),1,NA))

dem$School<-ifelse(dem$PCND_MedSchool %in% 'OTHER','International',
                  ifelse(!is.na(dem$PCND_MedSchool),'US',NA)


#dem_old$GOBA_Cert
# CLEAN OLD DATA----------------------------

# remove NA NPI or NA Year
dem_old <- dem_old[!is.na(dem_old$NPI_Match),]
dem_oldid <- as.data.frame(unique(dem_old$NPI_Match))
# there are 53533 in dem, but 52605 unique NPI, so 928 duplicates?

# only need to merge a subset of the demos
dem_oldsub <- dem_old[, c('NPI_Match', 'GOBA_Region', 'HG_Gender', 'HG_Age', 'HG_Medical.School.1...Year','MD')]
colnames(dem_oldsub) <- c('NPI', 'Region', 'Gender', 'Age', 'SchoolYear','MD')

dem_oldsubna <- dem_oldsub[(!is.na(dem_oldsub$Region)) & (dem_oldsub$Gender!='') & (!is.na(dem_oldsub$Age)) & (!is.na(dem_oldsub$SchoolYear)),]
# only 27058 when remove all missing


# CLEAN NEW DATA----------------------------

# remove NA NPI or NA Year
dem <- dem[!is.na(dem$NPI),]  # no missing npi number
demid <- as.data.frame(unique(dem$NPI))
# in this demo there are 12055 but 12047 unique NPI, so 8 duplicates? SHOULD LOOK INTO WHICH ONES TO BE REMOVING OF THE DUPS

# only need to merge a subset of the demos
demsub <- dem[, c('NPI', 'GOBA_Region', 'PCND_gender', 'Age', 'PCND_GradYr','MD')]
colnames(demsub) <- c('NPI', 'Region', 'Gender', 'Age', 'SchoolYear','MD')

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
demsub_m<-melt(demsub,id.vars=c('NPI','Region'),measure.vars=c('Male','SchoolYear','Age','MD'))
#demsub_c<-cast(unique(demsub_m[!is.na(demsub_m$value),]),NPI+Region~variable,value='value')

##Two unique School Years, 1987 & 1988 used most recent
#dem_oldsub[dem_oldsub$NPI==1255401733,'SchoolYear']<-1988

dem_oldsub$Male<-ifelse(dem_oldsub$Gender=='Male',1,0)
dem_oldsub_m<-melt(dem_oldsub,id.vars=c('NPI','Region'),measure.vars=c('Male','SchoolYear','Age','MD'))
#dem_oldsub_c<-cast(unique(dem_oldsub_m[!is.na(dem_oldsub_m$value),]),NPI+Region~variable,value='value')

demfinal_bind<-unique(rbind(demsub_m,dem_oldsub_m))
demfinal<-cast(unique(demfinal_bind[!is.na(demfinal_bind$value),]),NPI+Region~variable,value='value',fun.aggregate=max,na.rm=T)
# coalesce type merge

certs<-unique(dem_old[,c('NPI_Match','GOBA_Cert')])
certs$Specialist<-ifelse(certs$GOBA_Cert %in% 'Obstetrics and Gynecology',0,1)
demfinal<-merge(demfinal,certs[,c('NPI_Match','Specialist')],by.x='NPI',by.y='NPI_Match',all.x=T)
finites<-function(x){
  x[is.infinite(x)]<-NA
  return(x)
}
demfinal[,3:6]<-apply(demfinal[,3:6],2,finites)


length(unique(demfinal$NPI))
#53838
length(unique(Prescriber$NPI))
#19158

# MISSINGS: 1216 missing regions, 15908 missing gender, 20378 missing age, 20069 missing year
test <- demfinal[is.na(demfinal$Region),]
test <- demfinal[is.na(demfinal$Gender),]
test <- demfinal[is.na(demfinal$Age),]
test <- demfinal[is.na(demfinal$SchoolYear),]

############################################
# Change merging data 
#
##############################################
NPI.vals<-unique(c(unique(Prescriber$NPI),unique(PaySum$NPI)))
NPI.unique<-data.frame(NPI=NPI.vals,one=1)
length(unique(Prescriber$NPI))
dim(Prescriber)
length(unique(PaySum$NPI))
dim(PaySum)

# LOAD PAYMENT/PRESCRIBING DATA-----------------------------------
load_anti_inf  <- read.csv(file = 'Z:\\Pharma_Influence\\Data\\Updated_data_10_17_2019\\class_Anti_infective.csv', header = T)
load_antichol  <- read.csv(file = 'Z:\\Pharma_Influence\\Data\\Updated_data_10_17_2019\\class_Anticholinergics_for_overactive_bladder.csv', header = T)
load_antiviral <- read.csv(file = 'Z:\\Pharma_Influence\\Data\\Updated_data_10_17_2019\\class_Antiviral.csv', header = T)
load_bisph     <- read.csv(file = 'Z:\\Pharma_Influence\\Data\\Updated_data_10_17_2019\\class_Bisphosphonates.csv', header = T)
load_hormone   <- read.csv(file = 'Z:\\Pharma_Influence\\Data\\Updated_data_10_17_2019\\class_Hormone_therapy_single_ingredient_therapy.csv', header = T)
load_oral      <- read.csv(file = 'Z:\\Pharma_Influence\\Data\\Updated_data_10_17_2019\\class_Oral_Combined_Estrogen_and_Progestin_Products_for_Hormone_Therapy.csv', header = T)
load_transderm <- read.csv(file = 'Z:\\Pharma_Influence\\Data\\Updated_data_10_17_2019\\class_Transdermal_estrogen.csv', header = T)
load_vaginal   <- read.csv(file = 'Z:\\Pharma_Influence\\Data\\Updated_data_10_17_2019\\class_Vaginal_Estrogen_Hormone_Therapy.csv', header = T)


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
library(Hmisc)
library(table1)
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
pay_ind$PayInd<-pay_ind$Pay>0
pay_ind$SchoolDecade<-cut2(pay_ind$SchoolYear,c(1970,1980,1990,2000,2010,2020))
pay_ind$AgeDecade<-cut2(pay_ind$Age,c(30,40,50,60,70,80,120))
pay_ind_unique<-merge(pay_ind[!duplicated(pay_ind$NPI),],dem_old[,c('NPI_Match','GOBA_Cert')],by.x='NPI',by.y='NPI_Match',all.x=T)
############## Table 1 ###############



chisq.test(pay_ind_unique$Male  , pay_ind_unique$PayInd)$p.value
chisq.test(pay_ind_unique$SchoolDecade  , pay_ind_unique$PayInd)$p.value
chisq.test(pay_ind_unique$AgeDecade  , pay_ind_unique$PayInd)$p.value
chisq.test(pay_ind_unique$Region  , pay_ind_unique$PayInd)$p.value
chisq.test(pay_ind_unique$GOBA_Cert  , pay_ind_unique$PayInd)$p.value
chisq.test(pay_ind_unique$Specialist  , pay_ind_unique$PayInd)$p.value



#table1(~factor(Male) + factor(SchoolDecade) + factor(AgeDecade)+factor(MD) +factor(Region) +factor(GOBA_Cert)| PayInd,data=pay_ind_unique)


dat_anti_inf_nona <- dat_anti_inf[(!is.na(dat_anti_inf$Region)) & (dat_anti_inf$Male!='') & (!is.na(dat_anti_inf$Age)) & (!is.na(dat_anti_inf$SchoolYear)) & (!is.na(dat_anti_inf$Specialist)),]
dat_antichol_nona <- dat_antichol[(!is.na(dat_antichol$Region)) & (dat_antichol$Male!='') & (!is.na(dat_antichol$Age)) & (!is.na(dat_antichol$SchoolYear)) & (!is.na(dat_antichol$Specialist)),]
dat_antiviral_nona <- dat_antiviral[(!is.na(dat_antiviral$Region)) & (dat_antiviral$Male!='') & (!is.na(dat_antiviral$Age)) & (!is.na(dat_antiviral$SchoolYear)) & (!is.na(dat_antiviral$Specialist)),]
dat_bisp_nona <- dat_bisp[(!is.na(dat_bisp$Region)) & (dat_bisp$Male!='') & (!is.na(dat_bisp$Age)) & (!is.na(dat_bisp$SchoolYear)) & (!is.na(dat_bisp$Specialist)),]
dat_hormone_nona <- dat_hormone[(!is.na(dat_hormone$Region)) & (dat_hormone$Male!='') & (!is.na(dat_hormone$Age)) & (!is.na(dat_hormone$SchoolYear)) & (!is.na(dat_hormone$Specialist)),]
dat_oral_nona <- dat_oral[(!is.na(dat_oral$Region)) & (dat_oral$Male!='') & (!is.na(dat_oral$Age)) & (!is.na(dat_oral$SchoolYear)) & (!is.na(dat_oral$Specialist)),]
dat_transderm_nona <- dat_transderm[(!is.na(dat_transderm$Region)) & (dat_transderm$Male!='') & (!is.na(dat_transderm$Age)) & (!is.na(dat_transderm$SchoolYear)) & (!is.na(dat_transderm$Specialist)),]
dat_vaginal_nona <- dat_vaginal[(!is.na(dat_vaginal$Region)) & (dat_vaginal$Male!='') & (!is.na(dat_vaginal$Age)) & (!is.na(dat_vaginal$SchoolYear)) & (!is.na(dat_vaginal$Specialist)),]

nrow(dat_anti_inf_nona)/nrow(dat_anti_inf)    #60% of data left after complete case
nrow(dat_antichol_nona)/nrow(dat_antichol)    #63% left
nrow(dat_antiviral_nona)/nrow(dat_antiviral)  #63% left
nrow(dat_bisp_nona)/nrow(dat_bisp)            #62% 
nrow(dat_hormone_nona)/nrow(dat_hormone)      #65% 
nrow(dat_oral_nona)/nrow(dat_oral)            #67%
nrow(dat_transderm_nona)/nrow(dat_transderm)  #70%
nrow(dat_vaginal_nona)/nrow(dat_vaginal)      #63%

save.image('Z:\\Pharma_Influence\\Results\\results_john\\base_data.RData')
# RUN MODELS WITH ALL DEMOS--------------------------------------------------------------------
modeldat <- function(loaddata){
  
  model <<-MCMCglmm(data=loaddata, Pre~ as.factor(drug) + Cpay + Male + Age + Specialist+ as.factor(Region) , random = ~NPI2, 
                    family="poisson",burnin=5000,nitt=20000,thin=10) 
  
  return(model)
}


m1<-zeroinfl( Pre~ as.factor(drug) + Cpay + Male + Age + Specialist+ as.factor(Region),data=dat_anti_inf_nona)

modeldat <- function(loaddata){

m1<-MCMCglmm(Pre~ trait - 1 + at.level(trait,1):as.factor(drug) + at.level(trait, 1):Cpay + at.level(trait, 1):Male + at.level(trait, 1):Age + at.level(trait,1):Specialist + at.level(trait, 1):as.factor(Region), 
 rcov = ~idh(trait):units, data = loaddata,random=~idh(trait):NPI2,
  family = "zipoisson",
 verbose = FALSE,burnin=5000,nitt=20000,thin=10) 
return(m1)
}



model_anti_inf  <- modeldat(loaddata=dat_anti_inf_nona)
summary(model_anti_inf)  # gender sig
save(model_anti_inf,file = 'Z:\\Pharma_Influence\\Results\\results_john\\model_anti_inf_zip.RData')
model_antichol  <- modeldat(loaddata=dat_antichol_nona)
summary(model_antichol)  # nothing sig
save(model_antichol,file = 'Z:\\Pharma_Influence\\Results\\results_john\\model_antichol_zip.RData')
model_antiviral <- modeldat(loaddata=dat_antiviral_nona)
summary(model_antiviral) # gender close but nothing sig
save(model_antiviral,file = 'Z:\\Pharma_Influence\\Results\\results_john\\model_antiviral_zip.RData')
model_bisp      <- modeldat(loaddata=dat_bisp_nona)   
summary(model_bisp)      # gender sig (last time no results, needed rerun)
save(model_bisp,file = 'Z:\\Pharma_Influence\\Results\\results_john\\model_bisp_zip.RData')  
model_hormone   <- modeldat(loaddata=dat_hormone_nona)
summary(model_hormone)   # region sig
save(model_hormone,file = 'Z:\\Pharma_Influence\\Results\\results_john\\model_hormone_zip.RData')
model_oral      <- modeldat(loaddata=dat_oral_nona)
summary(model_oral)      # nothing sig (last time region sig, gender close)
save(model_oral,file = 'Z:\\Pharma_Influence\\Results\\results_john\\model_oral_zip.RData')
model_transderm <- modeldat(loaddata=dat_transderm_nona)
summary(model_transderm) # nothing sig (las time gender and age sig)
save(model_transderm,file = 'Z:\\Pharma_Influence\\Results\\results_john\\model_transderm_zip.RData')
model_vaginal   <- modeldat(loaddata=dat_vaginal_nona)
summary(model_vaginal)   # nothing sig
save(model_vaginal,file = 'Z:\\Pharma_Influence\\Results\\results_john\\model_vaginal_zip.RData')



# RESULTS SUMMARY-----------------------------
# all of the models had a warning about some effects not being estimable


load('Z:\\Pharma_Influence\\Results\\results_john\\model_anti_inf_zip.RData')
load('Z:\\Pharma_Influence\\Results\\results_john\\model_antichol_zip.RData')
load('Z:\\Pharma_Influence\\Results\\results_john\\model_antiviral_zip.RData')
load('Z:\\Pharma_Influence\\Results\\results_john\\model_bisp_zip.RData')  
load('Z:\\Pharma_Influence\\Results\\results_john\\model_hormone_zip.RData')
load('Z:\\Pharma_Influence\\Results\\results_john\\model_oral_zip.RData')
load('Z:\\Pharma_Influence\\Results\\results_john\\model_transderm_zip.RData')
load('Z:\\Pharma_Influence\\Results\\results_john\\model_vaginal_zip.RData')










make_pred_plot <- function(data, payN,lab,dataframe){
  
  save <- data.frame(summary(data)$solutions)                      # save coefficents
  num <- (1:dim(summary(data)$solutions)[1])[grepl(':Cpay',row.names(summary(data)$solutions))]-1       # only want to keep the drugs and cpay, not the four demos at the end
  numbers <- setNames(data.frame(t((save[3:num,1]))), row.names(save)[3:num])  #transpose to have column for each drug                        # get average age to use for average estimate
  meanage <- mean(dataframe$Age, na.rm=T)  
  numbersrep <- numbers[rep(seq_len(nrow(numbers)), each = payN),] + # repeat to get number of instances desired
                                    save[1,1]  +                             # add intercept
                                    save[num+3,1] + 
                                    save[num+2,1]*meanage             # add mean age

  numberspay <- save[num,1]*c(1:payN) #Cpay times payments
  numbersadd <- sweep(numbersrep, 1, numberspay, "+") #add to every drug column
  numbersexp <- exp(numbersadd)   # exponentiate
  numbersexp$pay <- c(1:payN)
  numbersexp.m<-melt(numbersexp,id.vars='pay')
  numbersexp.m$group<-lab
  return(numbersexp.m)
}


# bisph: Bisphosphonates---------------------------------
#summary(model)
pred_anti_inf <- make_pred_plot(model_anti_inf, 1000,lab='Anti-Infectives',dataframe=dat_anti_inf_nona)
pred_antichol <- make_pred_plot(model_antichol, 60000,lab='Anticholinergics',dataframe=dat_antichol_nona)
pred_antiviral <- make_pred_plot(model_antiviral, 10000,lab='Anti-Virals',dataframe=dat_antiviral_nona)
pred_bisp <- make_pred_plot(model_bisp, 1000000,lab='Anti-Bisphosphonates',dataframe=dat_bisp_nona)
pred_hormone <- make_pred_plot(model_hormone, 10000,lab='Hormone Therapy',dataframe=dat_hormone_nona)
pred_oral <- make_pred_plot(model_oral, 10000,lab='Oral Combined Estrogen and Progestin',dataframe=dat_oral_nona)
pred_transderm <- make_pred_plot(model_transderm, 10000,lab='Transdermal Estrogen',dataframe=dat_transderm_nona)
pred_vaginal <- make_pred_plot(model_vaginal, 600000,lab='Vaginal Estrogen Hormone Therapy',dataframe=dat_vaginal_nona)


total_plot_data<-rbind(pred_anti_inf,
#pred_antichol,
#pred_antiviral,
pred_bisp,
pred_hormone,
#pred_oral,
#pred_transderm,
pred_vaginal)


total_plot_data_sub1<-total_plot_data[total_plot_data$value<10000,]
total_plot_data_sub2<-total_plot_data[total_plot_data$value<1000,]
proper<-function(x) paste0(toupper(substr(x, 1, 1)), tolower(substring(x, 2)))

labels.drugs<-cast(total_plot_data_sub2,variable~.,value='value',max,an.rm=T)
total_plot_data_sub1$drug<-factor(proper(gsub('at.level\\(trait, 1\\)\\:as.factor\\(drug\\)','',total_plot_data_sub1$variable)))
names(labels.drugs)[2]<-'value'
labels.drugs<-merge(labels.drugs,total_plot_data_sub1,by=c('variable','value'))

total_plot_data_sub1<-total_plot_data_sub1[with(total_plot_data_sub1,order(group,pay)),]
library(lattice)
#pdf(file='Z:\\Pharma_Influence\\Results\\results_john\\figure_2.pdf',width=14,height=8)
#tiff(file='Z:\\Pharma_Influence\\Results\\results_john\\figure_2.tiff',width=14,height=8,units='in',res=800,compression='lzw')
xyplot(log10(value)~pay|group,groups=variable,data=total_plot_data_sub1,type='l',as.table=T,
  scales=list(x=list(tck=c(1,0),alt=1,relation='free'),col='black',
  y=list(tck=c(1,0),alt=1,at=c(0,1,2,3,4,5),lab=c(expression(10^0),expression(10^1),expression(10^2),expression(10^3),expression(10^4),expression(10^5)))),
  strip=strip.custom(bg='grey'),panel=panel.superpose,
  panel.groups=function(x,y,group.number,...){
    panel.xyplot(x,y,...)
    xt<-x[x==max(x)]
    yt<-y[x==max(x)]
    panel.text(xt, jitter(yt),labels= group.number,...)
  },

  xlab='Cumulative Pay ($)',ylab='Prescriptions',layout=c(2,2)
)


levels(total_plot_data_sub1$drug)[
xyplot(log10(value)~pay|group,groups=variable,data=total_plot_data_sub1,type='l',
  scales=list(x=list(tck=c(1,0),alt=1,relation='free'),y=list(tck=c(1,0),alt=1,rel='free')),
  #y=list(tck=c(1,0),alt=1,at=c(-20,-15,-10,-5,0,5),lab=c(expression(10^-20),expression(10^-15),expression(10^-10),expression(10^-5),expression(10^0),expression(10^5)))),
  strip=strip.custom(bg='grey'),xlim=c(0,100),
  xlab='Cumulative Pay ($)',ylab='Prescriptions',layout=c(2,2)
)

dev.off()


library(ggplot2)
library(scales) 
library(gridExtra)

p1 <- ggplot(total_plot_data_sub1[total_plot_data_sub1$group %in%'Anti-Bisphosphonates',], aes(x = pay, y = value,group=drug,color=drug)) + geom_line(size=1.2) +
     scale_y_log10(limits=c(1,10000),
              labels = trans_format("log10", math_format(10^.x))) +  
              xlab("Cumulative Payments ($)") +  ylab("Prescriptions") + labs(color='Anti-Bisphosphonates')  +
     theme_bw() + theme(legend.position = c(0.7, 0.3),legend.text = element_text(size=16, face="bold"),legend.title=element_text(face='bold',size=18), text = element_text(size=20))

p2 <- ggplot(total_plot_data_sub1[total_plot_data_sub1$group %in%'Anti-Infectives',], aes(x = pay, y = value,group=drug,color=drug)) + geom_line(size=1.2) +
     scale_y_log10(limits=c(1,10000),
              labels = trans_format("log10", math_format(10^.x))) +  
              xlab("Cumulative Payments ($)") +  ylab("Prescriptions") + labs(color='Anti-Infectives')  +
     theme_bw() + theme(legend.position = c(0.7, 0.3),legend.text = element_text(size=16, face="bold"),legend.title=element_text(face='bold',size=18), text = element_text(size=20))

p3 <- ggplot(total_plot_data_sub1[total_plot_data_sub1$group %in%'Hormone Therapy',], aes(x = pay, y = value,group=drug,color=drug)) + geom_line(size=1.2) +
     scale_y_log10(limits=c(1,10000),
              labels = trans_format("log10", math_format(10^.x))) +  
              xlab("Cumulative Payments ($)") +  ylab("Prescriptions") + labs(color='Hormone Therapy')  +
     theme_bw() + theme(legend.position = c(0.7, 0.3),legend.text = element_text(size=16, face="bold"),legend.title=element_text(face='bold',size=18), text = element_text(size=20))

p4 <- ggplot(total_plot_data_sub1[total_plot_data_sub1$group %in%'Vaginal Estrogen Hormone Therapy',], aes(x = pay, y = value,group=drug,color=drug)) + geom_line(size=1.2) +
     scale_y_log10(limits=c(1,10000),
              labels = trans_format("log10", math_format(10^.x))) +  
              xlab("Cumulative Payments ($)") +  ylab("Prescriptions") + labs(color='Vaginal Estrogen \n Hormone Therapy')  +
     theme_bw() + theme(legend.position = c(0.7, 0.3),legend.text = element_text(size=16, face="bold"),legend.title=element_text(face='bold',size=18), text = element_text(size=20))


pdf(file='Z:\\Pharma_Influence\\Results\\results_john\\figure_2.pdf',width=20,height=14)
print(grid.arrange(p1, p2, p3, p4,nrow = 2))
dev.off()

rbind(
data.frame(lab='Anti-Infectives',summary(model_anti_inf)$solutions) ,
data.frame(lab='Anticholinergics',summary(model_antichol)$solutions) ,
#data.frame(lab='Anti-Virals',summary(model_antiviral)$solutions) ,
data.frame(lab='Anti-Bisphosphonates',summary(model_bisp)$solutions) ,
data.frame(lab='Hormone Therapy',summary(model_hormone)$solutions) ,
data.frame(lab='Oral Combined Estrogen and Progestin',summary(model_oral)$solutions) ,
#data.frame(lab='Transdermal Estrogen',summary(model_transderm)$solutions) ,
data.frame(lab='Vaginal Estrogen Hormone Therapy',summary(model_vaginal)$solutions) 
)



prob.logit<-function(x){
  y<-exp(x)/(1+exp(x))
  return(y)
}

prob.logit((summary(model_anti_inf)$solutions)[2,1])
prob.logit((summary(model_antichol)$solutions)[2,1])
#prob.logit((summary(model_antiviral)$solutions)[2,1])
prob.logit((summary(model_bisp)$solutions)[2,1])
prob.logit((summary(model_hormone)$solutions)[2,1])
prob.logit((summary(model_oral)$solutions)[2,1])
#prob.logit((summary(model_transderm)$solutions)[2,1])
prob.logit((summary(model_vaginal)$solutions)[2,1])
