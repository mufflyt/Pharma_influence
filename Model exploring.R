# Created by: AF
# Date: 8/9/2019
# Purpose: model experiment for prescribing project


# Notes: some class data files do not have any payment numbers, all blank,
# and often when there is payment numbers the prescribing is blank

# Need the years, is this cumulative or just one year?

# if prescribing/payment info is missing does that mean it's zero or unknown? will make a big difference

# start with one file that has a less missing/0
class_bisph <- read.csv(file = 'Z:\\Pharma_Influence\\Results\\class_Bisphosphonates.csv', header = T)
names(class_bisph)
summary(class_bisph)


# histogram excluding the zeros for all the variables
myhists <- function(var){
  hist(subset(class_bisph[,var], class_bisph[,var] !=0), breaks = 100, main = var)
}
myhists(var = "Pay_atelvia")
myhists(var = "Pay_boniva") #all zeros or missings
myhists(var = "Pay_divigel")
myhists(var = "Pay_elestrin") #interesting two peaks
myhists(var = "Pay_fosamax") 
myhists(var = "Pay_prolia")
myhists(var = "Pay_risedronate")
myhists(var = "Pre_atelvia")
myhists(var = "Pre_boniva") 
myhists(var = "Pre_divigel")
myhists(var = "Pre_elestrin")
myhists(var = "Pre_fosamax")
myhists(var = "Pre_prolia") #peak at end
myhists(var = "Pre_risedronate")

# missing percents
colMeans(is.na(class_bisph))*100





# Play around with model 
# install.packages("MCMCglmm")
library(MCMCglmm)

# wide to long
dat <- reshape(class_bisph, direction = "long", idvar = "NPI", sep = '_', timevar = 'trait',
              varying = c("Pay_atelvia", "Pay_boniva", "Pay_divigel", "Pay_elestrin", "Pay_fosamax", "Pay_prolia", "Pay_risedronate",
                          "Pre_atelvia", "Pre_boniva", "Pre_divigel", "Pre_elestrin", "Pre_fosamax", "Pre_prolia", "Pre_risedronate"))

# Mutivariate error structures (i.e. using the term 'trait') are required for multinomial data with more than 2 categories,
# or zero-infalted/altered/hurdle models.


# don't have multiple years so not a multi response model yet, not sure if it will work
test <- MCMCglmm(data = dat, Pre~ 1 + drug + Pay, random = ~NPI, family = "ziPoisson")

test <- MCMCglmm(data = dat, Pre~ 1 + trait + Pay, random = ~us(trait):NPI, family = "ziPoisson")



