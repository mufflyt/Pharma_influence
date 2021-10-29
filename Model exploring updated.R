# Updated data where missings were changed to 0s





class_bisph <- read.csv(file = 'Z:\\Pharma_Influence\\Data\\Updated_data_10_9_2019\\class_Bisphosphonates.csv', header = T)
names(class_bisph)
summary(class_bisph)

class_bisph <- class_bisph[class_bisph$NPI != 0,]


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

# no more missing
colMeans(is.na(class_bisph))*100


library(MCMCglmm)

# wide to long
dat <- reshape(class_bisph, direction = "long", idvar = "NPI", sep = '_', timevar = 'drug',
               varying = c("Pay_atelvia", "Pay_boniva", "Pay_divigel", "Pay_elestrin", "Pay_fosamax", "Pay_prolia", "Pay_risedronate",
                           "Pre_atelvia", "Pre_boniva", "Pre_divigel", "Pre_elestrin", "Pre_fosamax", "Pre_prolia", "Pre_risedronate"))
dat$NPI2 <- as.factor(dat$NPI)


# Mutivariate error structures (i.e. using the term 'trait') are required for multinomial data with more than 2 categories,
# or zero-infalted/altered/hurdle models.


# don't have multiple years so not a multi response model yet?

# test <- MCMCglmm(data = dat, Pre~ 1 + drug + Pay, random = ~NPI, family = "zipoisson")

test <- MCMCglmm(data = dat, Pre~ 1 + drug + Pay, random = ~NPI2, rcov = ~us(trait):units, family = "zipoisson",
                burnin = 1000, nitt = 10000, thin = 20)
summary(test)
plot(test$VCV)

test2 <- MCMCglmm(data = dat, Pre~ 1 + drug + Pay, random = ~NPI2, rcov = ~us(trait):units, family = "zipoisson",
                 burnin = 4000, nitt = 10000, thin = 20)
summary(test2)
plot(test2$VCV)





