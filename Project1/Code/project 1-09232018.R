#project 1:How does the hard drug affect the HIV patients?
#author:Lingdi Zhang 
library(readr)
library(rjags)
library(mcmcse)

library(data.table) #read in the data online
library(mice) # for multiple imputation
library(VIM) # for visualizing missing data patterns
library(randomForest) # model of random forest
library(dplyr)  #for mutate and remove variables
library(rms)
library(caret)
library(gbm) # for model of GBM
library(car) 
library(pROC)

#read the raw data into R
HIV <- read_csv("C:/Repositories/bios6624-zhanglingdi68/Project1/DataRaw/hiv_6624_final.csv")

names(HIV)
summary(HIV)
head(HIV)
dim(HIV)
#HIV has 3935 subjects and 34 variables

## test HIV Serostatus, and all of the patients are HIV positive
paste("Are all participants in this study HIV Serostatus positive?")
paste(dim(HIV)[1] == sum(HIV$hivpos))
dim(HIV)
sum(HIV$hivpos)




#subset the years=2 for the dataset
dat1<- HIV[ which(HIV$years==2|HIV$years==0), ]
summary(dat1)
dim(dat1)
#dat1 has 1221 subjects and 34 variables

dat1$indices <- duplicated(dat1[ ,2]) | duplicated(dat1[ ,2], fromLast = TRUE)
dat2 <- subset(dat1, indices == "TRUE")

# describe the dataset
summary(dat2)
str(dat2)
dim(dat2)
head(dat2)
names(dat2)

#dat1 has 1012 subjects which has duplicated ID number and 35 variables.

#create categorical variables

dat2$hard_drugs<-as.factor(dat2$hard_drugs)
levels(dat2$hard_drugs) <- c("No", "Yes")

dat2$ART<-as.factor(dat2$ART)
levels(dat2$ART) <- c("No", "Yes")

dat2$HASHV <- ifelse(dat2$HASHV == 1, "No", dat2$HASHV)
dat2$HASHV <- ifelse(dat2$HASHV == 2, "Yes", dat2$HASHV)

dat2$CESD <- ifelse(dat2$CESD >= 16, "No", dat2$CESD)
dat2$CESD <- ifelse(dat2$CESD <= 16, "Yes", dat2$CESD)


dat2$income<-as.factor(dat2$income)
levels(dat2$income) <- c("Less than $10,000",
                             "10,000-19,999",
                             "20,000-29,999",
                             "30,000-39,999",
                             "40,000-49,999",
                             "50,000-59,999",
                             "60,000 or more",
                             "Do not wish to answer")



dat2$SMOKE<-as.factor(dat2$SMOKE)
levels(dat2$SMOKE) <- c("Never", "Former", "Current")

dat2$HBP <- factor(dat2$HBP)
levels(dat2$HBP) <- c("No", "Yes", "No, based on trajectory",
                          "Yes, based on trajectory", "Insufficient data")

dat2$DIAB <- factor(dat2$DIAB)
levels(dat2$DIAB) <- c("No", "Yes", "No, based on trajectory",
                           "Insufficient data")

dat2$LIV34 <- factor(dat2$LIV34)
levels(dat2$LIV34) <- c("No", "Yes", "Insufficient data")

dat2$KID <- factor(dat2$KID)
levels(dat2$KID) <- c("No", "Yes", "No, based on trajectory",
                          "Insufficient data")

dat2$FRP <- factor(dat2$FRP)
levels(dat2$FRP) <- c("No", "Yes", "Insufficient data")

dat2$FP <- factor(dat2$FP)
levels(dat2$FP) <- c("No", "Yes", "Insufficient data")

dat2$DYSLIP <- factor(dat2$DYSLIP)
levels(dat2$DYSLIP) <- c("No", "Yes", "No, based on trajectory",
                             "Yes, based on trajectory", "Insufficient data")

dat2$DKGRP <- factor(dat2$DKGRP)
levels(dat2$DKGRP) <- c("None", "1 to 3 drinks/week", 
                            "4 to 13 drinks/week", "More than 13 drinks/week")

dat2$HEROPIATE <- factor(dat2$HEROPIATE)
levels(dat2$HEROPIATE) <- c("Not specified", "No", "Yes")

dat2$IDU <- factor(dat2$IDU)
levels(dat2$IDU) <- c("No", "Yes")

dat2$ADH <- factor(dat2$ADH)
levels(dat2$ADH) <- c("100%", "95-99%", "75-94%", "<75%")

dat2$RACE <- factor(dat2$RACE)
levels(dat2$RACE) <- c("White, non-Hispanic",
                           "White, Hispanic",
                           "Black, non-Hispanic",
                           "Black, Hispanic",
                           "American Indian or Alaskan Native",
                           "Asian or Pacific Islander",
                           "Other",
                           "Other Hispanic")

dat2$EDUCBAS <- factor(dat2$EDUCBAS)
levels(dat2$EDUCBAS) <- c("8th grade or less",
                              "9, 10, or 11th grade",
                              "12th grade",
                              "At least one year college but no degree",
                              "Four years college / got degree",
                              "Some graduate work",
                              "Post-graduate degree")

dat2$hivpos <- factor(dat2$hivpos)
levels(dat2$hivpos) <- c("Yes")


dat2$everART <- factor(dat2$everART)
levels(dat2$everART) <- c("No", "Yes")

# summarize the categorical variables

cate_vs <- c("HASHV", "HASHF", "income", "HBP", "DIAB", "LIV34", "KID", 
             "FRP", "FP", "DYSLIP", "SMOKE", "DKGRP","HEROPIATE", "IDU", 
             "ADH", "RACE", "EDUCBAS", "CESD", "hivpos", "ART", "everART", "hard_drugs")
paste(cate_vs)

# summarize the continuous variables

conti_vs <- c("BMI", "TCHOL", "TRIG", "LDL", "age")  
paste(conti_vs)
paste("There are",length(cate_vs) + length(conti_vs), "covariates")

#there are 27 covariates

dat2 <- dat2[ , c(-1, -35)]


# get data into wide format
dat3 <- reshape(dat2, idvar = "newid", 
                    timevar = "years", direction = "wide")

dat2<- select(dat2, AGG_PHYS, AGG_MENT, VLOAD, LEU3N,
                   income, BMI, age,
                   RACE, ADH, SMOKE, hard_drugs)

dat4 <- select(dat3, AGG_PHYS.0, AGG_PHYS.2, AGG_MENT.0, AGG_MENT.2,
                   VLOAD.0, VLOAD.2, LEU3N.0, LEU3N.2, income.0, BMI.0, age.0,
                   RACE.0, ADH.2, SMOKE.0, hard_drugs.0, hard_drugs.2)

# seem to have a BMI of ~ 514, not physiologically possible (exclude)
# also want to exclude insufficient/improbable data ( = 999 or -1)

dat3 <- filter(dat2, !(BMI > 500 | BMI == -1))
dat4 <- filter(dat3, !(BMI.0 > 500 | BMI.0 == -1))

write.csv(dat3, "C:/Repositories/bios6624-zhanglingdi68/Project1/DataProcessed/hivclean_long.csv")
write.csv(dat4, "C:/Repositories/bios6624-zhanglingdi68/Project1/DataProcessed/hivclean_wide.csv")


#read in data set
setwd("C:/Repositories/bios6624-zhanglingdi68/Project1/DataProcessed")
hivclean.long <- read.csv("hivclean_long.csv")
HIVcw <- read.csv("hivclean_wide.csv")



#check the distribution of the outcomes
par(mfrow = c(2,4))
hist(HIVcw$VLOAD.2)

qqnorm(HIVcw$VLOAD.2,
       ylab="Sample Quantiles for vrial load")
qqline(HIVcw$VLOAD.2, 
       col="red")

hist(HIVcw$LEU3N.2)

qqnorm(HIVcw$LEU3N.2,
       ylab="Sample Quantiles for LEU3N.2")
qqline(HIVcw$LEU3N.2, 
       col="red")

hist(HIVcw$AGG_MENT.2)
qqnorm(HIVcw$AGG_MENT.2,
       ylab="Sample Quantiles for AGG_MENT.2")
qqline(hHIVcw$AGG_MENT.2, 
       col="red")

hist(HIVcw$AGG_PHYS.2)
qqnorm(HIVcw$AGG_PHYS.2,
       ylab="Sample Quantiles for AGG_PHYS.2")
qqline(HIVcw$AGG_PHYS.2, 
       col="red")


#skewed, log tranform the outcomes

HIVCLEAN<-mutate(HIVcw, lVLOAD.2 = log(VLOAD.2),lVLOAD.0 = log(VLOAD.0))


par(mfrow = c(2,2))

hist(HIVCLEAN$lVLOAD.2)

qqnorm(HIVCLEAN$lVLOAD.2,
       ylab="Sample Quantiles for vrial load")
qqline(HIVCLEAN$lVLOAD.2, col="red")




hist(HIVCLEAN$lVLOAD.0)

qqnorm(HIVCLEAN$lVLOAD.0,
       ylab="Sample Quantiles for vrial load")
qqline(HIVCLEAN$lVLOAD.0, col="red")

  
summary(HIVCLEAN)
      
# export the clean dataset for use in other analysis programs
write_csv(HIVCLEAN, "C:/Repositories/bios6624-zhanglingdi68/Project1/DataProcessed/HIVCLEAN.csv")


HIVclean <- read_csv("C:/Repositories/bios6624-zhanglingdi68/Project1/DataProcessed/HIVclean.csv")

#check the missing data
# check the missing data pattern and probability
md.pattern(HIVCLEAN)
#AGG_MENT and AGG_PHYS have 6 missing data;VLOAD and LEU3N have 19 missing data, 2.2% missing, which can be delegated.

aggr(HIVCLEAN, numbers=TRUE, sortVars=TRUE) # from VIM package

#complicated coding for missing pattern
mice_plot <- aggr(dat, col=c('navyblue','yellow'),
                  numbers=TRUE, sortVars=TRUE,
                  labels=names(HIVclean), cex.axis=.7,
                  gap=3, ylab=c("Missing data","Pattern"))

#delete the missing data, since there are only small portion for the missing data

newdata <- na.omit(HIVCLEAN)


#check the colinearity or correlation
data <- data.frame(data=newdata) 
cor(newdata,use ="pairwise.complete.obs")

# check the correlations of the numerical variables,
library (PerformanceAnalytics)
corr_data <- newdata[, c(2,3,4,5,6,7,8,9,11,12,18,19)]
chart.Correlation(corr_data, histogram=TRUE, pch=19)

library (PerformanceAnalytics)
corr_data <- HIVCLEAN[, c(2,3,4,5,6,7,8,9,11,12,18,19)]
chart.Correlation(corr_data, histogram=TRUE, pch=19)


#new variables for differences of the year 2 and year 0.
newdata$diff_MENT=newdata$AGG_MENT.2-newdata$AGG_MENT.0
newdata$diff_PHYS=newdata$AGG_PHYS.2-newdata$AGG_PHYS.0
newdata$diff_LEU3N=newdata$LEU3N.2-newdata$LEU3N.0
newdata$diff_VLOAD=newdata$VLOAD.2-newdata$VLOAD.0
newdata$diff_lVLOAD=newdata$lVLOAD.2-newdata$lVLOAD.0

label(newdata$diff_MENT)="Change in SF36 MCS score"
label(newdata$diff_PHYS)="Change in SF36 PCS score"
label(newdata$diff_LEU3N)="Change in # of CD4 positive cells"
label(newdata$diff_VLOAD)="Change in Standardized viral load (copies/ml)"
label(newdata$diff_lVLOAD)="Change in log Standardized viral load (copies/ml)"

var2<-c("diff_MENT","diff_PHYS","diff_LEU3N","diff_VLOAD","diff_lVLOAD")

#create the dataset with the changes of the outcomes
write_csv(newdata, "C:/Repositories/bios6624-zhanglingdi68/Project1/DataProcessed/DIFF.csv")


#Part two: model building

DIFF <- read_csv("C:/Repositories/bios6624-zhanglingdi68/Project1/DataProcessed/DIFF.csv")


#Model selection
#Mixed linear model

hdrugmodel <-glm(diff_MENT~hard_drugs.0, data=DIFF)

sink("C:/Repositories/bios6624-zhanglingdi68/Project1/Output/hdrugmodel.txt")
print(hdrugmodel)
summary.glm(hdrugmodel)
confint(hdrugmodel )
sink()
#AIC: 3604


hdrugmodelF <-glm(diff_MENT~hard_drugs.0+income.0+BMI.0+age.0+RACE.0+ADH.2+SMOKE.0, data=DIFF)
sink("C:/Repositories/bios6624-zhanglingdi68/Project1/Output/hdrugmodelF.txt")
print(hdrugmodelF)
summary.glm(hdrugmodelF)
confint(hdrugmodelF)
sink()
#AIC: 3625

# LRT test the significance of the more variables
anova(hdrugmodel, hdrugmodelF, test="LRT") 
#p-value: 0.6418

# diff_PHYS
PHYSmodel <-glm(diff_PHYS~hard_drugs.0, data=DIFF)

sink("C:/Repositories/bios6624-zhanglingdi68/Project1/Output/PHYSmodel.txt")
print(PHYSmodel)
summary.glm(PHYSmodel)
confint(PHYSmodel)
sink()
#AIC: 3243,p-value: 0.01407 *


PHYSmodelF <-glm(diff_PHYS~hard_drugs.0+income.0+BMI.0+age.0+RACE.0+ADH.2+SMOKE.0, data=DIFF)
sink("C:/Repositories/bios6624-zhanglingdi68/Project1/Output/PHYSmodelF.txt")
print(PHYSmodelF)
summary.glm(PHYSmodelF)
confint(PHYSmodelF)
sink()
#AIC: 3254,p-value:0.0324 *

# LRT test the significance of the more variables
anova(PHYSmodel, PHYSmodelF, test="LRT") 
#p-value: 0.114

#diff_lVLOAD

lVLOADmodel <-glm(diff_lVLOAD~hard_drugs.0, data=DIFF)

sink("C:/Repositories/bios6624-zhanglingdi68/Project1/Output/lVLOADmodel.txt")
print(lVLOADmodel)
summary.glm(lVLOADmodel)
confint(lVLOADmodel)
sink()
#AIC:2270, p-value: 0.961  


lVLOADmodelF<-glm(diff_lVLOAD~hard_drugs.0+income.0+BMI.0+age.0+RACE.0+ADH.2+SMOKE.0, data=DIFF)
sink("C:/Repositories/bios6624-zhanglingdi68/Project1/Output/lVLOADmodelF.txt")
print(lVLOADmodelF)
summary.glm(lVLOADmodelF)
confint(lVLOADmodelF)
sink()
#AIC:  2255;  p-value: 0.47045

# LRT test the significance of the more variables
anova(lVLOADmodel,lVLOADmodelF, test="LRT") 

#p-value: 3.345e-05 ***



#diff_LEU3N

LEU3Nmodel <-glm(diff_LEU3N~hard_drugs.0, data=DIFF)

sink("C:/Repositories/bios6624-zhanglingdi68/Project1/Output/LEU3Nmodel.txt")
print(LEU3Nmodel)
summary.glm(LEU3Nmodel)
confint(LEU3Nmodel)
sink()
#AIC: 6073, p-value: 5.43e-08 ***


LEU3NmodelF<-glm(diff_LEU3N~hard_drugs.0+income.0+BMI.0+age.0+RACE.0+ADH.2+SMOKE.0, data=DIFF)
sink("C:/Repositories/bios6624-zhanglingdi68/Project1/Output/LEU3NmodelF.txt")
print(LEU3NmodelF)
summary.glm(LEU3NmodelF)
confint(LEU3NmodelF)
sink()
#AIC: 6079;  p-value: 1.23e-07 ***

# LRT test the significance of the more variables
anova(LEU3Nmodel, LEU3NmodelF, test="LRT") 

#p-value: 0.03291 *














iter <- 10000 # number of draws

## This script will work for any linear model we wish to fit with JAGS
## In addition, we can use this script to generate models for the univariable interaction hypothesis

y <- c(dat1$LEU3N)
# change the formula to test the different models
X <- model.matrix(~ cavol + wt + age + bph + cappen + svi + grade6 + grade7, data = dat1) 
N <- nrow(X)
p <- ncol(X)

## Hyperparameters for prior distributions

a <- 2 # invgamma shape
b <- 1 # invgamma rate
m <- rep(0, p) # mvnorm mean
R <- matrix(0, p, p) # mvnorm covariance
diag(R) <- 0.001 # note that JAGS uses dispersion matrix (scalars)

# create data list to pass to JAGS
jags_dat <- list(y = y, X = X, N = N, p = p,
                 a = a, b = b, m = m, R = R)

# Initialize the sampler (it's a Gibbs sampler so no need for several adaptations)

mod <- jags.model("C:/Repositories/Bios6624ClassExamples/PSAExample/Code/jags/linMod.jags", data = jags_dat, n.adapt = 1000, n.chains = 2)

# Sample observations from the posterior distributions

samples <- coda.samples(mod, variable.names = c("beta", "sigma2"), n.iter = iter) # generates matrix of draws from posterior dist.
samples_dic <- dic.samples(mod, n.iter = iter, type = "pD") # determine DIC of model fit (requires two chains)

## Diagnostics

# determine if the support of the parameter space is explored well (more adaptations)
plot(samples)

# determine amount of autocorrelation between draws (thinning parameter)
par(mfrow = c(ceiling(ncol(as.matrix(samples))/2), 2), mar = rep(1, 4))
apply(as.matrix(samples), 2, acf) # again the Gibb's sampler should have almost no autocorrelation

## HPDI function

hpd <- function(x, alpha = 0.05){
  
  n = length(x)
  m = round(n * alpha)
  x = sort(x)
  y = x[(n - m + 1):n] - x[1:m]
  z = min(y)
  k = which(y == z)[1]
  c(x[k], x[n - m + k])
  
}


## Generate output table when the final model has been selected

draws <- as.matrix(samples)

out_mat <- matrix("", nrow = ncol(draws), ncol = 5)
out_mat[,1] <- c(colnames(X), "sigma2") # names
out_mat[,2] <- apply(draws, 2, function(x) round(mcse(x)$est, 3)) # batch mean
out_mat[,3] <- apply(draws, 2, function(x) round(mcse(x)$se, 3)) # MCSE
out_mat[,4] <- round(apply(draws, 2, sd), 3) # Std. Dev.
out_mat[,5] <- apply(draws, 2, function(x)
  paste("(", round(hpd(x)[1], 3), ", ", round(hpd(x)[2], 3), ")", sep = "")) # HPDI

colnames(out_mat) <- c("Parameter", "Estimate", "MCSE", "Std. Dev.", "95% HPDI")

write_csv(as.data.frame(out_mat), "C:/Repositories/Bios6624Class/PSAExample/Output/bayes-est.csv")





