# Set base directory
setwd("/Users/ilanatrumble/Repositories/bios6624-trumblei/Project0/")
# Import data
library(gdata)
data0=read.csv("DataRaw/Project0_dental_data.csv")

# Continuous variables
# mean,median,25th and 75th quartiles,min,max
sex_table
race_table
smoker_table
# MISSING DATA
#### MISSING DATA ###
# Person who has age missing:
subset(data0,is.na(age))
missing1year=subset(data0, is.na(pd1year))
# Descriptive statistics for people with 1 year missing
#Continous variables, 1 year missing
summary(missing1year$age)
summary(missing1year$attachbase)
summary(missing1year$pdbase)
summary(missing1year$sites)
# Categorical variables
# Categorical variables, 1 year missing
trtgroup_table1<-table(missing1year$trtgroup,exclude=NULL)
sex_table1<-table(missing1year$sex,exclude=NULL)
race_table1<-table(missing1year$race,exclude=NULL)
 sex_table1
race_table1
smoker_table1
#### CREATE TABLE 1 ###
library(tableone)
#Create a variable list which we want in Table 1
listVars<-c("age","attachbase","attach1year","pdbase","pd1year","sites",
        write.csv(tab1Export,file="Reports/table1.csv")
        print(table1,quote=T,nospaces=T,showAllLevels = T)
            #### SCATTER PLOTS ###
            # CLEAN DATA : WHERE NOT MISSING AT 1 YEAR
            datac=subset(data0, is.na(pd1year)==F)
          y <- predict(lm(pd1year ~ age, data = datac_age))
            x <- datac_age$age
            lines(y = y, x = x, col = "navy"
                  
                  # Histograms of pd1year and attach1year
            #### HISTOGRAMS to assess normality ###
            par(mfrow=c(1, 2)) 
            x=datac$attach1year
            h<-hist(x, breaks=10, xlab="AL at 1 year", 
                  
          yfit <- yfit*diff(h$mids[1:2])*length(x)
                    lines(xfit, yfit, col="blue", lwd=2)
                    # Box plots 
                    #### BOX PLOTS ###
                    boxplot(attach1year~trtgroup,data=data0,main="Attachment loss at 1 year",
                            xlab="Treatment Group",ylab="Attachment loss")
                    boxplot(pd1year~trtgroup,data=data0,main="PD at 1 year",
                            xlab="Treatment Group",ylab="Pocket Depth")
                    
                    # Fit linear model
                    #### LINEAR MODELS ###
                    dataclean<-na.omit(data0)
                    full1<-lm(attach1year ~ factor(trtgroup) + factor(sex) + factor(race) 
                                age + factor(smoker) + sites + attachbase, data=dataclean)
                    # Reduced linear model
                    # Reduced linear model for AL
                    step(full1,data=dataclean,direction="backward")
                    red1<-lm(attach1year ~ factor(trtgroup) + factor(smoker)
                               attachbase, data = dataclean)
                    summary(red1)
                    # AIC: most negative is best
                    # Reminder: for AIC, most negative is best
                    full2<-lm(pd1year ~ factor(trtgroup) + factor(sex) + factor(race)
                                age + factor(smoker) + sites + pdbase, data=dataclean)
                    
                    # Reduced linear model
                    # Reduced linear model for PD
 step(full2,data=dataclean,direction="backward")
 red2<-lm(pd1year ~ factor(trtgroup) + factor(sex) + pdbase, data = dataclean)
summary(red2)
                    