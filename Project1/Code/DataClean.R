library(Hmisc)
library(knitr)
library(devtools)
library(Table1)
#To get Table1, please visit 
#https://rdrr.io/github/palmercl/Table1/

#install_github("palmercl/Table1")

#But I made some revise to Table 1, therefore use Macro.R to mask it.
setwd("/Users/zhouwenru/Documents/UCD Courses/BIOS6624/PRO1/Program/")
source('Macro.R')

#Reading data
data<-read.csv("/Users/zhouwenru/Documents/UCD Courses/BIOS6624/PRO1/Data processed/hiv_6624_final.csv",header=T,na.strings = "NA")

#Make variables, some income and race has low number
data$income_ori<-data$income
data$income[data$income_ori==1]=1
data$income[data$income_ori==2]=2
data$income[data$income_ori==3]=3
data$income[data$income_ori==4]=4
data$income[data$income_ori==5|data$income_ori==6|data$income_ori==7]=5
data$income[data$income_ori==9]=9#unknown

data$RACE_ori<-data$RACE
table(data$income_ori)
data$RACE[data$RACE_ori==1]=1#white
data$RACE[data$RACE_ori==3]=3#black
#data$RACE[data$RACE_ori==5]=5#native
#data$RACE[data$RACE_ori==6]=6#asian
data$RACE[data$RACE_ori==7]=7#others
data$RACE[data$RACE_ori==2|data$RACE_ori==4|data$RACE_ori==8]=8#hispanic


#factorize;
data$ART<-as.factor(data$ART)
data$income<-as.factor(data$income)
data$SMOKE<-as.factor(data$SMOKE)
data$RACE<-as.factor(data$RACE)
data$hard_drugs<-as.factor(data$hard_drugs)
data$ADH<-as.factor(data$ADH)


#set lables;
label(data$newid)="Subject ID Number"
label(data$AGG_MENT)='SF36 MCS score'
label(data$AGG_PHYS)='SF36 PCS score'
label(data$years)="years since iniating ART"
label(data$hard_drugs)="Hard drug use since last visit"
label(data$LEU3N)="# of CD4 positive cells (helpers)"
label(data$ADH)="Adherence to meds taken since last visit"
label(data$ART)="Take ART at visit"
label(data$SMOKE)="Smoke"
label(data$income)="income"
label(data$BMI)="BMI(kg/m2)"
label(data$VLOAD)="Standardized viral load (copies/ml)"
label(data$age)="Age, years"
label(data$RACE)="Race"

#levels;
levels(data$ART)=c("No ART",'ART')
levels(data$ADH)=c("100%",'95-99%','75-94%','<75%')
levels(data$SMOKE)=c("Never smoked",'Former smoker','Current smoker')

levels(data$income)=c('Less than $10,000','10,000-19,999','20,000-29,999','30,000-39,999','40,000 or more','Do not wish to answer')
table(data$income_ori,data$income)

levels(data$RACE)=c('White, non-Hispanic','Black, non-Hispanic','Other','Hispanic')
table(data$RACE_ori,data$RACE)


#create a new subset with only year 0 and 2
data_easy<-subset(data[,c("newid","AGG_MENT","AGG_PHYS","LEU3N","VLOAD","ART","years","hard_drugs","age","BMI","ADH","RACE","income","SMOKE")],data$years==0|data$years==2)


#Check mising hard drug status at baseline
missing_hard_drugs_count<-dim(subset(data_easy,data_easy$years==0&is.na(data_easy$hard_drugs)==TRUE))[1]
missing_hard_drugs_rate<-dim(subset(data_easy,data_easy$years==0&is.na(data_easy$hard_drugs)==TRUE))[1]/dim(subset(data_easy,data_easy$years==0))[1]



#Create baseline data with hard drug status
data_base<-subset(data_easy,data_easy$years==0&is.na(data_easy$hard_drugs)==FALSE)
colnames(data_base)<-c("newid",paste(colnames(data_base)[2:length(colnames(data_base))],"_0",sep=""))
data_base$group[data_base$hard_drugs_0==1]=1
data_base$group[data_base$hard_drugs_0==0]=0
data_base$group<-as.factor(data_base$group)
levels(data_base$group)<-c("No hard drug use at baseline","Hard drug use at baseline")
table(data_base$group)

#factor years to compare as two groups
data_base$years_0<-as.factor(data_base$years_0)
label(data_base$group)="hard drug at baseline"

#Create year two data
data_three<-merge(data_base[,c("newid","group")],subset(data_easy,data_easy$years==2),by="newid")
colnames(data_three)<-c("newid","group",paste(colnames(data_three)[3:length(colnames(data_three))],"_2",sep=""))
table(data_three$group)

#factor years to compare as two groups
data_three$years_2<-as.factor(data_three$years_2)
label(data_three$hard_drugs)="Hard drug use since last visit"


#Make wide table for Baseline and Two years after
data_wide<-merge(data_base[,c("newid","AGG_MENT_0","AGG_PHYS_0","LEU3N_0","VLOAD_0","years_0","age_0","BMI_0","RACE_0","income_0","SMOKE_0","group")],data_three[,c("newid","AGG_MENT_2","AGG_PHYS_2","LEU3N_2","VLOAD_2","years_2","ADH_2")],by="newid")
label(data_wide$group)="hard drug at baseline"
