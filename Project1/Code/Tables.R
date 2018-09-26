source('DataClean.R')


#Combine unknowns as missing
#BMI missing and strange value
#strange, Improbable value, missing
#see values
data_wide$BMI_0[data_wide$BMI_0>500]
#see summary
unknown_bmi<-c(length(data_wide$BMI_0[data_wide$BMI_0>500]),length(data_wide$BMI_0[data_wide$BMI_0==-1]),length(data_wide$BMI_0[is.na(data_wide$BMI_0)==TRUE]))
unknown_bmi_count<-table(data_wide$BMI_0>500|data_wide$BMI_0==-1|is.na(data_wide$BMI_0)==TRUE)[2]
unknown_bmi_rate<-table(data_wide$BMI_0>500|data_wide$BMI_0==-1|is.na(data_wide$BMI_0)==TRUE)[2]/(table(data_wide$BMI_0>500|data_wide$BMI_0==-1|is.na(data_wide$BMI_0)==TRUE)[1]+table(data_wide$BMI_0>500|data_wide$BMI_0==-1|is.na(data_wide$BMI_0)==TRUE)[2])
data_wide$BMI_0[data_wide$BMI_0>500|data_wide$BMI_0==-1]<-NA
#min(data_wide$BMI_0,na.rm=TRUE)
#max(data_wide$BMI_0,na.rm=TRUE)

#income
unknown_income<-c(length(data_wide$income_0[data_wide$income_0=="Do not wish to answer"]),length(data_wide$income_0[is.na(data_wide$income_0)==TRUE]))
unknown_income_count<-table(data_wide$income_0=="Do not wish to answer"|is.na(data_wide$income_0)==TRUE)[2]
unknown_income_rate<-table(data_wide$income_0=="Do not wish to answer"|is.na(data_wide$income_0)==TRUE)[2]/(table(data_wide$income_0=="Do not wish to answer"|is.na(data_wide$income_0)==TRUE)[1]+table(data_wide$income_0=="Do not wish to answer"|is.na(data_wide$income_0)==TRUE)[2])
data_wide$income_0[data_wide$income_0=="Do not wish to answer"]<-NA

#table 0: Missing at baseline
var0<-c("AGG_MENT_0","AGG_PHYS_0","LEU3N_0","VLOAD_0","group","age_0","BMI_0","RACE_0","income_0","SMOKE_0")
label(data_wide$years_0)<-c("Baseline")
levels(data_wide$years_0)<-c("Baseline")
con_missing0<-missing_table(data_wide,var0,data_wide$years_0)

#table 0: Missing at 2
var2<-c("AGG_MENT_2","AGG_PHYS_2","LEU3N_2","VLOAD_2","ADH_2")
label(data_wide$years_2)<-c("Vist 2")
levels(data_wide$years_2)<-c("Visit 2")
con_missing2<-missing_table(data_wide,var2,data_wide$years_2)

#Those who have BMI missing, insufficient and improbable should not be included in this model
data_final<-subset(data_wide,is.na(data_wide$BMI)==FALSE&is.na(data_wide$income_0)==FALSE)

#table 1 at baseline and yr2
var1_0<-c("AGG_MENT_0","AGG_PHYS_0","LEU3N_0","VLOAD_0","age_0","BMI_0","RACE_0","income_0","SMOKE_0")
con_final_0<-final_table(data_final,var1_0,data_final$group,margin=2,single=F,ron=2)

var1_2<-c("AGG_MENT_2","AGG_PHYS_2","LEU3N_2","VLOAD_2","ADH_2")
con_final_2<-final_table(data_final,var1_2,data_final$group,margin=2,single=F,ron=2)

#table 2: difference from baseline to V2
data_final$diff_AGG_MENT=data_final$AGG_MENT_2-data_final$AGG_MENT_0
data_final$diff_AGG_PHYS=data_final$AGG_PHYS_2-data_final$AGG_PHYS_0
data_final$diff_LEU3N=data_final$LEU3N_2-data_final$LEU3N_0
data_final$diff_VLOAD=data_final$VLOAD_2-data_final$VLOAD_0
data_final$diff_L_VLOAD=log(data_final$VLOAD_2)-log(data_final$VLOAD_0)

label(data_final$diff_AGG_MENT)="Change in SF36 MCS score"
label(data_final$diff_AGG_PHYS)="Change in SF36 PCS score"
label(data_final$diff_LEU3N)="Change in # of CD4 positive cells"
label(data_final$diff_VLOAD)="Change in Standardized viral load (copies/ml)"
label(data_final$diff_L_VLOAD)="Change in log Standardized viral load (copies/ml)"

var2<-c("diff_AGG_MENT","diff_AGG_PHYS","diff_LEU3N","diff_VLOAD","diff_L_VLOAD")

con_diff<-final_table(data_final,var2,data_final$group,margin=2,single=F,ron=2)

