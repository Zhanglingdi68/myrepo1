source('Tables.R')
#fit a linear model to check residuals
png("/Users/zhouwenru/Documents/UCD Courses/BIOS6624/PRO1/Output/Des_diff_AGG_MENT.png",
    width = 1500, 
    height = 2000,
    res = 100, 
    units = "px")
par(mfrow = c(4,2))

#numeric
easy_p(data_final$diff_AGG_MENT,data_final$age_0)
easy_p(data_final$diff_AGG_MENT,data_final$BMI_0)

#categorical
easy_p(data_final$diff_AGG_MENT,data_final$group)
easy_p(data_final$diff_AGG_MENT,data_final$RACE_0)
easy_p(data_final$diff_AGG_MENT,data_final$income_0)
easy_p(data_final$diff_AGG_MENT,data_final$SMOKE_0)
easy_p(data_final$diff_AGG_MENT,data_final$ADH_2)
dev.off()

png("/Users/zhouwenru/Documents/UCD Courses/BIOS6624/PRO1/Output/Des_diff_AGG_PHYS.png",
    width = 1500, 
    height = 2000,
    res = 100, 
    units = "px")
par(mfrow = c(4,2))

#numeric
easy_p(data_final$diff_AGG_PHYS,data_final$age_0)
easy_p(data_final$diff_AGG_PHYS,data_final$BMI_0)

#categorical
easy_p(data_final$diff_AGG_PHYS,data_final$group)
easy_p(data_final$diff_AGG_PHYS,data_final$RACE_0)
easy_p(data_final$diff_AGG_PHYS,data_final$income_0)
easy_p(data_final$diff_AGG_PHYS,data_final$SMOKE_0)
easy_p(data_final$diff_AGG_PHYS,data_final$ADH_2)

dev.off()

png("/Users/zhouwenru/Documents/UCD Courses/BIOS6624/PRO1/Output/Des_diff_LEU3N.png",
    width = 1500, 
    height = 2000,
    res = 100, 
    units = "px")
par(mfrow = c(4,2))

#numeric
easy_p(data_final$diff_LEU3N,data_final$age_0)
easy_p(data_final$diff_LEU3N,data_final$BMI_0)

#categorical
easy_p(data_final$diff_LEU3N,data_final$group)
easy_p(data_final$diff_LEU3N,data_final$RACE_0)
easy_p(data_final$diff_LEU3N,data_final$income_0)
easy_p(data_final$diff_LEU3N,data_final$SMOKE_0)
easy_p(data_final$diff_LEU3N,data_final$ADH_2)

dev.off()



png("/Users/zhouwenru/Documents/UCD Courses/BIOS6624/PRO1/Output/Des_diff_VLOAD.png",
    width = 1500, 
    height = 2000,
    res = 100, 
    units = "px")
par(mfrow = c(4,2))

#numeric
easy_p(data_final$diff_VLOAD,data_final$group)
easy_p(data_final$diff_VLOAD,data_final$age_0)
easy_p(data_final$diff_VLOAD,data_final$BMI_0)

#categorical
easy_p(data_final$diff_VLOAD,data_final$RACE_0)
easy_p(data_final$diff_VLOAD,data_final$income_0)
easy_p(data_final$diff_VLOAD,data_final$SMOKE_0)
easy_p(data_final$diff_VLOAD,data_final$ADH_2)

dev.off()


png("/Users/zhouwenru/Documents/UCD Courses/BIOS6624/PRO1/Output/Des_diff_L_VLOAD.png",
    width = 1500, 
    height = 2000,
    res = 100, 
    units = "px")
par(mfrow = c(4,2))

#numeric
easy_p(data_final$diff_L_VLOAD,data_final$age_0)
easy_p(data_final$diff_L_VLOAD,data_final$BMI_0)

#categorical
easy_p(data_final$diff_L_VLOAD,data_final$group)
easy_p(data_final$diff_L_VLOAD,data_final$RACE_0)
easy_p(data_final$diff_L_VLOAD,data_final$income_0)
easy_p(data_final$diff_L_VLOAD,data_final$SMOKE_0)
easy_p(data_final$diff_L_VLOAD,data_final$ADH_2)

dev.off()