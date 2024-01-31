## load libraries into R
library(TwoSampleMR)
library(data.table)

## set our working directory
setwd("Y:/UKB/R Scripts and Analyses") 


##ICV

## Univariable MR - reading in our own exposure data frames
exp1 <- read.table('BMI_age_10_nov2.uvinput', header=T)
exp2 <- read.table('BMI_adult_nov2.uvinput', header=T)

names(exp1) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')
names(exp2) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')

## Multivariable MR - reading in our own exposure data frames
exp7 <- read.table('BMI_age_10_nov2.mvinput', header=T)
exp8 <- read.table('BMI_adult_nov2.mvinput', header=T)

names(exp7) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')
names(exp8) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')

## add exposure name
exp1$Phenotype <- 'age 10'
exp2$Phenotype <- 'adult'
exp7$Phenotype <- 'age 10'
exp8$Phenotype <- 'adult'

## formats datasets for TwoSampleMR
exp_dat1 <- format_data(exp1, type = "exposure")
exp_dat2 <- format_data(exp2, type = "exposure")
exp_dat7 <- format_data(exp7, type="exposure")
exp_dat8 <- format_data(exp8, type="exposure")
exp_dat <- rbind(exp_dat7, exp_dat8)

all_SNPs <- unique(c(exp_dat1$SNP,exp_dat2$SNP,exp_dat$SNP))

##new way to import outcome data due to funny bug in code in previous iteration
child_icv <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/child_icv.txt")
adult_icv <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/adult_icv.txt")
comb_icv <- rbind(child_icv, adult_icv)
names(comb_icv) <- c('chr','pos','SNP','other_allele','effect_allele','a1','eaf','test','obs','beta','se','l95','u95','tstat','pval','error')
comb_icv$Phenotype <- "icv"
out_dat <- format_data(comb_icv, type="outcome")

tryCatch(dat1 <- harmonise_data(
  exposure_dat = exp_dat1, 
  outcome_dat = out_dat, action=1
), error=function(e) NULL)

tryCatch(dat2 <- harmonise_data(
  exposure_dat = exp_dat2, 
  outcome_dat = out_dat, action=1
), error=function(e) NULL)

tryCatch(res1 <- mr(dat1, method_list = "mr_ivw"), error=function(e) NULL)
tryCatch(res2 <- mr(dat2, method_list = "mr_ivw"), error=function(e) NULL)


## MVMR
##tryCatch(out_dat_mvmr <- extract_outcome_data(exp_dat$SNP, try), error=function(e) NULL)

tryCatch(mvdat_mvmr <- mv_harmonise_data(
  exposure_dat = exp_dat,
  outcome_dat = out_dat, harmonise_strictness = 1
), error=function(e) NULL)


tryCatch(res_mv <- mv_multiple(mvdat_mvmr), error=function(e) NULL)

## write results to a file
res_mv1 <- res_mv$result[order(res_mv$result$exposure),]

tryCatch(all_res_early <- cbind(res1[1,-c(1,4,5)], res_mv1[2,]), error=function(e) NULL)
tryCatch(all_res_adult <- cbind(res2[1,-c(1,4,5)], res_mv1[1,]), error=function(e) NULL)
tryCatch(write.table(all_res_early, "Y:/UKB/R Scripts and Analyses/Output Files/icv_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)
tryCatch(write.table(all_res_adult, "Y:/UKB/R Scripts and Analyses/Output Files/icv_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)

##Clear everything
rm(list = ls())


##TOTAL VOLUME

## Univariable MR - reading in our own exposure data frames
exp1 <- read.table('BMI_age_10_nov2.uvinput', header=T)
exp2 <- read.table('BMI_adult_nov2.uvinput', header=T)

names(exp1) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')
names(exp2) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')

## Multivariable MR - reading in our own exposure data frames
exp7 <- read.table('BMI_age_10_nov2.mvinput', header=T)
exp8 <- read.table('BMI_adult_nov2.mvinput', header=T)

names(exp7) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')
names(exp8) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')

## add exposure name
exp1$Phenotype <- 'age 10'
exp2$Phenotype <- 'adult'
exp7$Phenotype <- 'age 10'
exp8$Phenotype <- 'adult'

## formats datasets for TwoSampleMR
exp_dat1 <- format_data(exp1, type = "exposure")
exp_dat2 <- format_data(exp2, type = "exposure")
exp_dat7 <- format_data(exp7, type="exposure")
exp_dat8 <- format_data(exp8, type="exposure")
exp_dat <- rbind(exp_dat7, exp_dat8)

all_SNPs <- unique(c(exp_dat1$SNP,exp_dat2$SNP,exp_dat$SNP))

##new way to import outcome data due to funny bug in code in previous iteration
child_tv <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/child_tv.txt")
adult_tv <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/adult_tv.txt")
comb_tv <- rbind(child_tv, adult_tv)
names(comb_tv) <- c('chr','pos','SNP','other_allele','effect_allele','a1','eaf','test','obs','beta','se','l95','u95','tstat','pval','error')
comb_tv$Phenotype <- "tv"
out_dat <- format_data(comb_tv, type="outcome")

tryCatch(dat1 <- harmonise_data(
  exposure_dat = exp_dat1, 
  outcome_dat = out_dat, action=1
), error=function(e) NULL)

tryCatch(dat2 <- harmonise_data(
  exposure_dat = exp_dat2, 
  outcome_dat = out_dat, action=1
), error=function(e) NULL)

tryCatch(res1 <- mr(dat1, method_list = "mr_ivw"), error=function(e) NULL)
tryCatch(res2 <- mr(dat2, method_list = "mr_ivw"), error=function(e) NULL)


## MVMR
##tryCatch(out_dat_mvmr <- extract_outcome_data(exp_dat$SNP, try), error=function(e) NULL)

tryCatch(mvdat_mvmr <- mv_harmonise_data(
  exposure_dat = exp_dat,
  outcome_dat = out_dat, harmonise_strictness = 1
), error=function(e) NULL)


tryCatch(res_mv <- mv_multiple(mvdat_mvmr), error=function(e) NULL)

## write results to a file
res_mv1 <- res_mv$result[order(res_mv$result$exposure),]

tryCatch(all_res_early <- cbind(res1[1,-c(1,4,5)], res_mv1[2,]), error=function(e) NULL)
tryCatch(all_res_adult <- cbind(res2[1,-c(1,4,5)], res_mv1[1,]), error=function(e) NULL)
tryCatch(write.table(all_res_early, "Y:/UKB/R Scripts and Analyses/Output Files/tv_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)
tryCatch(write.table(all_res_adult, "Y:/UKB/R Scripts and Analyses/Output Files/tv_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)

##Clear everything
rm(list = ls())


##GREY MATTER VOLUME

## Univariable MR - reading in our own exposure data frames
exp1 <- read.table('BMI_age_10_nov2.uvinput', header=T)
exp2 <- read.table('BMI_adult_nov2.uvinput', header=T)

names(exp1) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')
names(exp2) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')

## Multivariable MR - reading in our own exposure data frames
exp7 <- read.table('BMI_age_10_nov2.mvinput', header=T)
exp8 <- read.table('BMI_adult_nov2.mvinput', header=T)

names(exp7) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')
names(exp8) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')

## add exposure name
exp1$Phenotype <- 'age 10'
exp2$Phenotype <- 'adult'
exp7$Phenotype <- 'age 10'
exp8$Phenotype <- 'adult'

## formats datasets for TwoSampleMR
exp_dat1 <- format_data(exp1, type = "exposure")
exp_dat2 <- format_data(exp2, type = "exposure")
exp_dat7 <- format_data(exp7, type="exposure")
exp_dat8 <- format_data(exp8, type="exposure")
exp_dat <- rbind(exp_dat7, exp_dat8)

all_SNPs <- unique(c(exp_dat1$SNP,exp_dat2$SNP,exp_dat$SNP))

##new way to import outcome data due to funny bug in code in previous iteration
child_gm <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/child_gm.txt")
adult_gm <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/adult_gm.txt")
comb_gm <- rbind(child_gm, adult_gm)
names(comb_gm) <- c('chr','pos','SNP','other_allele','effect_allele','a1','eaf','test','obs','beta','se','l95','u95','tstat','pval','error')
comb_gm$Phenotype <- "gm"
out_dat <- format_data(comb_gm, type="outcome")

tryCatch(dat1 <- harmonise_data(
  exposure_dat = exp_dat1, 
  outcome_dat = out_dat, action=1
), error=function(e) NULL)

tryCatch(dat2 <- harmonise_data(
  exposure_dat = exp_dat2, 
  outcome_dat = out_dat, action=1
), error=function(e) NULL)

tryCatch(res1 <- mr(dat1, method_list = "mr_ivw"), error=function(e) NULL)
tryCatch(res2 <- mr(dat2, method_list = "mr_ivw"), error=function(e) NULL)


## MVMR
##tryCatch(out_dat_mvmr <- extract_outcome_data(exp_dat$SNP, try), error=function(e) NULL)

tryCatch(mvdat_mvmr <- mv_harmonise_data(
  exposure_dat = exp_dat,
  outcome_dat = out_dat, harmonise_strictness = 1
), error=function(e) NULL)


tryCatch(res_mv <- mv_multiple(mvdat_mvmr), error=function(e) NULL)

## write results to a file
res_mv1 <- res_mv$result[order(res_mv$result$exposure),]

tryCatch(all_res_early <- cbind(res1[1,-c(1,4,5)], res_mv1[2,]), error=function(e) NULL)
tryCatch(all_res_adult <- cbind(res2[1,-c(1,4,5)], res_mv1[1,]), error=function(e) NULL)
tryCatch(write.table(all_res_early, "Y:/UKB/R Scripts and Analyses/Output Files/gm_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)
tryCatch(write.table(all_res_adult, "Y:/UKB/R Scripts and Analyses/Output Files/gm_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)

##Clear everything
rm(list = ls())


##SURFACE AREA

## Univariable MR - reading in our own exposure data frames
exp1 <- read.table('BMI_age_10_nov2.uvinput', header=T)
exp2 <- read.table('BMI_adult_nov2.uvinput', header=T)

names(exp1) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')
names(exp2) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')

## Multivariable MR - reading in our own exposure data frames
exp7 <- read.table('BMI_age_10_nov2.mvinput', header=T)
exp8 <- read.table('BMI_adult_nov2.mvinput', header=T)

names(exp7) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')
names(exp8) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')

## add exposure name
exp1$Phenotype <- 'age 10'
exp2$Phenotype <- 'adult'
exp7$Phenotype <- 'age 10'
exp8$Phenotype <- 'adult'

## formats datasets for TwoSampleMR
exp_dat1 <- format_data(exp1, type = "exposure")
exp_dat2 <- format_data(exp2, type = "exposure")
exp_dat7 <- format_data(exp7, type="exposure")
exp_dat8 <- format_data(exp8, type="exposure")
exp_dat <- rbind(exp_dat7, exp_dat8)

all_SNPs <- unique(c(exp_dat1$SNP,exp_dat2$SNP,exp_dat$SNP))

##new way to import outcome data due to funny bug in code in previous iteration
child_area <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/child_area.txt")
adult_area <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/adult_area.txt")
comb_area <- rbind(child_area, adult_area)
names(comb_area) <- c('chr','pos','SNP','other_allele','effect_allele','a1','eaf','test','obs','beta','se','l95','u95','tstat','pval','error')
comb_area$Phenotype <- "area"
out_dat <- format_data(comb_area, type="outcome")

tryCatch(dat1 <- harmonise_data(
  exposure_dat = exp_dat1, 
  outcome_dat = out_dat, action=1
), error=function(e) NULL)

tryCatch(dat2 <- harmonise_data(
  exposure_dat = exp_dat2, 
  outcome_dat = out_dat, action=1
), error=function(e) NULL)

tryCatch(res1 <- mr(dat1, method_list = "mr_ivw"), error=function(e) NULL)
tryCatch(res2 <- mr(dat2, method_list = "mr_ivw"), error=function(e) NULL)


## MVMR
##tryCatch(out_dat_mvmr <- extract_outcome_data(exp_dat$SNP, try), error=function(e) NULL)

tryCatch(mvdat_mvmr <- mv_harmonise_data(
  exposure_dat = exp_dat,
  outcome_dat = out_dat, harmonise_strictness = 1
), error=function(e) NULL)


tryCatch(res_mv <- mv_multiple(mvdat_mvmr), error=function(e) NULL)

## write results to a file
res_mv1 <- res_mv$result[order(res_mv$result$exposure),]

tryCatch(all_res_early <- cbind(res1[1,-c(1,4,5)], res_mv1[2,]), error=function(e) NULL)
tryCatch(all_res_adult <- cbind(res2[1,-c(1,4,5)], res_mv1[1,]), error=function(e) NULL)
tryCatch(write.table(all_res_early, "Y:/UKB/R Scripts and Analyses/Output Files/area_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)
tryCatch(write.table(all_res_adult, "Y:/UKB/R Scripts and Analyses/Output Files/area_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)

##Clear everything
rm(list = ls())


##CORTICAL THICKNESS

## Univariable MR - reading in our own exposure data frames
exp1 <- read.table('BMI_age_10_nov2.uvinput', header=T)
exp2 <- read.table('BMI_adult_nov2.uvinput', header=T)

names(exp1) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')
names(exp2) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')

## Multivariable MR - reading in our own exposure data frames
exp7 <- read.table('BMI_age_10_nov2.mvinput', header=T)
exp8 <- read.table('BMI_adult_nov2.mvinput', header=T)

names(exp7) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')
names(exp8) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')

## add exposure name
exp1$Phenotype <- 'age 10'
exp2$Phenotype <- 'adult'
exp7$Phenotype <- 'age 10'
exp8$Phenotype <- 'adult'

## formats datasets for TwoSampleMR
exp_dat1 <- format_data(exp1, type = "exposure")
exp_dat2 <- format_data(exp2, type = "exposure")
exp_dat7 <- format_data(exp7, type="exposure")
exp_dat8 <- format_data(exp8, type="exposure")
exp_dat <- rbind(exp_dat7, exp_dat8)

all_SNPs <- unique(c(exp_dat1$SNP,exp_dat2$SNP,exp_dat$SNP))

##new way to import outcome data due to funny bug in code in previous iteration
child_thick <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/child_thick.txt")
adult_thick <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/adult_thick.txt")
comb_thick <- rbind(child_thick, adult_thick)
names(comb_thick) <- c('chr','pos','SNP','other_allele','effect_allele','a1','eaf','test','obs','beta','se','l95','u95','tstat','pval','error')
comb_thick$Phenotype <- "thick"
out_dat <- format_data(comb_thick, type="outcome")

tryCatch(dat1 <- harmonise_data(
  exposure_dat = exp_dat1, 
  outcome_dat = out_dat, action=1
), error=function(e) NULL)

tryCatch(dat2 <- harmonise_data(
  exposure_dat = exp_dat2, 
  outcome_dat = out_dat, action=1
), error=function(e) NULL)

tryCatch(res1 <- mr(dat1, method_list = "mr_ivw"), error=function(e) NULL)
tryCatch(res2 <- mr(dat2, method_list = "mr_ivw"), error=function(e) NULL)


## MVMR
##tryCatch(out_dat_mvmr <- extract_outcome_data(exp_dat$SNP, try), error=function(e) NULL)

tryCatch(mvdat_mvmr <- mv_harmonise_data(
  exposure_dat = exp_dat,
  outcome_dat = out_dat, harmonise_strictness = 1
), error=function(e) NULL)


tryCatch(res_mv <- mv_multiple(mvdat_mvmr), error=function(e) NULL)

## write results to a file
res_mv1 <- res_mv$result[order(res_mv$result$exposure),]

tryCatch(all_res_early <- cbind(res1[1,-c(1,4,5)], res_mv1[2,]), error=function(e) NULL)
tryCatch(all_res_adult <- cbind(res2[1,-c(1,4,5)], res_mv1[1,]), error=function(e) NULL)
tryCatch(write.table(all_res_early, "Y:/UKB/R Scripts and Analyses/Output Files/thick_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)
tryCatch(write.table(all_res_adult, "Y:/UKB/R Scripts and Analyses/Output Files/thick_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)

##Clear everything
rm(list = ls())


##WHITE MATTER VOLUME

## Univariable MR - reading in our own exposure data frames
exp1 <- read.table('BMI_age_10_nov2.uvinput', header=T)
exp2 <- read.table('BMI_adult_nov2.uvinput', header=T)

names(exp1) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')
names(exp2) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')

## Multivariable MR - reading in our own exposure data frames
exp7 <- read.table('BMI_age_10_nov2.mvinput', header=T)
exp8 <- read.table('BMI_adult_nov2.mvinput', header=T)

names(exp7) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')
names(exp8) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')

## add exposure name
exp1$Phenotype <- 'age 10'
exp2$Phenotype <- 'adult'
exp7$Phenotype <- 'age 10'
exp8$Phenotype <- 'adult'

## formats datasets for TwoSampleMR
exp_dat1 <- format_data(exp1, type = "exposure")
exp_dat2 <- format_data(exp2, type = "exposure")
exp_dat7 <- format_data(exp7, type="exposure")
exp_dat8 <- format_data(exp8, type="exposure")
exp_dat <- rbind(exp_dat7, exp_dat8)

all_SNPs <- unique(c(exp_dat1$SNP,exp_dat2$SNP,exp_dat$SNP))

##new way to import outcome data due to funny bug in code in previous iteration
child_wm <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/child_wm.txt")
adult_wm <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/adult_wm.txt")
comb_wm <- rbind(child_wm, adult_wm)
names(comb_wm) <- c('chr','pos','SNP','other_allele','effect_allele','a1','eaf','test','obs','beta','se','l95','u95','tstat','pval','error')
comb_wm$Phenotype <- "wm"
out_dat <- format_data(comb_wm, type="outcome")

tryCatch(dat1 <- harmonise_data(
  exposure_dat = exp_dat1, 
  outcome_dat = out_dat, action=1
), error=function(e) NULL)

tryCatch(dat2 <- harmonise_data(
  exposure_dat = exp_dat2, 
  outcome_dat = out_dat, action=1
), error=function(e) NULL)

tryCatch(res1 <- mr(dat1, method_list = "mr_ivw"), error=function(e) NULL)
tryCatch(res2 <- mr(dat2, method_list = "mr_ivw"), error=function(e) NULL)


## MVMR
##tryCatch(out_dat_mvmr <- extract_outcome_data(exp_dat$SNP, try), error=function(e) NULL)

tryCatch(mvdat_mvmr <- mv_harmonise_data(
  exposure_dat = exp_dat,
  outcome_dat = out_dat, harmonise_strictness = 1
), error=function(e) NULL)


tryCatch(res_mv <- mv_multiple(mvdat_mvmr), error=function(e) NULL)

## write results to a file
res_mv1 <- res_mv$result[order(res_mv$result$exposure),]

tryCatch(all_res_early <- cbind(res1[1,-c(1,4,5)], res_mv1[2,]), error=function(e) NULL)
tryCatch(all_res_adult <- cbind(res2[1,-c(1,4,5)], res_mv1[1,]), error=function(e) NULL)
tryCatch(write.table(all_res_early, "Y:/UKB/R Scripts and Analyses/Output Files/wm_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)
tryCatch(write.table(all_res_adult, "Y:/UKB/R Scripts and Analyses/Output Files/wm_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)

##Clear everything
rm(list = ls())


##LOG WMH

## Univariable MR - reading in our own exposure data frames
exp1 <- read.table('BMI_age_10_nov2.uvinput', header=T)
exp2 <- read.table('BMI_adult_nov2.uvinput', header=T)

names(exp1) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')
names(exp2) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')

## Multivariable MR - reading in our own exposure data frames
exp7 <- read.table('BMI_age_10_nov2.mvinput', header=T)
exp8 <- read.table('BMI_adult_nov2.mvinput', header=T)

names(exp7) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')
names(exp8) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')

## add exposure name
exp1$Phenotype <- 'age 10'
exp2$Phenotype <- 'adult'
exp7$Phenotype <- 'age 10'
exp8$Phenotype <- 'adult'

## formats datasets for TwoSampleMR
exp_dat1 <- format_data(exp1, type = "exposure")
exp_dat2 <- format_data(exp2, type = "exposure")
exp_dat7 <- format_data(exp7, type="exposure")
exp_dat8 <- format_data(exp8, type="exposure")
exp_dat <- rbind(exp_dat7, exp_dat8)

all_SNPs <- unique(c(exp_dat1$SNP,exp_dat2$SNP,exp_dat$SNP))

##new way to import outcome data due to funny bug in code in previous iteration
child_logwmh <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/child_logwmh.txt")
adult_logwmh <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/adult_logwmh.txt")
comb_logwmh <- rbind(child_logwmh, adult_logwmh)
names(comb_logwmh) <- c('chr','pos','SNP','other_allele','effect_allele','a1','eaf','test','obs','beta','se','l95','u95','tstat','pval','error')
comb_logwmh$Phenotype <- "logwmh"
out_dat <- format_data(comb_logwmh, type="outcome")

tryCatch(dat1 <- harmonise_data(
  exposure_dat = exp_dat1, 
  outcome_dat = out_dat, action=1
), error=function(e) NULL)

tryCatch(dat2 <- harmonise_data(
  exposure_dat = exp_dat2, 
  outcome_dat = out_dat, action=1
), error=function(e) NULL)

tryCatch(res1 <- mr(dat1, method_list = "mr_ivw"), error=function(e) NULL)
tryCatch(res2 <- mr(dat2, method_list = "mr_ivw"), error=function(e) NULL)


## MVMR
##tryCatch(out_dat_mvmr <- extract_outcome_data(exp_dat$SNP, try), error=function(e) NULL)

tryCatch(mvdat_mvmr <- mv_harmonise_data(
  exposure_dat = exp_dat,
  outcome_dat = out_dat, harmonise_strictness = 1
), error=function(e) NULL)


tryCatch(res_mv <- mv_multiple(mvdat_mvmr), error=function(e) NULL)

## write results to a file
res_mv1 <- res_mv$result[order(res_mv$result$exposure),]

tryCatch(all_res_early <- cbind(res1[1,-c(1,4,5)], res_mv1[2,]), error=function(e) NULL)
tryCatch(all_res_adult <- cbind(res2[1,-c(1,4,5)], res_mv1[1,]), error=function(e) NULL)
tryCatch(write.table(all_res_early, "Y:/UKB/R Scripts and Analyses/Output Files/logwmh_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)
tryCatch(write.table(all_res_adult, "Y:/UKB/R Scripts and Analyses/Output Files/logwmh_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)

##Clear everything
rm(list = ls())


##HIPPOCAMPUS

## Univariable MR - reading in our own exposure data frames
exp1 <- read.table('BMI_age_10_nov2.uvinput', header=T)
exp2 <- read.table('BMI_adult_nov2.uvinput', header=T)

names(exp1) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')
names(exp2) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')

## Multivariable MR - reading in our own exposure data frames
exp7 <- read.table('BMI_age_10_nov2.mvinput', header=T)
exp8 <- read.table('BMI_adult_nov2.mvinput', header=T)

names(exp7) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')
names(exp8) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')

## add exposure name
exp1$Phenotype <- 'age 10'
exp2$Phenotype <- 'adult'
exp7$Phenotype <- 'age 10'
exp8$Phenotype <- 'adult'

## formats datasets for TwoSampleMR
exp_dat1 <- format_data(exp1, type = "exposure")
exp_dat2 <- format_data(exp2, type = "exposure")
exp_dat7 <- format_data(exp7, type="exposure")
exp_dat8 <- format_data(exp8, type="exposure")
exp_dat <- rbind(exp_dat7, exp_dat8)

all_SNPs <- unique(c(exp_dat1$SNP,exp_dat2$SNP,exp_dat$SNP))

##new way to import outcome data due to funny bug in code in previous iteration
child_hippo <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/child_hippo.txt")
adult_hippo <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/adult_hippo.txt")
comb_hippo <- rbind(child_hippo, adult_hippo)
names(comb_hippo) <- c('chr','pos','SNP','other_allele','effect_allele','a1','eaf','test','obs','beta','se','l95','u95','tstat','pval','error')
comb_hippo$Phenotype <- "hippo"
out_dat <- format_data(comb_hippo, type="outcome")

tryCatch(dat1 <- harmonise_data(
  exposure_dat = exp_dat1, 
  outcome_dat = out_dat, action=1
), error=function(e) NULL)

tryCatch(dat2 <- harmonise_data(
  exposure_dat = exp_dat2, 
  outcome_dat = out_dat, action=1
), error=function(e) NULL)

tryCatch(res1 <- mr(dat1, method_list = "mr_ivw"), error=function(e) NULL)
tryCatch(res2 <- mr(dat2, method_list = "mr_ivw"), error=function(e) NULL)


## MVMR
##tryCatch(out_dat_mvmr <- extract_outcome_data(exp_dat$SNP, try), error=function(e) NULL)

tryCatch(mvdat_mvmr <- mv_harmonise_data(
  exposure_dat = exp_dat,
  outcome_dat = out_dat, harmonise_strictness = 1
), error=function(e) NULL)


tryCatch(res_mv <- mv_multiple(mvdat_mvmr), error=function(e) NULL)

## write results to a file
res_mv1 <- res_mv$result[order(res_mv$result$exposure),]

tryCatch(all_res_early <- cbind(res1[1,-c(1,4,5)], res_mv1[2,]), error=function(e) NULL)
tryCatch(all_res_adult <- cbind(res2[1,-c(1,4,5)], res_mv1[1,]), error=function(e) NULL)
tryCatch(write.table(all_res_early, "Y:/UKB/R Scripts and Analyses/Output Files/hippo_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)
tryCatch(write.table(all_res_adult, "Y:/UKB/R Scripts and Analyses/Output Files/hippo_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)

##Clear everything
rm(list = ls())


##AMYGDALA

## Univariable MR - reading in our own exposure data frames
exp1 <- read.table('BMI_age_10_nov2.uvinput', header=T)
exp2 <- read.table('BMI_adult_nov2.uvinput', header=T)

names(exp1) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')
names(exp2) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')

## Multivariable MR - reading in our own exposure data frames
exp7 <- read.table('BMI_age_10_nov2.mvinput', header=T)
exp8 <- read.table('BMI_adult_nov2.mvinput', header=T)

names(exp7) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')
names(exp8) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')

## add exposure name
exp1$Phenotype <- 'age 10'
exp2$Phenotype <- 'adult'
exp7$Phenotype <- 'age 10'
exp8$Phenotype <- 'adult'

## formats datasets for TwoSampleMR
exp_dat1 <- format_data(exp1, type = "exposure")
exp_dat2 <- format_data(exp2, type = "exposure")
exp_dat7 <- format_data(exp7, type="exposure")
exp_dat8 <- format_data(exp8, type="exposure")
exp_dat <- rbind(exp_dat7, exp_dat8)

all_SNPs <- unique(c(exp_dat1$SNP,exp_dat2$SNP,exp_dat$SNP))

##new way to import outcome data due to funny bug in code in previous iteration
child_amygdala <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/child_amygdala.txt")
adult_amygdala <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/adult_amygdala.txt")
comb_amygdala <- rbind(child_amygdala, adult_amygdala)
names(comb_amygdala) <- c('chr','pos','SNP','other_allele','effect_allele','a1','eaf','test','obs','beta','se','l95','u95','tstat','pval','error')
comb_amygdala$Phenotype <- "amygdala"
out_dat <- format_data(comb_amygdala, type="outcome")

tryCatch(dat1 <- harmonise_data(
  exposure_dat = exp_dat1, 
  outcome_dat = out_dat, action=1
), error=function(e) NULL)

tryCatch(dat2 <- harmonise_data(
  exposure_dat = exp_dat2, 
  outcome_dat = out_dat, action=1
), error=function(e) NULL)

tryCatch(res1 <- mr(dat1, method_list = "mr_ivw"), error=function(e) NULL)
tryCatch(res2 <- mr(dat2, method_list = "mr_ivw"), error=function(e) NULL)


## MVMR
##tryCatch(out_dat_mvmr <- extract_outcome_data(exp_dat$SNP, try), error=function(e) NULL)

tryCatch(mvdat_mvmr <- mv_harmonise_data(
  exposure_dat = exp_dat,
  outcome_dat = out_dat, harmonise_strictness = 1
), error=function(e) NULL)


tryCatch(res_mv <- mv_multiple(mvdat_mvmr), error=function(e) NULL)

## write results to a file
res_mv1 <- res_mv$result[order(res_mv$result$exposure),]

tryCatch(all_res_early <- cbind(res1[1,-c(1,4,5)], res_mv1[2,]), error=function(e) NULL)
tryCatch(all_res_adult <- cbind(res2[1,-c(1,4,5)], res_mv1[1,]), error=function(e) NULL)
tryCatch(write.table(all_res_early, "Y:/UKB/R Scripts and Analyses/Output Files/amygdala_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)
tryCatch(write.table(all_res_adult, "Y:/UKB/R Scripts and Analyses/Output Files/amygdala_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)

##Clear everything
rm(list = ls())


##THALAMUS

## Univariable MR - reading in our own exposure data frames
exp1 <- read.table('BMI_age_10_nov2.uvinput', header=T)
exp2 <- read.table('BMI_adult_nov2.uvinput', header=T)

names(exp1) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')
names(exp2) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')

## Multivariable MR - reading in our own exposure data frames
exp7 <- read.table('BMI_age_10_nov2.mvinput', header=T)
exp8 <- read.table('BMI_adult_nov2.mvinput', header=T)

names(exp7) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')
names(exp8) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')

## add exposure name
exp1$Phenotype <- 'age 10'
exp2$Phenotype <- 'adult'
exp7$Phenotype <- 'age 10'
exp8$Phenotype <- 'adult'

## formats datasets for TwoSampleMR
exp_dat1 <- format_data(exp1, type = "exposure")
exp_dat2 <- format_data(exp2, type = "exposure")
exp_dat7 <- format_data(exp7, type="exposure")
exp_dat8 <- format_data(exp8, type="exposure")
exp_dat <- rbind(exp_dat7, exp_dat8)

all_SNPs <- unique(c(exp_dat1$SNP,exp_dat2$SNP,exp_dat$SNP))

##new way to import outcome data due to funny bug in code in previous iteration
child_thalamus <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/child_thalamus.txt")
adult_thalamus <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/adult_thalamus.txt")
comb_thalamus <- rbind(child_thalamus, adult_thalamus)
names(comb_thalamus) <- c('chr','pos','SNP','other_allele','effect_allele','a1','eaf','test','obs','beta','se','l95','u95','tstat','pval','error')
comb_thalamus$Phenotype <- "thalamus"
out_dat <- format_data(comb_thalamus, type="outcome")

tryCatch(dat1 <- harmonise_data(
  exposure_dat = exp_dat1, 
  outcome_dat = out_dat, action=1
), error=function(e) NULL)

tryCatch(dat2 <- harmonise_data(
  exposure_dat = exp_dat2, 
  outcome_dat = out_dat, action=1
), error=function(e) NULL)

tryCatch(res1 <- mr(dat1, method_list = "mr_ivw"), error=function(e) NULL)
tryCatch(res2 <- mr(dat2, method_list = "mr_ivw"), error=function(e) NULL)


## MVMR
##tryCatch(out_dat_mvmr <- extract_outcome_data(exp_dat$SNP, try), error=function(e) NULL)

tryCatch(mvdat_mvmr <- mv_harmonise_data(
  exposure_dat = exp_dat,
  outcome_dat = out_dat, harmonise_strictness = 1
), error=function(e) NULL)


tryCatch(res_mv <- mv_multiple(mvdat_mvmr), error=function(e) NULL)

## write results to a file
res_mv1 <- res_mv$result[order(res_mv$result$exposure),]

tryCatch(all_res_early <- cbind(res1[1,-c(1,4,5)], res_mv1[2,]), error=function(e) NULL)
tryCatch(all_res_adult <- cbind(res2[1,-c(1,4,5)], res_mv1[1,]), error=function(e) NULL)
tryCatch(write.table(all_res_early, "Y:/UKB/R Scripts and Analyses/Output Files/thalamus_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)
tryCatch(write.table(all_res_adult, "Y:/UKB/R Scripts and Analyses/Output Files/thalamus_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)

##Clear everything
rm(list = ls())
