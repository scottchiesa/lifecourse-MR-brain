## load libraries into R
library(TwoSampleMR)
library(data.table)

## set our working directory
setwd("Y:/UKB/R Scripts and Analyses") 


##GREY MATTER VOLUME NORMALISED FOR TBV

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
child_gmb <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/child_gmb.txt")
adult_gmb <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/adult_gmb.txt")
comb_gmb <- rbind(child_gmb, adult_gmb)
names(comb_gmb) <- c('chr','pos','SNP','other_allele','effect_allele','a1','eaf','test','obs','beta','se','l95','u95','tstat','pval','error')
comb_gmb$Phenotype <- "gmb"
out_dat <- format_data(comb_gmb, type="outcome")

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
tryCatch(write.table(all_res_early, "Y:/UKB/R Scripts and Analyses/Output Files/gmb_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)
tryCatch(write.table(all_res_adult, "Y:/UKB/R Scripts and Analyses/Output Files/gmb_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)

##Clear everything
rm(list = ls())


##CORTICAL SURFACE AREA NORMALISED FOR TBV

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
child_areab <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/child_areab.txt")
adult_areab <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/adult_areab.txt")
comb_areab <- rbind(child_areab, adult_areab)
names(comb_areab) <- c('chr','pos','SNP','other_allele','effect_allele','a1','eaf','test','obs','beta','se','l95','u95','tstat','pval','error')
comb_areab$Phenotype <- "areab"
out_dat <- format_data(comb_areab, type="outcome")

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
tryCatch(write.table(all_res_early, "Y:/UKB/R Scripts and Analyses/Output Files/areab_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)
tryCatch(write.table(all_res_adult, "Y:/UKB/R Scripts and Analyses/Output Files/areab_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)

##Clear everything
rm(list = ls())


##CORTICAL THICKNESS (NOT NORMALISED AS NO JUSTIFICATION TO DO SO)

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


##WHITE MATTER VOLUME NORMALISED FOR TBV

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
child_wmb <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/child_wmb.txt")
adult_wmb <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/adult_wmb.txt")
comb_wmb <- rbind(child_wmb, adult_wmb)
names(comb_wmb) <- c('chr','pos','SNP','other_allele','effect_allele','a1','eaf','test','obs','beta','se','l95','u95','tstat','pval','error')
comb_wmb$Phenotype <- "wmb"
out_dat <- format_data(comb_wmb, type="outcome")

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
tryCatch(write.table(all_res_early, "Y:/UKB/R Scripts and Analyses/Output Files/wmb_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)
tryCatch(write.table(all_res_adult, "Y:/UKB/R Scripts and Analyses/Output Files/wmb_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)

##Clear everything
rm(list = ls())


##LOG WMH (NOT NORMALISED AS NO JUSTIFICATION TO DO SO)

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
comb_logwmh$Phenotype <- "hippo"
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


##HIPPOCAMPUS NORMALISED FOR TBV

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
child_hippob <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/child_hippob.txt")
adult_hippob <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/adult_hippob.txt")
comb_hippob <- rbind(child_hippob, adult_hippob)
names(comb_hippob) <- c('chr','pos','SNP','other_allele','effect_allele','a1','eaf','test','obs','beta','se','l95','u95','tstat','pval','error')
comb_hippob$Phenotype <- "hippob"
out_dat <- format_data(comb_hippob, type="outcome")

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
tryCatch(write.table(all_res_early, "Y:/UKB/R Scripts and Analyses/Output Files/hippob_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)
tryCatch(write.table(all_res_adult, "Y:/UKB/R Scripts and Analyses/Output Files/hippob_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)

##Clear everything
rm(list = ls())


##AMYGDALA NORMALISED FOR TBV

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
child_amygb <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/child_amygb.txt")
adult_amygb <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/adult_amygb.txt")
comb_amygb <- rbind(child_amygb, adult_amygb)
names(comb_amygb) <- c('chr','pos','SNP','other_allele','effect_allele','a1','eaf','test','obs','beta','se','l95','u95','tstat','pval','error')
comb_amygb$Phenotype <- "amygb"
out_dat <- format_data(comb_amygb, type="outcome")

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
tryCatch(write.table(all_res_early, "Y:/UKB/R Scripts and Analyses/Output Files/amygb_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)
tryCatch(write.table(all_res_adult, "Y:/UKB/R Scripts and Analyses/Output Files/amygb_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)

##Clear everything
rm(list = ls())


##THALAMUS NORMALISED FOR TBV

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
child_thalb <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/child_thalb.txt")
adult_thalb <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/adult_thalb.txt")
comb_thalb <- rbind(child_thalb, adult_thalb)
names(comb_thalb) <- c('chr','pos','SNP','other_allele','effect_allele','a1','eaf','test','obs','beta','se','l95','u95','tstat','pval','error')
comb_thalb$Phenotype <- "thalb"
out_dat <- format_data(comb_thalb, type="outcome")

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
tryCatch(write.table(all_res_early, "Y:/UKB/R Scripts and Analyses/Output Files/thalb_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)
tryCatch(write.table(all_res_adult, "Y:/UKB/R Scripts and Analyses/Output Files/thalb_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)

##Clear everything
rm(list = ls())
