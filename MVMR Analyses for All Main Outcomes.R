## load libraries into R
library(TwoSampleMR)
library(data.table)

## set our working directory
setwd("Y:/UKB/R Scripts and Analyses") 


##TOTAL VOLUME NORMALISED FOR ICV

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
child_tvn <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/child_tvn.txt")
adult_tvn <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/adult_tvn.txt")
comb_tvn <- rbind(child_tvn, adult_tvn)
names(comb_tvn) <- c('chr','pos','SNP','other_allele','effect_allele','a1','eaf','test','obs','beta','se','l95','u95','tstat','pval','error')
comb_tvn$Phenotype <- "tvn"
out_dat <- format_data(comb_tvn, type="outcome")

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
tryCatch(write.table(all_res_early, "Y:/UKB/R Scripts and Analyses/Output Files/tvn_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)
tryCatch(write.table(all_res_adult, "Y:/UKB/R Scripts and Analyses/Output Files/tvn_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)

##Clear everything
rm(list = ls())


##GREY MATTER VOLUME NORMALISED FOR ICV

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
child_gmn <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/child_gmn.txt")
adult_gmn <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/adult_gmn.txt")
comb_gmn <- rbind(child_gmn, adult_gmn)
names(comb_gmn) <- c('chr','pos','SNP','other_allele','effect_allele','a1','eaf','test','obs','beta','se','l95','u95','tstat','pval','error')
comb_gmn$Phenotype <- "gmn"
out_dat <- format_data(comb_gmn, type="outcome")

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
tryCatch(write.table(all_res_early, "Y:/UKB/R Scripts and Analyses/Output Files/gmn_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)
tryCatch(write.table(all_res_adult, "Y:/UKB/R Scripts and Analyses/Output Files/gmn_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)

##Clear everything
rm(list = ls())


##CORTICAL SURFACE AREA NORMALISED FOR ICV

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
child_arean <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/child_arean.txt")
adult_arean <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/adult_arean.txt")
comb_arean <- rbind(child_arean, adult_arean)
names(comb_arean) <- c('chr','pos','SNP','other_allele','effect_allele','a1','eaf','test','obs','beta','se','l95','u95','tstat','pval','error')
comb_arean$Phenotype <- "arean"
out_dat <- format_data(comb_arean, type="outcome")

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
tryCatch(write.table(all_res_early, "Y:/UKB/R Scripts and Analyses/Output Files/arean_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)
tryCatch(write.table(all_res_adult, "Y:/UKB/R Scripts and Analyses/Output Files/arean_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)

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


##WHITE MATTER VOLUME NORMALISED FOR ICV

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
child_wmn <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/child_wmn.txt")
adult_wmn <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/adult_wmn.txt")
comb_wmn <- rbind(child_wmn, adult_wmn)
names(comb_wmn) <- c('chr','pos','SNP','other_allele','effect_allele','a1','eaf','test','obs','beta','se','l95','u95','tstat','pval','error')
comb_wmn$Phenotype <- "wmn"
out_dat <- format_data(comb_wmn, type="outcome")

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
tryCatch(write.table(all_res_early, "Y:/UKB/R Scripts and Analyses/Output Files/wmn_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)
tryCatch(write.table(all_res_adult, "Y:/UKB/R Scripts and Analyses/Output Files/wmn_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)

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


##HIPPOCAMPUS NORMALISED FOR ICV

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
child_hippon <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/child_hippon.txt")
adult_hippon <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/adult_hippon.txt")
comb_hippon <- rbind(child_hippon, adult_hippon)
names(comb_hippon) <- c('chr','pos','SNP','other_allele','effect_allele','a1','eaf','test','obs','beta','se','l95','u95','tstat','pval','error')
comb_hippon$Phenotype <- "hippon"
out_dat <- format_data(comb_hippon, type="outcome")

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
tryCatch(write.table(all_res_early, "Y:/UKB/R Scripts and Analyses/Output Files/hippon_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)
tryCatch(write.table(all_res_adult, "Y:/UKB/R Scripts and Analyses/Output Files/hippon_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)

##Clear everything
rm(list = ls())


##AMYGDALA NORMALISED FOR ICV

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
child_amygn <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/child_amygn.txt")
adult_amygn <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/adult_amygn.txt")
comb_amygn <- rbind(child_amygn, adult_amygn)
names(comb_amygn) <- c('chr','pos','SNP','other_allele','effect_allele','a1','eaf','test','obs','beta','se','l95','u95','tstat','pval','error')
comb_amygn$Phenotype <- "amygn"
out_dat <- format_data(comb_amygn, type="outcome")

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
tryCatch(write.table(all_res_early, "Y:/UKB/R Scripts and Analyses/Output Files/amygn_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)
tryCatch(write.table(all_res_adult, "Y:/UKB/R Scripts and Analyses/Output Files/amygn_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)

##Clear everything
rm(list = ls())


##THALAMUS NORMALISED FOR ICV

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
child_thaln <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/child_thaln.txt")
adult_thaln <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/adult_thaln.txt")
comb_thaln <- rbind(child_thaln, adult_thaln)
names(comb_thaln) <- c('chr','pos','SNP','other_allele','effect_allele','a1','eaf','test','obs','beta','se','l95','u95','tstat','pval','error')
comb_thaln$Phenotype <- "thaln"
out_dat <- format_data(comb_thaln, type="outcome")

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
tryCatch(write.table(all_res_early, "Y:/UKB/R Scripts and Analyses/Output Files/thaln_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)
tryCatch(write.table(all_res_adult, "Y:/UKB/R Scripts and Analyses/Output Files/thaln_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)

##Clear everything
rm(list = ls())