## load libraries into R
library(TwoSampleMR)
library(data.table)

## set our working directory
setwd("Y:/UKB/R Scripts and Analyses/Birthweight Analysis") 


##TOTAL VOLUME NORMALISED FOR ICV

## Univariable MR - reading in our own exposure data frames
exp1 <- read.table('birthweight_nov2.uvinput', header=T)
exp2 <- read.table('BMI_age_10_bw_nov2.uvinput', header=T)

names(exp1) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')
names(exp2) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')

## Multivariable MR - reading in our own exposure data frames
exp7 <- read.table('birthweight_nov2.mvinput', header=T)
exp8 <- read.table('BMI_age_10_bw_nov2.mvinput', header=T)

names(exp7) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')
names(exp8) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')

## add exposure name
exp1$Phenotype <- 'birthweight'
exp2$Phenotype <- 'age 10'
exp7$Phenotype <- 'birthweight'
exp8$Phenotype <- 'age 10'

## formats datasets for TwoSampleMR
exp_dat1 <- format_data(exp1, type = "exposure")
exp_dat2 <- format_data(exp2, type = "exposure")
exp_dat7 <- format_data(exp7, type="exposure")
exp_dat8 <- format_data(exp8, type="exposure")
exp_dat <- rbind(exp_dat7, exp_dat8)

all_SNPs <- unique(c(exp_dat1$SNP,exp_dat2$SNP,exp_dat$SNP))

##new way to import outcome data due to funny bug in code in previous iteration
bw_tvn <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/bw_tvn.txt")
child_bw_tvn <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/child_bw_tvn.txt")
comb_bw_tvn <- rbind(bw_tvn, child_bw_tvn)
names(comb_bw_tvn) <- c('chr','pos','SNP','other_allele','effect_allele','a1','eaf','test','obs','beta','se','l95','u95','tstat','pval','error')
comb_bw_tvn$Phenotype <- "tvn"
out_dat <- format_data(comb_bw_tvn, type="outcome")

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
tryCatch(write.table(all_res_early, "Y:/UKB/R Scripts and Analyses/Output Files/tvn_bw_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)
tryCatch(write.table(all_res_adult, "Y:/UKB/R Scripts and Analyses/Output Files/tvn_bw_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)

##Clear everything
rm(list = ls())


##GREY MATTER VOLUME NORMALISED FOR ICV

## Univariable MR - reading in our own exposure data frames
exp1 <- read.table('birthweight_nov2.uvinput', header=T)
exp2 <- read.table('BMI_age_10_bw_nov2.uvinput', header=T)

names(exp1) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')
names(exp2) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')

## Multivariable MR - reading in our own exposure data frames
exp7 <- read.table('birthweight_nov2.mvinput', header=T)
exp8 <- read.table('BMI_age_10_bw_nov2.mvinput', header=T)

names(exp7) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')
names(exp8) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')

## add exposure name
exp1$Phenotype <- 'birthweight'
exp2$Phenotype <- 'age 10'
exp7$Phenotype <- 'birthweight'
exp8$Phenotype <- 'age 10'

## formats datasets for TwoSampleMR
exp_dat1 <- format_data(exp1, type = "exposure")
exp_dat2 <- format_data(exp2, type = "exposure")
exp_dat7 <- format_data(exp7, type="exposure")
exp_dat8 <- format_data(exp8, type="exposure")
exp_dat <- rbind(exp_dat7, exp_dat8)

all_SNPs <- unique(c(exp_dat1$SNP,exp_dat2$SNP,exp_dat$SNP))

##new way to import outcome data due to funny bug in code in previous iteration
bw_gmn <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/bw_gmn.txt")
child_bw_gmn <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/child_bw_gmn.txt")
comb_bw_gmn <- rbind(bw_gmn, child_bw_gmn)
names(comb_bw_gmn) <- c('chr','pos','SNP','other_allele','effect_allele','a1','eaf','test','obs','beta','se','l95','u95','tstat','pval','error')
comb_bw_gmn$Phenotype <- "gmn"
out_dat <- format_data(comb_bw_gmn, type="outcome")

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
tryCatch(write.table(all_res_early, "Y:/UKB/R Scripts and Analyses/Output Files/gmn_bw_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)
tryCatch(write.table(all_res_adult, "Y:/UKB/R Scripts and Analyses/Output Files/gmn_bw_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)

##Clear everything
rm(list = ls())


##CORTICAL SURFACE AREA NORMALISED FOR ICV

## Univariable MR - reading in our own exposure data frames
exp1 <- read.table('birthweight_nov2.uvinput', header=T)
exp2 <- read.table('BMI_age_10_bw_nov2.uvinput', header=T)

names(exp1) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')
names(exp2) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')

## Multivariable MR - reading in our own exposure data frames
exp7 <- read.table('birthweight_nov2.mvinput', header=T)
exp8 <- read.table('BMI_age_10_bw_nov2.mvinput', header=T)

names(exp7) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')
names(exp8) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')

## add exposure name
exp1$Phenotype <- 'birthweight'
exp2$Phenotype <- 'age 10'
exp7$Phenotype <- 'birthweight'
exp8$Phenotype <- 'age 10'

## formats datasets for TwoSampleMR
exp_dat1 <- format_data(exp1, type = "exposure")
exp_dat2 <- format_data(exp2, type = "exposure")
exp_dat7 <- format_data(exp7, type="exposure")
exp_dat8 <- format_data(exp8, type="exposure")
exp_dat <- rbind(exp_dat7, exp_dat8)

all_SNPs <- unique(c(exp_dat1$SNP,exp_dat2$SNP,exp_dat$SNP))

##new way to import outcome data due to funny bug in code in previous iteration
bw_arean <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/bw_arean.txt")
child_bw_arean <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/child_bw_arean.txt")
comb_bw_arean <- rbind(bw_arean, child_bw_arean)
names(comb_bw_arean) <- c('chr','pos','SNP','other_allele','effect_allele','a1','eaf','test','obs','beta','se','l95','u95','tstat','pval','error')
comb_bw_arean$Phenotype <- "arean"
out_dat <- format_data(comb_bw_arean, type="outcome")

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
tryCatch(write.table(all_res_early, "Y:/UKB/R Scripts and Analyses/Output Files/arean_bw_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)
tryCatch(write.table(all_res_adult, "Y:/UKB/R Scripts and Analyses/Output Files/arean_bw_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)

##Clear everything
rm(list = ls())


##CORTICAL THICKNESS (NOT NORMALISED AS NO JUSTIFICATION TO DO SO)

## Univariable MR - reading in our own exposure data frames
exp1 <- read.table('birthweight_nov2.uvinput', header=T)
exp2 <- read.table('BMI_age_10_bw_nov2.uvinput', header=T)

names(exp1) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')
names(exp2) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')

## Multivariable MR - reading in our own exposure data frames
exp7 <- read.table('birthweight_nov2.mvinput', header=T)
exp8 <- read.table('BMI_age_10_bw_nov2.mvinput', header=T)

names(exp7) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')
names(exp8) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')

## add exposure name
exp1$Phenotype <- 'birthweight'
exp2$Phenotype <- 'age 10'
exp7$Phenotype <- 'birthweight'
exp8$Phenotype <- 'age 10'

## formats datasets for TwoSampleMR
exp_dat1 <- format_data(exp1, type = "exposure")
exp_dat2 <- format_data(exp2, type = "exposure")
exp_dat7 <- format_data(exp7, type="exposure")
exp_dat8 <- format_data(exp8, type="exposure")
exp_dat <- rbind(exp_dat7, exp_dat8)

all_SNPs <- unique(c(exp_dat1$SNP,exp_dat2$SNP,exp_dat$SNP))

##new way to import outcome data due to funny bug in code in previous iteration
bw_thick <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/bw_thick.txt")
child_bw_thick <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/child_bw_thick.txt")
comb_thick <- rbind(bw_thick, child_bw_thick)
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
tryCatch(write.table(all_res_early, "Y:/UKB/R Scripts and Analyses/Output Files/thick_bw_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)
tryCatch(write.table(all_res_adult, "Y:/UKB/R Scripts and Analyses/Output Files/thick_bw_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)

##Clear everything
rm(list = ls())


##WHITE MATTER VOLUME NORMALISED FOR ICV

## Univariable MR - reading in our own exposure data frames
exp1 <- read.table('birthweight_nov2.uvinput', header=T)
exp2 <- read.table('BMI_age_10_bw_nov2.uvinput', header=T)

names(exp1) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')
names(exp2) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')

## Multivariable MR - reading in our own exposure data frames
exp7 <- read.table('birthweight_nov2.mvinput', header=T)
exp8 <- read.table('BMI_age_10_bw_nov2.mvinput', header=T)

names(exp7) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')
names(exp8) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')

## add exposure name
exp1$Phenotype <- 'birthweight'
exp2$Phenotype <- 'age 10'
exp7$Phenotype <- 'birthweight'
exp8$Phenotype <- 'age 10'

## formats datasets for TwoSampleMR
exp_dat1 <- format_data(exp1, type = "exposure")
exp_dat2 <- format_data(exp2, type = "exposure")
exp_dat7 <- format_data(exp7, type="exposure")
exp_dat8 <- format_data(exp8, type="exposure")
exp_dat <- rbind(exp_dat7, exp_dat8)

all_SNPs <- unique(c(exp_dat1$SNP,exp_dat2$SNP,exp_dat$SNP))

##new way to import outcome data due to funny bug in code in previous iteration
bw_wmn <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/bw_wmn.txt")
child_bw_wmn <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/child_bw_wmn.txt")
comb_wmn <- rbind(bw_wmn, child_bw_wmn)
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
tryCatch(write.table(all_res_early, "Y:/UKB/R Scripts and Analyses/Output Files/wmn_bw_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)
tryCatch(write.table(all_res_adult, "Y:/UKB/R Scripts and Analyses/Output Files/wmn_bw_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)

##Clear everything
rm(list = ls())


##LOG WMH (NOT NORMALISED AS NO JUSTIFICATION TO DO SO)

## Univariable MR - reading in our own exposure data frames
exp1 <- read.table('birthweight_nov2.uvinput', header=T)
exp2 <- read.table('BMI_age_10_bw_nov2.uvinput', header=T)

names(exp1) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')
names(exp2) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')

## Multivariable MR - reading in our own exposure data frames
exp7 <- read.table('birthweight_nov2.mvinput', header=T)
exp8 <- read.table('BMI_age_10_bw_nov2.mvinput', header=T)

names(exp7) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')
names(exp8) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')

## add exposure name
exp1$Phenotype <- 'birthweight'
exp2$Phenotype <- 'age 10'
exp7$Phenotype <- 'birthweight'
exp8$Phenotype <- 'age 10'

## formats datasets for TwoSampleMR
exp_dat1 <- format_data(exp1, type = "exposure")
exp_dat2 <- format_data(exp2, type = "exposure")
exp_dat7 <- format_data(exp7, type="exposure")
exp_dat8 <- format_data(exp8, type="exposure")
exp_dat <- rbind(exp_dat7, exp_dat8)

all_SNPs <- unique(c(exp_dat1$SNP,exp_dat2$SNP,exp_dat$SNP))

##new way to import outcome data due to funny bug in code in previous iteration
bw_logwmh <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/bw_logwmh.txt")
child_bw_logwmh <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/child_bw_logwmh.txt")
comb_logwmh <- rbind(bw_logwmh, child_bw_logwmh)
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
tryCatch(write.table(all_res_early, "Y:/UKB/R Scripts and Analyses/Output Files/logwmh_bw_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)
tryCatch(write.table(all_res_adult, "Y:/UKB/R Scripts and Analyses/Output Files/logwmh_bw_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)

##Clear everything
rm(list = ls())


##HIPPOCAMPUS NORMALISED FOR ICV

## Univariable MR - reading in our own exposure data frames
exp1 <- read.table('birthweight_nov2.uvinput', header=T)
exp2 <- read.table('BMI_age_10_bw_nov2.uvinput', header=T)

names(exp1) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')
names(exp2) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')

## Multivariable MR - reading in our own exposure data frames
exp7 <- read.table('birthweight_nov2.mvinput', header=T)
exp8 <- read.table('BMI_age_10_bw_nov2.mvinput', header=T)

names(exp7) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')
names(exp8) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')

## add exposure name
exp1$Phenotype <- 'birthweight'
exp2$Phenotype <- 'age 10'
exp7$Phenotype <- 'birthweight'
exp8$Phenotype <- 'age 10'

## formats datasets for TwoSampleMR
exp_dat1 <- format_data(exp1, type = "exposure")
exp_dat2 <- format_data(exp2, type = "exposure")
exp_dat7 <- format_data(exp7, type="exposure")
exp_dat8 <- format_data(exp8, type="exposure")
exp_dat <- rbind(exp_dat7, exp_dat8)

all_SNPs <- unique(c(exp_dat1$SNP,exp_dat2$SNP,exp_dat$SNP))

##new way to import outcome data due to funny bug in code in previous iteration
bw_hippon <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/bw_hippon.txt")
child_bw_hippon <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/child_bw_hippon.txt")
comb_hippon <- rbind(bw_hippon, child_bw_hippon)
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
tryCatch(write.table(all_res_early, "Y:/UKB/R Scripts and Analyses/Output Files/hippon_bw_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)
tryCatch(write.table(all_res_adult, "Y:/UKB/R Scripts and Analyses/Output Files/hippon_bw_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)

##Clear everything
rm(list = ls())


##AMYGDALA NORMALISED FOR ICV

## Univariable MR - reading in our own exposure data frames
exp1 <- read.table('birthweight_nov2.uvinput', header=T)
exp2 <- read.table('BMI_age_10_bw_nov2.uvinput', header=T)

names(exp1) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')
names(exp2) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')

## Multivariable MR - reading in our own exposure data frames
exp7 <- read.table('birthweight_nov2.mvinput', header=T)
exp8 <- read.table('BMI_age_10_bw_nov2.mvinput', header=T)

names(exp7) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')
names(exp8) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')

## add exposure name
exp1$Phenotype <- 'birthweight'
exp2$Phenotype <- 'age 10'
exp7$Phenotype <- 'birthweight'
exp8$Phenotype <- 'age 10'

## formats datasets for TwoSampleMR
exp_dat1 <- format_data(exp1, type = "exposure")
exp_dat2 <- format_data(exp2, type = "exposure")
exp_dat7 <- format_data(exp7, type="exposure")
exp_dat8 <- format_data(exp8, type="exposure")
exp_dat <- rbind(exp_dat7, exp_dat8)

all_SNPs <- unique(c(exp_dat1$SNP,exp_dat2$SNP,exp_dat$SNP))

##new way to import outcome data due to funny bug in code in previous iteration
bw_amygn <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/bw_amygn.txt")
child_bw_amygn <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/child_bw_amygn.txt")
comb_amygn <- rbind(bw_amygn, child_bw_amygn)
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
tryCatch(write.table(all_res_early, "Y:/UKB/R Scripts and Analyses/Output Files/amygn_bw_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)
tryCatch(write.table(all_res_adult, "Y:/UKB/R Scripts and Analyses/Output Files/amygn_bw_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)

##Clear everything
rm(list = ls())


##THALAMUS NORMALISED FOR ICV

## Univariable MR - reading in our own exposure data frames
exp1 <- read.table('birthweight_nov2.uvinput', header=T)
exp2 <- read.table('BMI_age_10_bw_nov2.uvinput', header=T)

names(exp1) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')
names(exp2) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')

## Multivariable MR - reading in our own exposure data frames
exp7 <- read.table('birthweight_nov2.mvinput', header=T)
exp8 <- read.table('BMI_age_10_bw_nov2.mvinput', header=T)

names(exp7) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')
names(exp8) <- c('SNP','CHR','BP','GENPOS','effect_allele','other_allele','eaf','info',
                 'chisq_lm','p','beta','se','chisq_bolt','pval')

## add exposure name
exp1$Phenotype <- 'birthweight'
exp2$Phenotype <- 'age 10'
exp7$Phenotype <- 'birthweight'
exp8$Phenotype <- 'age 10'

## formats datasets for TwoSampleMR
exp_dat1 <- format_data(exp1, type = "exposure")
exp_dat2 <- format_data(exp2, type = "exposure")
exp_dat7 <- format_data(exp7, type="exposure")
exp_dat8 <- format_data(exp8, type="exposure")
exp_dat <- rbind(exp_dat7, exp_dat8)

all_SNPs <- unique(c(exp_dat1$SNP,exp_dat2$SNP,exp_dat$SNP))

##new way to import outcome data due to funny bug in code in previous iteration
bw_thaln <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/bw_thaln.txt")
child_bw_thaln <- read.table("Y:/UKB/R Scripts and Analyses/Individual Data Files from Myriad/child_bw_thaln.txt")
comb_thaln <- rbind(bw_thaln, child_bw_thaln)
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
tryCatch(write.table(all_res_early, "Y:/UKB/R Scripts and Analyses/Output Files/thaln_bw_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)
tryCatch(write.table(all_res_adult, "Y:/UKB/R Scripts and Analyses/Output Files/thaln_bw_output_file.txt", row.names=F, quote=F, sep="\t", append=T, col.names=F), error=function(e) NULL)

##Clear everything
rm(list = ls())