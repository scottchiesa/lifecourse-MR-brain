
#Title: ABCD GRS analyses
#Author: Lydia Rader
#Date: Winter 2023
#Purpose: The purpose of this script is to run a genetic risk score analysis to see if an adiposity GRS predicts multiple brain structural features in the ABCD
#Models: The models all have the brain features as outcomes; the GRS (either continuous [model set 1] or ordinal [model set 2] based on expected percentiles) as the main predictor of interest; sex, age, and 10 PCs as covariates; and a random intercept for family number to control for relatedness in the ABCD

setwd("~/Documents/CURRENT/Collaborations/MR_UCL_0823")

library(data.table)
library(psych)
library(dplyr)
library(lme4)
library(lmerTest)
library(ggplot2)

pheno_geno <- read.csv('pheno_EUR_geno_020824.csv', header = T)

####1. Set of first models; Using a Continuous z-transformed GRS####

colnames(pheno_geno[,58:66]) # These are the indexed variables
image_phenos <- pheno_geno[,58:66]

# Initialize the list:
image_regressions_list <- list()  

#Run a loop to run each model
for (i in 58:66){       
  fit_regressions <- lmer(pheno_geno[,i] ~ z_score + z_sex + z_age + z_PC1 + z_PC2 + z_PC3 + z_PC4 + z_PC5 + z_PC6 + z_PC7 + z_PC8 + z_PC9 + z_PC10 + (1|rel_family_id), data=pheno_geno)
  image_regressions_list[[i]] <- data.frame(summary(fit_regressions)$coefficients[2,])
}

#Assign Imaging phenos as names to the list of results from regression: 
names(image_regressions_list) <- colnames(pheno_geno[, 58:66])

#pull out the regression coefficients for imaging phenos (main) and save each as a dataframe: 
pgs_estimate <- as.data.frame(unlist(sapply(image_regressions_list, "[[", c(1,1)))) 

#c(1,1) captures the 1st value from the row, which corresponds to the GRS estimate
colnames(pgs_estimate) <- "pgs_estimate"

#pull out the standard error coefficients for imaging phenos (main) and save each as a dataframe: 
standard_error <- as.data.frame(unlist(sapply(image_regressions_list, "[[", c(1,2)))) 

#c(1,1) captures the 1st value from the row, which corresponds to the estimate
colnames(standard_error) <- "standard_error"

#pull out the p-values for imaging phenos (main) and save each as a dataframe: 
#c(1,5) captures the 5th value from the row, which corresponds to p-value
pval_main <- as.data.frame(unlist(sapply(image_regressions_list, "[[", c(1,5)))) 

# Rename the column pval_main
colnames(pval_main) <- "pval_main"

#Combine the above dataframes into one:
df_regressions <- data.frame(pgs_estimate, standard_error, pval_main)
df_regressions$imaging_pheno <- names(image_regressions_list[1:9])
df <- df_regressions %>% relocate(imaging_pheno, .before = pgs_estimate)

#Calculate the 95% confidence intervals based on the standard errors
df$upper_ci <- df$pgs_estimate+(df$standard_error*1.96)
df$lower_ci <- df$pgs_estimate-(df$standard_error*1.96)
df <- df %>% relocate(pval_main, .after = lower_ci)

#write.csv(df, 'summary_icv_index_zPGS_021024.csv', row.names = F)


#create forest plot
ggplot(data=df, aes(y=imaging_pheno, x=pgs_estimate, xmin=lower_ci, xmax=upper_ci)) +
  xlim(-0.1,0.1) +
  geom_point() + 
  geom_errorbarh(height=.1) + 
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  theme_classic() +
  labs(title='Forest plot of continuous PGS estimates w ICV indexed', x='Continuous PGS estimate effect size', y = 'Imaging phenotype') 

# Clear workspace to begin second set of models
rm(df, df_regressions, fit_regressions, image_phenos, image_regressions_list, pgs_estimate, pval_main, standard_error)

####2. Second set of models; Using the ordinal GRS based on these percentiles 0-33;33-85;86+####

colnames(pheno_geno[,58:66]) # These are the indexed variables
image_phenos <- pheno_geno[,58:66]
# Initialize the list:
image_regressions_list <- list()  

for (i in 58:66){       
  fit_regressions <- lmer(pheno_geno[,i] ~ pgs_ord + z_sex + z_age + z_PC1 + z_PC2 + z_PC3 + z_PC4 + z_PC5 + z_PC6 + z_PC7 + z_PC8 + z_PC9 + z_PC10 + (1|rel_family_id), data=pheno_geno)
  image_regressions_list[[i]] <- data.frame(summary(fit_regressions)$coefficients[2,])
}

##Assign Imaging phenos as names to the list of results from regression: 
names(image_regressions_list) <- colnames(pheno_geno[, 58:66])

##pull out the regression coefficients for imaging phenos (main) and save each as a dataframe: 
pgs_estimate <- as.data.frame(unlist(sapply(image_regressions_list, "[[", c(1,1)))) 
#c(1,1) captures the 1st value from the row, which corresponds to the estimate
colnames(pgs_estimate) <- "pgs_estimate"

##pull out the standard error coefficients for imaging phenos (main) and save each as a dataframe: 
standard_error <- as.data.frame(unlist(sapply(image_regressions_list, "[[", c(1,2)))) 
#c(1,1) captures the 1st value from the row, which corresponds to the estimate
colnames(standard_error) <- "standard_error"

##pull out the p-values for imaging phenos (main) and save each as a dataframe: 
#c(1,5) captures the 5th value from the row, which corresponds to p-value
pval_main <- as.data.frame(unlist(sapply(image_regressions_list, "[[", c(1,5)))) 

# Rename the column pval_main
colnames(pval_main) <- "pval_main"

##Combine the above dataframes into one:
df_regressions <- data.frame(pgs_estimate, standard_error, pval_main)
df_regressions$imaging_pheno <- names(image_regressions_list[1:9])
df <- df_regressions %>% relocate(imaging_pheno, .before = pgs_estimate)

#Calculate the 95% confidence intervals based on the standard errors
df$upper_ci <- df$pgs_estimate+(df$standard_error*1.96)
df$lower_ci <- df$pgs_estimate-(df$standard_error*1.96)
df <- df %>% relocate(pval_main, .after = lower_ci)

#write.csv(df, 'summary_icv_index_pgsOrd_021024.csv', row.names = F)

####Forest plots####
#create forest plot
ggplot(data=df, aes(y=imaging_pheno, x=pgs_estimate, xmin=lower_ci, xmax=upper_ci)) +
  xlim(-0.1,0.1) +
  geom_point() + 
  geom_errorbarh(height=.1) + 
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  theme_classic() +
  labs(title='Forest plot of Ordinal PGS estimates w ICV indexed', x='Ordinal PGS estimate effect size', y = 'Imaging phenotype') 
