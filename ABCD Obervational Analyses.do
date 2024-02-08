***************************************************************************************************************************************************************
********************************************************ABCD Dataset Preparation*******************************************************************************
***************************************************************************************************************************************************************

//Importing and Merging Necessary Data//

//Age, sex, and anthropometrics//
import delimited "Y:\ABCD\Data\Package_1208597\abcd_ant01.txt", varnames(1) clear
drop in 1
keep subjectkey interview_age sex eventname anthroheightcalc anthroweight1lb anthroweight2lb anthro_weight1_hybrid_lb  anthro_waist_cm anthro_waist_hybrid_cm
save "Y:\ABCD\Data Files\ABCD Master", replace

//Child's Race, Parent's Education and Income//

import delimited "Y:\ABCD\Data\Package_1208597\pdem02.txt", varnames(1) clear
drop in 1
save "Y:\ABCD\Data Files\pdem02.dta", replace
use "Y:\ABCD\Data Files\ABCD Master.dta", clear
merge 1:1 subjectkey eventname using "Y:\ABCD\Data Files\pdem02.dta", keepusing(demo_race_a_p___10 - demo_race_a_p___99 demo_prnt_ed_v2 demo_prnt_empl_v2 demo_prnt_income_v2 demo_prtnr_ed_v2 demo_prtnr_empl_v2 demo_prtnr_income_v2 demo_comb_income_v2)
rename _merge _merge_race_edu_income
save "Y:\ABCD\Data Files\ABCD Master", replace

//Ethnicity//

import delimited "Y:\ABCD\Data\Package_1208597\multigrp_ethnic_id_meim01.txt", varnames(1) clear
drop in 1
save "Y:\ABCD\Data Files\multigrp_ethnic_id_meim01.dta", replace
use "Y:\ABCD\Data Files\ABCD Master.dta", clear
merge 1:1 subjectkey eventname using "Y:\ABCD\Data Files\multigrp_ethnic_id_meim01.dta", keepusing(meim_ethnic_id)
rename _merge _merge_ethnic
save "Y:\ABCD\Data Files\ABCD Master", replace

//blood pressure//

import delimited "Y:\ABCD\Data\Package_1208597\abcd_bp01.txt", varnames(1) clear
drop in 1
save "Y:\ABCD\Data Files\abcd_bp01.dta", replace
use "Y:\ABCD\Data Files\ABCD Master.dta", clear
merge 1:1 subjectkey eventname using "Y:\ABCD\Data Files\abcd_bp01.dta", keepusing(blood_pressure_test_sys_1 - blood_pressure_dia_mean)
rename _merge _merge_bp
save "Y:\ABCD\Data Files\ABCD Master", replace

//MRI models//
import delimited "Y:\ABCD\Data\Package_1208597\abcd_mri01.txt", varnames(1) clear
drop in 1
save "Y:\ABCD\Data Files\abcd_mri01.dta", replace
use "Y:\ABCD\Data Files\ABCD Master.dta", clear
merge 1:1 subjectkey eventname using "Y:\ABCD\Data Files\abcd_mri01.dta", keepusing(mri_info_manufacturer - mri_info_softwareversion)
rename _merge _merge_mri
save "Y:\ABCD\Data Files\ABCD Master", replace

//measures of brain structure//
import delimited "Y:\ABCD\Data\Package_1208597\abcd_smrip10201.txt", varnames(1) clear
drop in 1
save "Y:\ABCD\Data Files\abcd_smrip10201.dta", replace
use "Y:\ABCD\Data Files\ABCD Master.dta", clear
merge 1:1 subjectkey eventname using "Y:\ABCD\Data Files\abcd_smrip10201.dta", keepusing(smri_thick_cdk_meanlh smri_thick_cdk_meanrh smri_thick_cdk_mean smri_sulc_cdk_meanlh smri_sulc_cdk_meanrh smri_sulc_cdk_mean smri_area_cdk_totallh smri_area_cdk_totalrh smri_area_cdk_total smri_vol_cdk_totallh smri_vol_cdk_totalrh smri_vol_cdk_total smri_vol_scs_wholeb smri_vol_scs_tplh smri_vol_scs_tprh smri_vol_scs_hpuslh smri_vol_scs_hpusrh smri_vol_scs_amygdalalh smri_vol_scs_amygdalarh smri_vol_scs_wmhint smri_vol_scs_intracranialv)
rename _merge _merge_brain
save "Y:\ABCD\Data Files\ABCD Master", replace

//convert everything to numeric//
destring, replace

//Recoding certain string categorical variables//

replace eventname = "0" if eventname == "baseline_year_1_arm_1"
replace eventname = "1" if eventname == "1_year_follow_up_y_arm_1"
replace eventname = "2" if eventname == "2_year_follow_up_y_arm_1"
replace eventname = "3" if eventname == "3_year_follow_up_y_arm_1"
destring eventname, replace
label define eventname 0 "baseline_year_1_arm_1" 1 "1_year_follow_up_y_arm_1" 2 "2_year_follow_up_y_arm_1" 3 "3_year_follow_up_y_arm_1"
label values eventname eventname

replace mri_info_manufacturer = "1" if mri_info_manufacturer == "GE MEDICAL SYSTEMS"
replace mri_info_manufacturer = "2" if mri_info_manufacturer ==  "Philips Medical Systems"
replace mri_info_manufacturer = "3" if mri_info_manufacturer == "SIEMENS"
destring mri_info_manufacturer, replace
label define mri_info_manufacturer 1 "GE MEDICAL SYSTEMS" 2 "Philips Medical Systems" 3 "SIEMENS"
label values mri_info_manufacturer mri_info_manufacturer

//Renaming and creating new variables//

//Sex//

replace sex = "1" if sex == "M"
replace sex = "2" if sex == "F"
destring sex, replace
label define sex 1 "M" 2 "F", modify
label values sex sex

//Age//

gen age = interview_age/12

//Race//

gen race = .
replace race = 0 if demo_race_a_p___10 == 1
replace race = 1 if demo_race_a_p___11 == 1
replace race = 2 if demo_race_a_p___12 == 1
replace race = 3 if demo_race_a_p___13 == 1
replace race = 4 if demo_race_a_p___14 == 1
replace race = 5 if demo_race_a_p___15 == 1
replace race = 6 if demo_race_a_p___16 == 1
replace race = 7 if demo_race_a_p___17 == 1
replace race = 8 if demo_race_a_p___18 == 1
replace race = 9 if demo_race_a_p___19 == 1
replace race = 10 if demo_race_a_p___20 == 1
replace race = 11 if demo_race_a_p___21 == 1
replace race = 12 if demo_race_a_p___22 == 1
replace race = 13 if demo_race_a_p___23 == 1
replace race = 14 if demo_race_a_p___24 == 1
replace race = 15 if demo_race_a_p___25 == 1

recode race (0 = 0) (1 = 1) (2/15 = 2)
label define race 0 "White" 1 "Black" 2 "Other"
label values race race

//Anthropometrics//

gen weight = .
replace weight = (anthroweight1lb + anthroweight2lb)/2
replace weight = anthro_weight1_hybrid_lb if anthroweight1lb == . & anthroweight2lb == .
replace weight = weight*0.4536

gen height = anthroheightcalc*0.0254
replace height = . if height < 1

gen waist = . 
replace waist = anthro_waist_cm
replace anthro_waist_hybrid_cm = "." if anthro_waist_hybrid_cm == "2ol"
destring  anthro_waist_hybrid_cm, replace
replace waist = anthro_waist_hybrid_cm if anthro_waist_cm == .

gen bmi = weight/(height^2)	
drop if bmi > 40		//gets rid of three crazy values//
gen whr = waist/height

//Household Income//

gen income = .
replace income = demo_comb_income_v2
replace income = . if demo_comb_income_v2 > 10

//Parent's Education//

gen parent_ed = .
replace parent_ed = demo_prnt_ed_v2
replace parent_ed = . if parent_ed > 21

//Merge in race income and parent data that have been filled in for every time point in another dataset or later BMI stuff doesn't work//

drop race income parent_ed

merge 1:1 subjectkey eventname using "Y:\ABCD\Data Files\reshaped race income parent data.dta"

//Brain Structures//

clonevar total_vol = smri_vol_scs_wholeb
clonevar cort_thick = smri_thick_cdk_mean
clonevar sulc_depth = smri_sulc_cdk_mean
clonevar cort_area = smri_area_cdk_total
clonevar cort_vol = smri_vol_cdk_total
gen thalamus = smri_vol_scs_tplh + smri_vol_scs_tprh
gen hippocampus = smri_vol_scs_hpuslh + smri_vol_scs_hpusrh
gen amygdala = smri_vol_scs_amygdalalh + smri_vol_scs_amygdalarh
clonevar wmh = smri_vol_scs_wmhint
gen logwmh = log(wmh)
clonevar icv = smri_vol_scs_intracranialv
clonevar iq = nihtbx_totalcomp_agecorrected

zscore whr waist total_vol cort_vol cort_area cort_thick sulc_depth logwmh hippocampus amygdala thalamus

save "Y:\ABCD\Data Files\ABCD Master", replace


***************************************************************************************************************************************************************************************
**************************************************************************RUNNING OBSERVATIONAL ANALYSES*******************************************************************************
***************************************************************************************************************************************************************************************

//Keep only European ancestry at baseline//

keep if eventname == 0
keep if race == 0

//Generate normalised values//

ds total_vol cort_thick sulc_depth cort_area cort_vol logwmh thalamus hippocampus amygdala 
local outcomes `r(varlist)'
foreach outcome of local outcomes {
	gen `outcome'_n = `outcome'/icv
} 

//Create BMI z-scores using CDC Growth Charts//

egen bmi_CDC = zanthro(bmi,ba,US), xvar(age) gender(sex) gencode(male=1, female=2) ageunit(year) 		
gen obesity = 0
sum bmi_CDC, det
replace obesity = 1 if bmi_CDC >= 2.248584

//Create categories to match those identified in previous papers//

egen bmi_perc = cut(bmi_CDC), group(100)
gen bmi_cat = .
replace bmi_cat = 1 if bmi_perc <=33 
replace bmi_cat = 2 if bmi_perc >33 & bmi_perc <=84
replace bmi_cat = 3 if bmi_perc > 84 & bmi_perc < .

// Generate some descriptive stats//

sum age bmi bmi_CDC income parent_ed total_vol cort_area cort_thick sulc_depth cort_area wmh hippocampus amygdala thalamus icv if income !=. & bmi_CDC !=.
tab sex if income !=. & bmi_CDC !=.

//Run regressions//

zscore total_vol_n cort_thick_n sulc_depth_n cort_area_n cort_vol_n logwmh_n thalamus_n hippocampus_n amygdala_n

ds z_total_vol z_cort_thick z_sulc_depth z_cort_area z_cort_vol z_logwmh z_thalamus z_hippocampus z_amygdala z_total_vol_n z_cort_thick_n z_sulc_depth_n z_cort_area_n z_cort_vol_n z_logwmh_n z_thalamus_n z_hippocampus_n z_amygdala_n
local outcomes `r(varlist)'
foreach outcome of local outcomes {
	regress `outcome' bmi_cat age i.sex parent_ed income mri_info_manufacturer
} 