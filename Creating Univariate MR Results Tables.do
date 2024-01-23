*****************RUN THIS FILE AFTER USING STATA TO PERFORM ALL ANALYSES**********************
*****************WILL CREATE A NICE DATASET WITH ALL UNIVARIATE RESULTS***********************

//Normalised to ICV//

import excel "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\tvn_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\tvn_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\gmn_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\gmn_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\arean_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\arean_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thick_results.xlsx", sheet("Sheet1") firstrow clear		//not normalised//
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thick_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\wmn_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\wmn_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\logwmh_results.xlsx", sheet("Sheet1") firstrow clear	//not normalised//
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\logwmh_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thaln_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thaln_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\hippon_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\hippon_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\amygn_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\amygn_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\csfn_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\csfn_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\tvn_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\tvn_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\gmn_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\gmn_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\arean_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\arean_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thickn_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thickn_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\wmn_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\wmn_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\logwmhn_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\logwmhn_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thaln_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thaln_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\hippon_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\hippon_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\amygn_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\amygn_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\csfn_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\csfn_results.dta", replace

use "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\tvn_results.dta", clear
append using "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\gmn_results.dta", force
append using "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\arean_results.dta", force
append using "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thick_results.dta", force		//not normalised//
append using "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\wmn_results.dta", force
append using "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\logwmh_results.dta", force		//not normalised//
append using "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\hippon_results.dta", force
append using "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\amygn_results.dta", force
append using "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thaln_results.dta", force
append using "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\csfn_results.dta", force

gen Age = ""
replace Age = "Child"

append using "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\tvn_results.dta", force
append using "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\gmn_results.dta", force
append using "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\arean_results.dta", force
append using "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thick_results.dta", force		//not normalised//
append using "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\wmn_results.dta", force
append using "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\logwmh_results.dta", force		//not normalised//
append using "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\hippon_results.dta", force
append using "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\amygn_results.dta", force
append using "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thaln_results.dta", force
append using "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\csfn_results.dta", force

replace Age = "Adult" if Age == ""

order Age
 
save "Y:\UKB\MR\Data Files\Univariate Results.dta", replace


//Normalised to Brain Volume//

import excel "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\gmb_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\gmb_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\areab_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\areab_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thickb_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thickb_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\wmb_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\wmb_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\logwmhb_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\logwmhb_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thalb_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thalb_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\hippob_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\hippob_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\amygb_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\amygb_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\gmb_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\gmb_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\areab_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\areab_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thickb_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thickb_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\wmb_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\wmb_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\logwmhb_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\logwmhb_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thalb_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thalb_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\hippob_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\hippob_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\amygb_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\amygb_results.dta", replace

use "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\gmb_results.dta", clear
append using "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\areab_results.dta", force
append using "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thickb_results.dta", force
append using "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\wmb_results.dta", force
append using "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\logwmhb_results.dta", force
append using "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\hippob_results.dta", force
append using "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\amygb_results.dta", force
append using "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thalb_results.dta", force

gen Age = ""
replace Age = "Child"

append using "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\gmb_results.dta", force
append using "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\areab_results.dta", force
append using "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thickb_results.dta", force
append using "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\wmb_results.dta", force
append using "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\logwmhb_results.dta", force
append using "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\hippob_results.dta", force
append using "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\amygb_results.dta", force
append using "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thalb_results.dta", force

replace Age = "Adult" if Age == ""

order Age
 
save "Y:\UKB\MR\Data Files\Univariate Results Normalised to TBV.dta", replace


//Not normalised at all//

import excel "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\icv_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\icv_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\tv_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\tv_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\gm_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\gm_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\area_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\area_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thick_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thick_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\wm_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\wm_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\logwmh_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\logwmh_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thalamus_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thalamus_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\hippo_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\hippo_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\amygdala_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\amygdala_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\icv_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\icv_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\tv_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\tv_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\gm_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\gm_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\area_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\area_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thick_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thick_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\wm_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\wm_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\logwmh_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\logwmh_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thalamus_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thalamus_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\hippo_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\hippo_results.dta", replace
import excel "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\amygdala_results.xlsx", sheet("Sheet1") firstrow clear
destring Beta - p, replace
save "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\amygdala_results.dta", replace

use "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\icv_results.dta", clear
append using "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\tv_results.dta", force
append using "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\gm_results.dta", force
append using "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\area_results.dta", force
append using "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thick_results.dta", force
append using "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\wm_results.dta", force
append using "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\logwmh_results.dta", force
append using "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\hippo_results.dta", force
append using "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\amygdala_results.dta", force
append using "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thalamus_results.dta", force

gen Age = ""
replace Age = "Child"

append using "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\icv_results.dta", force
append using "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\tv_results.dta", force
append using "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\gm_results.dta", force
append using "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\area_results.dta", force
append using "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thick_results.dta", force
append using "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\wm_results.dta", force
append using "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\logwmh_results.dta", force
append using "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\hippo_results.dta", force
append using "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\amygdala_results.dta", force
append using "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thalamus_results.dta", force

replace Age = "Adult" if Age == ""

order Age
 
save "Y:\UKB\MR\Data Files\Univariate Results Not Normalised.dta", replace