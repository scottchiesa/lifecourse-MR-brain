****************RUN THIS FILE AFTER USING R TO PERFORM ALL OF THE ANALYSES******************
*****************WILL CREATE A NICE EXCEL TABLE WITH ALL MAIN MVMR RESULTS***********************

import delimited "Y:\UKB\R Scripts and Analyses\Output Files\tvn_output_file.txt", clear
save "Y:\UKB\R Scripts and Analyses\Output Files\tvn_output_file.dta", replace
import delimited "Y:\UKB\R Scripts and Analyses\Output Files\gmn_output_file.txt", clear
save "Y:\UKB\R Scripts and Analyses\Output Files\gmn_output_file.dta", replace
import delimited "Y:\UKB\R Scripts and Analyses\Output Files\arean_output_file.txt", clear
save "Y:\UKB\R Scripts and Analyses\Output Files\arean_output_file.dta", replace
import delimited "Y:\UKB\R Scripts and Analyses\Output Files\thick_output_file.txt", clear		//note that these aren't adjusted for icv//
save "Y:\UKB\R Scripts and Analyses\Output Files\thick_output_file.dta", replace
import delimited "Y:\UKB\R Scripts and Analyses\Output Files\wmn_output_file.txt", clear
save "Y:\UKB\R Scripts and Analyses\Output Files\wmn_output_file.dta", replace
import delimited "Y:\UKB\R Scripts and Analyses\Output Files\logwmh_output_file.txt", clear		//same for these//
save "Y:\UKB\R Scripts and Analyses\Output Files\logwmh_output_file.dta", replace
import delimited "Y:\UKB\R Scripts and Analyses\Output Files\amygn_output_file.txt", clear
save "Y:\UKB\R Scripts and Analyses\Output Files\amygn_output_file.dta", replace
import delimited "Y:\UKB\R Scripts and Analyses\Output Files\hippon_output_file.txt", clear
save "Y:\UKB\R Scripts and Analyses\Output Files\hippon_output_file.dta", replace
import delimited "Y:\UKB\R Scripts and Analyses\Output Files\thaln_output_file.txt", clear
save "Y:\UKB\R Scripts and Analyses\Output Files\thaln_output_file.dta", replace

use "Y:\UKB\R Scripts and Analyses\Output Files\tvn_output_file.dta", clear
append using "Y:\UKB\R Scripts and Analyses\Output Files\gmn_output_file.dta", force
append using "Y:\UKB\R Scripts and Analyses\Output Files\arean_output_file.dta", force
append using "Y:\UKB\R Scripts and Analyses\Output Files\thick_output_file.dta", force
append using "Y:\UKB\R Scripts and Analyses\Output Files\wmn_output_file.dta", force
append using "Y:\UKB\R Scripts and Analyses\Output Files\logwmh_output_file.dta", force
append using "Y:\UKB\R Scripts and Analyses\Output Files\amygn_output_file.dta", force
append using "Y:\UKB\R Scripts and Analyses\Output Files\hippon_output_file.dta", force
append using "Y:\UKB\R Scripts and Analyses\Output Files\thaln_output_file.dta", force

keep v2 v8 v11-v14
rename v2 outcome
rename v8 age
rename v11 n
rename v12 beta
rename v13 se
rename v14 p

save "Y:\UKB\R Scripts and Analyses\Output Files\MVMR Results.dta", replace

//Now do for adjusted for tbv

import delimited "Y:\UKB\R Scripts and Analyses\Output Files\gmb_output_file.txt", clear
save "Y:\UKB\R Scripts and Analyses\Output Files\gmb_output_file.dta", replace
import delimited "Y:\UKB\R Scripts and Analyses\Output Files\areab_output_file.txt", clear
save "Y:\UKB\R Scripts and Analyses\Output Files\areab_output_file.dta", replace
import delimited "Y:\UKB\R Scripts and Analyses\Output Files\thick_output_file.txt", clear		//note that these aren't adjusted for icv//
save "Y:\UKB\R Scripts and Analyses\Output Files\thick_output_file.dta", replace
import delimited "Y:\UKB\R Scripts and Analyses\Output Files\wmb_output_file.txt", clear
save "Y:\UKB\R Scripts and Analyses\Output Files\wmb_output_file.dta", replace
import delimited "Y:\UKB\R Scripts and Analyses\Output Files\logwmh_output_file.txt", clear		//same for these//
save "Y:\UKB\R Scripts and Analyses\Output Files\logwmh_output_file.dta", replace
import delimited "Y:\UKB\R Scripts and Analyses\Output Files\amygb_output_file.txt", clear
save "Y:\UKB\R Scripts and Analyses\Output Files\amygb_output_file.dta", replace
import delimited "Y:\UKB\R Scripts and Analyses\Output Files\hippob_output_file.txt", clear
save "Y:\UKB\R Scripts and Analyses\Output Files\hippob_output_file.dta", replace
import delimited "Y:\UKB\R Scripts and Analyses\Output Files\thalb_output_file.txt", clear
save "Y:\UKB\R Scripts and Analyses\Output Files\thalb_output_file.dta", replace

use "Y:\UKB\R Scripts and Analyses\Output Files\gmb_output_file.dta", clear
append using "Y:\UKB\R Scripts and Analyses\Output Files\areab_output_file.dta", force
append using "Y:\UKB\R Scripts and Analyses\Output Files\thick_output_file.dta", force
append using "Y:\UKB\R Scripts and Analyses\Output Files\wmb_output_file.dta", force
append using "Y:\UKB\R Scripts and Analyses\Output Files\logwmh_output_file.dta", force
append using "Y:\UKB\R Scripts and Analyses\Output Files\amygb_output_file.dta", force
append using "Y:\UKB\R Scripts and Analyses\Output Files\hippob_output_file.dta", force
append using "Y:\UKB\R Scripts and Analyses\Output Files\thalb_output_file.dta", force

keep v2 v8 v11-v14
rename v2 outcome
rename v8 age
rename v11 n
rename v12 beta
rename v13 se
rename v14 p

save "Y:\UKB\R Scripts and Analyses\Output Files\MVMR Results Adjusted for TBV.dta", replace

//Now do for unadjusted analyses//

import delimited "Y:\UKB\R Scripts and Analyses\Output Files\icv_output_file.txt", clear		//icv here too//
save "Y:\UKB\R Scripts and Analyses\Output Files\icv_output_file.dta", replace
import delimited "Y:\UKB\R Scripts and Analyses\Output Files\tv_output_file.txt", clear
save "Y:\UKB\R Scripts and Analyses\Output Files\tv_output_file.dta", replace
import delimited "Y:\UKB\R Scripts and Analyses\Output Files\gm_output_file.txt", clear
save "Y:\UKB\R Scripts and Analyses\Output Files\gm_output_file.dta", replace
import delimited "Y:\UKB\R Scripts and Analyses\Output Files\area_output_file.txt", clear
save "Y:\UKB\R Scripts and Analyses\Output Files\area_output_file.dta", replace
import delimited "Y:\UKB\R Scripts and Analyses\Output Files\thick_output_file.txt", clear		//duplicate of before//
save "Y:\UKB\R Scripts and Analyses\Output Files\thick_output_file.dta", replace
import delimited "Y:\UKB\R Scripts and Analyses\Output Files\wm_output_file.txt", clear
save "Y:\UKB\R Scripts and Analyses\Output Files\wm_output_file.dta", replace
import delimited "Y:\UKB\R Scripts and Analyses\Output Files\logwmh_output_file.txt", clear		//duplicate of before//
save "Y:\UKB\R Scripts and Analyses\Output Files\logwmh_output_file.dta", replace
import delimited "Y:\UKB\R Scripts and Analyses\Output Files\amygdala_output_file.txt", clear
save "Y:\UKB\R Scripts and Analyses\Output Files\amygdala_output_file.dta", replace
import delimited "Y:\UKB\R Scripts and Analyses\Output Files\hippo_output_file.txt", clear
save "Y:\UKB\R Scripts and Analyses\Output Files\hippo_output_file.dta", replace
import delimited "Y:\UKB\R Scripts and Analyses\Output Files\thalamus_output_file.txt", clear
save "Y:\UKB\R Scripts and Analyses\Output Files\thalamus_output_file.dta", replace

use "Y:\UKB\R Scripts and Analyses\Output Files\icv_output_file.dta", clear
append using "Y:\UKB\R Scripts and Analyses\Output Files\tv_output_file.dta", force
append using "Y:\UKB\R Scripts and Analyses\Output Files\gm_output_file.dta", force
append using "Y:\UKB\R Scripts and Analyses\Output Files\area_output_file.dta", force
append using "Y:\UKB\R Scripts and Analyses\Output Files\thick_output_file.dta", force
append using "Y:\UKB\R Scripts and Analyses\Output Files\wm_output_file.dta", force
append using "Y:\UKB\R Scripts and Analyses\Output Files\logwmh_output_file.dta", force
append using "Y:\UKB\R Scripts and Analyses\Output Files\amygdala_output_file.dta", force
append using "Y:\UKB\R Scripts and Analyses\Output Files\hippo_output_file.dta", force
append using "Y:\UKB\R Scripts and Analyses\Output Files\thalamus_output_file.dta", force

keep v2 v8 v11-v14
rename v2 outcome
rename v8 age
rename v11 n
rename v12 beta
rename v13 se
rename v14 p

save "Y:\UKB\R Scripts and Analyses\Output Files\Unadjusted MVMR Results.dta", replace

**********************************REPEAT FOR BIRTHWEIGHT ANALYSES**********************************


import delimited "Y:\UKB\R Scripts and Analyses\Output Files\tvn_bw_output_file.txt", clear
save "Y:\UKB\R Scripts and Analyses\Output Files\tvn_bw_output_file.dta", replace
import delimited "Y:\UKB\R Scripts and Analyses\Output Files\gmn_bw_output_file.txt", clear
save "Y:\UKB\R Scripts and Analyses\Output Files\gmn_bw_output_file.dta", replace
import delimited "Y:\UKB\R Scripts and Analyses\Output Files\arean_bw_output_file.txt", clear
save "Y:\UKB\R Scripts and Analyses\Output Files\arean_bw_output_file.dta", replace
import delimited "Y:\UKB\R Scripts and Analyses\Output Files\thick_bw_output_file.txt", clear
save "Y:\UKB\R Scripts and Analyses\Output Files\thick_bw_output_file.dta", replace
import delimited "Y:\UKB\R Scripts and Analyses\Output Files\wmn_bw_output_file.txt", clear
save "Y:\UKB\R Scripts and Analyses\Output Files\wmn_bw_output_file.dta", replace
import delimited "Y:\UKB\R Scripts and Analyses\Output Files\logwmh_bw_output_file.txt", clear
save "Y:\UKB\R Scripts and Analyses\Output Files\logwmh_bw_output_file.dta", replace
import delimited "Y:\UKB\R Scripts and Analyses\Output Files\amygn_bw_output_file.txt", clear
save "Y:\UKB\R Scripts and Analyses\Output Files\amygn_bw_output_file.dta", replace
import delimited "Y:\UKB\R Scripts and Analyses\Output Files\hippon_bw_output_file.txt", clear
save "Y:\UKB\R Scripts and Analyses\Output Files\hippon_bw_output_file.dta", replace
import delimited "Y:\UKB\R Scripts and Analyses\Output Files\thaln_bw_output_file.txt", clear
save "Y:\UKB\R Scripts and Analyses\Output Files\thaln_bw_output_file.dta", replace

use "Y:\UKB\R Scripts and Analyses\Output Files\tvn_bw_output_file.dta", clear
append using "Y:\UKB\R Scripts and Analyses\Output Files\gmn_bw_output_file.dta", force
append using "Y:\UKB\R Scripts and Analyses\Output Files\arean_bw_output_file.dta", force
append using "Y:\UKB\R Scripts and Analyses\Output Files\thick_bw_output_file.dta", force
append using "Y:\UKB\R Scripts and Analyses\Output Files\wmn_bw_output_file.dta", force
append using "Y:\UKB\R Scripts and Analyses\Output Files\logwmh_bw_output_file.dta", force
append using "Y:\UKB\R Scripts and Analyses\Output Files\amygn_bw_output_file.dta", force
append using "Y:\UKB\R Scripts and Analyses\Output Files\hippon_bw_output_file.dta", force
append using "Y:\UKB\R Scripts and Analyses\Output Files\thaln_bw_output_file.dta", force

keep v2 v8 v11-v14
rename v2 outcome
rename v8 age
rename v11 n
rename v12 beta
rename v13 se
rename v14 p

save "Y:\UKB\R Scripts and Analyses\Output Files\MVMR Birthweight Results.dta", replace