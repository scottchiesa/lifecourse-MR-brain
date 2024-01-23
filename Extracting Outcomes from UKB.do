* Script Edited by Scott Chiesa to select variables for BMI-COG project, July 2023

clear all
set maxvar 50000
set more off

use "/lustre/projects/MRCLHAgen/UKB/maindata_refreshed_061123/2023_11_maindata.dta"

keep ///
n_eid ///
n_25004_2_0 ///
n_25010_2_0 ///
n_25006_2_0 ///
n_25008_2_0 ///
n_25781_2_0 ///
n_25019_2_0 ///
n_25020_2_0 ///
n_25011_2_0 ///
n_25012_2_0 ///
n_25021_2_0 ///
n_25022_2_0 /// 
n_21001_0_0 ///
n_31_0_0 ///
n_34_0_0 ///
n_21022_0_0 ///
n_22189_0_0 ///
n_6138_0_0 ///

rename n_25004_2_0 csf
rename n_25010_2_0 tv
rename n_25006_2_0 gm
rename n_25008_2_0 wm
rename n_25781_2_0 wmh

rename n_21022_0_0 age
rename n_22189_0_0 ses
rename n_31_0_0 sex
rename n_21001_0_0 bmi
gen bmi_cat = .
replace bmi_cat = 1 if bmi <= 25
replace bmi_cat = 2 if bmi > 25 & bmi <= 31.7
replace bmi_cat = 3 if bmi > 31.7 & bmi < .

gen education = .
replace education = 1 if n_6138_0_0 == 1
replace education = 2 if n_6138_0_0 == 2
replace education = 2 if n_6138_0_0 == 5
replace education = 2 if n_6138_0_0 == 6
replace education = 3 if n_6138_0_0 == 3
replace education = 4 if n_6138_0_0 == 4

gen thalamus = (n_25011_2_0 + n_25012_2_0)/2
gen hippocampus = (n_25019_2_0 + n_25020_2_0)/2
gen amygdala = (n_25021_2_0 + n_25022_2_0)/2
gen logwmh = log(wmh)

drop n_25011_2_0 n_25012_2_0 n_25019_2_0 n_25020_2_0 n_25021_2_0 n_25022_2_0

merge 1:1 n_eid using "/lustre/projects/MRCLHAgen/UKB/NMR_brainIDPs_140623/2023_06_NMR_brainIDPs.dta", ///
keepusing(n_26521_2_0 n_26721_2_0 n_26822_2_0 n_26755_2_0 n_26856_2_0)
drop _merge

gen area = (n_26721_2_0 + n_26822_2_0)/2
gen thick = (n_26755_2_0 + n_26856_2_0)/2
rename n_26521_2_0 icv
drop if icv == .
drop n_26721_2_0 n_26755_2_0 n_26822_2_0 n_26856_2_0

gen tvn = tv/icv
gen gmn = gm/icv
gen arean = area/icv
gen thickn = thick/icv
gen wmn = wm/icv
gen wmhn = wmh/icv
gen logwmhn = logwmh/icv
gen thaln = thalamus/icv
gen hippon =  hippocampus/icv
gen amygn = amygdala/icv
gen csfn = csf/icv

gen gmb = gm/tv
gen areab = area/tv
gen thickb = thick/tv
gen wmb = wm/tv
gen wmhb = wmh/tv
gen logwmhb = logwmh/tv
gen thalb = thalamus/tv
gen hippob = hippocampus/tv
gen amygb = amygdala/tv
gen csfb = csf/tv

merge 1:1 n_eid using "/lustre/projects/MRCLHAgen/UKB/participants_to_exclude/participant_exclusions_w71702_20230425.dta"
drop if _m==3
drop _m

keep n_eid tv gm area thick wm wmh logwmh thalamus hippocampus amygdala icv csf tvn gmn arean thickn wmn wmhn logwmhn logwmhn2 thaln hippon amygn csfn gmb areab thickb wmb wmhb logwmhb thalb hippob amygb csfb bmi bmi_cat age sex ses education

clonevar fid = n_eid
rename n_eid iid
order fid iid age sex ses education  bmi bmi_cat tv gm area thick wm wmh logwmh thalamus hippocampus amygdala csf tvn gmn arean thickn wmn wmhn logwmhn thaln hippon amygn csfn gmb areab thickb wmb wmhb logwmhb thalb hippob amygb csfb icv

zscore tv gm area thick wm wmh logwmh thalamus hippocampus amygdala csf tvn gmn arean thickn wmn wmhn logwmhn logwmhn2 thaln hippon amygn csfn gmb areab thickb wmb wmhb logwmhb thalb hippob amygb csfb bmi icv

drop if z_tv == .
drop if z_gm == .
drop if z_area == .
drop if z_thick == . 
drop if z_wm == .
drop if z_wmh == .
drop if z_logwmh == .
drop if z_thalamus == .
drop if z_hippocampus == .
drop if z_amygdala == .
drop if z_csf == .
drop if z_tvn == .
drop if z_gmn == .
drop if z_arean == .
drop if z_thickn == .
drop if z_wmn == .
drop if z_wmhn == .
drop if z_logwmhn == .
drop if z_logwmhn2 == .
drop if z_thaln == .
drop if z_hippon == .
drop if z_amygn == .
drop if csfn == .
drop if z_gmb == .
drop if z_areab == . 
drop if z_thickb == .
drop if z_wmb == .
drop if z_wmhb == .
drop if z_logwmhb == .
drop if z_thalb == .
drop if z_hippob == .
drop if z_amygb == .
drop if csfb == .
drop if z_icv == .

preserve
keep fid iid z_icv
export delimited using UKB_icv.txt, delimiter(tab) replace
restore

preserve 
keep fid iid z_tv 
export delimited using UKB_tv.txt, delimiter(tab) replace
restore

preserve
keep fid iid z_gm
export delimited using UKB_gm.txt, delimiter(tab) replace
restore

preserve
keep fid iid z_area
export delimited using UKB_area.txt, delimiter(tab) replace
restore

preserve
keep fid iid z_thick
export delimited using UKB_thick.txt, delimiter(tab) replace
restore

preserve
keep fid iid z_wm 
export delimited using UKB_wm.txt, delimiter(tab) replace
restore

preserve
keep fid iid z_wmh
export delimited using UKB_wmh.txt, delimiter(tab) replace
restore

preserve
keep fid iid z_logwmh
export delimited using UKB_logwmh.txt, delimiter(tab) replace
restore

preserve
keep fid iid z_thalamus
export delimited using UKB_thalamus.txt, delimiter(tab) replace
restore

preserve
keep fid iid z_hippocampus
export delimited using UKB_hippocampus.txt, delimiter(tab) replace
restore

preserve
keep fid iid z_amygdala
export delimited using UKB_amygdala.txt, delimiter(tab) replace
restore

preserve
keep fid iid z_csf
export delimited using UKB_csf.txt, delimiter(tab) replace
restore

//icv adjusted from here//

preserve
keep fid iid z_tvn
export delimited using UKB_tvn.txt, delimiter(tab) replace
restore

preserve
keep fid iid z_gmn
export delimited using UKB_gmn.txt, delimiter(tab) replace
restore

preserve
keep fid iid z_arean
export delimited using UKB_arean.txt, delimiter(tab) replace
restore

preserve
keep fid iid z_thickn
export delimited using UKB_thickn.txt, delimiter(tab) replace
restore

preserve
keep fid iid z_wmn
export delimited using UKB_wmn.txt, delimiter(tab) replace
restore

preserve
keep fid iid z_amygn
export delimited using UKB_wmhn.txt, delimiter(tab) replace
restore

preserve
keep fid iid z_logwmhn
export delimited using UKB_logwmhn.txt, delimiter(tab) replace
restore

preserve
keep fid iid z_thaln
export delimited using UKB_thaln.txt, delimiter(tab) replace
restore

preserve
keep fid iid z_hippon
export delimited using UKB_hippon.txt, delimiter(tab) replace
restore

preserve
keep fid iid z_amygn
export delimited using UKB_amygn.txt, delimiter(tab) replace
restore

preserve
keep fid iid z_csfn
export delimited using UKB_csfn.txt, delimiter(tab) replace
restore

//tbv adjusted from here//

preserve
keep fid iid z_gmb
export delimited using UKB_gmb.txt, delimiter(tab) replace
restore

preserve
keep fid iid z_areab
export delimited using UKB_areab.txt, delimiter(tab) replace
restore

preserve
keep fid iid z_thickb
export delimited using UKB_thickb.txt, delimiter(tab) replace
restore

preserve
keep fid iid z_wmb
export delimited using UKB_wmb.txt, delimiter(tab) replace
restore

preserve
keep fid iid z_logwmhb
export delimited using UKB_logwmhb.txt, delimiter(tab) replace
restore

preserve
keep fid iid z_thalb
export delimited using UKB_thalb.txt, delimiter(tab) replace
restore

preserve
keep fid iid z_hippob
export delimited using UKB_hippob.txt, delimiter(tab) replace
restore

preserve
keep fid iid z_amygb
export delimited using UKB_amygb.txt, delimiter(tab) replace
restore

preserve
keep fid iid z_csfb
export delimited using UKB_csfb.txt, delimiter(tab) replace
restore

save "/home/rmgpstc/Scratch/UKB_MRBRAIN/UKB_Selected_Variables.dta", replace


//Observational BMI Results for Figure 1//

ds z_icv z_tv z_gm z_area z_thick z_wm z_wmh z_logwmh z_thalamus z_hippocampus z_amygdala z_csf z_tvn z_gmn z_arean z_thickn z_wmn z_wmhn z_logwmhn z_thaln z_hippon z_amygn z_csfn  z_gmb z_areab z_thickb z_wmb z_wmhb z_logwmhb z_thalb z_hippob z_amygb z_csfb
local outcomes `r(varlist)'
local i = 0
foreach outcome of local outcomes {
	local i = `i' + 1
	di "Outcome: `outcome'"
	eststo: regress `outcome' bmi_cat age i.sex ses i.education
	}
	esttab _all using obs_results.csv, replace cells(b(fmt(3)) se(fmt(3)) p(fmt(3))) gaps lines nostar

