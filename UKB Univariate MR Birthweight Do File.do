********************************************************************************************************************************************************************************************
*******************************************************************************UKB LIFECOURSE BODY SIZE - BRAIN PHENOTYPES MENDELIAN RANDOMIZATION********************************************
**************************************************************************************************SCOTT CHIESA 27/07/23*********************************************************************
********************************************************************************************************************************************************************************************

********************************************************************************************************************************************************************************************
***************************************************************************************Importing and Preparing all Summary Data for Analysis************************************************
********************************************************************************************************************************************************************************************

//Export all outcome phenotype .txt files from Myriad into local folder//

//First delete any pre-existing files in subfolders from previous analysis runs as otherwise the later parts won't work//

cd "Y:\UKB\MR\Data Files\Birthweight Outcomes"
local myfilelist : dir . files "*.dta"
foreach file of local myfilelist{
	erase "`file'"
}

cd "Y:\UKB\MR\Data Files\Child BW Outcomes"
local myfilelist : dir . files "*.dta"
foreach file of local myfilelist{
	erase "`file'"
}

//Now convert all summary data .txt files from Myriad into .dta files//	

cd "Y:\UKB\MR\Individual Outcome Files from Myriad for BW Analysis"	
local myfilelist : dir . files "*.txt"
foreach file of local myfilelist {
	di "`file'"		//display file name//
	import delimited "`file'", clear delimiter(tab) varnames(1)
	local outfile = subinstr("`file'", ".txt", "", .)
	save "`outfile'", replace
}

//Move all newly created birthweight and child bw .dta files to their respective folders for later//

mvfiles , infolder(".") outfolder("Y:\UKB\MR\Data Files\Birthweight Outcomes") match("bw_*.dta") erase
mvfiles , infolder(".") outfolder("Y:\UKB\MR\Data Files\Child BW Outcomes") match("child_bw_*.dta") erase


*****************************Preparing Exposures for Merging******************************************

//Add prefix to all birthweight exposure variables//

use "Y:\UKB\MR\Data Files\Birthweight Exposures\BMI_birthweight.dta", clear			
foreach var of varlist * {
	if "`var'" != "id" {
		rename `var' BMI_`var'
	}
	save "Y:\UKB\MR\Data Files\Birthweight Exposures\birthweight_exp.dta", replace
}

//Add prefix to all child bw exposure variables//

use "Y:\UKB\MR\Data Files\Child BW Exposures\BMI_child_bw.dta", clear
foreach var of varlist * {
	if "`var'" != "id" {
		rename `var' BMI_`var'
	}
	save "Y:\UKB\MR\Data Files\Child BW Exposures\child_bw_exp.dta", replace
}


*****************************Preparing Outcomes for Merging******************************************

//Add prefix to all birthweight outcome variables//

cd "Y:\UKB\MR\Data Files\Birthweight Outcomes"
local myfilelist : dir . files "*.dta"
foreach file of local myfilelist {
	use "`file'", clear
		foreach var of varlist * {
			if "`var'" != "id" {
				local newname = subinstr("`file'",".dta","",.)
		rename `var' `newname'_`var'	
		}
	}
save "`file'", replace
}

//Add prefix to all child bw outcome variables//

cd "Y:\UKB\MR\Data Files\Child BW Outcomes"
local myfilelist : dir . files "*.dta"
foreach file of local myfilelist {
	use "`file'", clear
		foreach var of varlist * {
			if "`var'" != "id" {
				local newname = subinstr("`file'",".dta","",.)
		rename `var' `newname'_`var'	
		}
	}
save "`file'", replace
}


************Merging Exposures and Outcomes for Both Birthweight and Child BW Effects******************

//Merge birthweight files//

cd "Y:\UKB\MR\Data Files\Birthweight Outcomes"
local myfilelist : dir . files "*.dta"
use "Y:\UKB\MR\Data Files\Birthweight Exposures\birthweight_exp.dta"
foreach file of local myfilelist {
	merge 1:1 id using `file'
	drop _merge
	save "Y:\UKB\MR\Data Files\Birthweight Outcomes\Birthweight Effects File.dta", replace
}

//Merge child bw files//

cd "Y:\UKB\MR\Data Files\Child BW Outcomes"
local myfilelist : dir . files "*.dta"
use "Y:\UKB\MR\Data Files\Child BW Exposures\child_bw_exp.dta"
foreach file of local myfilelist {
	merge 1:1 id using `file'
	drop _merge
	save "Y:\UKB\MR\Data Files\Child BW Outcomes\Child BW Effects File.dta", replace
}

//Drop missing rows in birthweight data which are messing up graphs//

use "Y:\UKB\MR\Data Files\Birthweight Outcomes\Birthweight Effects File.dta", clear
drop in 20
drop in 53
drop in 61
save "Y:\UKB\MR\Data Files\Birthweight Outcomes\Birthweight Effects File.dta", replace

/*Drop two missing rows in child bw data which is messing up graphs//

use "Y:\UKB\MR\Data Files\Child BW Outcomes\Child BW Effects File.dta", clear
drop in 171 
drop in 289
save "Y:\UKB\MR\Data Files\Adult Outcomes\Adult Effects File.dta", replace*/


********************************************************************************************************************************************************************************************
***********************************************************************************PERFORMING MENDELIAN RANDOMIZATIONS**********************************************************************
**************************************************************************UNIVARIATE ONLY - MVMR RUN IN SEPARATE R SCRIPT*******************************************************************
********************************************************************************************************************************************************************************************

**************************************************ALL OUTCOMES INDEXED TO ICV EXCEPT CORTICAL THICKNESS AND WHITE MATTER HYPERINTENSITIES***************************************************
*********************************************** EACH STARTS WITH BIRTHWEIGHT SNPS VS OUTCOMES THEN REPEATS WITH CHILDHOOD SNPS VS OUTCOMES**************************************************


***************************************************************************************BIRTHWEIGHT EFFECTS**********************************************************************************

log using "Y:\UKB\Updated_Birthweight_MR", replace

//TOTAL BRAIN VOLUME NORMALISED FOR INTRACRANIAL VOLUME//

use "Y:\UKB\MR\Data Files\Birthweight Outcomes\Birthweight Effects File.dta", clear

* IVW method (multiplicative random effects)
tab bw_tvn_beta
local snp_num = `r(N)'
mregger bw_tvn_beta BMI_beta [aw=1/(bw_tvn_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger bw_tvn_beta BMI_beta [aw=1/(bw_tvn_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian bw_tvn_beta bw_tvn_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal bw_tvn_beta bw_tvn_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\tvn_results", modify
	putexcel A1="Outcome" B1="No. SNPs" C1="Method" D1="Beta" E1="SE" F1="p"
	local x = `x'+1
	putexcel A`x' = "Total Brain Volume" 
	putexcel B`x' = `snp_num'
	putexcel C`x' = "`method_`estimator''"
	putexcel D`x' = "`beta_`estimator''"
	putexcel E`x' = "`se_`estimator''"
	putexcel F`x' = "`p_`estimator''"
}

* Is there evidence of heterogeneity in the genetic effects?
mregger bw_tvn_beta BMI_beta [aw=1/(bw_tvn_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger bw_tvn_beta BMI_beta [aw=1/(bw_tvn_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter bw_tvn_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\tvn_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest bw_tvn_beta bw_tvn_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\tvn_forest.png", replace

* Asymmetry
mrfunnel bw_tvn_beta bw_tvn_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\tvn_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = bw_tvn_beta/BMI_beta
gen wald_ratio_se = bw_tvn_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\tvn_loo.png", replace


//GREY MATTER VOLUME NORMALISED FOR INTRACRANIAL VOLUME//

use "Y:\UKB\MR\Data Files\Birthweight Outcomes\Birthweight Effects File.dta", clear

* IVW method (multiplicative random effects)
tab bw_gmn_beta
local snp_num = `r(N)'
mregger bw_gmn_beta BMI_beta [aw=1/(bw_gmn_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger bw_gmn_beta BMI_beta [aw=1/(bw_gmn_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian bw_gmn_beta bw_gmn_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal bw_gmn_beta bw_gmn_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\gmn_results", modify
	putexcel A1="Outcome" B1="No. SNPs" C1="Method" D1="Beta" E1="SE" F1="p"
	local x = `x'+1
	putexcel A`x' = "Grey Matter Volume" 
	putexcel B`x' = `snp_num'
	putexcel C`x' = "`method_`estimator''"
	putexcel D`x' = "`beta_`estimator''"
	putexcel E`x' = "`se_`estimator''"
	putexcel F`x' = "`p_`estimator''"
}

* Is there evidence of heterogeneity in the genetic effects?
mregger bw_gmn_beta BMI_beta [aw=1/(bw_gmn_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger bw_gmn_beta BMI_beta [aw=1/(bw_gmn_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter bw_gmn_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\gmn_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest bw_gmn_beta bw_gmn_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\gmn_forest.png", replace

* Asymmetry
mrfunnel bw_gmn_beta bw_gmn_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\gmn_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = bw_gmn_beta/BMI_beta
gen wald_ratio_se = bw_gmn_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\gmn_loo.png", replace


//CORTICAL SURFACE AREA NORMALISED FOR INTRACRANIAL VOLUME//

use "Y:\UKB\MR\Data Files\Birthweight Outcomes\Birthweight Effects File.dta", clear

* IVW method (multiplicative random effects)
tab bw_arean_beta
local snp_num = `r(N)'
mregger bw_arean_beta BMI_beta [aw=1/(bw_arean_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger bw_arean_beta BMI_beta [aw=1/(bw_arean_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian bw_arean_beta bw_arean_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal bw_arean_beta bw_arean_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\arean_results", modify
	putexcel A1="Outcome" B1="No. SNPs" C1="Method" D1="Beta" E1="SE" F1="p"
	local x = `x'+1
	putexcel A`x' = "Cortical Surface Area" 
	putexcel B`x' = `snp_num'
	putexcel C`x' = "`method_`estimator''"
	putexcel D`x' = "`beta_`estimator''"
	putexcel E`x' = "`se_`estimator''"
	putexcel F`x' = "`p_`estimator''"
}

* Is there evidence of heterogeneity in the genetic effects?
mregger bw_arean_beta BMI_beta [aw=1/(bw_arean_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger bw_arean_beta BMI_beta [aw=1/(bw_arean_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter bw_arean_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\arean_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest bw_arean_beta bw_arean_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\arean_forest.png", replace

* Asymmetry
mrfunnel bw_arean_beta bw_arean_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\arean_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = bw_arean_beta/BMI_beta
gen wald_ratio_se = bw_arean_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\arean_loo.png", replace


//CORTICAL THICKNESS NORMALISED FOR INTRACRANIAL VOLUME//

*****NOTE THAT THIS ISN'T NORMALISED AS NOT APPROPRIATE HERE*****

use "Y:\UKB\MR\Data Files\Birthweight Outcomes\Birthweight Effects File.dta", clear

* IVW method (multiplicative random effects)
tab bw_thick_beta
local snp_num = `r(N)'
mregger bw_thick_beta BMI_beta [aw=1/(bw_thick_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger bw_thick_beta BMI_beta [aw=1/(bw_thick_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian bw_thick_beta bw_thick_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal bw_thick_beta bw_thick_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\thick_results", modify
	putexcel A1="Outcome" B1="No. SNPs" C1="Method" D1="Beta" E1="SE" F1="p"
	local x = `x'+1
	putexcel A`x' = "Cortical Thickness" 
	putexcel B`x' = `snp_num'
	putexcel C`x' = "`method_`estimator''"
	putexcel D`x' = "`beta_`estimator''"
	putexcel E`x' = "`se_`estimator''"
	putexcel F`x' = "`p_`estimator''"
}

* Is there evidence of heterogeneity in the genetic effects?
mregger bw_thick_beta BMI_beta [aw=1/(bw_thick_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger bw_thick_beta BMI_beta [aw=1/(bw_thick_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter bw_thick_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\thick_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest bw_thick_beta bw_thick_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\thick_forest.png", replace

* Asymmetry
mrfunnel bw_thick_beta bw_thick_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\thick_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = bw_thick_beta/BMI_beta
gen wald_ratio_se = bw_thick_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\thick_loo.png", replace


//WHITE MATTER VOLUME NORMALISED FOR INTRACRANIAL VOLUME//

use "Y:\UKB\MR\Data Files\Birthweight Outcomes\Birthweight Effects File.dta", clear

* IVW method (multiplicative random effects)
tab bw_wmn_beta
local snp_num = `r(N)'
mregger bw_wmn_beta BMI_beta [aw=1/(bw_wmn_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger bw_wmn_beta BMI_beta [aw=1/(bw_wmn_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian bw_wmn_beta bw_wmn_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal bw_wmn_beta bw_wmn_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\wmn_results", modify
	putexcel A1="Outcome" B1="No. SNPs" C1="Method" D1="Beta" E1="SE" F1="p"
	local x = `x'+1
	putexcel A`x' = "White Matter Volume" 
	putexcel B`x' = `snp_num'
	putexcel C`x' = "`method_`estimator''"
	putexcel D`x' = "`beta_`estimator''"
	putexcel E`x' = "`se_`estimator''"
	putexcel F`x' = "`p_`estimator''"
}

* Is there evidence of heterogeneity in the genetic effects?
mregger bw_wmn_beta BMI_beta [aw=1/(bw_wmn_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger bw_wmn_beta BMI_beta [aw=1/(bw_wmn_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter bw_wmn_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\wmn_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest bw_wmn_beta bw_wmn_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\wmn_forest.png", replace

* Asymmetry
mrfunnel bw_wmn_beta bw_wmn_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\wmn_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = bw_wmn_beta/BMI_beta
gen wald_ratio_se = bw_wmn_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\wmn_loo.png", replace


//LOG-TRANSFORMED WHITE MATTER HYPERINTENSITIIES NORMALISED FOR INTRACRANIAL VOLUME//

*****NOTE THAT THIS ISN'T NORMALISED AS NOT APPROPRIATE HERE*****

use "Y:\UKB\MR\Data Files\Birthweight Outcomes\Birthweight Effects File.dta", clear

* IVW method (multiplicative random effects)
tab bw_logwmh_beta
local snp_num = `r(N)'
mregger bw_logwmh_beta BMI_beta [aw=1/(bw_logwmh_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger bw_logwmh_beta BMI_beta [aw=1/(bw_logwmh_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian bw_logwmh_beta bw_logwmh_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal bw_logwmh_beta bw_logwmh_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\logwmh_results", modify
	putexcel A1="Outcome" B1="No. SNPs" C1="Method" D1="Beta" E1="SE" F1="p"
	local x = `x'+1
	putexcel A`x' = "Log WMH" 
	putexcel B`x' = `snp_num'
	putexcel C`x' = "`method_`estimator''"
	putexcel D`x' = "`beta_`estimator''"
	putexcel E`x' = "`se_`estimator''"
	putexcel F`x' = "`p_`estimator''"
}

* Is there evidence of heterogeneity in the genetic effects?
mregger bw_logwmh_beta BMI_beta [aw=1/(bw_logwmh_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger bw_logwmh_beta BMI_beta [aw=1/(bw_logwmh_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter bw_logwmh_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\logwmh_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest bw_logwmh_beta bw_logwmh_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\logwmh_forest.png", replace

* Asymmetry
mrfunnel bw_logwmh_beta bw_logwmh_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\logwmh_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = bw_logwmh_beta/BMI_beta
gen wald_ratio_se = bw_logwmh_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\logwmh_loo.png", replace


//THALAMUS NORMALISED FOR INTRACRANIAL VOLUME//

use "Y:\UKB\MR\Data Files\Birthweight Outcomes\Birthweight Effects File.dta", clear

* IVW method (multiplicative random effects)
tab bw_thaln_beta
local snp_num = `r(N)'
mregger bw_thaln_beta BMI_beta [aw=1/(bw_thaln_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger bw_thaln_beta BMI_beta [aw=1/(bw_thaln_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian bw_thaln_beta bw_thaln_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal bw_thaln_beta bw_thaln_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\thaln_results", modify
	putexcel A1="Outcome" B1="No. SNPs" C1="Method" D1="Beta" E1="SE" F1="p"
	local x = `x'+1
	putexcel A`x' = "Thalamus" 
	putexcel B`x' = `snp_num'
	putexcel C`x' = "`method_`estimator''"
	putexcel D`x' = "`beta_`estimator''"
	putexcel E`x' = "`se_`estimator''"
	putexcel F`x' = "`p_`estimator''"
}

* Is there evidence of heterogeneity in the genetic effects?
mregger bw_thaln_beta BMI_beta [aw=1/(bw_thaln_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger bw_thaln_beta BMI_beta [aw=1/(bw_thaln_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter bw_thaln_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\thaln_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest bw_thaln_beta bw_thaln_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\thaln_forest.png", replace

* Asymmetry
mrfunnel bw_thaln_beta bw_thaln_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\thaln_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = bw_thaln_beta/BMI_beta
gen wald_ratio_se = bw_thaln_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\thaln_loo.png", replace


//HIPPOCAMPUS NORMALISED FOR INTRACRANIAL VOLUME//

use "Y:\UKB\MR\Data Files\Birthweight Outcomes\Birthweight Effects File.dta", clear

* IVW method (multiplicative random effects)
tab bw_hippon_beta
local snp_num = `r(N)'
mregger bw_hippon_beta BMI_beta [aw=1/(bw_hippon_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger bw_hippon_beta BMI_beta [aw=1/(bw_hippon_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian bw_hippon_beta bw_hippon_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal bw_hippon_beta bw_hippon_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\hippon_results", modify
	putexcel A1="Outcome" B1="No. SNPs" C1="Method" D1="Beta" E1="SE" F1="p"
	local x = `x'+1
	putexcel A`x' = "Hippocampus" 
	putexcel B`x' = `snp_num'
	putexcel C`x' = "`method_`estimator''"
	putexcel D`x' = "`beta_`estimator''"
	putexcel E`x' = "`se_`estimator''"
	putexcel F`x' = "`p_`estimator''"
}

* Is there evidence of heterogeneity in the genetic effects?
mregger bw_hippon_beta BMI_beta [aw=1/(bw_hippon_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger bw_hippon_beta BMI_beta [aw=1/(bw_hippon_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter bw_hippon_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\hippon_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest bw_hippon_beta bw_hippon_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\hippon_forest.png", replace

* Asymmetry
mrfunnel bw_hippon_beta bw_hippon_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\hippon_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = bw_hippon_beta/BMI_beta
gen wald_ratio_se = bw_hippon_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\hippon_loo.png", replace


//AMYGDALA NORMALISED FOR INTRACRANIAL VOLUME//

use "Y:\UKB\MR\Data Files\Birthweight Outcomes\Birthweight Effects File.dta", clear

* IVW method (multiplicative random effects)
tab bw_amygn_beta
local snp_num = `r(N)'
mregger bw_amygn_beta BMI_beta [aw=1/(bw_amygn_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger bw_amygn_beta BMI_beta [aw=1/(bw_amygn_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian bw_amygn_beta bw_amygn_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal bw_amygn_beta bw_amygn_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\amygn_results", modify
	putexcel A1="Outcome" B1="No. SNPs" C1="Method" D1="Beta" E1="SE" F1="p"
	local x = `x'+1
	putexcel A`x' = "Amygdala" 
	putexcel B`x' = `snp_num'
	putexcel C`x' = "`method_`estimator''"
	putexcel D`x' = "`beta_`estimator''"
	putexcel E`x' = "`se_`estimator''"
	putexcel F`x' = "`p_`estimator''"
}

* Is there evidence of heterogeneity in the genetic effects?
mregger bw_amygn_beta BMI_beta [aw=1/(bw_amygn_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger bw_amygn_beta BMI_beta [aw=1/(bw_amygn_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter bw_amygn_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\amygn_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest bw_amygn_beta bw_amygn_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\amygn_forest.png", replace

* Asymmetry
mrfunnel bw_amygn_beta bw_amygn_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\amygn_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = bw_amygn_beta/BMI_beta
gen wald_ratio_se = bw_amygn_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Birthweight Outcomes\Results Files\amygn_loo.png", replace




*********************************************************************************CHILDHOOD EFFECTS*******************************************************************************


//TOTAL BRAIN VOLUME NORMALISED FOR INTRACRANIAL VOLUME//

use "Y:\UKB\MR\Data Files\Child BW Outcomes\Child BW Effects File.dta", clear

* IVW method (multiplicative random effects)
tab child_bw_tvn_beta
local snp_num = `r(N)'
mregger child_bw_tvn_beta BMI_beta [aw=1/(child_bw_tvn_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger child_bw_tvn_beta BMI_beta [aw=1/(child_bw_tvn_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian child_bw_tvn_beta child_bw_tvn_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal child_bw_tvn_beta child_bw_tvn_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\tvn_results", modify
	putexcel A1="Outcome" B1="No. SNPs" C1="Method" D1="Beta" E1="SE" F1="p"
	local x = `x'+1
	putexcel A`x' = "Total Brain Volume" 
	putexcel B`x' = `snp_num'
	putexcel C`x' = "`method_`estimator''"
	putexcel D`x' = "`beta_`estimator''"
	putexcel E`x' = "`se_`estimator''"
	putexcel F`x' = "`p_`estimator''"
}

* Is there evidence of heterogeneity in the genetic effects?
mregger child_bw_tvn_beta BMI_beta [aw=1/(child_bw_tvn_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger child_bw_tvn_beta BMI_beta [aw=1/(child_bw_tvn_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter child_bw_tvn_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\tvn_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest child_bw_tvn_beta child_bw_tvn_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\tvn_forest.png", replace

* Asymmetry
mrfunnel child_bw_tvn_beta child_bw_tvn_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\tvn_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = child_bw_tvn_beta/BMI_beta
gen wald_ratio_se = child_bw_tvn_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\tvn_loo.png", replace


//GREY MATTER VOLUME NORMALISED FOR INTRACRAIAL VOLUME//

use "Y:\UKB\MR\Data Files\Child BW Outcomes\Child BW Effects File.dta", clear

* IVW method (multiplicative random effects)
tab child_bw_gmn_beta
local snp_num = `r(N)'
mregger child_bw_gmn_beta BMI_beta [aw=1/(child_bw_gmn_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger child_bw_gmn_beta BMI_beta [aw=1/(child_bw_gmn_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian child_bw_gmn_beta child_bw_gmn_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal child_bw_gmn_beta child_bw_gmn_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\gmn_results", modify
	putexcel A1="Outcome" B1="No. SNPs" C1="Method" D1="Beta" E1="SE" F1="p"
	local x = `x'+1
	putexcel A`x' = "Grey Matter Volume" 
	putexcel B`x' = `snp_num'
	putexcel C`x' = "`method_`estimator''"
	putexcel D`x' = "`beta_`estimator''"
	putexcel E`x' = "`se_`estimator''"
	putexcel F`x' = "`p_`estimator''"
}

* Is there evidence of heterogeneity in the genetic effects?
mregger child_bw_gmn_beta BMI_beta [aw=1/(child_bw_gmn_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger child_bw_gmn_beta BMI_beta [aw=1/(child_bw_gmn_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter child_bw_gmn_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\gmn_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest child_bw_gmn_beta child_bw_gmn_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\gmn_forest.png", replace

* Asymmetry
mrfunnel child_bw_gmn_beta child_bw_gmn_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\gmn_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = child_bw_gmn_beta/BMI_beta
gen wald_ratio_se = child_bw_gmn_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\gmn_loo.png", replace


//CORTICAL SURFACE AREA NORMALISED FOR INTRACRANIAL VOLUME//

use "Y:\UKB\MR\Data Files\Child BW Outcomes\Child BW Effects File.dta", clear

* IVW method (multiplicative random effects)
tab child_bw_arean_beta
local snp_num = `r(N)'
mregger child_bw_arean_beta BMI_beta [aw=1/(child_bw_arean_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger child_bw_arean_beta BMI_beta [aw=1/(child_bw_arean_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian child_bw_arean_beta child_bw_arean_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal child_bw_arean_beta child_bw_arean_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\arean_results", modify
	putexcel A1="Outcome" B1="No. SNPs" C1="Method" D1="Beta" E1="SE" F1="p"
	local x = `x'+1
	putexcel A`x' = "Cortical Surface Area" 
	putexcel B`x' = `snp_num'
	putexcel C`x' = "`method_`estimator''"
	putexcel D`x' = "`beta_`estimator''"
	putexcel E`x' = "`se_`estimator''"
	putexcel F`x' = "`p_`estimator''"
}

* Is there evidence of heterogeneity in the genetic effects?
mregger child_bw_arean_beta BMI_beta [aw=1/(child_bw_arean_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger child_bw_arean_beta BMI_beta [aw=1/(child_bw_arean_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter child_bw_arean_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\arean_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest child_bw_arean_beta child_bw_arean_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\arean_forest.png", replace

* Asymmetry
mrfunnel child_bw_arean_beta child_bw_arean_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\arean_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = child_bw_arean_beta/BMI_beta
gen wald_ratio_se = child_bw_arean_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\arean_loo.png", replace


//CORTICAL THICKNESS NORMALISED FOR INTRACRANIAL VOLUME//

*****NOTE THAT THIS ISN'T USED IN ANALYSIS AS NORMALISATION NOT APPROPRIATE HERE*****

use "Y:\UKB\MR\Data Files\Child BW Outcomes\Child BW Effects File.dta", clear

* IVW method (multiplicative random effects)
tab child_bw_thickn_beta
local snp_num = `r(N)'
mregger child_bw_thickn_beta BMI_beta [aw=1/(child_bw_thickn_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger child_bw_thickn_beta BMI_beta [aw=1/(child_bw_thickn_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian child_bw_thickn_beta child_bw_thickn_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal child_bw_thickn_beta child_bw_thickn_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\thickn_results", modify
	putexcel A1="Outcome" B1="No. SNPs" C1="Method" D1="Beta" E1="SE" F1="p"
	local x = `x'+1
	putexcel A`x' = "Cortical Thickness" 
	putexcel B`x' = `snp_num'
	putexcel C`x' = "`method_`estimator''"
	putexcel D`x' = "`beta_`estimator''"
	putexcel E`x' = "`se_`estimator''"
	putexcel F`x' = "`p_`estimator''"
}

* Is there evidence of heterogeneity in the genetic effects?
mregger child_bw_thickn_beta BMI_beta [aw=1/(child_bw_thickn_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger child_bw_thickn_beta BMI_beta [aw=1/(child_bw_thickn_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter child_bw_thickn_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\thickn_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest child_bw_thickn_beta child_bw_thickn_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\thickn_forest.png", replace

* Asymmetry
mrfunnel child_bw_thickn_beta child_bw_thickn_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\thickn_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = child_bw_thickn_beta/BMI_beta
gen wald_ratio_se = child_bw_thickn_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\thickn_loo.png", replace


//WHITE MATTER VOLUME NORMALISED FOR INTRACRANIAL VOLUME/

use "Y:\UKB\MR\Data Files\Child BW Outcomes\Child BW Effects File.dta", clear

* IVW method (multiplicative random effects)
tab child_bw_wmn_beta
local snp_num = `r(N)'
mregger child_bw_wmn_beta BMI_beta [aw=1/(child_bw_wmn_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger child_bw_wmn_beta BMI_beta [aw=1/(child_bw_wmn_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian child_bw_wmn_beta child_bw_wmn_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal child_bw_wmn_beta child_bw_wmn_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\wmn_results", modify
	putexcel A1="Outcome" B1="No. SNPs" C1="Method" D1="Beta" E1="SE" F1="p"
	local x = `x'+1
	putexcel A`x' = "White Matter Volume" 
	putexcel B`x' = `snp_num'
	putexcel C`x' = "`method_`estimator''"
	putexcel D`x' = "`beta_`estimator''"
	putexcel E`x' = "`se_`estimator''"
	putexcel F`x' = "`p_`estimator''"
}

* Is there evidence of heterogeneity in the genetic effects?
mregger child_bw_wmn_beta BMI_beta [aw=1/(child_bw_wmn_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger child_bw_wmn_beta BMI_beta [aw=1/(child_bw_wmn_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter child_bw_wmn_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\wmn_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest child_bw_wmn_beta child_bw_wmn_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\wmn_forest.png", replace

* Asymmetry
mrfunnel child_bw_wmn_beta child_bw_wmn_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\wmn_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = child_bw_wmn_beta/BMI_beta
gen wald_ratio_se = child_bw_wmn_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\wmn_loo.png", replace


//LOG-TRANSFORMED WHITE MATTER HYPERINTENSITIES NORMALISED FOR INTRACRANIAL VOLUME//

*****NOTE THAT THIS ISN'T USED IN ANALYSIS AS NORMALISATION NOT APPROPRIATE HERE*****

use "Y:\UKB\MR\Data Files\Child BW Outcomes\Child BW Effects File.dta", clear

* IVW method (multiplicative random effects)
tab child_bw_logwmhn_beta
local snp_num = `r(N)'
mregger child_bw_logwmhn_beta BMI_beta [aw=1/(child_bw_logwmhn_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger child_bw_logwmhn_beta BMI_beta [aw=1/(child_bw_logwmhn_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian child_bw_logwmhn_beta child_bw_logwmhn_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal child_bw_logwmhn_beta child_bw_logwmhn_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\logwmhn_results", modify
	putexcel A1="Outcome" B1="No. SNPs" C1="Method" D1="Beta" E1="SE" F1="p"
	local x = `x'+1
	putexcel A`x' = "Log WMH" 
	putexcel B`x' = `snp_num'
	putexcel C`x' = "`method_`estimator''"
	putexcel D`x' = "`beta_`estimator''"
	putexcel E`x' = "`se_`estimator''"
	putexcel F`x' = "`p_`estimator''"
}

* Is there evidence of heterogeneity in the genetic effects?
mregger child_bw_logwmhn_beta BMI_beta [aw=1/(child_bw_logwmhn_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger child_bw_logwmhn_beta BMI_beta [aw=1/(child_bw_logwmhn_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter child_bw_logwmhn_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\logwmhn_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest child_bw_logwmhn_beta child_bw_logwmhn_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\logwmhn_forest.png", replace

* Asymmetry
mrfunnel child_bw_logwmhn_beta child_bw_logwmhn_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\logwmhn_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = child_bw_logwmhn_beta/BMI_beta
gen wald_ratio_se = child_bw_logwmhn_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\logwmhn_loo.png", replace


//THALAMUS NORMALISED FOR INTRACRANIAL VOLUME//


use "Y:\UKB\MR\Data Files\Child BW Outcomes\Child BW Effects File.dta", clear

* IVW method (multiplicative random effects)
tab child_bw_thaln_beta
local snp_num = `r(N)'
mregger child_bw_thaln_beta BMI_beta [aw=1/(child_bw_thaln_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger child_bw_thaln_beta BMI_beta [aw=1/(child_bw_thaln_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian child_bw_thaln_beta child_bw_thaln_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal child_bw_thaln_beta child_bw_thaln_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\thaln_results", modify
	putexcel A1="Outcome" B1="No. SNPs" C1="Method" D1="Beta" E1="SE" F1="p"
	local x = `x'+1
	putexcel A`x' = "Thalamus" 
	putexcel B`x' = `snp_num'
	putexcel C`x' = "`method_`estimator''"
	putexcel D`x' = "`beta_`estimator''"
	putexcel E`x' = "`se_`estimator''"
	putexcel F`x' = "`p_`estimator''"
}

* Is there evidence of heterogeneity in the genetic effects?
mregger child_bw_thaln_beta BMI_beta [aw=1/(child_bw_thaln_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger child_bw_thaln_beta BMI_beta [aw=1/(child_bw_thaln_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter child_bw_thaln_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\thaln_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest child_bw_thaln_beta child_bw_thaln_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\thaln_forest.png", replace

* Asymmetry
mrfunnel child_bw_thaln_beta child_bw_thaln_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\thaln_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = child_bw_thaln_beta/BMI_beta
gen wald_ratio_se = child_bw_thaln_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\thaln_loo.png", replace


//HIPPOCAMPUS NORMALISED FOR INTRACRANIAL VOLUME//


use "Y:\UKB\MR\Data Files\Child BW Outcomes\Child BW Effects File.dta", clear

* IVW method (multiplicative random effects)
tab child_bw_hippon_beta
local snp_num = `r(N)'
mregger child_bw_hippon_beta BMI_beta [aw=1/(child_bw_hippon_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger child_bw_hippon_beta BMI_beta [aw=1/(child_bw_hippon_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian child_bw_hippon_beta child_bw_hippon_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal child_bw_hippon_beta child_bw_hippon_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\hippon_results", modify
	putexcel A1="Outcome" B1="No. SNPs" C1="Method" D1="Beta" E1="SE" F1="p"
	local x = `x'+1
	putexcel A`x' = "Hippocampus" 
	putexcel B`x' = `snp_num'
	putexcel C`x' = "`method_`estimator''"
	putexcel D`x' = "`beta_`estimator''"
	putexcel E`x' = "`se_`estimator''"
	putexcel F`x' = "`p_`estimator''"
}

* Is there evidence of heterogeneity in the genetic effects?
mregger child_bw_hippon_beta BMI_beta [aw=1/(child_bw_hippon_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger child_bw_hippon_beta BMI_beta [aw=1/(child_bw_hippon_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter child_bw_hippon_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\hippon_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest child_bw_hippon_beta child_bw_hippon_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\hippon_forest.png", replace

* Asymmetry
mrfunnel child_bw_hippon_beta child_bw_hippon_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\hippon_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = child_bw_hippon_beta/BMI_beta
gen wald_ratio_se = child_bw_hippon_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\hippon_loo.png", replace


//AMYGDALA NORMALISED FOR INTRACRANIAL VOLUME//


use "Y:\UKB\MR\Data Files\Child BW Outcomes\Child BW Effects File.dta", clear

* IVW method (multiplicative random effects)
tab child_bw_amygn_beta
local snp_num = `r(N)'
mregger child_bw_amygn_beta BMI_beta [aw=1/(child_bw_amygn_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger child_bw_amygn_beta BMI_beta [aw=1/(child_bw_amygn_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian child_bw_amygn_beta child_bw_amygn_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal child_bw_amygn_beta child_bw_amygn_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\amygn_results", modify
	putexcel A1="Outcome" B1="No. SNPs" C1="Method" D1="Beta" E1="SE" F1="p"
	local x = `x'+1
	putexcel A`x' = "Amygdala" 
	putexcel B`x' = `snp_num'
	putexcel C`x' = "`method_`estimator''"
	putexcel D`x' = "`beta_`estimator''"
	putexcel E`x' = "`se_`estimator''"
	putexcel F`x' = "`p_`estimator''"
}

* Is there evidence of heterogeneity in the genetic effects?
mregger child_bw_amygn_beta BMI_beta [aw=1/(child_bw_amygn_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger child_bw_amygn_beta BMI_beta [aw=1/(child_bw_amygn_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter child_bw_amygn_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\amygn_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest child_bw_amygn_beta child_bw_amygn_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\amygn_forest.png", replace

* Asymmetry
mrfunnel child_bw_amygn_beta child_bw_amygn_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\amygn_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = child_bw_amygn_beta/BMI_beta
gen wald_ratio_se = child_bw_amygn_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Child BW Outcomes\Results Files\amygn_loo.png", replace

log close
