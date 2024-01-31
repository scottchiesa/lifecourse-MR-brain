*****************************************************************************************************************************************************************************************
******************************************************************UKB LIFECOURSE BODY SIZE - BRAIN PHENOTYPES MENDELIAN RANDOMIZATION****************************************************
************************************************************************************SCOTT CHIESA 27/07/23********************************************************************************
*****************************************************************************************************************************************************************************************

*****************************************************************************************************************************************************************************************
**************************************************************************Importing and Preparing all Summary Data for Analysis**********************************************************
*****************************************************************************************************************************************************************************************


//Export all outcome phenotype .txt files from Myriad into local folder//

//First delete any pre-existing files in subfolders from previous analysis runs as otherwise the later parts won't work//

cd "Y:\UKB\MR\Data Files\Child Outcomes"
local myfilelist : dir . files "*.dta"
foreach file of local myfilelist{
	erase "`file'"
}

cd "Y:\UKB\MR\Data Files\Adult Outcomes"
local myfilelist : dir . files "*.dta"
foreach file of local myfilelist{
	erase "`file'"
}

//Now convert all summary data .txt files from Myriad into .dta files//	

cd "Y:\UKB\MR\Individual Outcome Files from Myriad"	
local myfilelist : dir . files "*.txt"
foreach file of local myfilelist {
	di "`file'"		//display file name//
	import delimited "`file'", clear delimiter(tab) varnames(1)
	local outfile = subinstr("`file'", ".txt", "", .)
	save "`outfile'", replace
}

//Move all newly created child and adult .dta files to their respective folders for later//

mvfiles , infolder(".") outfolder("Y:\UKB\MR\Data Files\Child Outcomes") match("child_*.dta") erase
mvfiles , infolder(".") outfolder("Y:\UKB\MR\Data Files\Adult Outcomes") match("adult_*.dta") erase


*****************************Preparing Exposures for Merging******************************************

//Add prefix to all child exposure variables//

use "Y:\UKB\MR\Data Files\Child Exposures\BMI_child.dta", clear			
foreach var of varlist * {
	if "`var'" != "id" {
		rename `var' BMI_`var'
	}
	save "Y:\UKB\MR\Data Files\Child Exposures\child_exp.dta", replace
}

//Add prefix to all adult exposure variables//

use "Y:\UKB\MR\Data Files\Adult Exposures\BMI_adult.dta", clear
foreach var of varlist * {
	if "`var'" != "id" {
		rename `var' BMI_`var'
	}
	save "Y:\UKB\MR\Data Files\Adult Exposures\adult_exp.dta", replace
}


*****************************Preparing Outcomes for Merging******************************************

//Add prefix to all child outcome variables//

cd "Y:\UKB\MR\Data Files\Child Outcomes"
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

//Add prefix to all adult outcome variables//

cd "Y:\UKB\MR\Data Files\Adult Outcomes"
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


************Merging Exposures and Outcomes for Both Child and Adult Effects******************

//Merge child files//

cd "Y:\UKB\MR\Data Files\Child Outcomes"
local myfilelist : dir . files "*.dta"
use "Y:\UKB\MR\Data Files\Child Exposures\child_exp.dta"
foreach file of local myfilelist {
	merge 1:1 id using `file'
	drop _merge
	save "Y:\UKB\MR\Data Files\Child Outcomes\Child Effects File.dta", replace
}

//Merge adult files//

cd "Y:\UKB\MR\Data Files\Adult Outcomes"
local myfilelist : dir . files "*.dta"
use "Y:\UKB\MR\Data Files\Adult Exposures\adult_exp.dta"
foreach file of local myfilelist {
	merge 1:1 id using `file'
	drop _merge
	save "Y:\UKB\MR\Data Files\Adult Outcomes\Adult Effects File.dta", replace
}

//Drop single missing row in childhood data which is messing up graphs//

use "Y:\UKB\MR\Data Files\Child Outcomes\Child Effects File.dta", clear
drop in 152
save "Y:\UKB\MR\Data Files\Child Outcomes\Child Effects File.dta", replace

//Drop two missing rows in adulthood data which is messing up graphs//

use "Y:\UKB\MR\Data Files\Adult Outcomes\Adult Effects File.dta", clear
drop in 171 
drop in 289
save "Y:\UKB\MR\Data Files\Adult Outcomes\Adult Effects File.dta", replace


********************************************************************************************************************************************************************************************
***********************************************************************************PERFORMING MENDELIAN RANDOMIZATIONS**********************************************************************
**************************************************************************UNIVARIATE ONLY - MVMR RUN IN SEPARATE R SCRIPT*******************************************************************
********************************************************************************************************************************************************************************************

*****************************************MOST OF THIS CODE IS JUST REPEATING THE SAME MEASURES BUT FOR DIFFERENT TRAITS INDEXED IN DIFFERENT WAYS*******************************************
*************************************ORDER OF ANALYSIS IS 1) BRAIN TRAITS INDEXED TO ICV, 2) NON-INDEXED BRAIN TRAITS, 3) BRAIN TRAITS INDEXED TO TBV***************************************
***************************************************** EACH STARTS WITH CHILDHOOD SNPS VS OUTCOMES THEN REPEATS WITH ADULT SNPS VS OUTCOMES**************************************************


****************************************************************************************CHILDHOOD EFFECTS***********************************************************************************

log using "Y:\UKB\Updated_MR", replace

//TOTAL BRAIN VOLUME NORMALISED FOR INTRACRANIAL VOLUME//

use "Y:\UKB\MR\Data Files\Child Outcomes\Child Effects File.dta", clear

* IVW method (multiplicative random effects)
tab child_tvn_beta
local snp_num = `r(N)'
mregger child_tvn_beta BMI_beta [aw=1/(child_tvn_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger child_tvn_beta BMI_beta [aw=1/(child_tvn_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian child_tvn_beta child_tvn_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal child_tvn_beta child_tvn_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\tvn_results", modify
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
mregger child_tvn_beta BMI_beta [aw=1/(child_tvn_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger child_tvn_beta BMI_beta [aw=1/(child_tvn_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter child_tvn_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\tvn_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest child_tvn_beta child_tvn_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\tvn_forest.png", replace

* Asymmetry
mrfunnel child_tvn_beta child_tvn_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\tvn_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = child_tvn_beta/BMI_beta
gen wald_ratio_se = child_tvn_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\tvn_loo.png", replace


//GREY MATTER VOLUME NORMALISED FOR INTRACRANIAL VOLUME//

use "Y:\UKB\MR\Data Files\Child Outcomes\Child Effects File.dta", clear

* IVW method (multiplicative random effects)
tab child_gmn_beta
local snp_num = `r(N)'
mregger child_gmn_beta BMI_beta [aw=1/(child_gmn_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger child_gmn_beta BMI_beta [aw=1/(child_gmn_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian child_gmn_beta child_gmn_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal child_gmn_beta child_gmn_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\gmn_results", modify
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
mregger child_gmn_beta BMI_beta [aw=1/(child_gmn_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger child_gmn_beta BMI_beta [aw=1/(child_gmn_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter child_gmn_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\gmn_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest child_gmn_beta child_gmn_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\gmn_forest.png", replace

* Asymmetry
mrfunnel child_gmn_beta child_gmn_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\gmn_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = child_gmn_beta/BMI_beta
gen wald_ratio_se = child_gmn_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\gmn_loo.png", replace


//CORTICAL SURFACE AREA NORMALISED FOR INTRACRANIAL VOLUME//

use "Y:\UKB\MR\Data Files\Child Outcomes\Child Effects File.dta", clear

* IVW method (multiplicative random effects)
tab child_arean_beta
local snp_num = `r(N)'
mregger child_arean_beta BMI_beta [aw=1/(child_arean_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger child_arean_beta BMI_beta [aw=1/(child_arean_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian child_arean_beta child_arean_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal child_arean_beta child_arean_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\arean_results", modify
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
mregger child_arean_beta BMI_beta [aw=1/(child_arean_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger child_arean_beta BMI_beta [aw=1/(child_arean_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter child_arean_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\arean_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest child_arean_beta child_arean_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\arean_forest.png", replace

* Asymmetry
mrfunnel child_arean_beta child_arean_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\arean_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = child_arean_beta/BMI_beta
gen wald_ratio_se = child_arean_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\arean_loo.png", replace


//CORTICAL THICKNESS NORMALISED FOR INTRACRANIAL VOLUME//

*****NOTE THAT THIS ISN'T USED IN ANALYSIS AS NORMALISATION NOT APPROPRIATE HERE*****

use "Y:\UKB\MR\Data Files\Child Outcomes\Child Effects File.dta", clear

* IVW method (multiplicative random effects)
tab child_thickn_beta
local snp_num = `r(N)'
mregger child_thickn_beta BMI_beta [aw=1/(child_thickn_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger child_thickn_beta BMI_beta [aw=1/(child_thickn_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian child_thickn_beta child_thickn_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal child_thickn_beta child_thickn_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thickn_results", modify
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
mregger child_thickn_beta BMI_beta [aw=1/(child_thickn_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger child_thickn_beta BMI_beta [aw=1/(child_thickn_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter child_thickn_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thickn_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest child_thickn_beta child_thickn_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thickn_forest.png", replace

* Asymmetry
mrfunnel child_thickn_beta child_thickn_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thickn_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = child_thickn_beta/BMI_beta
gen wald_ratio_se = child_thickn_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thickn_loo.png", replace

//WHITE MATTER VOLUME NORMALISED FOR INTRACRANIAL VOLUME//

use "Y:\UKB\MR\Data Files\Child Outcomes\Child Effects File.dta", clear

* IVW method (multiplicative random effects)
tab child_wmn_beta
local snp_num = `r(N)'
mregger child_wmn_beta BMI_beta [aw=1/(child_wmn_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger child_wmn_beta BMI_beta [aw=1/(child_wmn_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian child_wmn_beta child_wmn_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal child_wmn_beta child_wmn_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\wmn_results", modify
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
mregger child_wmn_beta BMI_beta [aw=1/(child_wmn_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger child_wmn_beta BMI_beta [aw=1/(child_wmn_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter child_wmn_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\wmn_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest child_wmn_beta child_wmn_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\wmn_forest.png", replace

* Asymmetry
mrfunnel child_wmn_beta child_wmn_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\wmn_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = child_wmn_beta/BMI_beta
gen wald_ratio_se = child_wmn_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\wmn_loo.png", replace


//LOG-TRANSFORMED WHITE MATTER HYPERINTENSITIIES NORMALISED FOR INTRACRANIAL VOLUME//

*****NOTE THAT THIS ISN'T USED IN ANALYSIS AS NORMALISATION NOT APPROPRIATE HERE*****

use "Y:\UKB\MR\Data Files\Child Outcomes\Child Effects File.dta", clear

* IVW method (multiplicative random effects)
tab child_logwmhn_beta
local snp_num = `r(N)'
mregger child_logwmhn_beta BMI_beta [aw=1/(child_logwmhn_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger child_logwmhn_beta BMI_beta [aw=1/(child_logwmhn_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian child_logwmhn_beta child_logwmhn_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal child_logwmhn_beta child_logwmhn_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\logwmhn_results", modify
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
mregger child_logwmhn_beta BMI_beta [aw=1/(child_logwmhn_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger child_logwmhn_beta BMI_beta [aw=1/(child_logwmhn_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter child_logwmhn_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\logwmhn_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest child_logwmhn_beta child_logwmhn_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\logwmhn_forest.png", replace

* Asymmetry
mrfunnel child_logwmhn_beta child_logwmhn_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\logwmhn_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = child_logwmhn_beta/BMI_beta
gen wald_ratio_se = child_logwmhn_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\logwmhn_loo.png", replace


//THALAMUS NORMALISED FOR INTRACRANIAL VOLUME//

use "Y:\UKB\MR\Data Files\Child Outcomes\Child Effects File.dta", clear

* IVW method (multiplicative random effects)
tab child_thaln_beta
local snp_num = `r(N)'
mregger child_thaln_beta BMI_beta [aw=1/(child_thaln_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger child_thaln_beta BMI_beta [aw=1/(child_thaln_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian child_thaln_beta child_thaln_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal child_thaln_beta child_thaln_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thaln_results", modify
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
mregger child_thaln_beta BMI_beta [aw=1/(child_thaln_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger child_thaln_beta BMI_beta [aw=1/(child_thaln_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter child_thaln_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thaln_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest child_thaln_beta child_thaln_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thaln_forest.png", replace

* Asymmetry
mrfunnel child_thaln_beta child_thaln_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thaln_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = child_thaln_beta/BMI_beta
gen wald_ratio_se = child_thaln_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thaln_loo.png", replace

//HIPPOCAMPUS NORMALISED FOR INTRACRANIAL VOLUME//

use "Y:\UKB\MR\Data Files\Child Outcomes\Child Effects File.dta", clear

* IVW method (multiplicative random effects)
tab child_hippon_beta
local snp_num = `r(N)'
mregger child_hippon_beta BMI_beta [aw=1/(child_hippon_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger child_hippon_beta BMI_beta [aw=1/(child_hippon_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian child_hippon_beta child_hippon_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal child_hippon_beta child_hippon_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\hippon_results", modify
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
mregger child_hippon_beta BMI_beta [aw=1/(child_hippon_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger child_hippon_beta BMI_beta [aw=1/(child_hippon_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter child_hippon_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\hippon_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest child_hippon_beta child_hippon_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\hippon_forest.png", replace

* Asymmetry
mrfunnel child_hippon_beta child_hippon_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\hippon_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = child_hippon_beta/BMI_beta
gen wald_ratio_se = child_hippon_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\hippon_loo.png", replace

//AMYGDALA NORMALISED FOR INTRACRANIAL VOLUME//

use "Y:\UKB\MR\Data Files\Child Outcomes\Child Effects File.dta", clear

* IVW method (multiplicative random effects)
tab child_amygn_beta
local snp_num = `r(N)'
mregger child_amygn_beta BMI_beta [aw=1/(child_amygn_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger child_amygn_beta BMI_beta [aw=1/(child_amygn_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian child_amygn_beta child_amygn_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal child_amygn_beta child_amygn_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\amygn_results", modify
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
mregger child_amygn_beta BMI_beta [aw=1/(child_amygn_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger child_amygn_beta BMI_beta [aw=1/(child_amygn_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter child_amygn_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\amygn_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest child_amygn_beta child_amygn_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\amygn_forest.png", replace

* Asymmetry
mrfunnel child_amygn_beta child_amygn_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\amygn_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = child_amygn_beta/BMI_beta
gen wald_ratio_se = child_amygn_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\amygn_loo.png", replace

//VENTRICULAR CSF NORMALISED FOR INTRACRANIAL VOLUME//

use "Y:\UKB\MR\Data Files\Child Outcomes\Child Effects File.dta", clear

* IVW method (multiplicative random effects)
tab child_csfn_beta
local snp_num = `r(N)'
mregger child_csfn_beta BMI_beta [aw=1/(child_csfn_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger child_csfn_beta BMI_beta [aw=1/(child_csfn_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian child_csfn_beta child_csfn_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal child_csfn_beta child_csfn_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\csfn_results", modify
	putexcel A1="Outcome" B1="No. SNPs" C1="Method" D1="Beta" E1="SE" F1="p"
	local x = `x'+1
	putexcel A`x' = "Ventricular CSF" 
	putexcel B`x' = `snp_num'
	putexcel C`x' = "`method_`estimator''"
	putexcel D`x' = "`beta_`estimator''"
	putexcel E`x' = "`se_`estimator''"
	putexcel F`x' = "`p_`estimator''"
}

* Is there evidence of heterogeneity in the genetic effects?
mregger child_csfn_beta BMI_beta [aw=1/(child_csfn_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger child_csfn_beta BMI_beta [aw=1/(child_csfn_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter child_csfn_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\csfn_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest child_csfn_beta child_csfn_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\csfn_forest.png", replace

* Asymmetry
mrfunnel child_csfn_beta child_csfn_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\csfn_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = child_csfn_beta/BMI_beta
gen wald_ratio_se = child_csfn_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\csfn_loo.png", replace




******************************************************************************ADULTHOOD EFFECTS***********************************************************************************



//TOTAL BRAIN VOLUME NORMALISED FOR INTRACRANIAL VOLUME//

use "Y:\UKB\MR\Data Files\Adult Outcomes\Adult Effects File.dta", clear

* IVW method (multiplicative random effects)
tab adult_tvn_beta
local snp_num = `r(N)'
mregger adult_tvn_beta BMI_beta [aw=1/(adult_tvn_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger adult_tvn_beta BMI_beta [aw=1/(adult_tvn_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian adult_tvn_beta adult_tvn_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal adult_tvn_beta adult_tvn_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\tvn_results", modify
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
mregger adult_tvn_beta BMI_beta [aw=1/(adult_tvn_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger adult_tvn_beta BMI_beta [aw=1/(adult_tvn_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter adult_tvn_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\tvn_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest adult_tvn_beta adult_tvn_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\tvn_forest.png", replace

* Asymmetry
mrfunnel adult_tvn_beta adult_tvn_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\tvn_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = adult_tvn_beta/BMI_beta
gen wald_ratio_se = adult_tvn_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\tvn_loo.png", replace


//GREY MATTER VOLUME NORMALISED FOR INTRACRAIAL VOLUME//

use "Y:\UKB\MR\Data Files\Adult Outcomes\Adult Effects File.dta", clear

* IVW method (multiplicative random effects)
tab adult_gmn_beta
local snp_num = `r(N)'
mregger adult_gmn_beta BMI_beta [aw=1/(adult_gmn_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger adult_gmn_beta BMI_beta [aw=1/(adult_gmn_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian adult_gmn_beta adult_gmn_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal adult_gmn_beta adult_gmn_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\gmn_results", modify
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
mregger adult_gmn_beta BMI_beta [aw=1/(adult_gmn_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger adult_gmn_beta BMI_beta [aw=1/(adult_gmn_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter adult_gmn_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\gmn_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest adult_gmn_beta adult_gmn_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\gmn_forest.png", replace

* Asymmetry
mrfunnel adult_gmn_beta adult_gmn_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\gmn_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = adult_gmn_beta/BMI_beta
gen wald_ratio_se = adult_gmn_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\gmn_loo.png", replace


//CORTICAL SURFACE AREA NORMALISED FOR INTRACRANIAL VOLUME//

use "Y:\UKB\MR\Data Files\Adult Outcomes\Adult Effects File.dta", clear

* IVW method (multiplicative random effects)
tab adult_arean_beta
local snp_num = `r(N)'
mregger adult_arean_beta BMI_beta [aw=1/(adult_arean_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger adult_arean_beta BMI_beta [aw=1/(adult_arean_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian adult_arean_beta adult_arean_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal adult_arean_beta adult_arean_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\arean_results", modify
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
mregger adult_arean_beta BMI_beta [aw=1/(adult_arean_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger adult_arean_beta BMI_beta [aw=1/(adult_arean_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter adult_arean_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\arean_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest adult_arean_beta adult_arean_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\arean_forest.png", replace

* Asymmetry
mrfunnel adult_arean_beta adult_arean_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\arean_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = adult_arean_beta/BMI_beta
gen wald_ratio_se = adult_arean_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\arean_loo.png", replace


//CORTICAL THICKNESS NORMALISED FOR INTRACRANIAL VOLUME//

*****NOTE THAT THIS ISN'T USED IN ANALYSIS AS NORMALISATION NOT APPROPRIATE HERE*****

use "Y:\UKB\MR\Data Files\Adult Outcomes\Adult Effects File.dta", clear

* IVW method (multiplicative random effects)
tab adult_thickn_beta
local snp_num = `r(N)'
mregger adult_thickn_beta BMI_beta [aw=1/(adult_thickn_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger adult_thickn_beta BMI_beta [aw=1/(adult_thickn_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian adult_thickn_beta adult_thickn_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal adult_thickn_beta adult_thickn_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thickn_results", modify
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
mregger adult_thickn_beta BMI_beta [aw=1/(adult_thickn_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger adult_thickn_beta BMI_beta [aw=1/(adult_thickn_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter adult_thickn_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thickn_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest adult_thickn_beta adult_thickn_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thickn_forest.png", replace

* Asymmetry
mrfunnel adult_thickn_beta adult_thickn_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thickn_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = adult_thickn_beta/BMI_beta
gen wald_ratio_se = adult_thickn_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thickn_loo.png", replace


//WHITE MATTER VOLUME NORMALISED FOR INTRACRANIAL VOLUME/

use "Y:\UKB\MR\Data Files\Adult Outcomes\Adult Effects File.dta", clear

* IVW method (multiplicative random effects)
tab adult_wmn_beta
local snp_num = `r(N)'
mregger adult_wmn_beta BMI_beta [aw=1/(adult_wmn_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger adult_wmn_beta BMI_beta [aw=1/(adult_wmn_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian adult_wmn_beta adult_wmn_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal adult_wmn_beta adult_wmn_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\wmn_results", modify
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
mregger adult_wmn_beta BMI_beta [aw=1/(adult_wmn_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger adult_wmn_beta BMI_beta [aw=1/(adult_wmn_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter adult_wmn_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\wmn_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest adult_wmn_beta adult_wmn_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\wmn_forest.png", replace

* Asymmetry
mrfunnel adult_wmn_beta adult_wmn_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\wmn_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = adult_wmn_beta/BMI_beta
gen wald_ratio_se = adult_wmn_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\wmn_loo.png", replace

//LOG-TRANSFORMED WHITE MATTER HYPERINTENSITIES NORMALISED FOR INTRACRANIAL VOLUME//

*****NOTE THAT THIS ISN'T USED IN ANALYSIS AS NORMALISATION NOT APPROPRIATE HERE*****

use "Y:\UKB\MR\Data Files\Adult Outcomes\Adult Effects File.dta", clear

* IVW method (multiplicative random effects)
tab adult_logwmhn_beta
local snp_num = `r(N)'
mregger adult_logwmhn_beta BMI_beta [aw=1/(adult_logwmhn_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger adult_logwmhn_beta BMI_beta [aw=1/(adult_logwmhn_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian adult_logwmhn_beta adult_logwmhn_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal adult_logwmhn_beta adult_logwmhn_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\logwmhn_results", modify
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
mregger adult_logwmhn_beta BMI_beta [aw=1/(adult_logwmhn_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger adult_logwmhn_beta BMI_beta [aw=1/(adult_logwmhn_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter adult_logwmhn_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\logwmhn_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest adult_logwmhn_beta adult_logwmhn_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\logwmhn_forest.png", replace

* Asymmetry
mrfunnel adult_logwmhn_beta adult_logwmhn_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\logwmhn_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = adult_logwmhn_beta/BMI_beta
gen wald_ratio_se = adult_logwmhn_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\logwmhn_loo.png", replace

//THALAMUS NORMALISED FOR INTRACRANIAL VOLUME//

use "Y:\UKB\MR\Data Files\Adult Outcomes\Adult Effects File.dta", clear

* IVW method (multiplicative random effects)
tab adult_thaln_beta
local snp_num = `r(N)'
mregger adult_thaln_beta BMI_beta [aw=1/(adult_thaln_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger adult_thaln_beta BMI_beta [aw=1/(adult_thaln_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian adult_thaln_beta adult_thaln_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal adult_thaln_beta adult_thaln_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thaln_results", modify
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
mregger adult_thaln_beta BMI_beta [aw=1/(adult_thaln_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger adult_thaln_beta BMI_beta [aw=1/(adult_thaln_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter adult_thaln_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thaln_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest adult_thaln_beta adult_thaln_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thaln_forest.png", replace

* Asymmetry
mrfunnel adult_thaln_beta adult_thaln_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thaln_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = adult_thaln_beta/BMI_beta
gen wald_ratio_se = adult_thaln_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thaln_loo.png", replace

//HIPPOCAMPUS NORMALISED FOR INTRACRANIAL VOLUME//

use "Y:\UKB\MR\Data Files\Adult Outcomes\Adult Effects File.dta", clear

* IVW method (multiplicative random effects)
tab adult_hippon_beta
local snp_num = `r(N)'
mregger adult_hippon_beta BMI_beta [aw=1/(adult_hippon_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger adult_hippon_beta BMI_beta [aw=1/(adult_hippon_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian adult_hippon_beta adult_hippon_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal adult_hippon_beta adult_hippon_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\hippon_results", modify
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
mregger adult_hippon_beta BMI_beta [aw=1/(adult_hippon_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger adult_hippon_beta BMI_beta [aw=1/(adult_hippon_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter adult_hippon_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\hippon_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest adult_hippon_beta adult_hippon_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\hippon_forest.png", replace

* Asymmetry
mrfunnel adult_hippon_beta adult_hippon_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\hippon_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = adult_hippon_beta/BMI_beta
gen wald_ratio_se = adult_hippon_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\hippon_loo.png", replace

//AMYGDALA NORMALISED FOR INTRACRANIAL VOLUME//

use "Y:\UKB\MR\Data Files\Adult Outcomes\Adult Effects File.dta", clear

* IVW method (multiplicative random effects)
tab adult_amygn_beta
local snp_num = `r(N)'
mregger adult_amygn_beta BMI_beta [aw=1/(adult_amygn_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger adult_amygn_beta BMI_beta [aw=1/(adult_amygn_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian adult_amygn_beta adult_amygn_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal adult_amygn_beta adult_amygn_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\amygn_results", modify
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
mregger adult_amygn_beta BMI_beta [aw=1/(adult_amygn_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger adult_amygn_beta BMI_beta [aw=1/(adult_amygn_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter adult_amygn_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\amygn_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest adult_amygn_beta adult_amygn_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\amygn_forest.png", replace

* Asymmetry
mrfunnel adult_amygn_beta adult_amygn_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\amygn_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = adult_amygn_beta/BMI_beta
gen wald_ratio_se = adult_amygn_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\amygn_loo.png", replace


//VENTRICULAR CSF NORMALISED FOR INTRACRANIAL VOLUME//

use "Y:\UKB\MR\Data Files\Adult Outcomes\Adult Effects File.dta", clear

* IVW method (multiplicative random effects)
tab adult_csfn_beta
local snp_num = `r(N)'
mregger adult_csfn_beta BMI_beta [aw=1/(adult_csfn_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger adult_csfn_beta BMI_beta [aw=1/(adult_csfn_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian adult_csfn_beta adult_csfn_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal adult_csfn_beta adult_csfn_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\csfn_results", modify
	putexcel A1="Outcome" B1="No. SNPs" C1="Method" D1="Beta" E1="SE" F1="p"
	local x = `x'+1
	putexcel A`x' = "Ventricular CSF" 
	putexcel B`x' = `snp_num'
	putexcel C`x' = "`method_`estimator''"
	putexcel D`x' = "`beta_`estimator''"
	putexcel E`x' = "`se_`estimator''"
	putexcel F`x' = "`p_`estimator''"
}

* Is there evidence of heterogeneity in the genetic effects?
mregger adult_csfn_beta BMI_beta [aw=1/(adult_csfn_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger adult_csfn_beta BMI_beta [aw=1/(adult_csfn_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter adult_csfn_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\csfn_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest adult_csfn_beta adult_csfn_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\csfn_forest.png", replace

* Asymmetry
mrfunnel adult_csfn_beta adult_csfn_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\csfn_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = adult_csfn_beta/BMI_beta
gen wald_ratio_se = adult_csfn_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\csfn_loo.png", replace




******************************************************************************************************************************************************************************************************************************************************************************SENSITIVITY ANALYSES USING DIFFERENT ADJUSTMENTS******************************************************************
************************************************************************************************************************************************************************************************

********************************************************************************ABSOLUTE BRAIN VOLUMES WITHOUT NORMALISATION********************************************************************


******************************************************************************************CHILDHOOD EFFECTS*************************************************************************************

//INTRACRANIAL VOLUME//

use "Y:\UKB\MR\Data Files\Child Outcomes\Child Effects File.dta", clear

* IVW method (multiplicative random effects)
tab child_icv_beta
local snp_num = `r(N)'
mregger child_icv_beta BMI_beta [aw=1/(child_icv_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger child_icv_beta BMI_beta [aw=1/(child_icv_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian child_icv_beta child_icv_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal child_icv_beta child_icv_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\icv_results", modify
	putexcel A1="Outcome" B1="No. SNPs" C1="Method" D1="Beta" E1="SE" F1="p"
	local x = `x'+1
	putexcel A`x' = "Intracranial Volume" 
	putexcel B`x' = `snp_num'
	putexcel C`x' = "`method_`estimator''"
	putexcel D`x' = "`beta_`estimator''"
	putexcel E`x' = "`se_`estimator''"
	putexcel F`x' = "`p_`estimator''"
}

* Is there evidence of heterogeneity in the genetic effects?
mregger child_icv_beta BMI_beta [aw=1/(child_icv_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger child_icv_beta BMI_beta [aw=1/(child_icv_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter child_icv_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\icv_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest child_icv_beta child_icv_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\icv_forest.png", replace

* Asymmetry
mrfunnel child_icv_beta child_icv_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\icv_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = child_icv_beta/BMI_beta
gen wald_ratio_se = child_icv_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\icv_loo.png", replace

//TOTAL BRAIN VOLUME //

use "Y:\UKB\MR\Data Files\Child Outcomes\Child Effects File.dta", clear

* IVW method (multiplicative random effects)
tab child_tv_beta
local snp_num = `r(N)'
mregger child_tv_beta BMI_beta [aw=1/(child_tv_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger child_tv_beta BMI_beta [aw=1/(child_tv_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian child_tv_beta child_tv_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal child_tv_beta child_tv_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\tv_results", modify
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
mregger child_tv_beta BMI_beta [aw=1/(child_tv_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger child_tv_beta BMI_beta [aw=1/(child_tv_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter child_tv_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\tv_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest child_tv_beta child_tv_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\tv_forest.png", replace

* Asymmetry
mrfunnel child_tv_beta child_tv_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\tv_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = child_tv_beta/BMI_beta
gen wald_ratio_se = child_tv_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\tv_loo.png", replace

//GREY MATTER VOLUME//

use "Y:\UKB\MR\Data Files\Child Outcomes\Child Effects File.dta", clear

* IVW method (multiplicative random effects)
tab child_gm_beta
local snp_num = `r(N)'
mregger child_gm_beta BMI_beta [aw=1/(child_gm_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger child_gm_beta BMI_beta [aw=1/(child_gm_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian child_gm_beta child_gm_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal child_gm_beta child_gm_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\gm_results", modify
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
mregger child_gm_beta BMI_beta [aw=1/(child_gm_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger child_gm_beta BMI_beta [aw=1/(child_gm_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter child_gm_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\gm_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest child_gm_beta child_gm_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\gm_forest.png", replace

* Asymmetry
mrfunnel child_gm_beta child_gm_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\gm_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = child_gm_beta/BMI_beta
gen wald_ratio_se = child_gm_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\gm_loo.png", replace


//CORTICAL SURFACE AREA//

use "Y:\UKB\MR\Data Files\Child Outcomes\Child Effects File.dta", clear

* IVW method (multiplicative random effects)
tab child_area_beta
local snp_num = `r(N)'
mregger child_area_beta BMI_beta [aw=1/(child_area_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger child_area_beta BMI_beta [aw=1/(child_area_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian child_area_beta child_area_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal child_area_beta child_area_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\area_results", modify
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
mregger child_area_beta BMI_beta [aw=1/(child_area_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger child_area_beta BMI_beta [aw=1/(child_area_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter child_area_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\area_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest child_area_beta child_area_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\area_forest.png", replace

* Asymmetry
mrfunnel child_area_beta child_area_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\area_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = child_area_beta/BMI_beta
gen wald_ratio_se = child_area_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\area_loo.png", replace


//CORTICAL THICKNESS//

use "Y:\UKB\MR\Data Files\Child Outcomes\Child Effects File.dta", clear

* IVW method (multiplicative random effects)
tab child_thick_beta
local snp_num = `r(N)'
mregger child_thick_beta BMI_beta [aw=1/(child_thick_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger child_thick_beta BMI_beta [aw=1/(child_thick_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian child_thick_beta child_thick_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal child_thick_beta child_thick_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thick_results", modify
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
mregger child_thick_beta BMI_beta [aw=1/(child_thick_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger child_thick_beta BMI_beta [aw=1/(child_thick_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter child_thick_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thick_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest child_thick_beta child_thickb_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thick_forest.png", replace

* Asymmetry
mrfunnel child_thick_beta child_thickb_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thick_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = child_thick_beta/BMI_beta
gen wald_ratio_se = child_thick_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thick_loo.png", replace


//WHITE MATTER VOLUME//

use "Y:\UKB\MR\Data Files\Child Outcomes\Child Effects File.dta", clear

* IVW method (multiplicative random effects)
tab child_wm_beta
local snp_num = `r(N)'
mregger child_wm_beta BMI_beta [aw=1/(child_wm_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger child_wm_beta BMI_beta [aw=1/(child_wm_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian child_wm_beta child_wm_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal child_wm_beta child_wm_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\wm_results", modify
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
mregger child_wm_beta BMI_beta [aw=1/(child_wm_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger child_wm_beta BMI_beta [aw=1/(child_wm_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter child_wm_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\wm_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest child_wm_beta child_wm_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\wm_forest.png", replace

* Asymmetry
mrfunnel child_wm_beta child_wm_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\wm_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = child_wm_beta/BMI_beta
gen wald_ratio_se = child_wm_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\wm_loo.png", replace


//LOG-TRANSFORMED WHITE MATTER HYPERINTENSITIES//

use "Y:\UKB\MR\Data Files\Child Outcomes\Child Effects File.dta", clear

* IVW method (multiplicative random effects)
tab child_logwmh_beta
local snp_num = `r(N)'
mregger child_logwmh_beta BMI_beta [aw=1/(child_logwmh_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger child_logwmh_beta BMI_beta [aw=1/(child_logwmh_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian child_logwmh_beta child_logwmh_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal child_logwmh_beta child_logwmh_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\logwmh_results", modify
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
mregger child_logwmh_beta BMI_beta [aw=1/(child_logwmh_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger child_logwmh_beta BMI_beta [aw=1/(child_logwmh_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter child_logwmh_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\logwmh_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest child_logwmh_beta child_logwmh_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\logwmh_forest.png", replace

* Asymmetry
mrfunnel child_logwmh_beta child_logwmh_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\logwmh_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = child_logwmh_beta/BMI_beta
gen wald_ratio_se = child_logwmh_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\logwmh_loo.png", replace


//THALAMUS//

use "Y:\UKB\MR\Data Files\Child Outcomes\Child Effects File.dta", clear

* IVW method (multiplicative random effects)
tab child_thalamus_beta
local snp_num = `r(N)'
mregger child_thalamus_beta BMI_beta [aw=1/(child_thalamus_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger child_thalamus_beta BMI_beta [aw=1/(child_thalamus_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian child_thalamus_beta child_thalamus_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal child_thalamus_beta child_thalamus_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thalamus_results", modify
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
mregger child_thalamus_beta BMI_beta [aw=1/(child_thalamus_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger child_thalamus_beta BMI_beta [aw=1/(child_thalamus_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter child_thalamus_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thalamus_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest child_thalamus_beta child_thalamus_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thalamus_forest.png", replace

* Asymmetry
mrfunnel child_thalamus_beta child_thalamus_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thalamus_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = child_thalamus_beta/BMI_beta
gen wald_ratio_se = child_thalamus_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thalamus_loo.png", replace


//HIPPOCAMPUS//

use "Y:\UKB\MR\Data Files\Child Outcomes\Child Effects File.dta", clear

* IVW method (multiplicative random effects)
tab child_hippo_beta
local snp_num = `r(N)'
mregger child_hippo_beta BMI_beta [aw=1/(child_hippo_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger child_hippo_beta BMI_beta [aw=1/(child_hippo_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian child_hippo_beta child_hippo_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal child_hippo_beta child_hippo_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\hippo_results", modify
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
mregger child_hippo_beta BMI_beta [aw=1/(child_hippo_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger child_hippo_beta BMI_beta [aw=1/(child_hippo_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter child_hippo_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\hippo_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest child_hippo_beta child_hippob_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\hippo_forest.png", replace

* Asymmetry
mrfunnel child_hippo_beta child_hippo_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\hippo_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = child_hippo_beta/BMI_beta
gen wald_ratio_se = child_hippo_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\hippo_loo.png", replace


//AMYGDALA//

use "Y:\UKB\MR\Data Files\Child Outcomes\Child Effects File.dta", clear

* IVW method (multiplicative random effects)
tab child_amygdala_beta
local snp_num = `r(N)'
mregger child_amygdala_beta BMI_beta [aw=1/(child_amygdala_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger child_amygdala_beta BMI_beta [aw=1/(child_amygdala_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian child_amygdala_beta child_amygdala_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal child_amygdala_beta child_amygdala_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\amygdala_results", modify
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
mregger child_amygdala_beta BMI_beta [aw=1/(child_amygdala_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger child_amygdala_beta BMI_beta [aw=1/(child_amygdala_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter child_amygdala_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\amygdala_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest child_amygdala_beta child_amygdala_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\amygdala_forest.png", replace

* Asymmetry
mrfunnel child_amygdala_beta child_amygdala_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\amygdala_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = child_amygdala_beta/BMI_beta
gen wald_ratio_se = child_amygdala_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\amygdala_loo.png", replace


//VENTRICULAR CSF//

use "Y:\UKB\MR\Data Files\Child Outcomes\Child Effects File.dta", clear

* IVW method (multiplicative random effects)
tab child_csf_beta
local snp_num = `r(N)'
mregger child_csf_beta BMI_beta [aw=1/(child_csf_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger child_csf_beta BMI_beta [aw=1/(child_csf_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian child_csf_beta child_csf_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal child_csf_beta child_csf_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\csf_results", modify
	putexcel A1="Outcome" B1="No. SNPs" C1="Method" D1="Beta" E1="SE" F1="p"
	local x = `x'+1
	putexcel A`x' = "Ventricular CSF" 
	putexcel B`x' = `snp_num'
	putexcel C`x' = "`method_`estimator''"
	putexcel D`x' = "`beta_`estimator''"
	putexcel E`x' = "`se_`estimator''"
	putexcel F`x' = "`p_`estimator''"
}

* Is there evidence of heterogeneity in the genetic effects?
mregger child_csf_beta BMI_beta [aw=1/(child_csf_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger child_csf_beta BMI_beta [aw=1/(child_csf_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter child_csf_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\csf_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest child_csf_beta child_csf_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\csf_forest.png", replace

* Asymmetry
mrfunnel child_csf_beta child_csf_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\csf_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = child_csf_beta/BMI_beta
gen wald_ratio_se = child_csf_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\csf_loo.png", replace


**********************************************************************************ADULTHOOD EFFECTS****************************************************************************************

//INTRACRANIAL VOLUME//

use "Y:\UKB\MR\Data Files\Adult Outcomes\Adult Effects File.dta", clear

* IVW method (multiplicative random effects)
tab adult_icv_beta
local snp_num = `r(N)'
mregger adult_icv_beta BMI_beta [aw=1/(adult_icv_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger adult_icv_beta BMI_beta [aw=1/(adult_icv_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian adult_icv_beta adult_icv_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal adult_icv_beta adult_icv_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\icv_results", modify
	putexcel A1="Outcome" B1="No. SNPs" C1="Method" D1="Beta" E1="SE" F1="p"
	local x = `x'+1
	putexcel A`x' = "Intracranial Volume" 
	putexcel B`x' = `snp_num'
	putexcel C`x' = "`method_`estimator''"
	putexcel D`x' = "`beta_`estimator''"
	putexcel E`x' = "`se_`estimator''"
	putexcel F`x' = "`p_`estimator''"
}

* Is there evidence of heterogeneity in the genetic effects?
mregger adult_icv_beta BMI_beta [aw=1/(adult_icv_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger adult_icv_beta BMI_beta [aw=1/(adult_icv_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter adult_icv_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\icv_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest adult_icv_beta adult_icv_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\icv_forest.png", replace

* Asymmetry
mrfunnel adult_icv_beta adult_icv_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\icv_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = adult_icv_beta/BMI_beta
gen wald_ratio_se = adult_icv_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\icv_loo.png", replace


//TOTAL VOLUME//

use "Y:\UKB\MR\Data Files\Adult Outcomes\Adult Effects File.dta", clear

* IVW method (multiplicative random effects)
tab adult_tv_beta
local snp_num = `r(N)'
mregger adult_tv_beta BMI_beta [aw=1/(adult_tv_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger adult_tv_beta BMI_beta [aw=1/(adult_tv_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian adult_tv_beta adult_tv_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal adult_tv_beta adult_tv_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\tv_results", modify
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
mregger adult_tv_beta BMI_beta [aw=1/(adult_tv_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger adult_tv_beta BMI_beta [aw=1/(adult_tv_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter adult_tv_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\tv_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest adult_tv_beta adult_tv_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\tv_forest.png", replace

* Asymmetry
mrfunnel adult_tv_beta adult_tv_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\tv_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = adult_tv_beta/BMI_beta
gen wald_ratio_se = adult_tv_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\tv_loo.png", replace


//GREY MATTER VOLUME//

use "Y:\UKB\MR\Data Files\Adult Outcomes\Adult Effects File.dta", clear

* IVW method (multiplicative random effects)
tab adult_gm_beta
local snp_num = `r(N)'
mregger adult_gm_beta BMI_beta [aw=1/(adult_gm_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger adult_gm_beta BMI_beta [aw=1/(adult_gm_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian adult_gm_beta adult_gm_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal adult_gm_beta adult_gm_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\gm_results", modify
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
mregger adult_gm_beta BMI_beta [aw=1/(adult_gm_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger adult_gm_beta BMI_beta [aw=1/(adult_gm_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter adult_gm_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\gm_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest adult_gm_beta adult_gm_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\gm_forest.png", replace

* Asymmetry
mrfunnel adult_gm_beta adult_gm_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\gm_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = adult_gm_beta/BMI_beta
gen wald_ratio_se = adult_gm_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\gm_loo.png", replace


//CORTICAL SURFACE AREA NORMALISED FOR TBV//

use "Y:\UKB\MR\Data Files\Adult Outcomes\Adult Effects File.dta", clear

* IVW method (multiplicative random effects)
tab adult_area_beta
local snp_num = `r(N)'
mregger adult_area_beta BMI_beta [aw=1/(adult_area_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger adult_area_beta BMI_beta [aw=1/(adult_area_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian adult_area_beta adult_area_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal adult_area_beta adult_area_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\area_results", modify
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
mregger adult_area_beta BMI_beta [aw=1/(adult_area_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger adult_area_beta BMI_beta [aw=1/(adult_area_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter adult_area_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\area_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest adult_area_beta adult_area_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\area_forest.png", replace

* Asymmetry
mrfunnel adult_area_beta adult_area_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\area_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = adult_area_beta/BMI_beta
gen wald_ratio_se = adult_area_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\area_loo.png", replace


//CORTICAL THICKNESS//

use "Y:\UKB\MR\Data Files\Adult Outcomes\Adult Effects File.dta", clear

* IVW method (multiplicative random effects)
tab adult_thick_beta
local snp_num = `r(N)'
mregger adult_thick_beta BMI_beta [aw=1/(adult_thick_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger adult_thick_beta BMI_beta [aw=1/(adult_thick_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian adult_thick_beta adult_thick_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal adult_thick_beta adult_thick_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thick_results", modify
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
mregger adult_thick_beta BMI_beta [aw=1/(adult_thick_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger adult_thick_beta BMI_beta [aw=1/(adult_thick_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter adult_thick_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thick_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest adult_thick_beta adult_thick_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thick_forest.png", replace

* Asymmetry
mrfunnel adult_thick_beta adult_thick_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thick_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = adult_thick_beta/BMI_beta
gen wald_ratio_se = adult_thick_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thick_loo.png", replace


//WHITE MATTER VOLUME//

use "Y:\UKB\MR\Data Files\Adult Outcomes\Adult Effects File.dta", clear

* IVW method (multiplicative random effects)
tab adult_wm_beta
local snp_num = `r(N)'
mregger adult_wm_beta BMI_beta [aw=1/(adult_wm_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger adult_wm_beta BMI_beta [aw=1/(adult_wm_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian adult_wm_beta adult_wm_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal adult_wm_beta adult_wm_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\wm_results", modify
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
mregger adult_wm_beta BMI_beta [aw=1/(adult_wm_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger adult_wm_beta BMI_beta [aw=1/(adult_wm_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter adult_wm_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\wm_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest adult_wm_beta adult_wm_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\wm_forest.png", replace

* Asymmetry
mrfunnel adult_wm_beta adult_wm_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\wm_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = adult_wm_beta/BMI_beta
gen wald_ratio_se = adult_wm_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\wm_loo.png", replace


//LOG-TRANSFORMED WHITE MATTER HYPERINTENSITIES//

use "Y:\UKB\MR\Data Files\Adult Outcomes\Adult Effects File.dta", clear

* IVW method (multiplicative random effects)
tab adult_logwmh_beta
local snp_num = `r(N)'
mregger adult_logwmh_beta BMI_beta [aw=1/(adult_logwmh_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger adult_logwmh_beta BMI_beta [aw=1/(adult_logwmh_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian adult_logwmh_beta adult_logwmh_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal adult_logwmh_beta adult_logwmh_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\logwmh_results", modify
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
mregger adult_logwmh_beta BMI_beta [aw=1/(adult_logwmh_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger adult_logwmh_beta BMI_beta [aw=1/(adult_logwmh_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter adult_logwmh_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\logwmh_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest adult_logwmh_beta adult_logwmh_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\logwmh_forest.png", replace

* Asymmetry
mrfunnel adult_logwmh_beta adult_logwmh_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\logwmh_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = adult_logwmh_beta/BMI_beta
gen wald_ratio_se = adult_logwmh_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\logwmh_loo.png", replace


//THALAMUS//

use "Y:\UKB\MR\Data Files\Adult Outcomes\Adult Effects File.dta", clear

* IVW method (multiplicative random effects)
tab adult_thalamus_beta
local snp_num = `r(N)'
mregger adult_thalamus_beta BMI_beta [aw=1/(adult_thalamus_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger adult_thalamus_beta BMI_beta [aw=1/(adult_thalamus_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian adult_thalamus_beta adult_thalamus_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal adult_thalamus_beta adult_thalamus_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thalamus_results", modify
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
mregger adult_thalamus_beta BMI_beta [aw=1/(adult_thalamus_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger adult_thalamus_beta BMI_beta [aw=1/(adult_thalamus_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter adult_thalamus_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thalamus_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest adult_thalamus_beta adult_thalamus_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thalamus_forest.png", replace

* Asymmetry
mrfunnel adult_thalamus_beta adult_thalamus_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thalamus_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = adult_thalamus_beta/BMI_beta
gen wald_ratio_se = adult_thalamus_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thalamus_loo.png", replace


//HIPPOCAMPUS//

use "Y:\UKB\MR\Data Files\Adult Outcomes\Adult Effects File.dta", clear

* IVW method (multiplicative random effects)
tab adult_hippo_beta
local snp_num = `r(N)'
mregger adult_hippo_beta BMI_beta [aw=1/(adult_hippo_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger adult_hippo_beta BMI_beta [aw=1/(adult_hippo_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian adult_hippo_beta adult_hippo_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal adult_hippo_beta adult_hippo_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\hippo_results", modify
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
mregger adult_hippo_beta BMI_beta [aw=1/(adult_hippo_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger adult_hippo_beta BMI_beta [aw=1/(adult_hippo_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter adult_hippo_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\hippo_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest adult_hippo_beta adult_hippo_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\hippo_forest.png", replace

* Asymmetry
mrfunnel adult_hippo_beta adult_hippo_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\hippo_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = adult_hippo_beta/BMI_beta
gen wald_ratio_se = adult_hippo_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\hippo_loo.png", replace


//AMYGDALA//

use "Y:\UKB\MR\Data Files\Adult Outcomes\Adult Effects File.dta", clear

* IVW method (multiplicative random effects)
tab adult_amygdala_beta
local snp_num = `r(N)'
mregger adult_amygdala_beta BMI_beta [aw=1/(adult_amygdala_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger adult_amygdala_beta BMI_beta [aw=1/(adult_amygdala_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian adult_amygdala_beta adult_amygdala_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal adult_amygdala_beta adult_amygdala_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\amygdala_results", modify
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
mregger adult_amygdala_beta BMI_beta [aw=1/(adult_amygdala_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger adult_amygdala_beta BMI_beta [aw=1/(adult_amygdala_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter adult_amygdala_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\amygdala_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest adult_amygdala_beta adult_amygdala_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\amygdala_forest.png", replace

* Asymmetry
mrfunnel adult_amygdala_beta adult_amygdala_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\amygdala_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = adult_amygdala_beta/BMI_beta
gen wald_ratio_se = adult_amygdala_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\amygdala_loo.png", replace


//VENTRICULAR CSF//

use "Y:\UKB\MR\Data Files\Adult Outcomes\Adult Effects File.dta", clear

* IVW method (multiplicative random effects)
tab adult_csf_beta
local snp_num = `r(N)'
mregger adult_csf_beta BMI_beta [aw=1/(adult_csf_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger adult_csf_beta BMI_beta [aw=1/(adult_csf_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian adult_csf_beta adult_csf_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal adult_csf_beta adult_csf_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\csfn_results", modify
	putexcel A1="Outcome" B1="No. SNPs" C1="Method" D1="Beta" E1="SE" F1="p"
	local x = `x'+1
	putexcel A`x' = "Ventricular CSF" 
	putexcel B`x' = `snp_num'
	putexcel C`x' = "`method_`estimator''"
	putexcel D`x' = "`beta_`estimator''"
	putexcel E`x' = "`se_`estimator''"
	putexcel F`x' = "`p_`estimator''"
}

* Is there evidence of heterogeneity in the genetic effects?
mregger adult_csf_beta BMI_beta [aw=1/(adult_csf_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger adult_csf_beta BMI_beta [aw=1/(adult_csf_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter adult_csf_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\csf_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest adult_csf_beta adult_csf_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\csf_forest.png", replace

* Asymmetry
mrfunnel adult_csf_beta adult_csf_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\csf_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = adult_csf_beta/BMI_beta
gen wald_ratio_se = adult_csf_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\csf_loo.png", replace



*************************************************************************VARIOUS BRAIN VOLUMES NORMALISED FOR TBV RATHER THAN ICV******************************************************



***************************************************************************************CHILDHOOD EFFECTS*******************************************************************************


//GREY MATTER VOLUME NORMALISED FOR TBV//

use "Y:\UKB\MR\Data Files\Child Outcomes\Child Effects File.dta", clear

* IVW method (multiplicative random effects)
tab child_gmb_beta
local snp_num = `r(N)'
mregger child_gmb_beta BMI_beta [aw=1/(child_gmb_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger child_gmb_beta BMI_beta [aw=1/(child_gmb_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian child_gmb_beta child_gmb_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal child_gmb_beta child_gmb_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\gmb_results", modify
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
mregger child_gmb_beta BMI_beta [aw=1/(child_gmb_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger child_gmb_beta BMI_beta [aw=1/(child_gmb_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter child_gmb_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\gmb_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest child_gmb_beta child_gmb_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\gmb_forest.png", replace

* Asymmetry
mrfunnel child_gmb_beta child_gmb_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\gmb_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = child_gmb_beta/BMI_beta
gen wald_ratio_se = child_gmb_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\gmb_loo.png", replace


//CORTICAL SURFACE AREA NORMALISED FOR TBV//

use "Y:\UKB\MR\Data Files\Child Outcomes\Child Effects File.dta", clear

* IVW method (multiplicative random effects)
tab child_areab_beta
local snp_num = `r(N)'
mregger child_areab_beta BMI_beta [aw=1/(child_areab_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger child_areab_beta BMI_beta [aw=1/(child_areab_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian child_areab_beta child_areab_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal child_areab_beta child_areab_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\areab_results", modify
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
mregger child_areab_beta BMI_beta [aw=1/(child_areab_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger child_areab_beta BMI_beta [aw=1/(child_areab_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter child_areab_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\areab_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest child_areab_beta child_areab_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\areab_forest.png", replace

* Asymmetry
mrfunnel child_areab_beta child_areab_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\areab_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = child_areab_beta/BMI_beta
gen wald_ratio_se = child_areab_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\areab_loo.png", replace


//CORTICAL THICKNESS NORMALISED FOR TBV//

use "Y:\UKB\MR\Data Files\Child Outcomes\Child Effects File.dta", clear

* IVW method (multiplicative random effects)
tab child_thickb_beta
local snp_num = `r(N)'
mregger child_thickb_beta BMI_beta [aw=1/(child_thickb_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger child_thickb_beta BMI_beta [aw=1/(child_thickb_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian child_thickb_beta child_thickb_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal child_thickb_beta child_thickb_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thickb_results", modify
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
mregger child_thickb_beta BMI_beta [aw=1/(child_thickb_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger child_thickb_beta BMI_beta [aw=1/(child_thickb_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter child_thickb_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thickb_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest child_thickb_beta child_thickb_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thickb_forest.png", replace

* Asymmetry
mrfunnel child_thickb_beta child_thickb_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thickb_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = child_thickb_beta/BMI_beta
gen wald_ratio_se = child_thickb_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thickb_loo.png", replace


//WHITE MATTER VOLUME NORMALISED FOR TBV//

use "Y:\UKB\MR\Data Files\Child Outcomes\Child Effects File.dta", clear

* IVW method (multiplicative random effects)
tab child_wmb_beta
local snp_num = `r(N)'
mregger child_wmb_beta BMI_beta [aw=1/(child_wmb_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger child_wmb_beta BMI_beta [aw=1/(child_wmb_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian child_wmb_beta child_wmb_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal child_wmb_beta child_wmb_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\wmb_results", modify
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
mregger child_wmb_beta BMI_beta [aw=1/(child_wmb_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger child_wmb_beta BMI_beta [aw=1/(child_wmb_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter child_wmb_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\wmb_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest child_wmb_beta child_wmb_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\wmb_forest.png", replace

* Asymmetry
mrfunnel child_wmb_beta child_wmb_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\wmb_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = child_wmb_beta/BMI_beta
gen wald_ratio_se = child_wmb_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\wmb_loo.png", replace


//LOG-TRANSFORMED WHITE MATTER HYPERINTENSITIES NORMALISED FOR TBV//

use "Y:\UKB\MR\Data Files\Child Outcomes\Child Effects File.dta", clear

* IVW method (multiplicative random effects)
tab child_logwmhb_beta
local snp_num = `r(N)'
mregger child_logwmhb_beta BMI_beta [aw=1/(child_logwmhb_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger child_logwmhb_beta BMI_beta [aw=1/(child_logwmhb_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian child_logwmhb_beta child_logwmhb_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal child_logwmhb_beta child_logwmhb_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\logwmhb_results", modify
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
mregger child_logwmhb_beta BMI_beta [aw=1/(child_logwmhb_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger child_logwmhb_beta BMI_beta [aw=1/(child_logwmhb_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter child_logwmhb_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\logwmhb_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest child_logwmhb_beta child_logwmhb_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\logwmhb_forest.png", replace

* Asymmetry
mrfunnel child_logwmhb_beta child_logwmhb_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\logwmhb_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = child_logwmhb_beta/BMI_beta
gen wald_ratio_se = child_logwmhb_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\logwmhb_loo.png", replace


//THALAMUS NORMALISED FOR TBV//

use "Y:\UKB\MR\Data Files\Child Outcomes\Child Effects File.dta", clear

* IVW method (multiplicative random effects)
tab child_thalb_beta
local snp_num = `r(N)'
mregger child_thalb_beta BMI_beta [aw=1/(child_thalb_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger child_thalb_beta BMI_beta [aw=1/(child_thalb_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian child_thalb_beta child_thalb_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal child_thalb_beta child_thalb_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thalb_results", modify
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
mregger child_thalb_beta BMI_beta [aw=1/(child_thalb_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger child_thalb_beta BMI_beta [aw=1/(child_thalb_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter child_thalb_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thalb_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest child_thalb_beta child_thalb_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thalb_forest.png", replace

* Asymmetry
mrfunnel child_thalb_beta child_thalb_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thalb_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = child_thalb_beta/BMI_beta
gen wald_ratio_se = child_thalb_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\thalb_loo.png", replace


//HIPPOCAMPUS NORMALISED FOR TBV//

use "Y:\UKB\MR\Data Files\Child Outcomes\Child Effects File.dta", clear

* IVW method (multiplicative random effects)
tab child_hippob_beta
local snp_num = `r(N)'
mregger child_hippob_beta BMI_beta [aw=1/(child_hippob_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger child_hippob_beta BMI_beta [aw=1/(child_hippob_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian child_hippob_beta child_hippob_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal child_hippob_beta child_hippob_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\hippob_results", modify
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
mregger child_hippob_beta BMI_beta [aw=1/(child_hippob_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger child_hippob_beta BMI_beta [aw=1/(child_hippob_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter child_hippob_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\hippob_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest child_hippob_beta child_hippob_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\hippob_forest.png", replace

* Asymmetry
mrfunnel child_hippob_beta child_hippob_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\hippob_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = child_hippob_beta/BMI_beta
gen wald_ratio_se = child_hippob_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\hippob_loo.png", replace


//AMYGDALA NORMALISED FOR TBV//

use "Y:\UKB\MR\Data Files\Child Outcomes\Child Effects File.dta", clear

* IVW method (multiplicative random effects)
tab child_amygb_beta
local snp_num = `r(N)'
mregger child_amygb_beta BMI_beta [aw=1/(child_amygb_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger child_amygb_beta BMI_beta [aw=1/(child_amygb_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian child_amygb_beta child_amygb_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal child_amygb_beta child_amygb_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\amygb_results", modify
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
mregger child_amygb_beta BMI_beta [aw=1/(child_amygb_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger child_amygb_beta BMI_beta [aw=1/(child_amygb_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter child_amygb_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\amygb_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest child_amygb_beta child_amygb_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\amygb_forest.png", replace

* Asymmetry
mrfunnel child_amygb_beta child_amygb_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\amygb_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = child_amygb_beta/BMI_beta
gen wald_ratio_se = child_amygb_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Child Outcomes\Results Files\amygb_loo.png", replace


**********************************************************************************ADULTHOOD EFFECTS***************************************************************************************


//GREY MATTER VOLUME NORMALISED FOR TBV//

use "Y:\UKB\MR\Data Files\Adult Outcomes\Adult Effects File.dta", clear

* IVW method (multiplicative random effects)
tab adult_gmb_beta
local snp_num = `r(N)'
mregger adult_gmb_beta BMI_beta [aw=1/(adult_gmb_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger adult_gmb_beta BMI_beta [aw=1/(adult_gmb_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian adult_gmb_beta adult_gmb_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal adult_gmb_beta adult_gmb_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\gmb_results", modify
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
mregger adult_gmb_beta BMI_beta [aw=1/(adult_gmb_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger adult_gmb_beta BMI_beta [aw=1/(adult_gmb_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter adult_gmb_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\gmb_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest adult_gmb_beta adult_gmb_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\gmb_forest.png", replace

* Asymmetry
mrfunnel adult_gmb_beta adult_gmb_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\gmb_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = adult_gmb_beta/BMI_beta
gen wald_ratio_se = adult_gmb_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\gmb_loo.png", replace


//CORTICAL SURFACE AREA NORMALISED FOR TBV//

use "Y:\UKB\MR\Data Files\Adult Outcomes\Adult Effects File.dta", clear

* IVW method (multiplicative random effects)
tab adult_areab_beta
local snp_num = `r(N)'
mregger adult_areab_beta BMI_beta [aw=1/(adult_areab_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger adult_areab_beta BMI_beta [aw=1/(adult_areab_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian adult_areab_beta adult_areab_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal adult_areab_beta adult_areab_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\areab_results", modify
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
mregger adult_areab_beta BMI_beta [aw=1/(adult_areab_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger adult_areab_beta BMI_beta [aw=1/(adult_areab_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter adult_areab_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\areab_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest adult_areab_beta adult_areab_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\areab_forest.png", replace

* Asymmetry
mrfunnel adult_areab_beta adult_areab_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\areab_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = adult_areab_beta/BMI_beta
gen wald_ratio_se = adult_areab_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\areab_loo.png", replace


//CORTICAL THICKNESS NORMALISED FOR TBV//

use "Y:\UKB\MR\Data Files\Adult Outcomes\Adult Effects File.dta", clear

* IVW method (multiplicative random effects)
tab adult_thickb_beta
local snp_num = `r(N)'
mregger adult_thickb_beta BMI_beta [aw=1/(adult_thickb_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger adult_thickb_beta BMI_beta [aw=1/(adult_thickb_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian adult_thickb_beta adult_thickb_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal adult_thickb_beta adult_thickb_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thickb_results", modify
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
mregger adult_thickb_beta BMI_beta [aw=1/(adult_thickb_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger adult_thickb_beta BMI_beta [aw=1/(adult_thickb_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter adult_thickb_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thickb_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest adult_thickb_beta adult_thickb_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thickb_forest.png", replace

* Asymmetry
mrfunnel adult_thickb_beta adult_thickb_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thickb_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = adult_thickb_beta/BMI_beta
gen wald_ratio_se = adult_thickb_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thickb_loo.png", replace


//WHITE MATTER VOLUME NORMALISED FOR TBV//

use "Y:\UKB\MR\Data Files\Adult Outcomes\Adult Effects File.dta", clear

* IVW method (multiplicative random effects)
tab adult_wmb_beta
local snp_num = `r(N)'
mregger adult_wmb_beta BMI_beta [aw=1/(adult_wmb_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger adult_wmb_beta BMI_beta [aw=1/(adult_wmb_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian adult_wmb_beta adult_wmb_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal adult_wmb_beta adult_wmb_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\wmb_results", modify
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
mregger adult_wmb_beta BMI_beta [aw=1/(adult_wmb_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger adult_wmb_beta BMI_beta [aw=1/(adult_wmb_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter adult_wmb_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\wmb_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest adult_wmb_beta adult_wmb_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\wmb_forest.png", replace

* Asymmetry
mrfunnel adult_wmb_beta adult_wmb_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\wmb_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = adult_wmb_beta/BMI_beta
gen wald_ratio_se = adult_wmb_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\wmb_loo.png", replace


//LOG-TRANSFORMED WHITE MATTER HYPERINTENSITIES NORMALISED FOR TBV//

use "Y:\UKB\MR\Data Files\Adult Outcomes\Adult Effects File.dta", clear

* IVW method (multiplicative random effects)
tab adult_logwmhb_beta
local snp_num = `r(N)'
mregger adult_logwmhb_beta BMI_beta [aw=1/(adult_logwmhb_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger adult_logwmhb_beta BMI_beta [aw=1/(adult_logwmhb_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian adult_logwmhb_beta adult_logwmhb_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal adult_logwmhb_beta adult_logwmhb_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\logwmhb_results", modify
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
mregger adult_logwmhb_beta BMI_beta [aw=1/(adult_logwmhb_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger adult_logwmhb_beta BMI_beta [aw=1/(adult_logwmhb_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter adult_logwmhb_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\logwmhb_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest adult_logwmhb_beta adult_logwmhb_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\logwmhb_forest.png", replace

* Asymmetry
mrfunnel adult_logwmhb_beta adult_logwmhb_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\logwmhb_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = adult_logwmhb_beta/BMI_beta
gen wald_ratio_se = adult_logwmhb_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\logwmhb_loo.png", replace


//THALAMUS NORMALISED FOR TBV//

use "Y:\UKB\MR\Data Files\Adult Outcomes\Adult Effects File.dta", clear

* IVW method (multiplicative random effects)
tab adult_thalb_beta
local snp_num = `r(N)'
mregger adult_thalb_beta BMI_beta [aw=1/(adult_thalb_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger adult_thalb_beta BMI_beta [aw=1/(adult_thalb_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian adult_thalb_beta adult_thalb_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal adult_thalb_beta adult_thalb_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thalb_results", modify
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
mregger adult_thalb_beta BMI_beta [aw=1/(adult_thalb_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger adult_thalb_beta BMI_beta [aw=1/(adult_thalb_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter adult_thalb_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thalb_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest adult_thalb_beta adult_thalb_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thalb_forest.png", replace

* Asymmetry
mrfunnel adult_thalb_beta adult_thalb_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thalb_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = adult_thalb_beta/BMI_beta
gen wald_ratio_se = adult_thalb_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\thalb_loo.png", replace


//HIPPOCAMPUS NORMALISED FOR TBV//

use "Y:\UKB\MR\Data Files\Adult Outcomes\Adult Effects File.dta", clear

* IVW method (multiplicative random effects)
tab adult_hippob_beta
local snp_num = `r(N)'
mregger adult_hippob_beta BMI_beta [aw=1/(adult_hippob_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger adult_hippob_beta BMI_beta [aw=1/(adult_hippob_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian adult_hippob_beta adult_hippob_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal adult_hippob_beta adult_hippob_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\hippob_results", modify
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
mregger adult_hippob_beta BMI_beta [aw=1/(adult_hippob_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger adult_hippob_beta BMI_beta [aw=1/(adult_hippob_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter adult_hippob_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\hippob_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest adult_hippob_beta adult_hippob_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\hippob_forest.png", replace

* Asymmetry
mrfunnel adult_hippob_beta adult_hippob_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\hippob_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = adult_hippob_beta/BMI_beta
gen wald_ratio_se = adult_hippob_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\hippob_loo.png", replace


//AMYGDALA NORMALISED FOR TBV//

use "Y:\UKB\MR\Data Files\Adult Outcomes\Adult Effects File.dta", clear

* IVW method (multiplicative random effects)
tab adult_amygb_beta
local snp_num = `r(N)'
mregger adult_amygb_beta BMI_beta [aw=1/(adult_amygb_se^2)], ivw 
lincom BMI_beta
local beta_ivw = `r(estimate)'
local se_ivw = `r(se)'
local p_ivw = `r(p)'
local method_ivw = "Inverse variance weighted"

* MR-Egger method (multiplicative random effects)
mregger adult_amygb_beta BMI_beta [aw=1/(adult_amygb_se^2)] 
lincom slope
local beta_mregger_slope = `r(estimate)'
local se_mregger_slope = `r(se)'
local p_mregger_slope = `r(p)'
local method_mregger_slope = "MR Egger"

* Weighted median estimator
mrmedian adult_amygb_beta adult_amygb_se BMI_beta BMI_se, weighted seed(85674)
lincom beta
local beta_median = `r(estimate)'
local se_median = `r(se)'
local p_median = `r(p)'
local method_median = "Weighted median"

* Weighted mode estimator
mrmodal adult_amygb_beta adult_amygb_se BMI_beta BMI_se, weighted seed(85674) 
lincom beta
local beta_mode = `r(estimate)'
local se_mode = `r(se)'
local p_mode = `r(p)'
local method_mode = "Weighted mode"

* Export to an Excel spreadsheet
local x=1
foreach estimator in ivw mregger_slope median mode {
	putexcel set "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\amygb_results", modify
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
mregger adult_amygb_beta BMI_beta [aw=1/(adult_amygb_se^2)] , ivw fe heterogi

* Is there evidence of horizontal pleiotropy?
mregger adult_amygb_beta BMI_beta [aw=1/(adult_amygb_se^2)]  
lincom _cons
local mregger_int = `r(estimate)'

* Scatter plot of SNP-Outcome vs SNP-Exposure effects
twoway scatter adult_amygb_beta BMI_beta || function y = (`beta_ivw')*x, range(0.1 -0.055) lcol(ltblue) || ///
function y = (`beta_mregger_slope')*x + (`mregger_int'), lcol(blue) range(0.1 -0.055) || function y = (`beta_median')*x, lcol(mint) range(0.1 -0.055) || ///
function y = (`beta_mode')*x, lcol(green) range(0.1 -0.055) legend(order(2 "IVW" 3 "MR-Egger" 4 "Weighted median" 5 "Weighted mode")) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\amygb_scatter.png", replace							

* Forest plot of Wald ratios 
mrforest adult_amygb_beta adult_amygb_se BMI_beta BMI_se, ivid(id) effect(Wald ratio) ividlabel(SNP) models(2) xlabel(-1,-.5,0,.5,1) nonote ivwopts(fe) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white)) 
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\amygb_forest.png", replace

* Asymmetry
mrfunnel adult_amygb_beta adult_amygb_se BMI_beta BMI_se, metric(invse) plotregion(fcolor(white) lcolor(white)) graphregion(fcolor(white) lcolor(white))
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\amygb_funnel.png", replace

* Perform leave-one-out analysis to check for influential individual SNPs 
set more off
gen wald_ratio = adult_amygb_beta/BMI_beta
gen wald_ratio_se = adult_amygb_se/abs(BMI_beta)  
metaninf wald_ratio wald_ratio_se, label(namevar=id)
graph export "Y:\UKB\MR\Data Files\Adult Outcomes\Results Files\amygb_loo.png", replace


log close


