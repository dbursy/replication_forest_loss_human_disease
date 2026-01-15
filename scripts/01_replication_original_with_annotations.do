/********************************************************************
Title: Replication of Does Forest Loss Increase Human Disease? Evidence from Nigeria — Commented Version
Author: Dominik Bursy, Jerico Fiestas-Flores
Date: January 2026

Description:
This do-file replicates the results from:
Berazneva, Julia, and Tanya S. Byker. 2017. "Does Forest Loss Increase Human Disease? Evidence from Nigeria." American Economic Review 107 (5): 516–21.

No substantive changes to the original analysis are made. Code for Figure 2 has been added.

Original authors' code: https://doi.org/10.3886/E113539V1
********************************************************************/

clear all
set maxvar 32000
set matsize 10000

use "../resources/113539-V1/Data_programs_readme/dta.dta", clear

*------------------------------------------------------------
* Identify LGAs containing at least one cluster per DHS wave(DHS 2008 and 2013)
*------------------------------------------------------------

foreach i of numlist 2008 2013 {
	gen _`i'_flag=DHSYEAR==`i'
	bys LGA: egen _`i'=max(_`i'_flag)
	drop _`i'_flag
}

gen _2008_only=_2008==1 & _2013==0
gen _2013_only=_2008==0 & _2013==1
gen _2008_2013=_2008==1 & _2013==1

gen LGA_overlap=.
replace LGA_overlap=1 if _2008_only==1
replace LGA_overlap=2 if _2013_only==1
replace LGA_overlap=3 if _2008_2013


label define LGAcat 1 "2008" 2 "2013" 3 "2008 & 2013" , replace
label values LGA_overlap LGAcat

*------------------------------------------------------------
* Setup  GIS variables (Lagged forest loss and luminosity, initial forest cover and soil characteristics)	
*------------------------------------------------------------

forvalues i=0/12{
gen forest_loss_`i'=.
gen luminosity_chg_`i'=.
}

gen treecover_2000=.

foreach i of numlist 2008 2013 {
	forvalues j=2001/2013 {
	local lag=`i'-`j'
	di `lag'
	if `lag'<0 {
	di "neg"
	}
	else{
	replace forest_loss_`lag'=loss_`j'_pt_`i'MEAN if DHSYEAR==`i'
	}
	}
}


foreach yr of numlist 2006/2008 {
local lag=2008-`yr'
local prev=`yr'-1
di `lag'
replace luminosity_chg_`lag'=f16`yr'_cluster_2008MEAN-f16`prev'_cluster_2008MEAN if DHSYEAR==2008
}

foreach yr of numlist 2011/2013 {
local lag=2013-`yr'
local prev=`yr'-1
di `lag'
replace luminosity_chg_`lag'=f18`yr'_cluster_2013MEAN-f18`prev'_cluster_2013MEAN if DHSYEAR==2013
}

foreach i of numlist 2008 2013 {
	replace treecover_2000=treecover_`i'MEAN if DHSYEAR==`i'
}

order forest_loss_0 forest_loss_1 forest_loss_2 forest_loss_3 forest_loss_4 forest_loss_5 forest_loss_6 forest_loss_7 luminosity_chg_0  luminosity_chg_1 luminosity_chg_2 

*------------------------------------------------------------		
* Weight De-Normalization
* When pooling multiple waves of DHS data it is recommended to de-normalize weights using populaiton estimates
* We follow this procedure based on the following guidelines: http://userforum.dhsprogram.com/index.php?t=getfile&id=4
* Population estimates are based on UN Population estimates for women aged 15-49 by year in Nigeria: https://esa.un.org/unpd/wpp/Download/Standard/ASCII/
*------------------------------------------------------------

cap drop DN_wgt stratid wave_psu
gen DN_wgt = .
replace DN_wgt=(v005/1000000)*(34442230/33385) if DHSYEAR==2008
replace DN_wgt=(v005/1000000)*(39172540/38185) if DHSYEAR==2013
* Updated Weight De-Normalization
// replace DN_wgt=(v005/1000000)*(35882027/28647) if DHSYEAR==2008
// replace DN_wgt=(v005/1000000)*(41018918/31225) if DHSYEAR==2013

egen stratid=group(DHSYEAR v022)
egen wave_psu=group(DHSYEAR v021)
svyset wave_psu [pweight=DN_wgt], strata(stratid) singleunit(certainty) 

/********************************************************************
 Table 1: The Impact of Forest Loss on Disease Incidence
 - Outcomes: fever, diarrhea, cough
********************************************************************/

*------------------------------------------------------------
* Requires: ssc install estout
*------------------------------------------------------------

cap which esttab
if _rc ssc install estout, replace

*------------------------------------------------------------
* Run models and store estimates
*------------------------------------------------------------

eststo clear

local outcomes "fever diarrhea cough"

foreach y of local outcomes {

    * Spec 1: baseline
    eststo `y'_1: svy: reg `y' forest_loss_0-forest_loss_3 if LGA_overlap==3

    * Spec 2: controls
    eststo `y'_2: svy: reg `y' forest_loss_0-forest_loss_3 ///
        i.DHSYEAR i.LGA i.month#i.region#i.DHSYEAR i.month ///
        treecover_2000 no_HH_members no_kids_under_5 time_to_water ///
        head_HH_age HH_head_edu_years rural toilet poorest firewood floor ///
        age age_resp edu_years no_child_total no_child_living livewith ///
        christian muslim yoruba igbo hausa resp_slept_net married resp_works pregnant ///
        if LGA_overlap==3

    * Spec 3: controls + extra covariates
    eststo `y'_3: svy: reg `y' forest_loss_0-forest_loss_3 ///
        i.DHSYEAR i.LGA i.month#i.region#i.DHSYEAR i.month ///
        treecover_2000 cec_ave_pt ph_ave_pt occ_ave_pt altitude ///
        luminosity_chg_0-luminosity_chg_2 ///
        no_HH_members no_kids_under_5 time_to_water ///
        head_HH_age HH_head_edu_years rural toilet poorest firewood floor ///
        age age_resp edu_years no_child_total no_child_living livewith ///
        christian muslim yoruba igbo hausa resp_slept_net married resp_works pregnant ///
        if LGA_overlap==3
}

*------------------------------------------------------------
* Create and export table
*------------------------------------------------------------

esttab fever_1 fever_2 fever_3 diarrhea_1 diarrhea_2 diarrhea_3 cough_1 cough_2 cough_3 ///
	using "../output/tables/Berazneva_Byker_Table_1.rtf", replace ///
    keep(forest_loss_0 forest_loss_1 forest_loss_2 forest_loss_3) ///
    order(forest_loss_0 forest_loss_1 forest_loss_2 forest_loss_3) ///
    cells( ///
        b(star fmt(%9.3f)) se(par fmt(%9.3f)) p(par([ ]) fmt(%9.3f)) ///
    ) ///
    collabels(none) ///
    mlabels(none) ///
    stats(N r2, labels("N" "R2") fmt(%9.0fc %9.3f)) ///
    compress

/********************************************************************
 Figure 2: Cumulative impact of forest loss (Spec 3)
********************************************************************/

*------------------------------------------------------------
* Run Spec 3 and store
*------------------------------------------------------------

local outcome "fever diarrhea cough"

foreach y of local outcome {
    svy: reg `y' forest_loss_0-forest_loss_3 ///
        i.DHSYEAR i.LGA i.month#i.region#i.DHSYEAR i.month ///
        treecover_2000 cec_ave_pt ph_ave_pt occ_ave_pt altitude ///
        luminosity_chg_0-luminosity_chg_2 ///
        no_HH_members no_kids_under_5 time_to_water ///
        head_HH_age HH_head_edu_years rural toilet poorest firewood floor ///
        age age_resp edu_years no_child_total no_child_living livewith ///
        christian muslim yoruba igbo hausa resp_slept_net married resp_works pregnant ///
        if LGA_overlap==3

    estimates store spec3_`y'
}

*------------------------------------------------------------
* Build cumulative effects dataset 
*------------------------------------------------------------

tempfile cumdata
tempname handle

postfile `handle' str10 outcome byte year double beta ll95 ul95 using `cumdata', replace

foreach y in fever diarrhea cough {

    estimates restore spec3_`y'

    forvalues t = 0/3 {

        * Build cumulative expression: b0, b0+b1, b0+b1+b2, b0+b1+b2+b3
        local expr "forest_loss_0"
        if `t' >= 1 local expr "`expr' + forest_loss_1"
        if `t' >= 2 local expr "`expr' + forest_loss_2"
        if `t' >= 3 local expr "`expr' + forest_loss_3"

        lincom `expr', level(95)

        * lincom returns r(estimate), r(lb), r(ub)
        post `handle' ("`y'") (`t') (r(estimate)) (r(lb)) (r(ub))
    }
}

postclose `handle'
use `cumdata', clear

*------------------------------------------------------------
* Panel labels/order
*------------------------------------------------------------

set scheme s1mono

gen byte out = .
replace out = 1 if outcome == "fever"
replace out = 2 if outcome == "diarrhea"
replace out = 3 if outcome == "cough"

label define outlbl 1 "Panel A. Malaria (fever)" ///
                    2 "Panel B. Diarrhea" ///
                    3 "Panel C. Respiratory (cough)", replace
label values out outlbl

label var year "Years since forest loss"
label var beta "Cumulative impact of forest loss"

*------------------------------------------------------------
* Common graph options 
*------------------------------------------------------------

local common ///
    xtitle("Years since forest loss", size(medsmall)) ///
    ytitle("Cumulative impact of forest loss", size(medsmall)) ///
    xlabel(0(1)3, nogrid labsize(small)) ///
    ylabel(-6(2)6, angle(horizontal) nogrid labsize(small)) ///
    yscale(range(-6 6)) ///
    legend(off) ///
    graphregion(color(white)) ///
    plotregion(style(none))

local layers ///
    (function y=0, range(0 3) lcolor(gs10) lwidth(vthin)) ///
    (line ll95 year, lpattern(shortdash_dot) lcolor(gs8) lwidth(medthick)) ///
    (line ul95 year, lpattern(shortdash_dot) lcolor(gs8) lwidth(medthick)) ///
    (line beta year, lcolor(black) lwidth(medthick)) ///
    (scatter beta year, mcolor(gs8) msymbol(circle) msize(vsmall))

*------------------------------------------------------------
* Panel A: Fever
*------------------------------------------------------------

twoway ///
    `layers' ///
    if outcome=="fever", ///
    title("Panel A. Malaria (fever)", size(medsmall)) ///
    `common'

graph save gA_fever.gph, replace

*------------------------------------------------------------
* Panel B: Diarrhea
*------------------------------------------------------------

twoway ///
    `layers' ///
    if outcome=="diarrhea", ///
    title("Panel B. Diarrhea", size(medsmall)) ///
    `common'

graph save gB_diarrhea.gph, replace

*------------------------------------------------------------
* Panel C: Cough
*------------------------------------------------------------

twoway ///
    `layers' ///
    if outcome=="cough", ///
    title("Panel C. Respiratory (cough)", size(medsmall)) ///
    `common'

graph save gC_cough.gph, replace

*------------------------------------------------------------
* Combine graphs and export 
*------------------------------------------------------------

graph combine gA_fever.gph gB_diarrhea.gph gC_cough.gph, ///
    col(3) ///
    graphregion(color(white))

graph export "../output/figures/Berazneva_Byker_Figure_2.png", replace width(3000)

* Remove intermediate graph files before regenerating figures
foreach f in gA_fever gB_diarrhea gC_cough {
    capture erase "`f'.gph"
}

/********************************************************************
End of the File 
********************************************************************/
