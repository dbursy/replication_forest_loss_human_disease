**********************************************************************
* Full Replication Do-File 											 *
* Bursy and Fiestas-Flores 						                     *
* A comment on Does Forest Loss Increase Human Disease?              *
* Evidence from Nigeria                    					         *
*                                            						 *
* October 2025 							                             *
*                                          						     *
**********************************************************************

**********************************************************************
* Analysis of the Original Dataset 
**********************************************************************

clear all

use "../resources/113539-V1/Data_programs_readme/dta.dta", clear
 
drop if missing(LGA)

******************************
* Identify LGAs containing at least one cluster per DHS wave(DHS 2008 and 2013)
****************************** 

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

******************************
* Setup  GIS variables (Lagged forest loss and luminosity, initial forest cover and soil characteristics)	
****************************** 

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

******************************
* ORIGINAL DENORMALIZATION
* 
* When pooling multiple waves of DHS data it is recommended to de-normalize weights using populaiton estimates
* We follow this procedure based on the following guidelines: http://userforum.dhsprogram.com/index.php?t=getfile&id=4
* Population estimates are based on UN Population estimates for women aged 15-49 by year in Nigeria: https://esa.un.org/unpd/wpp/Download/Standard/ASCII/
******************************

cap drop DN_wgt stratid wave_psu
gen DN_wgt = .
replace DN_wgt=(v005/1000000)*(34442230/33385) if DHSYEAR==2008 
// 35882027 in 2008 UN 2025 / v005 for 2008 is 28647
replace DN_wgt=(v005/1000000)*(39172540/38185) if DHSYEAR==2013 
// 41018918 in 2013 UN 2025 / v005 for 2013 is 31225

egen stratid=group(DHSYEAR v022)
egen wave_psu=group(DHSYEAR v021)
svyset wave_psu [pweight=DN_wgt], strata(stratid) singleunit(certainty) 

******************************
* Alternative DENORMALIZATION
*
* Number come from the UN DATA Portal Population Division
* https://population.un.org/dataportal/data/indicators/41/locations/566/start/2008/end/2013/table/pivotbylocation
******************************

*** Robustness 1: Alt. Denormalization with original cleaned raw data ***
gen DN_wgt_c = . // De normalized weight
replace DN_wgt_c=(v005/1000000)*(35882027/28647) if DHSYEAR==2008 // 35882027 in 2008 UN 2025 / v005 for 2008 is 28647
replace DN_wgt_c=(v005/1000000)*(41018918/31225) if DHSYEAR==2013 // 41018918 in 2013 UN 2025 / v005 for 2013 is 31225

gen var = DN_wgt - DN_wgt_c // The corrected weight is larger, that is, larger population estimates (bigger standard errors, possibly)

egen stratid_c=group(DHSYEAR v022)
egen wave_psu_c=group(DHSYEAR v021)
svyset wave_psu [pweight=DN_wgt_c], strata(stratid_c) singleunit(certainty) 

******************************
* Regressions from Full Specification of Table 1 
* 
* Column 1 
******************************
 
local outcome "fever diarrhea cough"

foreach y of local outcome {   

svy: reg `y' forest_loss_0-forest_loss_3 i.DHSYEAR i.LGA i.month#i.region#i.DHSYEAR i.month treecover_2000 cec_ave_pt ph_ave_pt occ_ave_pt altitude luminosity_chg_0-luminosity_chg_2 no_HH_members no_kids_under_5 time_to_water head_HH_age HH_head_edu_years rural toilet poorest firewood floor age  age_resp edu_years no_child_total no_child_living livewith christian muslim yoruba igbo hausa resp_slept_net married resp_works pregnant if LGA_overlap==3 
estimates  store table_column_1_`y' // added

}

estimates dir // added

******************************
* Regressions Avoiding Overspecification
* 
* Column 4
******************************

local outcome "fever diarrhea cough"

foreach y of local outcome {   

* Dropping "extra" geographical controls
svy: reg `y' forest_loss_0-forest_loss_3 i.DHSYEAR i.LGA i.month#i.region treecover_2000 altitude no_HH_members no_kids_under_5 time_to_water head_HH_age HH_head_edu_years rural toilet poorest firewood floor age  age_resp edu_years no_child_total no_child_living livewith christian muslim yoruba igbo hausa resp_slept_net married resp_works pregnant if LGA_overlap==3
estimates  store table_column_4_`y' // added // reducing the interaction (drop i.month#i.region#i.DHSYEAR) Drop soil variables as colinear with def (cec_ave_pt ph_ave_pt occ_ave_pt soil?) drop limunisty as is already colinear with deforestation ( luminosity_chg_0-luminosity_chg_2)

* Dropping "extra" sociodemographic controls
svy: reg `y' forest_loss_0-forest_loss_3 i.DHSYEAR i.LGA i.month#i.region  treecover_2000 cec_ave_pt ph_ave_pt occ_ave_pt altitude luminosity_chg_0-luminosity_chg_2 no_HH_members no_kids_under_5 time_to_water HH_head_edu_years rural poorest  age livewith christian muslim yoruba igbo hausa  married  pregnant if LGA_overlap==3
estimates  store table_column_4_`y' // added // focus only on the kid age (dropping head_HH_age, age_resp) focusing only on HH education (dropping edu_years resp_works resp_slept_net) focus on children living there (dropping no_child_total no_child_living) ///
* avoid correlation with other SES (dropping toilet firewood floor)

* Dropping "extra" geographical and sociodemographic controls
svy: reg `y' forest_loss_0-forest_loss_3 i.DHSYEAR i.LGA i.month#i.region  treecover_2000 altitude no_HH_members no_kids_under_5 time_to_water HH_head_edu_years rural poorest  age livewith christian muslim yoruba igbo hausa  married  pregnant if LGA_overlap==3
estimates  store table_column_4_`y' // added // Drop all previous ones

}

estimates dir // added

******************************
* Save estimates
******************************

/*
local outcomes fever diarrhea cough

foreach y of local outcomes {
    foreach c in 1 4 {
        estimates use table_column_`c'_`y'
        estimates save "../output/table_column_`c'_`y'.ster", replace
    }
}
*/

**********************************************************************
* Analysis of the Replication Dataset – Survey Wave 2008 to 2013     *
**********************************************************************

use "../output/dta_replication_full.dta", clear

******************************
* Identify LGAs containing at least one cluster per DHS wave(DHS 2008 and 2013)
****************************** 

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

******************************
* Setup  GIS variables (Lagged forest loss and luminosity, initial forest cover and soil characteristics)	
****************************** 

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

******************************
* ORIGINAL DENORMALIZATION
* 
* When pooling multiple waves of DHS data it is recommended to de-normalize weights using populaiton estimates
* We follow this procedure based on the following guidelines: http://userforum.dhsprogram.com/index.php?t=getfile&id=4
* Population estimates are based on UN Population estimates for women aged 15-49 by year in Nigeria: https://esa.un.org/unpd/wpp/Download/Standard/ASCII/
******************************

cap drop DN_wgt stratid wave_psu
gen DN_wgt = .
replace DN_wgt=(v005/1000000)*(34442230/33385) if DHSYEAR==2008 
// 35882027 in 2008 UN 2025 / v005 for 2008 is 28647
replace DN_wgt=(v005/1000000)*(39172540/38185) if DHSYEAR==2013 
// 41018918 in 2013 UN 2025 / v005 for 2013 is 31225

egen stratid=group(DHSYEAR v022)
egen wave_psu=group(DHSYEAR v021)
svyset wave_psu [pweight=DN_wgt], strata(stratid) singleunit(certainty) 

******************************
* Alternative DENORMALIZATION
*
* Number come from the UN DATA Portal Population Division
* https://population.un.org/dataportal/data/indicators/41/locations/566/start/2008/end/2013/table/pivotbylocation
******************************

*** Robustness 1: Alt. Denormalization with original cleaned raw data ***
gen DN_wgt_c = . // De normalized weight
replace DN_wgt_c=(v005/1000000)*(35882027/28647) if DHSYEAR==2008 // 35882027 in 2008 UN 2025 / v005 for 2008 is 28647
replace DN_wgt_c=(v005/1000000)*(41018918/31225) if DHSYEAR==2013 // 41018918 in 2013 UN 2025 / v005 for 2013 is 31225

gen var = DN_wgt - DN_wgt_c // The corrected weight is larger, that is, larger population estimates (bigger standard errors, possibly)

egen stratid_c=group(DHSYEAR v022)
egen wave_psu_c=group(DHSYEAR v021)
svyset wave_psu [pweight=DN_wgt_c], strata(stratid_c) singleunit(certainty) 

******************************
* Regressions from Full Specification of Table 1 
* 
* Column 2 
******************************
 
local outcome "fever diarrhea cough"

foreach y of local outcome {   

svy: reg `y' forest_loss_0-forest_loss_3 i.DHSYEAR i.LGA i.month#i.region#i.DHSYEAR i.month treecover_2000 cec_ave_pt ph_ave_pt occ_ave_pt altitude luminosity_chg_0-luminosity_chg_2 no_HH_members no_kids_under_5 time_to_water head_HH_age HH_head_edu_years rural toilet poorest firewood floor age  age_resp edu_years no_child_total no_child_living livewith christian muslim yoruba igbo hausa resp_slept_net married resp_works pregnant if LGA_overlap==3 
estimates  store table_column_2_`y' // added

}

estimates dir // added

******************************
* Regressions Avoiding Overspecification
* 
* Column 5 
******************************

local outcome "fever diarrhea cough"

foreach y of local outcome {   

* Dropping "extra" geographical controls
svy: reg `y' forest_loss_0-forest_loss_3 i.DHSYEAR i.LGA i.month#i.region treecover_2000 altitude no_HH_members no_kids_under_5 time_to_water head_HH_age HH_head_edu_years rural toilet poorest firewood floor age  age_resp edu_years no_child_total no_child_living livewith christian muslim yoruba igbo hausa resp_slept_net married resp_works pregnant if LGA_overlap==3
estimates  store table_column_5_`y' // added // reducing the interaction (drop i.month#i.region#i.DHSYEAR) Drop soil variables as colinear with def (cec_ave_pt ph_ave_pt occ_ave_pt soil?) drop limunisty as is already colinear with deforestation ( luminosity_chg_0-luminosity_chg_2)

* Dropping "extra" sociodemographic controls
svy: reg `y' forest_loss_0-forest_loss_3 i.DHSYEAR i.LGA i.month#i.region  treecover_2000 cec_ave_pt ph_ave_pt occ_ave_pt altitude luminosity_chg_0-luminosity_chg_2 no_HH_members no_kids_under_5 time_to_water HH_head_edu_years rural poorest  age livewith christian muslim yoruba igbo hausa  married  pregnant if LGA_overlap==3
estimates  store table_column_5_`y' // added // focus only on the kid age (dropping head_HH_age, age_resp) focusing only on HH education (dropping edu_years resp_works resp_slept_net) focus on children living there (dropping no_child_total no_child_living) ///
* avoid correlation with other SES (dropping toilet firewood floor)

* Dropping "extra" geographical and sociodemographic controls
svy: reg `y' forest_loss_0-forest_loss_3 i.DHSYEAR i.LGA i.month#i.region  treecover_2000 altitude no_HH_members no_kids_under_5 time_to_water HH_head_edu_years rural poorest  age livewith christian muslim yoruba igbo hausa  married  pregnant if LGA_overlap==3
estimates  store table_column_5_`y' // added // Drop all previous ones

}

estimates dir // added

******************************
* Save estimates
******************************

/*
local outcomes fever diarrhea cough

foreach y of local outcomes {
    foreach c in 2 5 {
        estimates use table_column_`c'_`y'
        estimates save "../output/table_column_`c'_`y'.ster", replace
    }
}
*/

**********************************************************************
* Analysis of the Replication Dataset – Survey Wave 2008 to 2018     *
**********************************************************************

use "../output/dta_replication_full_extended.dta", clear

******************************
* Identify LGAs containing at least one cluster per DHS wave(DHS 2008 and 2013)
****************************** 


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


* Add new dummy to include 2018
foreach i of numlist 2008 2013 2018 { // added
	gen _`i'_flag=DHSYEAR==`i'
	bys LGA: egen _`i'_2=max(_`i'_flag)
	drop _`i'_flag
}


gen _2008_only_2=_2008_2==1 & _2013_2==0 & _2018_2==0 // added
gen _2013_only_2=_2008_2==0 & _2013_2==1 & _2018_2==0 // added
gen _2018_only_2=_2008_2==0 & _2013_2==0 & _2018_2==1 // added
gen _2008_2013_2=_2008_2==1 & _2013_2==1 & _2018_2==0 // added
gen _2008_2018_2=_2008_2==1 & _2013_2==0 & _2018_2==1 // added
gen _2013_2018_2=_2008_2==0 & _2013_2==1 & _2018_2==1 // added
gen _ALL	  =_2008_2==1  & _2013_2==1 & _2018_2==1 // added

gen LGA_overlap_2 = .
replace LGA_overlap_2 = 1 if _2008_only_2 == 1
replace LGA_overlap_2 = 2 if _2013_only_2 == 1
replace LGA_overlap_2 = 3 if _2018_only_2 == 1
replace LGA_overlap_2 = 4 if _2008_2013_2 == 1
replace LGA_overlap_2 = 5 if _2008_2018_2 == 1
replace LGA_overlap_2 = 6 if _2013_2018_2 == 1
replace LGA_overlap_2 = 7 if _ALL == 1

* Label the categories
label define LGAcat_2 1 "2008 only" 2 "2013 only" 3 "2018 only" 4 "2008 & 2013" 5 "2008 & 2018" 6 "2013 & 2018" 7 "2008, 2013 & 2018", replace
label values LGA_overlap_2 LGAcat_2

******************************
* Setup  GIS variables (Lagged forest loss and luminosity, initial forest cover and soil characteristics)	
****************************** 

forvalues i=0/17{ // added
	gen forest_loss_`i'=.
	gen luminosity_chg_`i'=.
}

gen treecover_2000=.

foreach i of numlist 2008 2013 2018 { 
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
	
	foreach yr of numlist 3/17 { // Eliminat surpluss
	drop luminosity_chg_`yr'
	}
	
	foreach i of numlist 2008 2013 2018 {  // added
		replace treecover_2000=treecover_`i'MEAN if DHSYEAR==`i'
	}

order forest_loss_0 forest_loss_1 forest_loss_2 forest_loss_3 forest_loss_4 forest_loss_5 forest_loss_6 forest_loss_7 forest_loss_8 forest_loss_9 forest_loss_10 forest_loss_11 forest_loss_12 forest_loss_13 forest_loss_14 forest_loss_15 forest_loss_16 forest_loss_17 luminosity_chg_0  luminosity_chg_1 luminosity_chg_2

******************************
* ORIGINAL DENORMALIZATION
* 
* When pooling multiple waves of DHS data it is recommended to de-normalize weights using populaiton estimates
* We follow this procedure based on the following guidelines: http://userforum.dhsprogram.com/index.php?t=getfile&id=4
* Population estimates are based on UN Population estimates for women aged 15-49 by year in Nigeria: https://esa.un.org/unpd/wpp/Download/Standard/ASCII/
******************************

cap drop DN_wgt stratid wave_psu
gen DN_wgt = .
replace DN_wgt=(v005/1000000)*(34442230/33385) if DHSYEAR==2008 // 35882027 in 2008 UN 2025 / v005 for 2008 is 28647
replace DN_wgt=(v005/1000000)*(39172540/38185) if DHSYEAR==2013 // 41018918 in 2013 UN 2025 / v005 for 2013 is 31225

egen stratid=group(DHSYEAR v022)
egen wave_psu=group(DHSYEAR v021)
svyset wave_psu [pweight=DN_wgt], strata(stratid) singleunit(certainty) 

******************************
* Alternative DENORMALIZATION
*
* Number come from the UN DATA Portal Population Division
* https://population.un.org/dataportal/data/indicators/41/locations/566/start/2008/end/2013/table/pivotbylocation
****************************** 

gen DN_wgt_c = . // De normalized weight
replace DN_wgt_c=(v005/1000000)*(35882027/28647) if DHSYEAR==2008 
replace DN_wgt_c=(v005/1000000)*(41018918/31225) if DHSYEAR==2013 
replace DN_wgt_c=(v005/1000000)*(47146386/33775) if DHSYEAR==2018 
 **** used this for female at reprodutive age: https://population.un.org/dataportal/ -- FRA 2008 - 2018

gen var = DN_wgt - DN_wgt_c // The corrected weight is larger, that is, larger population estimates (bigger standard errors, possibly)

egen stratid_c=group(DHSYEAR v022)
egen wave_psu_c=group(DHSYEAR v021)
svyset wave_psu [pweight=DN_wgt_c], strata(stratid_c) singleunit(certainty) 

******************************
* Regressions from Full Specification of Table 1 
* 
* Column 3 
******************************


local outcome "fever diarrhea cough"

foreach y of local outcome {   

svy: reg `y' forest_loss_0-forest_loss_3 if  LGA_overlap_2==7  
estimates store table_column_3_`y' // added

svy: reg `y' forest_loss_0-forest_loss_3 i.DHSYEAR i.LGA i.month#i.region#i.DHSYEAR i.month treecover_2000 no_HH_members no_kids_under_5 time_to_water head_HH_age HH_head_edu_years rural toilet poorest firewood floor age  age_resp edu_years no_child_total no_child_living livewith christian muslim yoruba igbo hausa resp_slept_net married resp_works pregnant if LGA_overlap_2==7
estimates  store table_column_3_`y' // added

svy: reg `y' forest_loss_0-forest_loss_3 i.DHSYEAR i.LGA i.month#i.region#i.DHSYEAR i.month treecover_2000 cec_ave_pt ph_ave_pt occ_ave_pt altitude luminosity_chg_0-luminosity_chg_2 no_HH_members no_kids_under_5 time_to_water head_HH_age HH_head_edu_years rural toilet poorest firewood floor age  age_resp edu_years no_child_total no_child_living livewith christian muslim yoruba igbo hausa resp_slept_net married resp_works pregnant if LGA_overlap_2==7 
estimates  store table_column_3_`y' // added


}

estimates dir

******************************
* Regressions Avoiding Overspecification
* 
* Column 6
******************************


local outcome "fever diarrhea cough"

foreach y of local outcome {   

* Dropping "extra" geographical controls
svy: reg `y' forest_loss_0-forest_loss_3 i.DHSYEAR i.LGA i.month#i.region treecover_2000 altitude no_HH_members no_kids_under_5 time_to_water head_HH_age HH_head_edu_years rural toilet poorest firewood floor age  age_resp edu_years no_child_total no_child_living livewith christian muslim yoruba igbo hausa resp_slept_net married resp_works pregnant if LGA_overlap_2==7 
estimates  store table_column_6_`y' // added // reducing the interaction (drop i.month#i.region#i.DHSYEAR) Drop soil variables as colinear with def (cec_ave_pt ph_ave_pt occ_ave_pt soil?) drop limunisty as is already colinear with deforestation ( luminosity_chg_0-luminosity_chg_2)

* Dropping "extra" sociodemographic controls
svy: reg `y' forest_loss_0-forest_loss_3 i.DHSYEAR i.LGA i.month#i.region  treecover_2000 cec_ave_pt ph_ave_pt occ_ave_pt altitude luminosity_chg_0-luminosity_chg_2 no_HH_members no_kids_under_5 time_to_water HH_head_edu_years rural poorest  age livewith christian muslim yoruba igbo hausa  married  pregnant if LGA_overlap_2==7 
estimates  store table_column_6_`y' // added // focus only on the kid age (dropping head_HH_age, age_resp) focusing only on HH education (dropping edu_years resp_works resp_slept_net) focus on children living there (dropping no_child_total no_child_living) ///
* avoid correlation with other SES (dropping toilet firewood floor)

* Dropping "extra" geographical and sociodemographic controls
svy: reg `y' forest_loss_0-forest_loss_3 i.DHSYEAR i.LGA i.month#i.region  treecover_2000 altitude no_HH_members no_kids_under_5 time_to_water HH_head_edu_years rural poorest  age livewith christian muslim yoruba igbo hausa  married  pregnant if LGA_overlap_2==7 
estimates  store table_column_6_`y' // added // Drop all previous ones

}

estimates dir

******************************
* Save estimates
******************************

/*
local outcomes fever diarrhea cough

foreach y of local outcomes {
    foreach c in 3 6 {
        estimates use table_column_`c'_`y'
        estimates save "../output/table_column_`c'_`y'.ster", replace
    }
}
*/

**********************************************************************
* Export Results											         *
**********************************************************************

* clear all

******************************
* Load estimates
******************************

/*
local outcomes fever diarrhea cough

foreach y of local outcomes {
    foreach c in 1 2 3 4 5 6 {
        estimates use "../output/table_column_`c'_`y'.ster"
    }
}
*/

* estimates dir

******************************
* Export estimates
******************************

estimates dir

* change to . or , depending on the system and country
cd "../output/tables"

* Fever regressions combined
estout table_column_1_fever table_column_2_fever table_column_3_fever ///
       table_column_4_fever table_column_5_fever table_column_6_fever ///
       using "table_fever.xls", ///
       cells(b(star fmt(%7.3f)) se(par fmt(%7.3f))) ///
       stats(N r2, fmt(%17.3f)) ///
       keep(forest_* _cons) ///
       replace

******************************

estout table_column_1_fever table_column_2_fever table_column_3_fever ///
       table_column_4_fever table_column_5_fever table_column_6_fever ///
       using "table_fever.tex", ///
       style(tex) ///
       cells(b(star fmt(%7.3f)) se(par fmt(%7.3f))) ///
       stats(N r2, fmt(%17.3f) labels("Observations" "R-squared")) ///
       keep(forest_* _cons) ///
       replace ///
       label ///
       prehead("\begin{table}[htbp]\centering" ///
               "\caption{Fever regressions combined}" ///
               "\begin{tabular}{lcccccc}" ///
               "\toprule" ///
               "& \multicolumn{3}{c}{Full Specification} & \multicolumn{3}{c}{Overspecification Correction} \\" ///
               "\cmidrule(lr){2-4} \cmidrule(lr){5-7}" ///
               "Malaria (Fever) & (1) & (2) & (3) & (4) & (5) & (6) \\") ///
       postfoot("\bottomrule" ///
                "\end{tabular}" ///
                "\end{table}")

estout table_column_1_diarrhea table_column_2_diarrhea table_column_3_diarrhea ///
       table_column_4_diarrhea table_column_5_diarrhea table_column_6_diarrhea ///
       using "table_diarrhea.tex", ///
       style(tex) ///
       cells(b(star fmt(%7.3f)) se(par fmt(%7.3f))) ///
       stats(N r2, fmt(%17.3f) labels("Observations" "R-squared")) ///
       keep(forest_* _cons) ///
       replace ///
       label ///
       prehead("\begin{table}[htbp]\centering" ///
               "\caption{Fever regressions combined}" ///
               "\begin{tabular}{lcccccc}" ///
               "\toprule" ///
               "& \multicolumn{3}{c}{Full Specification} & \multicolumn{3}{c}{Overspecification Correction} \\" ///
               "\cmidrule(lr){2-4} \cmidrule(lr){5-7}" ///
               "Malaria (Fever) & (1) & (2) & (3) & (4) & (5) & (6) \\") ///
       postfoot("\bottomrule" ///
                "\end{tabular}" ///
                "\end{table}")

estout table_column_1_cough table_column_2_cough table_column_3_cough ///
       table_column_4_cough table_column_5_cough table_column_6_cough ///
       using "table_cough.tex", ///
       style(tex) ///
       cells(b(star fmt(%7.3f)) se(par fmt(%7.3f))) ///
       stats(N r2, fmt(%17.3f) labels("Observations" "R-squared")) ///
       keep(forest_* _cons) ///
       replace ///
       label ///
       prehead("\begin{table}[htbp]\centering" ///
               "\caption{Fever regressions combined}" ///
               "\begin{tabular}{lcccccc}" ///
               "\toprule" ///
               "& \multicolumn{3}{c}{Full Specification} & \multicolumn{3}{c}{Overspecification Correction} \\" ///
               "\cmidrule(lr){2-4} \cmidrule(lr){5-7}" ///
               "Respiratory (Cough) & (1) & (2) & (3) & (4) & (5) & (6) \\") ///
       postfoot("\bottomrule" ///
                "\end{tabular}" ///
                "\end{table}")

**********************************************************************
* End of the File 												     *
**********************************************************************
