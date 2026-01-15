/********************************************************************
Title: Replication of Does Forest Loss Increase Human Disease? Evidence from Nigeria — Using Replicated and Extended Data
Author: Dominik Bursy, Jerico Fiestas-Flores
Date: January 2026

Description:
This do-file replicates and extents the main results of:
Berazneva, Julia, and Tanya S. Byker. 2017. "Does Forest Loss Increase Human Disease? Evidence from Nigeria." American Economic Review 107 (5): 516–21.

The code is based on the authors' original replication script and uses both the original and an extended dataset.

Original authors' code: https://doi.org/10.3886/E113539V1
********************************************************************/

/********************************************************************
Analysis of the Original Dataset 
********************************************************************/

clear all

use "../resources/113539-V1/Data_programs_readme/dta.dta", clear
 
drop if missing(LGA)

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
* 
* When pooling multiple waves of DHS data it is recommended to de-normalize weights using populaiton estimates
* We follow this procedure based on the following guidelines: http://userforum.dhsprogram.com/index.php?t=getfile&id=4
* Population estimates are based on UN Population estimates for women aged 15-49 by year in Nigeria: https://esa.un.org/unpd/wpp/Download/Standard/ASCII/
*------------------------------------------------------------

cap drop DN_wgt stratid wave_psu
gen DN_wgt = .
replace DN_wgt=(v005/1000000)*(34442230/33385) if DHSYEAR==2008 
// 35882027 in 2008 UN 2025 / v005 for 2008 is 28647
replace DN_wgt=(v005/1000000)*(39172540/38185) if DHSYEAR==2013 
// 41018918 in 2013 UN 2025 / v005 for 2013 is 31225

egen stratid=group(DHSYEAR v022)
egen wave_psu=group(DHSYEAR v021)
svyset wave_psu [pweight=DN_wgt], strata(stratid) singleunit(certainty) 

*------------------------------------------------------------
* Updated Weight De-Normalization
*
* Numbers come from the UN DATA Portal Population Division
* https://population.un.org/dataportal/data/indicators/41/locations/566/start/2008/end/2013/table/pivotbylocation
*------------------------------------------------------------

*Robustness 1: Updated De-Normalization with original cleaned raw data
gen DN_wgt_c = . // De normalized weight
replace DN_wgt_c=(v005/1000000)*(35882027/28647) if DHSYEAR==2008 
// 35882027 in 2008 UN 2025 / v005 for 2008 is 28647
replace DN_wgt_c=(v005/1000000)*(41018918/31225) if DHSYEAR==2013 
// 41018918 in 2013 UN 2025 / v005 for 2013 is 31225

gen DN_wgt_comparison = DN_wgt - DN_wgt_c 
// The updated weight is larger, that is, larger population estimates (bigger standard errors, possibly)

egen stratid_c=group(DHSYEAR v022)
egen wave_psu_c=group(DHSYEAR v021)
svyset wave_psu [pweight=DN_wgt_c], strata(stratid_c) singleunit(certainty) 

*------------------------------------------------------------
* Full Specification 3 of Berazneva and Byker Table 1 
* 
* Replication Report Column 1 
*------------------------------------------------------------
 
local outcome "fever diarrhea cough"

foreach y of local outcome {   

svy: reg `y' forest_loss_0-forest_loss_3 i.DHSYEAR i.LGA i.month#i.region#i.DHSYEAR i.month treecover_2000 cec_ave_pt ph_ave_pt occ_ave_pt altitude luminosity_chg_0-luminosity_chg_2 no_HH_members no_kids_under_5 time_to_water head_HH_age HH_head_edu_years rural toilet poorest firewood floor age  age_resp edu_years no_child_total no_child_living livewith christian muslim yoruba igbo hausa resp_slept_net married resp_works pregnant if LGA_overlap==3 
estimates  store table_column_1_`y' // added

}

estimates dir // added

*------------------------------------------------------------
* Reduced Specification 3 of Berazneva and Byker Table 1 
* 
* Replication Report Column 4
*------------------------------------------------------------

local outcome "fever diarrhea cough"

foreach y of local outcome {   

* Dropping "extra" geographical controls
svy: reg `y' forest_loss_0-forest_loss_3 i.DHSYEAR i.LGA i.month#i.region treecover_2000 altitude no_HH_members no_kids_under_5 time_to_water head_HH_age HH_head_edu_years rural toilet poorest firewood floor age  age_resp edu_years no_child_total no_child_living livewith christian muslim yoruba igbo hausa resp_slept_net married resp_works pregnant if LGA_overlap==3
estimates  store table_column_4_`y' // added // reducing the interaction (drop i.month#i.region#i.DHSYEAR) Drop soil variables as colinear with def (cec_ave_pt ph_ave_pt occ_ave_pt soil?) drop limunisty as is potentially colinear with deforestation ( luminosity_chg_0-luminosity_chg_2)

* Dropping "extra" sociodemographic controls
svy: reg `y' forest_loss_0-forest_loss_3 i.DHSYEAR i.LGA i.month#i.region  treecover_2000 cec_ave_pt ph_ave_pt occ_ave_pt altitude luminosity_chg_0-luminosity_chg_2 no_HH_members no_kids_under_5 time_to_water HH_head_edu_years rural poorest  age livewith christian muslim yoruba igbo hausa  married  pregnant if LGA_overlap==3
estimates  store table_column_4_`y' // added // focus only on the kid age (dropping head_HH_age, age_resp) focusing only on HH education (dropping edu_years resp_works resp_slept_net) focus on children living there (dropping no_child_total no_child_living) ///
* avoid correlation with other SES (dropping toilet firewood floor)

* Dropping "extra" geographical and sociodemographic controls
svy: reg `y' forest_loss_0-forest_loss_3 i.DHSYEAR i.LGA i.month#i.region  treecover_2000 altitude no_HH_members no_kids_under_5 time_to_water HH_head_edu_years rural poorest  age livewith christian muslim yoruba igbo hausa  married  pregnant if LGA_overlap==3
estimates  store table_column_4_`y' // added // Drop all previous ones

}

estimates dir // added

*------------------------------------------------------------
* Save estimates (Optional)
*------------------------------------------------------------

/*
local outcomes fever diarrhea cough

foreach y of local outcomes {
    foreach c in 1 4 {
        estimates use table_column_`c'_`y'
        estimates save "../output/table_column_`c'_`y'.ster", replace
    }
}
*/

/********************************************************************
* Analysis of the Replicated Raw Dataset – Survey Wave 2008 to 2013
********************************************************************/

use "../output/dta_replication_full.dta", clear

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
* 
* When pooling multiple waves of DHS data it is recommended to de-normalize weights using populaiton estimates
* We follow this procedure based on the following guidelines: http://userforum.dhsprogram.com/index.php?t=getfile&id=4
* Population estimates are based on UN Population estimates for women aged 15-49 by year in Nigeria: https://esa.un.org/unpd/wpp/Download/Standard/ASCII/
*------------------------------------------------------------

cap drop DN_wgt stratid wave_psu
gen DN_wgt = .
replace DN_wgt=(v005/1000000)*(34442230/33385) if DHSYEAR==2008 
// 35882027 in 2008 UN 2025 / v005 for 2008 is 28647
replace DN_wgt=(v005/1000000)*(39172540/38185) if DHSYEAR==2013 
// 41018918 in 2013 UN 2025 / v005 for 2013 is 31225

egen stratid=group(DHSYEAR v022)
egen wave_psu=group(DHSYEAR v021)
svyset wave_psu [pweight=DN_wgt], strata(stratid) singleunit(certainty) 

*------------------------------------------------------------
* Updated Weight De-Normalization
*
* Number come from the UN DATA Portal Population Division
* https://population.un.org/dataportal/data/indicators/41/locations/566/start/2008/end/2013/table/pivotbylocation
*------------------------------------------------------------

* Robustness 1: Updated De-Normalization with original cleaned raw data ***
gen DN_wgt_c = . // De normalized weight
replace DN_wgt_c=(v005/1000000)*(35882027/28647) if DHSYEAR==2008 
// 35882027 in 2008 UN 2025 / v005 for 2008 is 28647
replace DN_wgt_c=(v005/1000000)*(41018918/31225) if DHSYEAR==2013 
// 41018918 in 2013 UN 2025 / v005 for 2013 is 31225

gen var = DN_wgt - DN_wgt_c 
// The updated weight is larger, that is, larger population estimates (bigger standard errors, possibly)

egen stratid_c=group(DHSYEAR v022)
egen wave_psu_c=group(DHSYEAR v021)
svyset wave_psu [pweight=DN_wgt_c], strata(stratid_c) singleunit(certainty) 

*------------------------------------------------------------
* Full Specification 3 of Berazneva and Byker Table 1 
* 
* Replication Report Column 2
*------------------------------------------------------------

local outcome "fever diarrhea cough"

foreach y of local outcome {   

svy: reg `y' forest_loss_0-forest_loss_3 i.DHSYEAR i.LGA i.month#i.region#i.DHSYEAR i.month treecover_2000 cec_ave_pt ph_ave_pt occ_ave_pt altitude luminosity_chg_0-luminosity_chg_2 no_HH_members no_kids_under_5 time_to_water head_HH_age HH_head_edu_years rural toilet poorest firewood floor age  age_resp edu_years no_child_total no_child_living livewith christian muslim yoruba igbo hausa resp_slept_net married resp_works pregnant if LGA_overlap==3 
estimates  store table_column_2_`y' // added

}

estimates dir // added

*------------------------------------------------------------
* Reduced Specification 3 of Berazneva and Byker Table 1 
* 
* Replication Report Column 5
*------------------------------------------------------------

local outcome "fever diarrhea cough"

foreach y of local outcome {   

* Dropping "extra" geographical controls
svy: reg `y' forest_loss_0-forest_loss_3 i.DHSYEAR i.LGA i.month#i.region treecover_2000 altitude no_HH_members no_kids_under_5 time_to_water head_HH_age HH_head_edu_years rural toilet poorest firewood floor age  age_resp edu_years no_child_total no_child_living livewith christian muslim yoruba igbo hausa resp_slept_net married resp_works pregnant if LGA_overlap==3
estimates  store table_column_5_`y' // added // reducing the interaction (drop i.month#i.region#i.DHSYEAR) Drop soil variables as colinear with def (cec_ave_pt ph_ave_pt occ_ave_pt soil?) drop limunisty as is potentially colinear with deforestation ( luminosity_chg_0-luminosity_chg_2)

* Dropping "extra" sociodemographic controls
svy: reg `y' forest_loss_0-forest_loss_3 i.DHSYEAR i.LGA i.month#i.region  treecover_2000 cec_ave_pt ph_ave_pt occ_ave_pt altitude luminosity_chg_0-luminosity_chg_2 no_HH_members no_kids_under_5 time_to_water HH_head_edu_years rural poorest  age livewith christian muslim yoruba igbo hausa  married  pregnant if LGA_overlap==3
estimates  store table_column_5_`y' // added // focus only on the kid age (dropping head_HH_age, age_resp) focusing only on HH education (dropping edu_years resp_works resp_slept_net) focus on children living there (dropping no_child_total no_child_living) ///
* avoid correlation with other SES (dropping toilet firewood floor)

* Dropping "extra" geographical and sociodemographic controls
svy: reg `y' forest_loss_0-forest_loss_3 i.DHSYEAR i.LGA i.month#i.region  treecover_2000 altitude no_HH_members no_kids_under_5 time_to_water HH_head_edu_years rural poorest  age livewith christian muslim yoruba igbo hausa  married  pregnant if LGA_overlap==3
estimates  store table_column_5_`y' // added // Drop all previous ones

}

estimates dir // added

*------------------------------------------------------------
* Save estimates (Optional)
*------------------------------------------------------------

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
* Analysis of the Replucated and Extended Raw Dataset – Survey Wave 2008 to 2018     *
**********************************************************************

use "../output/dta_replication_full_extended.dta", clear

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

*------------------------------------------------------------
* Setup  GIS variables (Lagged forest loss and luminosity, initial forest cover and soil characteristics)	
*------------------------------------------------------------ 

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

*------------------------------------------------------------
* Weight De-Normalization
* 
* When pooling multiple waves of DHS data it is recommended to de-normalize weights using populaiton estimates
* We follow this procedure based on the following guidelines: http://userforum.dhsprogram.com/index.php?t=getfile&id=4
* Population estimates are based on UN Population estimates for women aged 15-49 by year in Nigeria: https://esa.un.org/unpd/wpp/Download/Standard/ASCII/
*------------------------------------------------------------

cap drop DN_wgt stratid wave_psu
gen DN_wgt = .
replace DN_wgt=(v005/1000000)*(34442230/33385) if DHSYEAR==2008 
// 35882027 in 2008 UN 2025 / v005 for 2008 is 28647
replace DN_wgt=(v005/1000000)*(39172540/38185) if DHSYEAR==2013 
// 41018918 in 2013 UN 2025 / v005 for 2013 is 31225

egen stratid=group(DHSYEAR v022)
egen wave_psu=group(DHSYEAR v021)
svyset wave_psu [pweight=DN_wgt], strata(stratid) singleunit(certainty) 

*------------------------------------------------------------
* Alternative Weight De-Normalization
*
* Number come from the UN DATA Portal Population Division
* https://population.un.org/dataportal/data/indicators/41/locations/566/start/2008/end/2013/table/pivotbylocation
*------------------------------------------------------------ 

gen DN_wgt_c = . // De-Normalized weight
replace DN_wgt_c=(v005/1000000)*(35882027/28647) if DHSYEAR==2008 
replace DN_wgt_c=(v005/1000000)*(41018918/31225) if DHSYEAR==2013 
replace DN_wgt_c=(v005/1000000)*(47146386/33775) if DHSYEAR==2018 
* use this for female at reprodutive age:
* https://population.un.org/dataportal/ -- FRA 2008 - 2018

gen var = DN_wgt - DN_wgt_c 
// The corrected weight is larger, that is, larger population estimates (bigger standard errors, possibly)

egen stratid_c=group(DHSYEAR v022)
egen wave_psu_c=group(DHSYEAR v021)
svyset wave_psu [pweight=DN_wgt_c], strata(stratid_c) singleunit(certainty) 

*------------------------------------------------------------
* Full Specification 3 of Berazneva and Byker Table 1 
* 
* Replication Report Column 3
*------------------------------------------------------------

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

*------------------------------------------------------------
* Reduced Specification 3 of Berazneva and Byker Table 1 
* 
* Replication Report Column 6
*------------------------------------------------------------

local outcome "fever diarrhea cough"

foreach y of local outcome {   

* Dropping "extra" geographical controls
svy: reg `y' forest_loss_0-forest_loss_3 i.DHSYEAR i.LGA i.month#i.region treecover_2000 altitude no_HH_members no_kids_under_5 time_to_water head_HH_age HH_head_edu_years rural toilet poorest firewood floor age  age_resp edu_years no_child_total no_child_living livewith christian muslim yoruba igbo hausa resp_slept_net married resp_works pregnant if LGA_overlap_2==7 
estimates  store table_column_6_`y' // added // reducing the interaction (drop i.month#i.region#i.DHSYEAR) Drop soil variables as colinear with def (cec_ave_pt ph_ave_pt occ_ave_pt soil?) drop limunisty as is potentially colinear with deforestation ( luminosity_chg_0-luminosity_chg_2)

* Dropping "extra" sociodemographic controls
svy: reg `y' forest_loss_0-forest_loss_3 i.DHSYEAR i.LGA i.month#i.region  treecover_2000 cec_ave_pt ph_ave_pt occ_ave_pt altitude luminosity_chg_0-luminosity_chg_2 no_HH_members no_kids_under_5 time_to_water HH_head_edu_years rural poorest  age livewith christian muslim yoruba igbo hausa  married  pregnant if LGA_overlap_2==7 
estimates  store table_column_6_`y' // added // focus only on the kid age (dropping head_HH_age, age_resp) focusing only on HH education (dropping edu_years resp_works resp_slept_net) focus on children living there (dropping no_child_total no_child_living) ///
* avoid correlation with other SES (dropping toilet firewood floor)

* Dropping "extra" geographical and sociodemographic controls
svy: reg `y' forest_loss_0-forest_loss_3 i.DHSYEAR i.LGA i.month#i.region  treecover_2000 altitude no_HH_members no_kids_under_5 time_to_water HH_head_edu_years rural poorest  age livewith christian muslim yoruba igbo hausa  married  pregnant if LGA_overlap_2==7 
estimates  store table_column_6_`y' // added // Drop all previous ones

}

*------------------------------------------------------------
* Save estimates (Optional)
*------------------------------------------------------------

/*
local outcomes fever diarrhea cough

foreach y of local outcomes {
    foreach c in 3 6 {
        estimates use table_column_`c'_`y'
        estimates save "../output/table_column_`c'_`y'.ster", replace
    }
}
*/

/********************************************************************
Export Results
********************************************************************/

// clear all

*------------------------------------------------------------
* Load estimates
*------------------------------------------------------------

/*
local outcomes fever diarrhea cough

foreach y of local outcomes {
    foreach c in 1 2 3 4 5 6 {
        estimates use "../output/table_column_`c'_`y'.ster"
    }
}
*/

*------------------------------------------------------------
* Export estimates
*------------------------------------------------------------

estimates dir

* change to . or , depending on the system and country
cd "../output/tables"

esttab table_column_1_fever table_column_2_fever table_column_3_fever ///
    table_column_4_fever table_column_5_fever table_column_6_fever ///
	using "Berazneva_Byker_Table_Fever.rtf", replace ///
    keep(forest_loss_0 forest_loss_1 forest_loss_2 forest_loss_3 _cons) ///
    order(forest_loss_0 forest_loss_1 forest_loss_2 forest_loss_3 _cons) ///
    cells( ///
        b(star fmt(%9.3f)) se(par fmt(%9.3f)) p(par([ ]) fmt(%9.3f)) ///
    ) ///
	collabels(none) ///
    mlabels(none) ///
    stats(N r2, labels("N" "R2") fmt(%9.0fc %9.3f)) ///
	compress
	

esttab table_column_1_diarrhea table_column_2_diarrhea /// 
	table_column_3_diarrhea table_column_4_diarrhea ///
	table_column_5_diarrhea table_column_6_diarrhea ///
	using "Berazneva_Byker_Table_Diarrhea.rtf", replace ///
    keep(forest_loss_0 forest_loss_1 forest_loss_2 forest_loss_3) ///
    order(forest_loss_0 forest_loss_1 forest_loss_2 forest_loss_3) ///
    cells( ///
        b(star fmt(%9.3f)) se(par fmt(%9.3f)) p(par([ ]) fmt(%9.3f)) ///
    ) ///
    collabels(none) ///
    mlabels(none) ///
    stats(N r2, labels("N" "R2") fmt(%9.0fc %9.3f)) ///
    compress

esttab table_column_1_cough table_column_2_cough table_column_3_cough ///
       table_column_4_cough table_column_5_cough table_column_6_cough ///
	using "Berazneva_Byker_Table_Cough.rtf", replace ///
    keep(forest_loss_0 forest_loss_1 forest_loss_2 forest_loss_3) ///
    order(forest_loss_0 forest_loss_1 forest_loss_2 forest_loss_3) ///
    cells( ///
        b(star fmt(%9.3f)) se(par fmt(%9.3f)) p(par([ ]) fmt(%9.3f)) ///
    ) ///
    collabels(none) ///
    mlabels(none) ///
    stats(N r2, labels("N" "R2") fmt(%9.0fc %9.3f)) ///
    compress

*------------------------------------------------------------
* Export estimates (Optional Latex)
*------------------------------------------------------------

/*
esttab table_column_1_fever table_column_2_fever table_column_3_fever ///
    table_column_4_fever table_column_5_fever table_column_6_fever ///
	using "Berazneva_Byker_Table_Fever.tex", replace ///
    keep(forest_loss_0 forest_loss_1 forest_loss_2 forest_loss_3 _cons) ///
    order(forest_loss_0 forest_loss_1 forest_loss_2 forest_loss_3 _cons) ///
	coeflabels( ///
    forest_loss_0 "\; This year" ///
    forest_loss_1 "\; 1 year ago" ///
    forest_loss_2 "\; 2 years ago" ///
    forest_loss_3 "\; 3 years ago" ///
    _cons         "Constant" ///
	) ///
    cells( ///
        b(star fmt(%9.3f)) se(par fmt(%9.3f)) p(par([ ]) fmt(%9.3f)) ///
    ) ///
	collabels(none) ///
    mlabels(none) ///
    stats(N r2, labels("Observations" "R-squared") fmt(%9.0fc %9.3f)) ///
	compress
	

esttab table_column_1_diarrhea table_column_2_diarrhea /// 
	table_column_3_diarrhea table_column_4_diarrhea ///
	table_column_5_diarrhea table_column_6_diarrhea ///
	using "Berazneva_Byker_Table_Diarrhea.tex", replace ///
    keep(forest_loss_0 forest_loss_1 forest_loss_2 forest_loss_3 _cons) ///
    order(forest_loss_0 forest_loss_1 forest_loss_2 forest_loss_3 _cons) ///
	coeflabels( ///
    forest_loss_0 "\; This year" ///
    forest_loss_1 "\; 1 year ago" ///
    forest_loss_2 "\; 2 years ago" ///
    forest_loss_3 "\; 3 years ago" ///
    _cons         "Constant" ///
	) ///
    cells( ///
        b(star fmt(%9.3f)) se(par fmt(%9.3f)) p(par([ ]) fmt(%9.3f)) ///
    ) ///
    collabels(none) ///
    mlabels(none) ///
    stats(N r2, labels("Observations" "R-squared") fmt(%9.0fc %9.3f)) ///
    compress

esttab table_column_1_cough table_column_2_cough table_column_3_cough ///
       table_column_4_cough table_column_5_cough table_column_6_cough ///
	using "Berazneva_Byker_Table_Cough.tex", replace ///
    keep(forest_loss_0 forest_loss_1 forest_loss_2 forest_loss_3 _cons) ///
    order(forest_loss_0 forest_loss_1 forest_loss_2 forest_loss_3 _cons) ///
	coeflabels( ///
    forest_loss_0 "\; This year" ///
    forest_loss_1 "\; 1 year ago" ///
    forest_loss_2 "\; 2 years ago" ///
    forest_loss_3 "\; 3 years ago" ///
    _cons         "Constant" ///
	) ///
    cells( ///
        b(star fmt(%9.3f)) se(par fmt(%9.3f)) p(par([ ]) fmt(%9.3f)) ///
    ) ///
    collabels(none) ///
    mlabels(none) ///
    stats(N r2, labels("Observations" "R-squared") fmt(%9.0fc %9.3f)) ///
    compress
*/

/********************************************************************
Figure 2: Cumulative impact of forest loss 
Extended Dataset Full Specification 3

********************************************************************/

*------------------------------------------------------------
* Build cumulative effects dataset 
*------------------------------------------------------------

tempfile cumdata
tempname handle

postfile `handle' str10 outcome byte year double beta ll95 ul95 using `cumdata', replace

foreach y in fever diarrhea cough {

    estimates restore table_column_3_`y'

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

graph export "../figures/Berazneva_Byker_Figure_2_Updated.png", replace width(3000)

* Remove intermediate graph files before regenerating figures
foreach f in gA_fever gB_diarrhea gC_cough {
    capture erase "`f'.gph"
}

/********************************************************************
End of the File 
********************************************************************/
