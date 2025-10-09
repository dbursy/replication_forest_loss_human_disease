**********************************************
* Replication do file
* Bursy and Fiestas-Flores                   *
* A comment onf Does Forest Loss Increase 	 *
* Human Disease?   							 *
* Evidence from Nigeria                      *              
*                                            *
* September 2025                             *
*                                            *
**********************************************
clear all
* set maxvar 32000 /
* set matsize 10000

* Using data from the original source. Raw data was not provided in the replication file
 use "C:\LocalData\fiestas-flores\Seafile\Personal\Replication Games\Replication\Data\dta_replication_full_extended.dta", clear	

 
** Identify LGAs containing at least one cluster per DHS wave(DHS 2008 and 2013)
 
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
	

** Setup  GIS variables (Lagged forest loss and luminosity, initial forest cover and soil characteristics)	
	
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
		
tabstat treecover_2008MEAN	treecover_2013MEAN, stats(N mean sd) 	// Added for control
* When pooling multiple waves of DHS data it is recommended to de-normalize weights using populaiton estimates
* We follow this procedure based on the following guidelines: http://userforum.dhsprogram.com/index.php?t=getfile&id=4
* Population estimates are based on UN Population estimates for women aged 15-49 by year in Nigeria: https://esa.un.org/unpd/wpp/Download/Standard/ASCII/

* sort DHSYEAR
* by DHSYEAR: sum v005 if age_resp>14 & age_resp<50

* ORIGINAL DENORMALIZATION
cap drop DN_wgt stratid wave_psu
gen DN_wgt = .
replace DN_wgt=(v005/1000000)*(34442230/33385) if DHSYEAR==2008 // 35882027 in 2008 UN 2025 / v005 for 2008 is 28647
replace DN_wgt=(v005/1000000)*(39172540/38185) if DHSYEAR==2013 // 41018918 in 2013 UN 2025 / v005 for 2013 is 31225

egen stratid=group(DHSYEAR v022)
egen wave_psu=group(DHSYEAR v021)
svyset wave_psu [pweight=DN_wgt], strata(stratid) singleunit(certainty) 


** Regressions to replica Table 1: 
 
local outcome "fever diarrhea cough"

 foreach y of local outcome {   

svy: reg `y' forest_loss_0-forest_loss_3 if  LGA_overlap==3  
estimates store m1_`y' // added

svy: reg `y' forest_loss_0-forest_loss_3 i.DHSYEAR i.LGA i.month#i.region#i.DHSYEAR i.month treecover_2000 no_HH_members no_kids_under_5 time_to_water head_HH_age HH_head_edu_years rural toilet poorest firewood floor age  age_resp edu_years no_child_total no_child_living livewith christian muslim yoruba igbo hausa resp_slept_net married resp_works pregnant if LGA_overlap==3 
estimates  store m2_`y' // added

svy: reg `y' forest_loss_0-forest_loss_3 i.DHSYEAR i.LGA i.month#i.region#i.DHSYEAR i.month treecover_2000 cec_ave_pt ph_ave_pt occ_ave_pt altitude luminosity_chg_0-luminosity_chg_2 no_HH_members no_kids_under_5 time_to_water head_HH_age HH_head_edu_years rural toilet poorest firewood floor age  age_resp edu_years no_child_total no_child_living livewith christian muslim yoruba igbo hausa resp_slept_net married resp_works pregnant if LGA_overlap==3 
estimates  store m3_`y' // added


}

estimates dir // added

**********************************
*** Computational Replication ****
**********************************

* change to . or , depending on the system and country
cd "C:\LocalData\fiestas-flores\Seafile\Personal\Replication Games\Replication\Tables"
estout m1_fever m2_fever m3_fever m1_diarrhea m2_diarrhea m3_diarrhea m1_cough m2_cough m3_cough using "raw.xls", ///
cells(b(fmt(%7,3f)) se(abs fmt(%7,3f)) t(abs fmt(%7,3f))) stats(N r2, fmt(%17,3f)) keep(forest_* _cons) replace


*** Conclusion: It does not replicate with the original raw data

**********************************
*** Roburstness ***
**********************************

*** Robustness 1: Alt. Denormalization with original cleaned raw data ***
gen DN_wgt_c = . // De normalized weight
replace DN_wgt_c=(v005/1000000)*(35882027/28647) if DHSYEAR==2008 // 35882027 in 2008 UN 2025 / v005 for 2008 is 28647
replace DN_wgt_c=(v005/1000000)*(41018918/31225) if DHSYEAR==2013 // 41018918 in 2013 UN 2025 / v005 for 2013 is 31225

gen var = DN_wgt - DN_wgt_c // The corrected weight is larger, that is, larger population estimates (bigger standard errors, possibly)

egen stratid_c=group(DHSYEAR v022)
egen wave_psu_c=group(DHSYEAR v021)
svyset wave_psu [pweight=DN_wgt_c], strata(stratid_c) singleunit(certainty) 

		
local outcome "fever diarrhea cough"

 foreach y of local outcome {   

svy: reg `y' forest_loss_0-forest_loss_3 if  LGA_overlap==3  
estimates store m1_`y' // added

svy: reg `y' forest_loss_0-forest_loss_3 i.DHSYEAR i.LGA i.month#i.region#i.DHSYEAR i.month treecover_2000 no_HH_members no_kids_under_5 time_to_water head_HH_age HH_head_edu_years rural toilet poorest firewood floor age  age_resp edu_years no_child_total no_child_living livewith christian muslim yoruba igbo hausa resp_slept_net married resp_works pregnant if LGA_overlap==3 
estimates  store m2_`y' // added

svy: reg `y' forest_loss_0-forest_loss_3 i.DHSYEAR i.LGA i.month#i.region#i.DHSYEAR i.month treecover_2000 cec_ave_pt ph_ave_pt occ_ave_pt altitude luminosity_chg_0-luminosity_chg_2 no_HH_members no_kids_under_5 time_to_water head_HH_age HH_head_edu_years rural toilet poorest firewood floor age  age_resp edu_years no_child_total no_child_living livewith christian muslim yoruba igbo hausa resp_slept_net married resp_works pregnant if LGA_overlap==3 
estimates  store m3_`y' // added


}

estimates dir
	
* change to . or , depending on the system and country
cd "C:\LocalData\fiestas-flores\Seafile\Personal\Replication Games\Replication\Tables"
estout m1_fever m2_fever m3_fever m1_diarrhea m2_diarrhea m3_diarrhea m1_cough m2_cough m3_cough using "raw_altdn.xls", ///
cells(b(fmt(%7,3f)) se(abs fmt(%7,3f)) t(abs fmt(%7,3f))) stats(N r2, fmt(%17,3f)) keep(forest_* _cons) replace

*** Conclusion: It does not replicate with the original raw data using the correct de normalization


*** Robustness 4: Avoid potential overspecification problems ***

local outcome "fever diarrhea cough"

 foreach y of local outcome {   

* Dropping "extra" geographical controls
svy: reg `y' forest_loss_0-forest_loss_3 i.DHSYEAR i.LGA i.month#i.region treecover_2000 altitude no_HH_members no_kids_under_5 time_to_water head_HH_age HH_head_edu_years rural toilet poorest firewood floor age  age_resp edu_years no_child_total no_child_living livewith christian muslim yoruba igbo hausa resp_slept_net married resp_works pregnant if LGA_overlap==3 
estimates  store m3_1_`y' // added // reducing the interaction (drop i.month#i.region#i.DHSYEAR) Drop soil variables as colinear with def (cec_ave_pt ph_ave_pt occ_ave_pt soil?) drop limunisty as is already colinear with deforestation ( luminosity_chg_0-luminosity_chg_2)

* Dropping "extra" sociodemographic controls
svy: reg `y' forest_loss_0-forest_loss_3 i.DHSYEAR i.LGA i.month#i.region  treecover_2000 cec_ave_pt ph_ave_pt occ_ave_pt altitude luminosity_chg_0-luminosity_chg_2 no_HH_members no_kids_under_5 time_to_water HH_head_edu_years rural poorest  age livewith christian muslim yoruba igbo hausa  married  pregnant if LGA_overlap==3
estimates  store m3_2_`y' // added // focus only on the kid age (dropping head_HH_age, age_resp) focusing only on HH education (dropping edu_years resp_works resp_slept_net) focus on children living there (dropping no_child_total no_child_living) ///
* avoid correlation with other SES (dropping toilet firewood floor)

* Dropping "extra" geographical and sociodemographic controls
svy: reg `y' forest_loss_0-forest_loss_3 i.DHSYEAR i.LGA i.month#i.region  treecover_2000 altitude no_HH_members no_kids_under_5 time_to_water HH_head_edu_years rural poorest  age livewith christian muslim yoruba igbo hausa  married  pregnant if LGA_overlap==3 
estimates  store m3_3_`y' // added // Drop all previous ones

}

estimates dir

* change to . or , depending on the system and country
cd "C:\LocalData\fiestas-flores\Seafile\Personal\Replication Games\Replication\Tables"
estout m3_1_fever m3_2_fever m3_3_fever m3_1_diarrhea m3_2_diarrhea m3_3_diarrhea m3_1_cough m3_2_cough m3_3_cough using "raw_over.xls", ///
cells(b(fmt(%7,3f)) se(abs fmt(%7,3f)) t(abs fmt(%7,3f))) stats(N r2, fmt(%17,3f)) keep(forest_* _cons) replace



*** Conclusion: It does not replicate with the original raw data using the correct de normalization with a more parsimonious model




