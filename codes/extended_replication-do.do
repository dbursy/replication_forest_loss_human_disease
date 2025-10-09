**********************************************
* Replication do file (2)
* Bursy and Fiestas-Flores                   *
* A comment onf Does Forest Loss Increase 	 *
* Human Disease?   							 *
* Evidence from Nigeria                      *              
*                                            *
* September 2025                             *
*                                            *
**********************************************
clear all
* set maxvar 32000 
* set matsize 10000

**********************************
*** Roburstness ***
**********************************

*** Robustness 2: Alt. Denormalization with original cleaned raw data with 2018 wave ***

* Using data from the original source. Raw data was not provided in the replication file
 use "C:\LocalData\fiestas-flores\Seafile\Personal\Replication Games\Replication\Data\dta_replication_full_extended.dta", clear	

 ** Identify LGAs containing at least one cluster per DHS wave(DHS 2008, 2013 and 2018)
 
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

	


** Setup  GIS variables (Lagged forest loss and luminosity, initial forest cover and soil characteristics)	
	
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

	
order forest_loss_0 forest_loss_1 forest_loss_2 forest_loss_3 forest_loss_4 forest_loss_5 forest_loss_6 forest_loss_7  ///
forest_loss_8 forest_loss_9 forest_loss_10 forest_loss_11 forest_loss_12 ///  // Added
forest_loss_13 forest_loss_14 forest_loss_15 forest_loss_16 forest_loss_17 /// 
luminosity_chg_0  luminosity_chg_1 luminosity_chg_2

tabstat treecover_2008MEAN	treecover_2013MEAN treecover_2018MEAN, stats(N mean sd) 	// Added for control

	
* sort DHSYEAR
* by DHSYEAR: sum v005 if age_resp>14 & age_resp<50
	
* CORRECTED DENORMALIZATION
gen DN_wgt_c = . // De normalized weight
replace DN_wgt_c=(v005/1000000)*(35882027/28647) if DHSYEAR==2008 
replace DN_wgt_c=(v005/1000000)*(41018918/31225) if DHSYEAR==2013 
replace DN_wgt_c=(v005/1000000)*(47146386/33775) if DHSYEAR==2018 
 **** used this for female at reprodutive age: https://population.un.org/dataportal/ -- FRA 2008 - 2018

gen var = DN_wgt - DN_wgt_c // The corrected weight is larger, that is, larger population estimates (bigger standard errors, possibly)

egen stratid_c=group(DHSYEAR v022)
egen wave_psu_c=group(DHSYEAR v021)
svyset wave_psu [pweight=DN_wgt_c], strata(stratid_c) singleunit(certainty) 


** Regressions: 
 
local outcome "fever diarrhea cough"

 foreach y of local outcome {   

svy: reg `y' forest_loss_0-forest_loss_3 if  LGA_overlap_2==7  
estimates store m1_`y' // added

svy: reg `y' forest_loss_0-forest_loss_3 i.DHSYEAR i.LGA i.month#i.region#i.DHSYEAR i.month treecover_2000 no_HH_members no_kids_under_5 time_to_water head_HH_age HH_head_edu_years rural toilet poorest firewood floor age  age_resp edu_years no_child_total no_child_living livewith christian muslim yoruba igbo hausa resp_slept_net married resp_works pregnant if LGA_overlap_2==7
estimates  store m2_`y' // added

svy: reg `y' forest_loss_0-forest_loss_3 i.DHSYEAR i.LGA i.month#i.region#i.DHSYEAR i.month treecover_2000 cec_ave_pt ph_ave_pt occ_ave_pt altitude luminosity_chg_0-luminosity_chg_2 no_HH_members no_kids_under_5 time_to_water head_HH_age HH_head_edu_years rural toilet poorest firewood floor age  age_resp edu_years no_child_total no_child_living livewith christian muslim yoruba igbo hausa resp_slept_net married resp_works pregnant if LGA_overlap_2==7 
estimates  store m3_`y' // added


}

estimates dir
estout m1_fever m2_fever m3_fever , cells(b(fmt(%7,3f)) se(abs fmt(%7,3f)) t(abs fmt(%7,3f))) stats(N, fmt(%17,3f)) keep(forest_* _cons)

* change to . or , depending on the system and country
cd "C:\LocalData\fiestas-flores\Seafile\Personal\Replication Games\Replication\Tables"
estout m1_fever m2_fever m3_fever m1_diarrhea m2_diarrhea m3_diarrhea m1_cough m2_cough m3_cough using "raw_08_18_altdn.xls", ///
cells(b(fmt(%7,3f)) se(abs fmt(%7,3f)) t(abs fmt(%7,3f))) stats(N r2, fmt(%17,3f)) keep(forest_* _cons) replace


*** Robustness 5: Avoid potential overspecification problems ***

** Regressions: 
 
pwcorr forest_loss_0-forest_loss_3 cec_ave_pt ph_ave_pt occ_ave_pt  ///
luminosity_chg_0-luminosity_chg_2 altitude , star(0.05) 

pwcorr rural poorest HH_head_edu_years head_HH_age resp_slept_net resp_works, star(0.05)
pwcorr rural poorest luminosity_chg_0-luminosity_chg_2 age age_resp edu_years no_child_total no_child_living resp_slept_net resp_works, star(0.05) 
 
local outcome "fever diarrhea cough"

 foreach y of local outcome {   

* Dropping "extra" geographical controls
svy: reg `y' forest_loss_0-forest_loss_3 i.DHSYEAR i.LGA i.month#i.region treecover_2000 altitude no_HH_members no_kids_under_5 time_to_water head_HH_age HH_head_edu_years rural toilet poorest firewood floor age  age_resp edu_years no_child_total no_child_living livewith christian muslim yoruba igbo hausa resp_slept_net married resp_works pregnant if LGA_overlap_2==7 
estimates  store m3_1_`y' // added // reducing the interaction (drop i.month#i.region#i.DHSYEAR) Drop soil variables as colinear with def (cec_ave_pt ph_ave_pt occ_ave_pt soil?) drop limunisty as is already colinear with deforestation ( luminosity_chg_0-luminosity_chg_2)

* Dropping "extra" sociodemographic controls
svy: reg `y' forest_loss_0-forest_loss_3 i.DHSYEAR i.LGA i.month#i.region  treecover_2000 cec_ave_pt ph_ave_pt occ_ave_pt altitude luminosity_chg_0-luminosity_chg_2 no_HH_members no_kids_under_5 time_to_water HH_head_edu_years rural poorest  age livewith christian muslim yoruba igbo hausa  married  pregnant if LGA_overlap_2==7 
estimates  store m3_2_`y' // added // focus only on the kid age (dropping head_HH_age, age_resp) focusing only on HH education (dropping edu_years resp_works resp_slept_net) focus on children living there (dropping no_child_total no_child_living) ///
* avoid correlation with other SES (dropping toilet firewood floor)

* Dropping "extra" geographical and sociodemographic controls
svy: reg `y' forest_loss_0-forest_loss_3 i.DHSYEAR i.LGA i.month#i.region  treecover_2000 altitude no_HH_members no_kids_under_5 time_to_water HH_head_edu_years rural poorest  age livewith christian muslim yoruba igbo hausa  married  pregnant if LGA_overlap_2==7 
estimates  store m3_3_`y' // added // Drop all previous ones

}

estimates dir

* change to . or , depending on the system and country
cd "C:\LocalData\fiestas-flores\Seafile\Personal\Replication Games\Replication\Tables"
estout m3_1_fever m3_2_fever m3_3_fever m3_1_diarrhea m3_2_diarrhea m3_3_diarrhea m3_1_cough m3_2_cough m3_3_cough using "raw_08_18_over.xls", ///
cells(b(fmt(%7,3f)) se(abs fmt(%7,3f)) t(abs fmt(%7,3f))) stats(N r2, fmt(%17,3f)) keep(forest_* _cons) replace

*** Conclusion: It does not replicate with the original raw data using the correct de normalization with a more parsimonious model






