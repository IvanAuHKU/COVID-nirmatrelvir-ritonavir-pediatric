cd "D:\Data\COVID"
capture mkdir "Pediatric Paxlovid"
cd "Pediatric Paxlovid"

**# ***** Notes
*** Characteristics variables
* pseudo_key: anonymized pseudo-id of patients
* date_index: index date, defined as date of SARS-CoV-2 infection diagnosis or symptom onset, whichever occurred earlier
* date_paxlovid: date of nirmatrelvir/ritonavir initiation
* group: 1 for nirmatrelvir/ritonavir users (initiated within 5 day from index date); and 0 refers to non-users (beyond 5 days from index date)
* age: Age of patients (years)
* sex: Sex of patients (1 for male; 0 for female)
* symptomatic: Symptomatic presentation
* asthma: Asthma
* cancer: Cancer
* cvd: Cardiac disease
* lung: Lung disease
* mental: Mental disease
* neuro_cate: Neurologic disease
* obesity: Obesity
* dm: Diabetes mellitus
* disability: Disabilities
* ADHD: ADHD
* autism: Autism
* immuno_hist: Immunocompromised
* healthcare: Healthcare utilization
* vaccine_status: COVID-19 vaccination status of patients (1: Not fully vaccinated, 2: Fully vaccinated but not boosted, 3: Boosted)

*** Study outcomes
* admission: Hospitalization
* death: All-cause mortality
* composite_ip: In-hospital disease progression
* death_ip: In-hospital death
* intubation: Invasive mechanical ventilation
* ICU: Intensive care unit admission
* outcome_encephalitis: Encephalitis
* outcome_encephalopathy: Encephalopathy
* outcome_seizure: Seizure
* outcome_stroke: Stroke
* outcome_gbs: Guillain-Barre syndrome
* outcome_cardio: Cardiovascular complications
* outcome_carditis: Myocarditis/pericarditis
* outcome_arrhythmias: Arrhythmias
* outcome_cardiogenic_shock: Cardiogenic shock
* outcome_thrombosis: Thrombosis
* outcome_thromboembolism: Thromboembolism
* outcome_thromboembolic: Thromboembolic event
* outcome_thrombocytopenia: (Idiopathic) Thrombocytopenia
* outcome_MISC: Multisystem inflammatory syndrome in children (MIS-C)
* ali: Acute liver injury
* aki: Acute renal failure
* ARDS: Acute respiratory distress syndrome

*** Datasets
* "covid_patient_list_extracted.dta": list of patients eligible analysis, together with the characteristics and dates of outcome
* "dx_extracted.dta": list of all diagnosis codes

cls
**# Table 1 - Baseline characteristics of original cohort
* N/mean, %/SD
qui forvalues k = 1/1 {
use "covid paxlovid pediatric", clear
noi di _newline "overall" _newline "var" _col(30) "N" _col(45) "mean" _col(60) "sd"
sum age
noi di "age" _col(30) r(N) _col(45) r(mean) _col(60) r(sd)
noi di "var" _col(30) "N" _col(45) "%"
qui foreach var in symptomatic asthma cancer cvd lung mental neuro_cate obesity dm disability ADHD autism immuno_hist healthcare {
	count
	scalar N_all = r(N)
	count if `var' == 1
	scalar N = r(N)
	noi di "`var'=1" _col(30) N _col(45) N/N_all*100
}
foreach var in age_gp gender vaccine_status ym_covid_cate {
	if inlist("`var'", "age_gp", "vaccine_status", "ym_covid_cate") {
	forvalues j = 1/3 {
		count
		scalar N_all = r(N)
		count if `var' == `j'
		scalar N = r(N)
		noi di "`var'=`j'" _col(30) N _col(45) N/N_all*100
	}
	}
	if inlist("`var'", "gender") {
	forvalues j = 1(-1)0 {
		count
		scalar N_all = r(N)
		count if `var' == `j'
		scalar N = r(N)
		noi di "`var'=`j'" _col(30) N _col(45) N/N_all*100
	}
	}
}
}
*
********************************************************************************
cls
**# Table 2 - Main analysis - Target trial emulation
capture mkdir "main"
capture mkdir "main/cloned"
capture mkdir "main/split"
capture mkdir "main/Treatment"
capture mkdir "main/NoTreatment"
capture mkdir "main/emulated"
capture mkdir "main/KM ipcw_stab"

* Target trial emulation
qui forvalues k = 1/1 {
forvalues j = 0/500 {
*** Prepare dataset for analysis
	noi di "main" _col(15) "bs`j'" _col(30) "setup" 
	* Setup for trial emulation
	if !fileexists("main/split/split bs`j'.dta") {
	use "covid paxlovid pediatric.dta", clear
	if `j' > 0 {
		set seed `j'
		bsample
	}
	replace date_paxlovid = . if date_paxlovid >= date_admission
	* organize / rename / generate variables
	gen date_last_fu = min(date_death, mdy(02,12,2023), date_baseline + 28)
	gen date_event = date_admission
	gen event = inrange(date_event, date_baseline, date_last_fu)
	gen fup_obs = min(date_event-date_baseline, 28) if event == 1
	replace fup_obs = min(date_last_fu-date_baseline, 28) if event == 0
	tab fup_obs
	gen time_to_treatment = date_paxlovid - date_baseline
	tab time_to_treatment
	gen treatment = inrange(date_paxlovid - date_baseline, 0, 5)
	* keep necessary variables
	keep pseudo_key fup_obs event time_to_treatment treatment group
	gen bs = `j'	
	
*** Cloning & censoring
	noi di "main" _col(15) "bs`j'" _col(30) "cloning" 
	* Prepare dataset for analysis
	stset fup_obs, failure(event)

	* Arm A: no treatment within 5 days (control: non-exposed group)
	gen outcomeA = _d // _d = `event'
	gen fupA = _t // _t = follow up time

	/// if the patient received treatment within 6 months:
	/// 1. no event outcome, since the patient survived till censoring (treatment)
	replace outcomeA = 0 if treatment==1 & time_to_treatment <=5 
	/// 2. follow up is censored at treatment
	replace fupA = time_to_treatment if treatment==1 & time_to_treatment <=5

	* Arm B: treatment within 5 days (treated: exposed group)
	gen outcomeB = _d 
	gen fupB = _t 

	/// if the patient survived the first 5 days and did not receive treatment within 5 days:
	/// 1. no event outcome if the patient survived the first 6 months
	replace outcomeB = 0 if (treatment==0 & _t>5) | (treatment==1 & time_to_treatment >5 & time_to_treatment !=.)
	/// 2. follow up is censored at 6 months
	replace fupB = 5 if (treatment==0 & _t>5) | (treatment==1 & time_to_treatment >5 & time_to_treatment != .)

	** append clones 
	preserve
		drop outcomeB fupB
		rename outcomeA outcome
		rename fupA fup
		gen arm = "NoTreatment"
		tempfile a
		save "`a'", replace
	restore
		drop outcomeA fupA
		rename outcomeB outcome
		rename fupB fup
		gen arm = "Treatment"	
		cap append using "`a'"

	// Weight models

	sort _all
	gen NewID = _n

	** add 1 day to 0-survivors
	replace fup= 1 if fup==0

	** Weight model: define survival time and event indicator	
	* treatment Arm
	** Case 1: they do not deviate at time of treatment, but are not at risk of deviating any more
	gen wm_fup = time_to_treatment if arm == "Treatment" & time_to_treatment<=5 & time_to_treatment!=. & treatment == 1 
	gen wm_outcome = 0 if arm == "Treatment" & time_to_treatment<=5 & time_to_treatment!=. & treatment == 1 

	** Case 2: they deviate at 5 days
	replace wm_fup = 5 if arm == "Treatment" & ((treatment == 0 & fup >= 5) | (time_to_treatment>5 & treatment == 1))
	replace wm_outcome = 1 if arm == "Treatment" & ((treatment == 0 & fup >= 5) | (time_to_treatment>5 & treatment == 1)) 

	** Case 3: they do not deviate, but we need to keep their survival as observed and censor them as we do not know what happens afterwards
	replace wm_fup = fup if arm == "Treatment" & treatment == 0 & fup < 5
	replace wm_outcome = 0 if arm == "Treatment" & treatment == 0 & fup < 5
	** add 1 days to 0-survivors
	replace wm_fup= 1 if arm == "Treatment" & wm_fup==0 

	* No treatment Arm
	** Case 1: they do not deviate at time of treatment, but are not at risk of deviating any more
	replace wm_fup = time_to_treatment if arm == "NoTreatment" & time_to_treatment<=5 & treatment == 1 
	replace wm_outcome = 1 if arm == "NoTreatment" & time_to_treatment<=5 & treatment == 1 

	** Case 2: they deviate at 5 days
	replace wm_fup = 5 if arm == "NoTreatment" & ((treatment == 0 & fup >= 5) | (time_to_treatment>5 & treatment == 1)) 
	replace wm_outcome = 0 if arm == "NoTreatment" & ((treatment == 0 & fup >= 5) | (time_to_treatment>5 & treatment == 1)) 

	** Case 3: they do not deviate, but we need to keep their survival as observed and censor them as we do not know what happens afterwards
	replace wm_fup = fup if arm == "NoTreatment" & treatment == 0 & fup < 5 
	replace wm_outcome = 0 if arm == "NoTreatment" & treatment == 0 & fup < 5

	** add 1 days to 0-survivors
	replace wm_fup= 1 if arm == "NoTreatment" & wm_fup == 0

	gen arm_value = 1 if arm == "Treatment"
	replace arm_value = 0 if arm == "NoTreatment"
	
	drop arm
	order pseudo_key treatment arm_value
	
	if `j' == 0 {
		save "main/cloned/cloned bs0.dta", replace
	}
	
*** Split times
	noi di "main" _col(15) "bs`j'" _col(30) "split"

	** times of event
	stset fup, fail(outcome) id(NewID)
	stsplit, at(failures)

	gen tstart = _t0
	gen tstop = _t

	** times of censoring
	gen TrialEmul_cens = 1-outcome
	stset fup, fail(TrialEmul_cens) id(NewID)
	stsplit, at(failures)

	replace tstart = _t0 if tstart<_t0 & _t0 != . & _t != .
	replace tstop = _t if tstop>_t & _t0 != . & _t != .

	order pseudo_key arm_value tstart tstop
	sort _all
	compress
	save "main/split/split bs`j'.dta", replace

	}
	
*** IPCW - Treatment arm
	noi di "main" _col(15) "bs`j'" _col(30) "Treatment"
	if !fileexists("main/Treatment/Treatment bs`j'.dta") {
	use "main/split/split bs`j'.dta", clear
	keep if arm_value == 1
	replace outcome = . if TrialEmul_cens == .
	merge m:1 pseudo_key using "covid pediatric paxlovid characteristics.dta", keep(3) nogen

	* adapt records to the long format
	sort NewID tstart
	bysort NewID (tstart): replace wm_outcome = 0 if _n!=_N

	* Weight model:
	stset tstop, origin(time tstart) failure(wm_outcome) id(NewID)
	logit wm_outcome age gender vaccine_status_2 vaccine_status_3 ym_covid_cate_2 ym_covid_cate_3 symptomatic lung immuno_hist healthcare ///
	i.symptomatic##i.ym_covid_cate_3 i.symptomatic##i.gender i.symptomatic##i.healthcare if tstart == 4
	predict pr_cens if e(sample), pr
	replace pr_cens = 0 if pr_cens == .
	gen pr_uncens = 1-pr_cens
	bysort arm_value pseudo_key (tstart): replace pr_uncens=pr_uncens*pr_uncens[_n-1] if _n!=1
	gen weight = 1/pr_uncens
	replace weight = 1 if weight == .
	gen invalid = 0
	keep pseudo_key arm_value tstart tstop event fup_obs time_to_treatment treatment bs outcome fup NewID wm_fup wm_outcome TrialEmul_cens weight invalid
	compress
	save "main/Treatment/Treatment bs`j'.dta", replace
	}
	
*** IPCW - Control arm
	noi di "main" _col(15) "bs`j'" _col(30) "NoTreatment"
	if !fileexists("main/NoTreatment/Notreatment bs`j'.dta") {
	use "main/split/split bs`j'.dta", clear
	keep if arm_value == 0
	replace outcome = . if TrialEmul_cens == .
	merge m:1 pseudo_key using "covid pediatric paxlovid characteristics.dta", keep(3) nogen

	* adapt records to the long format
	sort NewID tstart
	bysort NewID: replace wm_outcome = 0 if _n!=_N

	* Weight model:
	stset tstop, origin(time tstart) failure(wm_outcome) id(NewID)
	xi: logit wm_outcome tstart age gender vaccine_status_2 vaccine_status_3 ym_covid_cate_2 ym_covid_cate_3 symptomatic immuno_hist healthcare, cluster(NewID)
	predict pr_cens if e(sample), pr
	replace pr_cens = 0 if pr_cens == .
	gen pr_uncens = 1-pr_cens
	bysort arm_value pseudo_key (tstart): replace pr_uncens=pr_uncens*pr_uncens[_n-1] if _n!=1
	gen weight = 1/pr_uncens
	replace weight = 1 if weight == .
	gen invalid = 0
	keep pseudo_key arm_value tstart tstop event fup_obs time_to_treatment treatment bs outcome fup NewID wm_fup wm_outcome TrialEmul_cens weight invalid
	compress
	save "main/NoTreatment/Notreatment bs`j'.dta", replace
	}
	
*** Combine & Generate weights
	noi di "main" _col(15) "bs`j'" _col(30) "IPCW"
	if !fileexists("main/emulated/emulated bs`j'.dta") {
	use "main/Treatment/Treatment bs`j'.dta", clear
	append using "main/NoTreatment/Notreatment bs`j'.dta"

	// create a new ID variable for each clone in each arm
	tostring NewID, gen(NewID_str)
	replace NewID_str = "0000000" + NewID_str if length(NewID_str)== 1
	replace NewID_str = "000000" + NewID_str if length(NewID_str)== 2
	replace NewID_str = "00000" + NewID_str if length(NewID_str)== 3
	replace NewID_str = "0000" + NewID_str if length(NewID_str)== 4
	replace NewID_str = "000" + NewID_str if length(NewID_str)== 5
	replace NewID_str = "00" + NewID_str if length(NewID_str)== 6
	replace NewID_str = "0" + NewID_str if length(NewID_str)== 7
	
	gen Anal_ID = "1" + NewID_str if arm_value == 1
	replace Anal_ID = "2" + NewID_str if arm_value == 0

	replace weight = 1 if wm_outcome == 1 & tstop == 5

	keep pseudo_key tstart tstop fup_obs event time_to_treatment treatment bs outcome fup wm_fup wm_outcome TrialEmul_cens weight Anal_ID arm_value invalid
	destring Anal_ID, replace

	* generate denominator
	rename weight _ipcw
	gen prob_denom = 1/_ipcw
	* generate numerator
	bysort arm_value tstart : gen N_group = _N
	bysort arm_value tstart : egen N_censored = sum(wm_outcome)
	gen prob_uncensored = (N_group-N_censored)/N_group
	gen prob_num = prob_uncensored
	sort arm_value pseudo_key tstart
	by arm_value pseudo_key : replace prob_num=prob_num*prob_num[_n-1] if _n!=1
	* generate stabilized IPCW
	gen _ipcw_stab = prob_num / prob_denom
	compress
	save "main/emulated/emulated bs`j'.dta", replace
	
	if `j' > 1 {
		capture erase "main/split/split bs`j'.dta"
		capture erase "main/Treatment/Treatment bs`j'.dta"
		capture erase "main/NoTreatment/NoTreatment bs`j'.dta"
	}
	}
	
*** Generate KM estimate
	noi di "main" _col(15) "bs`j'" _col(30) "KM"
	if !fileexists("main/KM ipcw_stab/KM ipcw_stab bs`j'.dta") {
	use "main/emulated/emulated bs`j'.dta", clear
	count if invalid == 1
	if r(N) == 0 {
		stset tstop [pweight = _ipcw_stab], enter(time tstart) failure(outcome)
		sts generate KM_s_w = s if arm_value == 1
		sts generate KM_ns_w = s if arm_value == 0
	}
	else {
		gen KM_s_w = .
		gen KM_ns_w = .
	}
	collapse (firstnm) KM_s_w KM_ns_w, by(fup bs invalid)
	save "main/KM ipcw_stab/KM ipcw_stab bs`j'.dta", replace
	
	if `j' > 1 {
		capture erase "main/emulated/emulated bs`j'.dta"
	}
	}
}
}
*
* Finalize bootstrap datasets
clear
qui forvalues k = 1/1 {
forvalues j = 0/500 {
	capture append using "main/KM ipcw_stab/KM ipcw_stab bs`j'.dta"
	*erase "main/KM/KM bs`j'.dta"
}
	gen hazard_s_w = 1 - KM_s_w
	gen hazard_ns_w = 1 - KM_ns_w
	gen odds_s_w = hazard_s_w/(1-hazard_s_w)
	gen odds_ns_w = hazard_ns_w/(1-hazard_ns_w)
	gen RR_w = hazard_s_w/hazard_ns_w
	gen diff_w = hazard_s_w - hazard_ns_w
	gen OR_w = odds_s_w / odds_ns_w
	compress
	save "main/KM ipcw_stab bs_all.dta", replace
}
*
cls
* KM estimate
qui forvalues k = 1/1 {
	use "main/KM ipcw_stab bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	
	sum hazard_s_w if bs == 0, d
	scalar hazard_s_mean = r(mean)
	centile hazard_s_w if bs > 0, centile(2.5 97.5)
	scalar hazard_s_cil = r(c_1)
	scalar hazard_s_ciu = r(c_2)
	
	sum hazard_ns_w if bs == 0, d
	scalar hazard_ns_mean = r(mean)
	centile hazard_ns_w if bs > 0, centile(2.5 97.5)
	scalar hazard_ns_cil = r(c_1)
	scalar hazard_ns_ciu = r(c_2)
	
	noi di "main" _col(10) hazard_s_mean _col(25) hazard_s_cil _col(40) hazard_s_ciu _col(55) ///
	hazard_ns_mean _col(70) hazard_s_cil _col(85) hazard_s_ciu
}
*
* Absolute risk reduction
qui forvalues k = 1/1 {
	use "main/KM ipcw_stab bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	
	sum diff_w if bs == 0, d
	scalar bs_mean = r(mean)
	
	centile diff_w if bs > 0, centile(2.5 97.5)
	scalar bs_cil = r(c_1)
	scalar bs_ciu = r(c_2)
	noi di "main" _col(10) -bs_mean _col(25) -bs_ciu _col(40) -bs_cil
}
*
* Relative risk
qui forvalues k = 1/1 {
	use "main/KM ipcw_stab bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	
	sum RR_w if bs == 0, d
	scalar bs_mean = r(mean)
	centile RR_w if bs > 0, centile(2.5 97.5)
	scalar bs_cil = r(c_1)
	scalar bs_ciu = r(c_2)
	noi di "main" _col(10) bs_mean _col(25) bs_cil _col(40) bs_ciu
}
*
* N
qui forvalues k = 1/1 {
	use "main/cloned/cloned bs0.dta", clear
	keep if (arm_value == 1 & group == 1) | (arm_value == 0 & group == 0)
	count if group == 1
	scalar N_1 = r(N)
	count if group == 0
	scalar N_0 = r(N)
	count if group == 1 & outcome == 1
	scalar n_e_1 = r(N)
	count if group == 0 & outcome == 1
	scalar n_e_0 = r(N)
	noi di "main" _col(15) N_1 _col(30) N_0 _col(45) n_e_1 _col(60) n_e_0
}
*
qui forvalues k = 1/1 {
	if `k' == 1 {
		noi di "main" _col(15) "N_1" _col(30) "N_0" _col(45) "n_e_1" _col(60) "n_e_0"
	}
	use "main/cloned/cloned bs0.dta", clear
	count if arm_value == 1
	scalar N_1 = r(N)
	count if arm_value == 1 & outcome == 1
	scalar n_e_1 = r(N)
	count if arm_value == 0
	scalar N_0 = r(N)
	count if arm_value == 0 & outcome == 1
	scalar n_e_0 = r(N)
	
	noi di "main" _col(15) N_1 _col(30) N_0 _col(45) n_e_1 _col(60) n_e_0
}
*

cls
**# Table 1 - day5 characteristics
forvalues t = 5/5 {
use "main/emulated/emulated bs0.dta", clear
merge m:1 pseudo_key using "covid paxlovid pediatric.dta", keep(3) nogen
keep if tstart == `t'
compress
save "main/emulated bs0 day`t'", replace
}
*
cls
* After IPCW
* N/mean, %/SD
qui forvalues t = 5/5 {
	noi di _newline "day `t'"
qui forvalues k = 1(-1)0 {
use "main/emulated bs0 day`t'", clear
noi di _newline "arm_value=`k'" _newline "var" _col(30) "N" _col(45) "mean" _col(60) "sd"
sum age [w=_ipcw_stab] if arm_value == `k'
noi di "age" _col(30) r(N) _col(45) r(mean) _col(60) r(sd)
noi di "var" _col(30) "N" _col(45) "%"
gen any = 0
foreach var in asthma cancer cvd lung mental neuro_cate obesity dm disability ADHD autism immuno_hist {
	replace any = 1 if `var' == 1
}
gen any_sum = 0
foreach var in asthma cancer cvd lung mental neuro_cate obesity dm disability ADHD autism immuno_hist {
	replace any_sum = any_sum + 1 if `var' == 1
}
qui foreach var in symptomatic asthma cancer cvd lung mental neuro_cate obesity dm disability ADHD autism immuno_hist healthcare {
	sum `var' if arm_value == `k' [w=_ipcw_stab]
	noi di "`var'=1" _col(30) r(N)*r(mean) _col(45) r(mean)*100
}
foreach var in age_gp gender vaccine_status ym_covid_cate {
	if inlist("`var'", "age_gp", "vaccine_status", "ym_covid_cate") {
	forvalues j = 1/3 {
		gen `var'`j' = `var' == `j'
		sum `var'`j' if arm_value == `k' [w=_ipcw_stab]
		noi di "`var'=`j'" _col(30) r(N)*r(mean) _col(45) r(mean)*100
	}
	}
	if inlist("`var'", "gender") {
	forvalues j = 1(-1)0 {
		gen `var'`j' = `var' == `j'
		sum `var'`j' if arm_value == `k' [w=_ipcw_stab]
		noi di "`var'=`j'" _col(30) r(N)*r(mean) _col(45) r(mean)*100
	}
	}
}
}
*
noi di
* SMD
qui forvalues k = 1/1 {
use "main/emulated bs0 day`t'", clear

noi di "var" _col(20) "smd"
* continuous
sum age if arm_value == 1 [w=_ipcw_stab]
local mean_1 = r(mean)
local sd_1 = r(sd)
sum age if arm_value == 0 [w=_ipcw_stab]
local mean_0 = r(mean)
local sd_0 = r(sd)
stddiffi `mean_1' `sd_1' `mean_0' `sd_0'
scalar SMD = abs(r(std_diff))
noi di "age" _col(20) SMD 

foreach var in symptomatic asthma cancer cvd lung mental neuro_cate obesity dm disability ADHD autism immuno_hist healthcare age_gp gender vaccine_status ym_covid_cate {
	tab `var' arm_value
	if r(r) == 2 {
		tab `var' if arm_value == 1 [aw=_ipcw_stab], matcell(temp)
		local m11 = int(temp[1,1])
		local m21 = int(temp[2,1])
		tab `var' if arm_value == 0 [aw=_ipcw_stab], matcell(temp)
		local m12 = int(temp[1,1])
		local m22 = int(temp[2,1])
		capture stddiffi `m11' `m12' \ `m21' `m22'
	}
	if r(r) == 3 {
		tab `var' if arm_value == 1 [aw=_ipcw_stab], matcell(temp)
		local m11 = int(temp[1,1])
		local m21 = int(temp[2,1])
		local m31 = int(temp[3,1])
		tab `var' if arm_value == 0 [aw=_ipcw_stab], matcell(temp)
		local m12 = int(temp[1,1])
		local m22 = int(temp[2,1])
		local m32 = int(temp[3,1])
		capture stddiffi `m11' `m12' \ `m21' `m22' \ `m31' `m32'
	}
	scalar SMD = abs(r(std_diff))
	noi di "`var'" _col(20) SMD 
}
}
}
*

********************************************************************************
cls
**# Supp T3
*** Supplementary Table 3. Baseline characteristics of nirmatrelvir/ritonavir group and control group before IPCW
* N/mean, %/SD
qui forvalues k = 1(-1)0 {
use "covid pediatric paxlovid characteristics", clear
noi di _newline "group=`k'" _newline "var" _col(30) "N" _col(45) "mean" _col(60) "sd"
sum age if group == `k'
noi di "age" _col(30) r(N) _col(45) r(mean) _col(60) r(sd)
noi di "var" _col(30) "N" _col(45) "%"
qui foreach var in symptomatic  asthma cancer cvd lung mental neuro_cate obesity dm disability ADHD autism immuno_hist healthcare {
	count if group == `k'
	scalar N_all = r(N)
	count if `var' == 1 & group == `k'
	scalar N = r(N)
	noi di "`var'=1" _col(30) N _col(45) N/N_all*100
}
foreach var in age_gp gender vaccine_status ym_covid_cate {
	if inlist("`var'", "vaccine_status", "ym_covid_cate") {
	forvalues j = 1/3 {
		count if group == `k'
		scalar N_all = r(N)
		count if `var' == `j' & group == `k'
		scalar N = r(N)
		noi di "`var'=`j'" _col(30) N _col(45) N/N_all*100
	}
	}
	if inlist("`var'", "age_gp") {
	forvalues j = 1/4 {
		count if group == `k'
		scalar N_all = r(N)
		count if `var' == `j' & group == `k'
		scalar N = r(N)
		noi di "`var'=`j'" _col(30) N _col(45) N/N_all*100
	}
	}
	if inlist("`var'", "gender") {
	forvalues j = 1(-1)0 {
		count if group == `k'
		scalar N_all = r(N)
		count if `var' == `j' & group == `k'
		scalar N = r(N)
		noi di "`var'=`j'" _col(30) N _col(45) N/N_all*100
	}
	}
}
}
*
* SMD
qui forvalues k = 1/1 {
use "covid pediatric paxlovid characteristics", clear
noi di "var" _col(20) "smd"
stddiff age, by(group) abs cohensd
scalar SMD = abs(r(stddiff)[1,1])
noi di "age" _col(20) SMD
gen any = 0
foreach var in asthma cancer cvd lung mental neuro_cate obesity dm disability ADHD autism immuno_hist {
	replace any = 1 if `var' == 1
}
foreach var in symptomatic  asthma cancer cvd lung mental neuro_cate obesity dm disability ADHD autism immuno_hist healthcare age_gp gender vaccine_status ym_covid_cate {
	tab `var' group
	if r(r) > 1 { 
		stddiff i.`var', by(group) abs cohensd
		scalar SMD = abs(r(stddiff)[1,1])
		noisily di "`var'" _col(20) SMD
	}
	else {
		noisily di "`var'" _col(20) .
	}
}
}
*

**# Supp T5
*** Supplementary Table 5. Principal diagnosis of hospitalization among nirmatrelvir/ritonavir users and controls
use "main/cloned/cloned subgp_1 bs0.dta", clear
keep if outcome == 1
gen arm_value_1 = arm_value == 1
gen arm_value_0 = arm_value == 0
collapse (max) arm_value_*, by(pseudo_key)
save "arm_value outcome", replace

use "dx_extracted", clear
merge m:1 pseudo_key using "date_admission_list", keep(3) nogen
merge m:1 pseudo_key using "arm_value outcome", nogen
*keep if group < .
keep if reference_date_id >= date_admission
* if covid was recorded in as primary or secondary, take covid
keep if ps == "P" | icd9_cd == "519.8"
gen covid = inlist(icd9_cd, "519.8")
bysort pseudo_key (reference_date_id covid icd9_cd) : keep if _n == _N

gen icd9_cd_cate_0 = inlist(icd9_cd, "519.8")
icd9 gen icd9_cd_cate_1 = icd9_cd, r(001/139 139*)
icd9 gen icd9_cd_cate_2 = icd9_cd, r(140/239 239*)
icd9 gen icd9_cd_cate_3 = icd9_cd, r(240/279 279*)
icd9 gen icd9_cd_cate_4 = icd9_cd, r(280/289 289*)
icd9 gen icd9_cd_cate_5 = icd9_cd, r(290/319 319*)
icd9 gen icd9_cd_cate_6 = icd9_cd, r(320/389 389*)
icd9 gen icd9_cd_cate_7 = icd9_cd, r(390/459 459*)
icd9 gen icd9_cd_cate_8 = icd9_cd, r(460/519 519*)
icd9 gen icd9_cd_cate_9 = icd9_cd, r(520/579 579*)
icd9 gen icd9_cd_cate_10 = icd9_cd, r(580/629 629*)
icd9 gen icd9_cd_cate_11 = icd9_cd, r(630/679 679*)
icd9 gen icd9_cd_cate_12 = icd9_cd, r(680/709 709*)
icd9 gen icd9_cd_cate_13 = icd9_cd, r(710/739 739*)
icd9 gen icd9_cd_cate_14 = icd9_cd, r(740/759 759*)
icd9 gen icd9_cd_cate_15 = icd9_cd, r(760/779 779*)
icd9 gen icd9_cd_cate_16 = icd9_cd, r(780/799 799*)
icd9 gen icd9_cd_cate_17 = icd9_cd, r(800/999 999*)
gen icd9_cd_cate_18 = substr(icd9_cd, 1, 1) == "V"
gen icd9_cd_cate = .
forvalues j = 0/18 {
	replace icd9_cd_cate = `j' if icd9_cd_cate_`j' == 1 & icd9_cd_cate == .
}
*
qui forvalues k = 0/18 {
	count if icd9_cd_cate == `k' & arm_value_1 == 1
	scalar n_1 = r(N)
	count if icd9_cd_cate == `k' & arm_value_0 == 1
	scalar n_0 = r(N)
	count if icd9_cd_cate == `k'
	scalar n_all = r(N)
	noi di `k' _col(10) n_0 _col(20) n_1 _col(30) n_all
}
*
********************************************************************************
cls
**# Supp T7
*** Supplementary Table 7. Sensitivity analyses on the 28-day all-cause hospitalization outcome 
**# Sen - Truncating IPCW at the 1st and 99th percentiles
capture mkdir "sen_trim"
capture mkdir "sen_trim/emulated"
capture mkdir "sen_trim/KM ipcw_stab_trim"

* Target trial emulation
qui forvalues k = 1/1 {
forvalues j = 0/500 {
*** Combine & Generate weights
	noi di "sen_trim" _col(15) "bs`j'" _col(30) "IPCW"
	if !fileexists("sen_trim/emulated/emulated bs`j'.dta") {
	use "main/emulated/emulated bs`j'.dta", clear
	* generate trimmed stabilized IPCW
	gen _ipcw_stab_trim = _ipcw_stab
	forvalues t = 0/27 {
		forvalues i = 0/1 {
			sum _ipcw_stab_trim if arm_value == `i' & tstart == `t', d
			replace _ipcw_stab_trim = r(p1) if arm_value == `i' & tstart == `t' & _ipcw_stab_trim < r(p1)
			replace _ipcw_stab_trim = r(p99) if arm_value == `i' & tstart == `t' & _ipcw_stab_trim > r(p99)
		}
	}
	*
	compress
	save "sen_trim/emulated/emulated bs`j'.dta", replace
	}
	
*** Generate KM estimate
	noi di "sen_trim" _col(15) "bs`j'" _col(30) "KM"
	if !fileexists("sen_trim/KM ipcw_stab_trim/KM ipcw_stab_trim bs`j'.dta") {
	use "sen_trim/emulated/emulated bs`j'.dta", clear
	count if invalid == 1
	if r(N) == 0 {
		stset tstop [pweight = _ipcw_stab_trim], enter(time tstart) failure(outcome)
		sts generate KM_s_w = s if arm_value == 1
		sts generate KM_ns_w = s if arm_value == 0
	}
	else {
		gen KM_s_w = .
		gen KM_ns_w = .
	}
	collapse (firstnm) KM_s_w KM_ns_w, by(fup bs invalid)
	save "sen_trim/KM ipcw_stab_trim/KM ipcw_stab_trim bs`j'.dta", replace
	
	if `j' > 1 {
		capture erase "sen_trim/emulated/emulated bs`j'.dta"
	}
	}
}
}
*
* Finalize bootstrap datasets
clear
qui forvalues k = 1/1 {
forvalues j = 0/500 {
	capture append using "sen_trim/KM ipcw_stab_trim/KM ipcw_stab_trim bs`j'.dta"
}
	gen hazard_s_w = 1 - KM_s_w
	gen hazard_ns_w = 1 - KM_ns_w
	gen odds_s_w = hazard_s_w/(1-hazard_s_w)
	gen odds_ns_w = hazard_ns_w/(1-hazard_ns_w)
	gen RR_w = hazard_s_w/hazard_ns_w
	gen diff_w = hazard_s_w - hazard_ns_w
	gen OR_w = odds_s_w / odds_ns_w
	compress
	save "sen_trim/KM ipcw_stab_trim bs_all.dta", replace
}
*
cls
* KM estimate
qui forvalues k = 1/1 {
	use "sen_trim/KM ipcw_stab_trim bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	
	sum hazard_s_w if bs == 0, d
	scalar hazard_s_mean = r(mean)
	centile hazard_s_w if bs > 0, centile(2.5 97.5)
	scalar hazard_s_cil = r(c_1)
	scalar hazard_s_ciu = r(c_2)
	
	sum hazard_ns_w if bs == 0, d
	scalar hazard_ns_mean = r(mean)
	centile hazard_ns_w if bs > 0, centile(2.5 97.5)
	scalar hazard_ns_cil = r(c_1)
	scalar hazard_ns_ciu = r(c_2)
	
	noi di "sen_trim" _col(10) hazard_s_mean _col(25) hazard_s_cil _col(40) hazard_s_ciu _col(55) ///
	hazard_ns_mean _col(70) hazard_s_cil _col(85) hazard_s_ciu
}
*
* Absolute risk reduction
qui forvalues k = 1/1 {
	use "sen_trim/KM ipcw_stab_trim bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	
	sum diff_w if bs == 0, d
	scalar bs_mean = r(mean)
	
	centile diff_w if bs > 0, centile(2.5 97.5)
	scalar bs_cil = r(c_1)
	scalar bs_ciu = r(c_2)
	noi di "sen_trim" _col(10) -bs_mean _col(25) -bs_ciu _col(40) -bs_cil
}
*
* Relative risk
qui forvalues k = 1/1 {
	use "sen_trim/KM ipcw_stab_trim bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	
	sum RR_w if bs == 0, d
	scalar bs_mean = r(mean)
	centile RR_w if bs > 0, centile(2.5 97.5)
	scalar bs_cil = r(c_1)
	scalar bs_ciu = r(c_2)
	noi di "sen_trim" _col(10) bs_mean _col(25) bs_cil _col(40) bs_ciu
}
*

**# Sen - Not restricted to 28 days
capture mkdir "sen_extend"
capture mkdir "sen_extend/cloned"
capture mkdir "sen_extend/split"
capture mkdir "sen_extend/Treatment"
capture mkdir "sen_extend/NoTreatment"
capture mkdir "sen_extend/emulated"
capture mkdir "sen_extend/KM ipcw_stab"

* Target trial emulation
qui forvalues k = 1/1 {
forvalues j = 0/500 {
*** Prepare dataset for analysis
	noi di "sen_extend" _col(15) "bs`j'" _col(30) "setup" 
	* Setup for trial emulation
	if !fileexists("sen_extend/KM ipcw_stab/KM ipcw_stab bs`j'.dta") {
	use "covid paxlovid pediatric.dta", clear
	if `j' > 0 {
		set seed `j'
		bsample
	}
	replace date_paxlovid = . if date_paxlovid >= date_admission
	* organize / rename / generate variables
	gen date_last_fu = min(date_death, mdy(02,12,2023))
	gen date_event = date_admission
	gen event = inrange(date_event, date_baseline, date_last_fu)
	gen fup_obs = date_event-date_baseline if event == 1
	replace fup_obs = date_last_fu-date_baseline if event == 0
	tab fup_obs
	gen time_to_treatment = date_paxlovid - date_baseline
	tab time_to_treatment
	gen treatment = inrange(date_paxlovid - date_baseline, 0, 5)
	* keep necessary variables
	keep pseudo_key fup_obs event time_to_treatment treatment group
	gen bs = `j'	
	
*** Cloning & censoring
	noi di "sen_extend" _col(15) "bs`j'" _col(30) "cloning" 
	* Prepare dataset for analysis
	stset fup_obs, failure(event)

	* Arm A: no treatment within 5 days (control: non-exposed group)
	gen outcomeA = _d // _d = `event'
	gen fupA = _t // _t = follow up time

	/// if the patient received treatment within 6 months:
	/// 1. no event outcome, since the patient survived till censoring (treatment)
	replace outcomeA = 0 if treatment==1 & time_to_treatment <=5 
	/// 2. follow up is censored at treatment
	replace fupA = time_to_treatment if treatment==1 & time_to_treatment <=5

	* Arm B: treatment within 5 days (treated: exposed group)
	gen outcomeB = _d 
	gen fupB = _t 

	/// if the patient survived the first 5 days and did not receive treatment within 5 days:
	/// 1. no event outcome if the patient survived the first 6 months
	replace outcomeB = 0 if (treatment==0 & _t>5) | (treatment==1 & time_to_treatment >5 & time_to_treatment !=.)
	/// 2. follow up is censored at 6 months
	replace fupB = 5 if (treatment==0 & _t>5) | (treatment==1 & time_to_treatment >5 & time_to_treatment != .)

	** append clones 
	preserve
		drop outcomeB fupB
		rename outcomeA outcome
		rename fupA fup
		gen arm = "NoTreatment"
		tempfile a
		save "`a'", replace
	restore
		drop outcomeA fupA
		rename outcomeB outcome
		rename fupB fup
		gen arm = "Treatment"	
		cap append using "`a'"

	// Weight models

	sort _all
	gen NewID = _n

	** add 1 day to 0-survivors
	replace fup= 1 if fup==0

	** Weight model: define survival time and event indicator	
	* treatment Arm
	** Case 1: they do not deviate at time of treatment, but are not at risk of deviating any more
	gen wm_fup = time_to_treatment if arm == "Treatment" & time_to_treatment<=5 & time_to_treatment!=. & treatment == 1 
	gen wm_outcome = 0 if arm == "Treatment" & time_to_treatment<=5 & time_to_treatment!=. & treatment == 1 

	** Case 2: they deviate at 5 days
	replace wm_fup = 5 if arm == "Treatment" & ((treatment == 0 & fup >= 5) | (time_to_treatment>5 & treatment == 1))
	replace wm_outcome = 1 if arm == "Treatment" & ((treatment == 0 & fup >= 5) | (time_to_treatment>5 & treatment == 1)) 

	** Case 3: they do not deviate, but we need to keep their survival as observed and censor them as we do not know what happens afterwards
	replace wm_fup = fup if arm == "Treatment" & treatment == 0 & fup < 5
	replace wm_outcome = 0 if arm == "Treatment" & treatment == 0 & fup < 5
	** add 1 days to 0-survivors
	replace wm_fup= 1 if arm == "Treatment" & wm_fup==0 

	* No treatment Arm
	** Case 1: they do not deviate at time of treatment, but are not at risk of deviating any more
	replace wm_fup = time_to_treatment if arm == "NoTreatment" & time_to_treatment<=5 & treatment == 1 
	replace wm_outcome = 1 if arm == "NoTreatment" & time_to_treatment<=5 & treatment == 1 

	** Case 2: they deviate at 5 days
	replace wm_fup = 5 if arm == "NoTreatment" & ((treatment == 0 & fup >= 5) | (time_to_treatment>5 & treatment == 1)) 
	replace wm_outcome = 0 if arm == "NoTreatment" & ((treatment == 0 & fup >= 5) | (time_to_treatment>5 & treatment == 1)) 

	** Case 3: they do not deviate, but we need to keep their survival as observed and censor them as we do not know what happens afterwards
	replace wm_fup = fup if arm == "NoTreatment" & treatment == 0 & fup < 5 
	replace wm_outcome = 0 if arm == "NoTreatment" & treatment == 0 & fup < 5

	** add 1 days to 0-survivors
	replace wm_fup= 1 if arm == "NoTreatment" & wm_fup == 0

	gen arm_value = 1 if arm == "Treatment"
	replace arm_value = 0 if arm == "NoTreatment"
	
	drop arm
	order pseudo_key treatment arm_value
	
	if `j' == 0 {
		save "sen_extend/cloned/cloned bs0.dta", replace
	}
	
*** Split times
	noi di "sen_extend" _col(15) "bs`j'" _col(30) "split"

	** times of event
	stset fup, fail(outcome) id(NewID)
	stsplit, at(failures)

	gen tstart = _t0
	gen tstop = _t

	** times of censoring
	gen TrialEmul_cens = 1-outcome
	stset fup, fail(TrialEmul_cens) id(NewID)
	stsplit, at(failures)

	replace tstart = _t0 if tstart<_t0 & _t0 != . & _t != .
	replace tstop = _t if tstop>_t & _t0 != . & _t != .

	order pseudo_key arm_value tstart tstop
	sort _all
	compress
	save "sen_extend/split/split bs`j'.dta", replace

	}
	
*** IPCW - Treatment arm
	noi di "sen_extend" _col(15) "bs`j'" _col(30) "Treatment"
	if !fileexists("sen_extend/KM ipcw_stab/KM ipcw_stab bs`j'.dta") {
	use "sen_extend/split/split bs`j'.dta", clear
	keep if arm_value == 1
	replace outcome = . if TrialEmul_cens == .
	merge m:1 pseudo_key using "covid pediatric paxlovid characteristics.dta", keep(3) nogen

	* adapt records to the long format
	sort NewID tstart
	bysort NewID (tstart): replace wm_outcome = 0 if _n!=_N

	* Weight model:
	stset tstop, origin(time tstart) failure(wm_outcome) id(NewID)
	logit wm_outcome age gender vaccine_status_2 vaccine_status_3 ym_covid_cate_2 ym_covid_cate_3 symptomatic lung immuno_hist healthcare ///
	i.symptomatic##i.ym_covid_cate_3 i.symptomatic##i.gender i.symptomatic##i.healthcare if tstart == 4
	predict pr_cens if e(sample), pr
	replace pr_cens = 0 if pr_cens == .
	gen pr_uncens = 1-pr_cens
	bysort arm_value pseudo_key (tstart): replace pr_uncens=pr_uncens*pr_uncens[_n-1] if _n!=1
	gen weight = 1/pr_uncens
	replace weight = 1 if weight == .
	gen invalid = 0
	keep pseudo_key arm_value tstart tstop event fup_obs time_to_treatment treatment bs outcome fup NewID wm_fup wm_outcome TrialEmul_cens weight invalid
	compress
	save "sen_extend/Treatment/Treatment bs`j'.dta", replace
	}
	
*** IPCW - Control arm
	noi di "sen_extend" _col(15) "bs`j'" _col(30) "NoTreatment"
	if !fileexists("sen_extend/KM ipcw_stab/KM ipcw_stab bs`j'.dta") {
	use "sen_extend/split/split bs`j'.dta", clear
	keep if arm_value == 0
	replace outcome = . if TrialEmul_cens == .
	merge m:1 pseudo_key using "covid pediatric paxlovid characteristics.dta", keep(3) nogen

	* adapt records to the long format
	sort NewID tstart
	bysort NewID: replace wm_outcome = 0 if _n!=_N

	* Weight model:
	stset tstop, origin(time tstart) failure(wm_outcome) id(NewID)
	xi: logit wm_outcome tstart age gender vaccine_status_2 vaccine_status_3 ym_covid_cate_2 ym_covid_cate_3 symptomatic immuno_hist healthcare, cluster(NewID)
	predict pr_cens if e(sample), pr
	replace pr_cens = 0 if pr_cens == .
	gen pr_uncens = 1-pr_cens
	bysort arm_value pseudo_key (tstart): replace pr_uncens=pr_uncens*pr_uncens[_n-1] if _n!=1
	gen weight = 1/pr_uncens
	replace weight = 1 if weight == .
	gen invalid = 0
	keep pseudo_key arm_value tstart tstop event fup_obs time_to_treatment treatment bs outcome fup NewID wm_fup wm_outcome TrialEmul_cens weight invalid
	compress
	save "sen_extend/NoTreatment/Notreatment bs`j'.dta", replace
	}
	
*** Combine & Generate weights
	noi di "sen_extend" _col(15) "bs`j'" _col(30) "IPCW"
	if !fileexists("sen_extend/KM ipcw_stab/KM ipcw_stab bs`j'.dta") {
	use "sen_extend/Treatment/Treatment bs`j'.dta", clear
	append using "sen_extend/NoTreatment/Notreatment bs`j'.dta"

	// create a new ID variable for each clone in each arm
	tostring NewID, gen(NewID_str)
	replace NewID_str = "0000000" + NewID_str if length(NewID_str)== 1
	replace NewID_str = "000000" + NewID_str if length(NewID_str)== 2
	replace NewID_str = "00000" + NewID_str if length(NewID_str)== 3
	replace NewID_str = "0000" + NewID_str if length(NewID_str)== 4
	replace NewID_str = "000" + NewID_str if length(NewID_str)== 5
	replace NewID_str = "00" + NewID_str if length(NewID_str)== 6
	replace NewID_str = "0" + NewID_str if length(NewID_str)== 7
	
	gen Anal_ID = "1" + NewID_str if arm_value == 1
	replace Anal_ID = "2" + NewID_str if arm_value == 0

	replace weight = 1 if wm_outcome == 1 & tstop == 5

	keep pseudo_key tstart tstop fup_obs event time_to_treatment treatment bs outcome fup wm_fup wm_outcome TrialEmul_cens weight Anal_ID arm_value invalid
	destring Anal_ID, replace

	* generate denominator
	rename weight _ipcw
	gen prob_denom = 1/_ipcw
	* generate numerator
	bysort arm_value tstart : gen N_group = _N
	bysort arm_value tstart : egen N_censored = sum(wm_outcome)
	gen prob_uncensored = (N_group-N_censored)/N_group
	gen prob_num = prob_uncensored
	sort arm_value pseudo_key tstart
	by arm_value pseudo_key : replace prob_num=prob_num*prob_num[_n-1] if _n!=1
	* generate stabilized IPCW
	gen _ipcw_stab = prob_num / prob_denom
	compress
	save "sen_extend/emulated/emulated bs`j'.dta", replace
	
	if `j' > 1 {
		capture erase "sen_extend/split/split bs`j'.dta"
		capture erase "sen_extend/Treatment/Treatment bs`j'.dta"
		capture erase "sen_extend/NoTreatment/NoTreatment bs`j'.dta"
	}
	
	}
	
*** Generate KM estimate
	noi di "sen_extend" _col(15) "bs`j'" _col(30) "KM"
	if !fileexists("sen_extend/KM ipcw_stab/KM ipcw_stab bs`j'.dta") {
	use "sen_extend/emulated/emulated bs`j'.dta", clear
	count if invalid == 1
	if r(N) == 0 {
		stset tstop [pweight = _ipcw_stab], enter(time tstart) failure(outcome)
		sts generate KM_s_w = s if arm_value == 1
		sts generate KM_ns_w = s if arm_value == 0
	}
	else {
		gen KM_s_w = .
		gen KM_ns_w = .
	}
	collapse (firstnm) KM_s_w KM_ns_w, by(fup bs invalid)
	save "sen_extend/KM ipcw_stab/KM ipcw_stab bs`j'.dta", replace
	if `j' > 1 {
		capture erase "sen_extend/emulated/emulated bs`j'.dta"
	}
	}
}
}
*
* Finalize bootstrap datasets
clear
qui forvalues k = 1/1 {
forvalues j = 0/500 {
	capture append using "sen_extend/KM ipcw_stab/KM ipcw_stab bs`j'.dta"
	*erase "sen_extend/KM/KM bs`j'.dta"
}
	gen hazard_s_w = 1 - KM_s_w
	gen hazard_ns_w = 1 - KM_ns_w
	gen odds_s_w = hazard_s_w/(1-hazard_s_w)
	gen odds_ns_w = hazard_ns_w/(1-hazard_ns_w)
	gen RR_w = hazard_s_w/hazard_ns_w
	gen diff_w = hazard_s_w - hazard_ns_w
	gen OR_w = odds_s_w / odds_ns_w
	compress
	save "sen_extend/KM ipcw_stab bs_all.dta", replace
}
*
cls
* KM estimate
qui forvalues k = 1/1 {
	use "sen_extend/KM ipcw_stab bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	
	sum hazard_s_w if bs == 0, d
	scalar hazard_s_mean = r(mean)
	centile hazard_s_w if bs > 0, centile(2.5 97.5)
	scalar hazard_s_cil = r(c_1)
	scalar hazard_s_ciu = r(c_2)
	
	sum hazard_ns_w if bs == 0, d
	scalar hazard_ns_mean = r(mean)
	centile hazard_ns_w if bs > 0, centile(2.5 97.5)
	scalar hazard_ns_cil = r(c_1)
	scalar hazard_ns_ciu = r(c_2)
	
	noi di "sen_extend" _col(10) hazard_s_mean _col(25) hazard_s_cil _col(40) hazard_s_ciu _col(55) ///
	hazard_ns_mean _col(70) hazard_s_cil _col(85) hazard_s_ciu
}
*
* Absolute risk reduction
qui forvalues k = 1/1 {
	use "sen_extend/KM ipcw_stab bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	
	sum diff_w if bs == 0, d
	scalar bs_mean = r(mean)
	
	centile diff_w if bs > 0, centile(2.5 97.5)
	scalar bs_cil = r(c_1)
	scalar bs_ciu = r(c_2)
	noi di "sen_extend" _col(10) -bs_mean _col(25) -bs_ciu _col(40) -bs_cil
}
*
* Relative risk
qui forvalues k = 1/1 {
	use "sen_extend/KM ipcw_stab bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	
	sum RR_w if bs == 0, d
	scalar bs_mean = r(mean)
	centile RR_w if bs > 0, centile(2.5 97.5)
	scalar bs_cil = r(c_1)
	scalar bs_ciu = r(c_2)
	noi di "sen_extend" _col(10) bs_mean _col(25) bs_cil _col(40) bs_ciu
}
*
* N
qui forvalues k = 1/1 {
	use "sen_extend/cloned/cloned bs0.dta", clear
	keep if (arm_value == 1 & group == 1) | (arm_value == 0 & group == 0)
	count if group == 1
	scalar N_1 = r(N)
	count if group == 0
	scalar N_0 = r(N)
	count if group == 1 & outcome == 1
	scalar n_e_1 = r(N)
	count if group == 0 & outcome == 1
	scalar n_e_0 = r(N)
	noi di "sen_extend" _col(15) N_1 _col(30) N_0 _col(45) n_e_1 _col(60) n_e_0
}
*
qui forvalues k = 1/1 {
	if `k' == 1 {
		noi di "sen" _col(15) "N_1" _col(30) "N_0" _col(45) "n_e_1" _col(60) "n_e_0"
	}
	use "sen_extend/cloned/cloned bs0.dta", clear
	count if arm_value == 1
	scalar N_1 = r(N)
	count if arm_value == 1 & outcome == 1
	scalar n_e_1 = r(N)
	count if arm_value == 0
	scalar N_0 = r(N)
	count if arm_value == 0 & outcome == 1
	scalar n_e_0 = r(N)
	
	noi di "sen_extend" _col(15) N_1 _col(30) N_0 _col(45) n_e_1 _col(60) n_e_0
}
*

**# Supp T8
*** Supplementary Table 8. Cumulative incidences of secondary study outcomes at 28 days in nirmatrelvir/ritonavir and control groups
use "covid paxlovid pediatric.dta", clear
keep if group < .
merge 1:1 pseudo_key using "date_covid_hosp", keep(1 3) nogen
merge 1:1 pseudo_key using "date_death_ip", keep(1 3) nogen
merge 1:1 pseudo_key using "date_admission", keep(1 3) nogen
merge 1:1 pseudo_key using "date_intubation", keep(1 3) nogen
merge 1:1 pseudo_key using "date_ICU.dta", keep(1 3) nogen
merge 1:1 pseudo_key using "date_ae.dta", keep(1 3) nogen
merge 1:1 pseudo_key using "date_outcome.dta", keep(1 3) nogen
merge 1:1 pseudo_key using "date_ali.dta", keep(1 3) nogen
merge 1:1 pseudo_key using "date_aki.dta", keep(1 3) nogen
merge 1:1 pseudo_key using "date_ARDS.dta", keep(1 3) nogen
merge 1:1 pseudo_key using "date_ctv_30.dta", keep(1 3) nogen
gen date_last_fu = min(date_death, mdy(02,12,2023), date_baseline + 28)

foreach var in death admission covid_hosp ae death_ip intubation ICU outcome_encephalitis outcome_encephalopathy outcome_seizure outcome_stroke outcome_gbs outcome_carditis outcome_myocarditis outcome_pericarditis outcome_arrhythmias outcome_cardiogenic_shock outcome_thrombosis outcome_thromboembolism outcome_thromboembolic outcome_thrombocytopenia outcome_MISC ali aki ARDS ctv_30 {
	gen `var' = inrange(date_`var', date_baseline, date_last_fu) if date_`var' > date_baseline
	gen `var'_dur = date_`var' - date_baseline if `var' == 1
	replace `var'_dur = date_last_fu - date_baseline if `var' == 0
	gen `var'_dur_r = `var'_dur if `var' == 1
}
gen date_composite = .
foreach var in death admission {
	replace date_composite = min(date_composite, date_`var') if `var' == 1
}
*
gen date_composite_ip = .
foreach var in death_ip intubation ICU {
	replace date_composite_ip = min(date_composite_ip, date_`var') if `var' == 1
}
*
gen date_outcome_cardio = .
foreach var in outcome_myocarditis outcome_pericarditis outcome_arrhythmias outcome_cardiogenic_shock {
	replace date_outcome_cardio = min(date_outcome_cardio, date_`var') if `var' == 1
}
*
foreach var in composite composite_ip outcome_cardio {
	gen `var' = inrange(date_`var', date_baseline, date_last_fu) if date_`var' > date_baseline
	gen `var'_dur = date_`var' - date_baseline if `var' == 1
	replace `var'_dur = date_last_fu - date_baseline if `var' == 0
	gen `var'_dur_r = `var'_dur if `var' == 1
}
*
format date* %td
compress
save "covid paxlovid pediatric secondary", replace

cls
**# Secondary outcome - Target trial emulation
qui foreach event in covid_hosp composite_ip ICU ali aki ARDS {
capture mkdir "`event'"
capture mkdir "`event'/prepare"
capture mkdir "`event'/cloned"
capture mkdir "`event'/split"
capture mkdir "`event'/Treatment"
capture mkdir "`event'/NoTreatment"
capture mkdir "`event'/emulated"
capture mkdir "`event'/KM ipcw_stab"

forvalues j = 0/500 {
*** Prepare dataset for analysis
	noi di "`event'" _col(30) "bs`j'" _col(45) "setup" 
	* Setup for trial emulation
	if !fileexists("`event'/KM ipcw_stab/KM ipcw_stab subgp_1 bs`j'.dta") {
	use "covid paxlovid pediatric secondary.dta", clear
	drop date_last_fu
	if `j' > 0 {
		set seed `j'
		bsample
	}
	replace date_paxlovid = . if date_paxlovid >= date_admission
	* organize / rename / generate variables
	gen date_last_fu = min(date_death, mdy(02,12,2023), date_baseline + 28)
	gen date_event = date_`event'
	gen event = inrange(date_event, date_baseline, date_last_fu)
	gen fup_obs = min(date_event-date_baseline, 28) if event == 1
	replace fup_obs = min(date_last_fu-date_baseline, 28) if event == 0
	gen time_to_treatment = date_paxlovid - date_baseline
	tab time_to_treatment
	gen treatment = inrange(date_paxlovid - date_baseline, 0, 5)
	* keep necessary variables
	keep pseudo_key fup_obs event time_to_treatment treatment group
	gen bs = `j'	
	
*** Cloning & censoring
	noi di "`event'" _col(30) "bs`j'" _col(45) "cloning" 
	* Prepare dataset for analysis
	stset fup_obs, failure(event)

	* Arm A: no treatment within 5 days (control: non-exposed group)
	gen outcomeA = _d // _d = `event'
	gen fupA = _t // _t = follow up time

	/// if the patient received treatment within 6 months:
	/// 1. no event outcome, since the patient survived till censoring (treatment)
	replace outcomeA = 0 if treatment==1 & time_to_treatment <=5 
	/// 2. follow up is censored at treatment
	replace fupA = time_to_treatment if treatment==1 & time_to_treatment <=5

	* Arm B: treatment within 5 days (treated: exposed group)
	gen outcomeB = _d 
	gen fupB = _t 

	/// if the patient survived the first 5 days and did not receive treatment within 5 days:
	/// 1. no event outcome if the patient survived the first 6 months
	replace outcomeB = 0 if (treatment==0 & _t>5) | (treatment==1 & time_to_treatment >5 & time_to_treatment !=.)
	/// 2. follow up is censored at 6 months
	replace fupB = 5 if (treatment==0 & _t>5) | (treatment==1 & time_to_treatment >5 & time_to_treatment != .)

	** append clones 
	preserve
		drop outcomeB fupB
		rename outcomeA outcome
		rename fupA fup
		gen arm = "NoTreatment"
		tempfile a
		save "`a'", replace
	restore
		drop outcomeA fupA
		rename outcomeB outcome
		rename fupB fup
		gen arm = "Treatment"	
		cap append using "`a'"

	// Weight models

	sort _all
	gen NewID = _n

	** add 1 day to 0-survivors
	replace fup= 1 if fup==0

	** Weight model: define survival time and event indicator	
	* treatment Arm
	** Case 1: they do not deviate at time of treatment, but are not at risk of deviating any more
	gen wm_fup = time_to_treatment if arm == "Treatment" & time_to_treatment<=5 & time_to_treatment!=. & treatment == 1 
	gen wm_outcome = 0 if arm == "Treatment" & time_to_treatment<=5 & time_to_treatment!=. & treatment == 1 

	** Case 2: they deviate at 5 days
	replace wm_fup = 5 if arm == "Treatment" & ((treatment == 0 & fup >= 5) | (time_to_treatment>5 & treatment == 1))
	replace wm_outcome = 1 if arm == "Treatment" & ((treatment == 0 & fup >= 5) | (time_to_treatment>5 & treatment == 1)) 

	** Case 3: they do not deviate, but we need to keep their survival as observed and censor them as we do not know what happens afterwards
	replace wm_fup = fup if arm == "Treatment" & treatment == 0 & fup < 5
	replace wm_outcome = 0 if arm == "Treatment" & treatment == 0 & fup < 5
	** add 1 days to 0-survivors
	replace wm_fup= 1 if arm == "Treatment" & wm_fup==0 

	* No treatment Arm
	** Case 1: they do not deviate at time of treatment, but are not at risk of deviating any more
	replace wm_fup = time_to_treatment if arm == "NoTreatment" & time_to_treatment<=5 & treatment == 1 
	replace wm_outcome = 1 if arm == "NoTreatment" & time_to_treatment<=5 & treatment == 1 

	** Case 2: they deviate at 5 days
	replace wm_fup = 5 if arm == "NoTreatment" & ((treatment == 0 & fup >= 5) | (time_to_treatment>5 & treatment == 1)) 
	replace wm_outcome = 0 if arm == "NoTreatment" & ((treatment == 0 & fup >= 5) | (time_to_treatment>5 & treatment == 1)) 

	** Case 3: they do not deviate, but we need to keep their survival as observed and censor them as we do not know what happens afterwards
	replace wm_fup = fup if arm == "NoTreatment" & treatment == 0 & fup < 5 
	replace wm_outcome = 0 if arm == "NoTreatment" & treatment == 0 & fup < 5

	** add 1 days to 0-survivors
	replace wm_fup= 1 if arm == "NoTreatment" & wm_fup == 0

	gen arm_value = 1 if arm == "Treatment"
	replace arm_value = 0 if arm == "NoTreatment"
	
	drop arm
	order pseudo_key treatment arm_value
	
	if `j' == 0 {
		save "`event'/cloned/cloned subgp_1 bs0.dta", replace
	}
	
*** Split times
	noi di "`event'" _col(30) "bs`j'" _col(45) "split" 

	** times of event
	stset fup, fail(outcome) id(NewID)
	stsplit, at(failures)

	gen tstart = _t0
	gen tstop = _t

	** times of censoring
	gen TrialEmul_cens = 1-outcome
	stset fup, fail(TrialEmul_cens) id(NewID)
	stsplit, at(failures)

	replace tstart = _t0 if tstart<_t0 & _t0 != . & _t != .
	replace tstop = _t if tstop>_t & _t0 != . & _t != .

	order pseudo_key arm_value tstart tstop
	sort _all
	compress
	save "`event'/split/split subgp_1 bs`j'.dta", replace
	}
	
*** IPCW - Treatment arm
	noi di "`event'" _col(30) "bs`j'" _col(45) "Treatment" 
	if !fileexists("`event'/KM ipcw_stab/KM ipcw_stab subgp_1 bs`j'.dta") {
	use "`event'/split/split subgp_1 bs`j'.dta", clear
	keep if arm_value == 1
	replace outcome = . if TrialEmul_cens == .
	merge m:1 pseudo_key using "covid pediatric paxlovid characteristics.dta", keep(3) nogen

	* adapt records to the long format
	sort NewID tstart
	bysort NewID (tstart): replace wm_outcome = 0 if _n!=_N

	* Weight model:
	stset tstop, origin(time tstart) failure(wm_outcome) id(NewID)
	logit wm_outcome age gender vaccine_status_2 vaccine_status_3 ym_covid_cate_2 ym_covid_cate_3 symptomatic lung immuno_hist healthcare ///
	i.symptomatic##i.ym_covid_cate_3 i.symptomatic##i.gender i.symptomatic##i.healthcare if tstop == 5
	predict pr_cens if e(sample), pr
	replace pr_cens = 0 if pr_cens == .
	gen pr_uncens = 1-pr_cens
	bysort arm_value pseudo_key (tstart): replace pr_uncens=pr_uncens*pr_uncens[_n-1] if _n!=1
	gen weight = 1/pr_uncens
	replace weight = 1 if weight == .
	gen invalid = 0
	keep pseudo_key arm_value tstart tstop event fup_obs time_to_treatment treatment bs outcome fup NewID wm_fup wm_outcome TrialEmul_cens weight invalid
	compress
	save "`event'/Treatment/Treatment subgp_1 bs`j'.dta", replace
	}
	
*** IPCW - Control arm
	noi di "`event'" _col(30) "bs`j'" _col(45) "NoTreatment" 
	if !fileexists("`event'/KM ipcw_stab/KM ipcw_stab subgp_1 bs`j'.dta") {
	use "`event'/split/split subgp_1 bs`j'.dta", clear
	keep if arm_value == 0
	replace outcome = . if TrialEmul_cens == .
	merge m:1 pseudo_key using "covid pediatric paxlovid characteristics.dta", keep(3) nogen

	* adapt records to the long format
	sort NewID tstart
	bysort NewID: replace wm_outcome = 0 if _n!=_N

	* Weight model:
	stset tstop, origin(time tstart) failure(wm_outcome) id(NewID)
	
	capture xi: logit wm_outcome tstart age gender vaccine_status_2 vaccine_status_3 ym_covid_cate_2 ym_covid_cate_3 symptomatic immuno_hist healthcare, cluster(NewID)
	if _rc == 0 {
	predict pr_cens if e(sample), pr
	replace pr_cens = 0 if pr_cens == .
	gen pr_uncens = 1-pr_cens
	bysort arm_value pseudo_key (tstart): replace pr_uncens=pr_uncens*pr_uncens[_n-1] if _n!=1
	gen weight = 1/pr_uncens
	replace weight = 1 if weight == .
	gen invalid = 0
	}
	if _rc != 0 {
		gen invalid = 1
		gen weight = .
	}
	keep pseudo_key arm_value tstart tstop event fup_obs time_to_treatment treatment bs outcome fup NewID wm_fup wm_outcome TrialEmul_cens weight invalid
	compress
	save "`event'/NoTreatment/Notreatment subgp_1 bs`j'.dta", replace
	}
	
*** Combine & Generate weights
	noi di "`event'" _col(30) "bs`j'" _col(45) "IPCW" 
	if !fileexists("`event'/KM ipcw_stab/KM ipcw_stab subgp_1 bs`j'.dta") {
	use "`event'/Treatment/Treatment subgp_1 bs`j'.dta", clear
	append using "`event'/NoTreatment/Notreatment subgp_1 bs`j'.dta"

	// create a new ID variable for each clone in each arm
	tostring NewID, gen(NewID_str)
	replace NewID_str = "0000000" + NewID_str if length(NewID_str)== 1
	replace NewID_str = "000000" + NewID_str if length(NewID_str)== 2
	replace NewID_str = "00000" + NewID_str if length(NewID_str)== 3
	replace NewID_str = "0000" + NewID_str if length(NewID_str)== 4
	replace NewID_str = "000" + NewID_str if length(NewID_str)== 5
	replace NewID_str = "00" + NewID_str if length(NewID_str)== 6
	replace NewID_str = "0" + NewID_str if length(NewID_str)== 7
	
	gen Anal_ID = "1" + NewID_str if arm_value == 1
	replace Anal_ID = "2" + NewID_str if arm_value == 0

	replace weight = 1 if wm_outcome == 1 & tstop == 5

	keep pseudo_key tstart tstop fup_obs event time_to_treatment treatment bs outcome fup wm_fup wm_outcome TrialEmul_cens weight Anal_ID arm_value invalid
	destring Anal_ID, replace

	* generate denominator
	rename weight _ipcw
	gen prob_denom = 1/_ipcw
	* generate numerator
	bysort arm_value tstart : gen N_group = _N
	bysort arm_value tstart : egen N_censored = sum(wm_outcome)
	gen prob_uncensored = (N_group-N_censored)/N_group
	gen prob_num = prob_uncensored
	sort arm_value pseudo_key tstart
	by arm_value pseudo_key : replace prob_num=prob_num*prob_num[_n-1] if _n!=1
	* generate stabilized IPCW
	gen _ipcw_stab = prob_num / prob_denom
	compress
	save "`event'/emulated/emulated subgp_1 bs`j'.dta", replace
	
	if `j' > 1 {
		capture erase "`event'/split/split subgp_1 bs`j'.dta"
		capture erase "`event'/Treatment/Treatment subgp_1 bs`j'.dta"
		capture erase "`event'/NoTreatment/NoTreatment subgp_1 bs`j'.dta"
	}
	}
	
*** Generate KM estimate
	noi di "`event'" _col(30) "bs`j'" _col(45) "KM" 
	if !fileexists("`event'/KM ipcw_stab/KM ipcw_stab subgp_1 bs`j'.dta") {
	use "`event'/emulated/emulated subgp_1 bs`j'.dta", clear
	count if invalid == 1
	if r(N) == 0 {
		stset tstop [pweight = _ipcw_stab], enter(time tstart) failure(outcome)
		sts generate KM_s_w = s if arm_value == 1
		sts generate KM_ns_w = s if arm_value == 0
	}
	else {
		gen KM_s_w = .
		gen KM_ns_w = .
	}
	collapse (firstnm) KM_s_w KM_ns_w, by(fup bs invalid)
	save "`event'/KM ipcw_stab/KM ipcw_stab subgp_1 bs`j'.dta", replace
	
	if `j' > 1 {
		capture erase "`event'/emulated/emulated subgp_1 bs`j'.dta"
	}
	}
}
}
*
cls
* Finalize bootstrap datasets
qui foreach event in death covid_hosp composite_ip death_ip intubation ICU outcome_MISC ali aki ARDS {
clear
forvalues j = 0/500 {
	capture append using "`event'/KM ipcw_stab/KM ipcw_stab subgp_1 bs`j'.dta"
	*erase "main/KM/KM bs`j'.dta"
}
	gen hazard_s_w = 1 - KM_s_w
	gen hazard_ns_w = 1 - KM_ns_w
	gen odds_s_w = hazard_s_w/(1-hazard_s_w)
	gen odds_ns_w = hazard_ns_w/(1-hazard_ns_w)
	gen RR_w = hazard_s_w/hazard_ns_w
	gen diff_w = hazard_s_w - hazard_ns_w
	gen OR_w = odds_s_w / odds_ns_w
	compress
	count if bs > 0 & fup == 28 & invalid == 0
	noi di "`event'" _col(30) r(N)
	save "`event'/KM ipcw_stab subgp_1 bs_all.dta", replace
}
*
cls
* KM estimate
qui foreach event in death covid_hosp composite_ip death_ip intubation ICU outcome_MISC ali aki ARDS {
	use "`event'/KM ipcw_stab subgp_1 bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	
	sum hazard_s_w if bs == 0, d
	scalar hazard_s_mean = r(mean)
	centile hazard_s_w if bs > 0, centile(2.5 97.5)
	scalar hazard_s_cil = r(c_1)
	scalar hazard_s_ciu = r(c_2)
	
	sum hazard_ns_w if bs == 0, d
	scalar hazard_ns_mean = r(mean)
	centile hazard_ns_w if bs > 0, centile(2.5 97.5)
	scalar hazard_ns_cil = r(c_1)
	scalar hazard_ns_ciu = r(c_2)
	
	noi di "`event'" _col(30) hazard_s_mean _col(45) hazard_s_cil _col(60) hazard_s_ciu _col(75) ///
	hazard_ns_mean _col(90) hazard_s_cil _col(105) hazard_s_ciu
}
*
* Absolute risk reduction
qui foreach event in death covid_hosp composite_ip death_ip intubation ICU outcome_MISC ali aki ARDS {
	use "`event'/KM ipcw_stab subgp_1 bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	
	sum diff_w if bs == 0, d
	scalar bs_mean = r(mean)
	
	centile diff_w if bs > 0, centile(2.5 97.5)
	scalar bs_cil = r(c_1)
	scalar bs_ciu = r(c_2)
	noi di "`event'" _col(30) -bs_mean _col(45) -bs_ciu _col(60) -bs_cil
}
*
* Relative risk
qui foreach event in death covid_hosp composite_ip death_ip intubation ICU outcome_MISC ali aki ARDS {
	use "`event'/KM ipcw_stab subgp_1 bs_all.dta", clear
	
	keep if KM_s_w < . & KM_ns_w < .
	drop if invalid == 1
	bysort bs (fup) : keep if _n == _N
	
	sum RR_w if bs == 0, d
	scalar bs_mean = r(mean)
	centile RR_w if bs > 0, centile(2.5 97.5)
	scalar bs_cil = r(c_1)
	scalar bs_ciu = r(c_2)
	noi di "`event'" _col(30) bs_mean _col(45) bs_cil _col(60) bs_ciu
}
*
* N
qui foreach event in death covid_hosp composite_ip death_ip intubation ICU outcome_MISC ali aki ARDS {
	use "`event'/cloned/cloned subgp_1 bs0.dta", clear
	keep if (arm_value == 1 & group == 1) | (arm_value == 0 & group == 0)
	count if group == 1
	scalar N_1 = r(N)
	count if group == 0
	scalar N_0 = r(N)
	count if group == 1 & outcome == 1
	scalar n_e_1 = r(N)
	count if group == 0 & outcome == 1
	scalar n_e_0 = r(N)
	noi di "`event'" _col(30) N_1 _col(45) N_0 _col(60) n_e_1 _col(75) n_e_0
}
*
qui foreach event in death covid_hosp composite_ip death_ip intubation ICU outcome_MISC ali aki ARDS {
	use "`event'/cloned/cloned subgp_1 bs0.dta", clear
	count if arm_value == 1
	scalar N_1 = r(N)
	count if arm_value == 1 & outcome == 1
	scalar n_e_1 = r(N)
	count if arm_value == 0
	scalar N_0 = r(N)
	count if arm_value == 0 & outcome == 1
	scalar n_e_0 = r(N)
	
	noi di "`event'" _col(30) N_1 _col(45) N_0 _col(60) n_e_1 _col(75) n_e_0
}
*

********************************************************************************
**# Figure - KM curves
clear
clear
use "main/KM ipcw_stab subgp_1 bs_all.dta", clear
tsset bs fup
tsfill
sort bs fup
foreach var of varlist KM_s_w-OR_w {
	replace `var' = `var'[_n-1] if `var' == .
}
*
sum fup, d
forvalues t = 1/`r(max)' {
foreach var in KM hazard {
foreach gp in s ns {
	capture gen `var'_`gp'_w_l = .
	capture gen `var'_`gp'_w_u = .
	_pctile `var'_`gp'_w if fup == `t' & bs > 0, p(2.5 97.5)
	replace `var'_`gp'_w_l = r(r1) if fup == `t'
	replace `var'_`gp'_w_u = r(r2) if fup == `t'
}
}
}
*
keep if bs == 0
insobs 1
foreach var of varlist _all {
	replace `var' = 0 if `var' == .
}

sort bs fup
rename fup tstart
compress
save "KM_CI", replace

use "main/emulated/emulated subgp_1 bs0.dta", clear
merge m:1 tstart using "KM_CI", keep(1 3)
stset tstop [pweight = _ipcw_stab], enter(time tstart) failure(outcome)
sts graph, failure by(arm_value) xlabel(0(7)28) ytitle("%", orientation(horizontal)) ///
ylabel(0 "0" .005 "0.5" .01 "1.0", format(%5.1f)) tmax(28) ttitle("Days") ///
plot1opts(lcolor(maroon)lpattern(dash)) plot2opts(lcolor(navy)) ///
title("Hospitalization") ///
risktable(0(7)28, order(2 "Nirmatrelvir/ritonavir" 1 "Control")) ///
legend(order(2 "Nirmatrelvir/ritonavir" 1 "Control") size(*0.8) pos(5) ring(0) col(2))
addplot : rarea hazard_s_w_l hazard_s_w_u _t if arm_value == 1, sort connect(step) color(navy%30) ylabel(0 "0" .005 "0.5" .01 "1.0", format(%5.1f)) xlabel(0(7)28)
addplot : rarea hazard_ns_w_l hazard_ns_w_u _t if arm_value == 0, sort connect(step) color(maroon%30) ylabel(0 "0" .005 "0.5" .01 "1.0", format(%5.1f)) xlabel(0(7)28) legend(order(2 "Nirmatrelvir/ritonavir" 1 "Control") size(*0.8) pos(5) ring(0) col(2))

********************************************************************************
**# Supp Fig 1
*** Supplementary Figure 1. Distribution of logarithm of the inverse probability of censoring weights during treatment initiation period (day 1-5) in nirmatrelvir/ritonavir and control groups
use "main/emulated/emulated subgp_1 bs0.dta", clear
keep if inrange(tstart,0,4)
replace tstart = tstart + 1
gen arm_value_0 = 1 - arm_value
gen ipcw_log = log10(_ipcw)
gen ipcw_stab_log = log10(_ipcw_stab)
keep *_log arm_value_0 tstart

graph box ipcw_log, over(tstart, relabel(1 "Day 1" 2 "Day 2" 3 "Day 3" 4 "Day 4" 5 "Day 5")) over(arm_value_0, relabel(1 "Nirmatrelvir/ritonavir" 2 "Control")) ytitle("Logarithm of the IPCW")
graph save "Graph" "ipcw_log.gph", replace
graph export "ipcw_log.png", as(png) name("Graph") replace

