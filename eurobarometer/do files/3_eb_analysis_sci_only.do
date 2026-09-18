use eb_clean,clear

eststo clear

* Model 1: Country SCI only
mixed qb7_comfort_index ///
    c.z_country_sci ///
	c.age i.female i.educ5 i.bills3 i.urban3 i.foreign_born ///
    c.log_gdp_pc c.log_pop c.z_internet ///
    [pweight=w1_norm] ///
    || country_id:

estat icc
estadd scalar ICC_country = r(icc2)

levelsof country_id if e(sample), local(countries)
local ncountry : word count `countries'
estadd scalar Countries = `ncountry'

eststo m1


*============================================================*
* Manager
*============================================================*


mixed qb7_1_r ///
    c.z_country_sci ///
    c.age i.female i.educ5 i.bills3 i.urban3 i.foreign_born ///
    c.log_gdp_pc c.log_pop c.z_internet ///
    [pweight=w1_norm] ///
    || country_id: ///

estat icc
estadd scalar ICC_country = r(icc2)

levelsof country_id if e(sample), local(countries)
local ncountry : word count `countries'
estadd scalar Countries = `ncountry'

local p_manager = 2*normal(-abs(_b[z_country_sci]/_se[z_country_sci]))

eststo m2


*============================================================*
* Colleague
*============================================================*

mixed qb7_2_r ///
    c.z_country_sci ///
    c.age i.female i.educ5 i.bills3 i.urban3 i.foreign_born ///
    c.log_gdp_pc c.log_pop c.z_internet ///
    [pweight=w1_norm] ///
    || country_id:

estat icc
estadd scalar ICC_country = r(icc2)

levelsof country_id if e(sample), local(countries)
local ncountry : word count `countries'
estadd scalar Countries = `ncountry'
local p_colleague = 2*normal(-abs(_b[z_country_sci]/_se[z_country_sci]))

eststo m3


*============================================================*
* Neighbour
*============================================================*

mixed qb7_3_r ///
    c.z_country_sci ///
    c.age i.female i.educ5 i.bills3 i.urban3 i.foreign_born ///
    c.log_gdp_pc c.log_pop c.z_internet ///
    [pweight=w1_norm] ///
    || country_id:

estat icc
estadd scalar ICC_country = r(icc2)

levelsof country_id if e(sample), local(countries)
local ncountry : word count `countries'
estadd scalar Countries = `ncountry'
local p_neighbour = 2*normal(-abs(_b[z_country_sci]/_se[z_country_sci]))

eststo m4


*============================================================*
* Doctor
*============================================================*

mixed qb7_4_r ///
    c.z_country_sci ///
    c.age i.female i.educ5 i.bills3 i.urban3 i.foreign_born ///
    c.log_gdp_pc c.log_pop c.z_internet ///
    [pweight=w1_norm] ///
    || country_id:

estat icc
estadd scalar ICC_country = r(icc2)

levelsof country_id if e(sample), local(countries)
local ncountry : word count `countries'
estadd scalar Countries = `ncountry'
local p_doctor = 2*normal(-abs(_b[z_country_sci]/_se[z_country_sci]))

eststo m5


*============================================================*
* Family
*============================================================*

mixed qb7_5_r ///
    c.z_country_sci ///
    c.age i.female i.educ5 i.bills3 i.urban3 i.foreign_born ///
    c.log_gdp_pc c.log_pop c.z_internet ///
    [pweight=w1_norm] ///
    || country_id:

estat icc
estadd scalar ICC_country = r(icc2)

levelsof country_id if e(sample), local(countries)
local ncountry : word count `countries'
estadd scalar Countries = `ncountry'
local p_family = 2*normal(-abs(_b[z_country_sci]/_se[z_country_sci]))

eststo m6


*============================================================*
* Friends
*============================================================*

mixed qb7_6_r ///
    c.z_country_sci ///
    c.age i.female i.educ5 i.bills3 i.urban3 i.foreign_born ///
    c.log_gdp_pc c.log_pop c.z_internet ///
    [pweight=w1_norm] ///
    || country_id:

estat icc
estadd scalar ICC_country = r(icc2)

levelsof country_id if e(sample), local(countries)
local ncountry : word count `countries'
estadd scalar Countries = `ncountry'
local p_friends = 2*normal(-abs(_b[z_country_sci]/_se[z_country_sci]))

eststo m7


*============================================================*
* Holm adjustment
*============================================================*

preserve

clear
set obs 6

gen outcome = ""
gen p_unadjusted = .
gen order = _n

replace outcome = "Manager"   in 1
replace outcome = "Colleague" in 2
replace outcome = "Neighbour" in 3
replace outcome = "Doctor"    in 4
replace outcome = "Family"    in 5
replace outcome = "Friends"   in 6

replace p_unadjusted = `p_manager'   in 1
replace p_unadjusted = `p_colleague' in 2
replace p_unadjusted = `p_neighbour' in 3
replace p_unadjusted = `p_doctor'    in 4
replace p_unadjusted = `p_family'    in 5
replace p_unadjusted = `p_friends'   in 6

* Sort p-values from smallest to largest
sort p_unadjusted

gen rank = _n
gen p_holm = p_unadjusted * (6 - rank + 1)

* Holm adjusted p-values must be monotonic
forvalues i = 2/6 {
    replace p_holm = max(p_holm, p_holm[_n-1]) if _n == `i'
}

* Cap at 1
replace p_holm = min(p_holm, 1)

sort outcome

list outcome p_unadjusted p_holm, noobs sep(0)

restore

*------------------------------------------------------------
* Export table
*------------------------------------------------------------
esttab m1 m2 m3 m4 m5 m6 m7  using "eurob_mixed_effect_sci_only.rtf", replace ///
    label ///
    mtitles("pooled" "Manager" "Colleague" "Neighbour" "Doctor" "Family" "Friends") ///
    b(3) se(3) ///
    star(* 0.10 ** 0.05 *** 0.01) ///
	    coeflabels( ///
        z_country_sci "Global SCI" ///
        age "Age" ///
        1.female "Female" ///
        1.foreign_born "Foreign born" ///
        1.urban3 "Urban" ///
        log_gdp_pc "GDP per capita (log)" ///
        log_pop "Population (log)" ///
        z_internet "Internet use (z)" ///
    ) ///
    stats(N Countries ICC_country, ///
        labels("N" "Countries" "Country ICC") ///
        fmt(0 0 3)) ///
    title("Migrant comfort index") ///
    compress
	
