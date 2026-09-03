	
use eb_clean,clear

eststo clear


* Model 1: Contact only
mixed qb7_comfort_index ///
    c.z_contact ///
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


* Model 2: Country SCI only
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

eststo m2


* Model 3: Contact + Country SCI
mixed qb7_comfort_index ///
    c.z_contact ///
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

eststo m3


* Model 4: Contact × Country SCI
mixed qb7_comfort_index ///
    c.z_contact##c.z_country_sci ///
	    c.age i.female i.educ5 i.bills3 i.urban3 i.foreign_born ///
    c.log_gdp_pc c.log_pop c.z_internet ///
    [pweight=w1_norm] ///
    || country_id:

estat icc
estadd scalar ICC_country = r(icc2)

levelsof country_id if e(sample), local(countries)
local ncountry : word count `countries'
estadd scalar Countries = `ncountry'

eststo m4

esttab m1 m2 m3 m4 using "eurob_mix_effect_full_model.rtf", replace ///
    label ///
    mtitles("Contact" "Country SCI" "Contact + SCI" "SCI × Contact") ///
    b(3) se(3) ///
    star(* 0.10 ** 0.05 *** 0.01) ///
    stats(N Countries ICC_country, ///
        labels("N" "Countries" "Country ICC") ///
        fmt(0 0 3)) ///
    title("Migrant comfort index") ///
    compress
