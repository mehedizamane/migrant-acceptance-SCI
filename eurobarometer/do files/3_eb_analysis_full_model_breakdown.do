use eb_clean,clear

***breakdown by groups:


eststo clear

*============================================================*
* Manager
*============================================================*

mixed qb7_1_r ///
    c.z_contact##c.z_country_sci ///
    c.age i.female i.educ5 i.bills3 i.urban3 i.foreign_born ///
    c.log_gdp_pc c.log_pop c.z_internet ///
    [pweight=w1_norm] ///
    || country_id: ///

estat icc
estadd scalar ICC_country = r(icc2)

levelsof country_id if e(sample), local(countries)
local ncountry : word count `countries'
estadd scalar Countries = `ncountry'

eststo M1


*============================================================*
* Colleague
*============================================================*

mixed qb7_2_r ///
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

eststo M2


*============================================================*
* Neighbour
*============================================================*

mixed qb7_3_r ///
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

eststo M3


*============================================================*
* Doctor
*============================================================*

mixed qb7_4_r ///
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

eststo M4


*============================================================*
* Family
*============================================================*

mixed qb7_5_r ///
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

eststo M5


*============================================================*
* Friends
*============================================================*

mixed qb7_6_r ///
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

eststo M6




*============================================================*
* Export to Word (RTF)
*============================================================*

esttab M1 M2 M3 M4 M5 M6 using "eurob_mixed_breakdown.rtf", ///
    replace ///
    label ///
    b(3) se(3) ///
    star(* 0.10 ** 0.05 *** 0.01) ///
    mtitles("Manager" "Colleague" "Neighbour" "Doctor" "Family" "Friends") ///
    coeflabels( ///
        z_contact "Individual contact (z)" ///
        z_country_sci "Country SCI (z)" ///
        c.z_contact#c.z_country_sci "Contact × Country SCI" ///
        age "Age" ///
        1.female "Female" ///
        1.foreign_born "Foreign born" ///
        1.urban3 "Urban" ///
        log_gdp_pc "GDP per capita (log)" ///
        log_pop "Population (log)" ///
        z_internet "Internet use (z)" ///
    ) ///
    stats(N Countries ICC_country ll AIC BIC, ///
        labels("Observations" ///
        "Countries" ///
        "Country ICC" ///
        "Log likelihood" ///
        "AIC" ///
        "BIC"))
		