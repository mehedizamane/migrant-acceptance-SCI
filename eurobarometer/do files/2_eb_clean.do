use eb_merged,clear


* Comfort level score
foreach i of numlist 1/6 {
    gen qb7_`i'_r = 5 - qb7_`i' if inrange(qb7_`i',1,4)   // reverse: 4->1 ... 1->4
}


*N= 3179 dropped 
drop if missing(qb7_1_r, qb7_2_r, qb7_3_r, qb7_4_r, qb7_5_r, qb7_6_r)


egen qb7_comfort_index = rowmean(qb7_1_r qb7_2_r qb7_3_r qb7_4_r qb7_5_r qb7_6_r)
label var qb7_comfort_index "Migrant Acceptance"


*foreign born, N=31 missing 
gen native_born = (qb18_1 == 1) if inrange(qb18_1,1,5)

gen foreign_born = 1 - native_born



*education, N=462 missing 
gen educ5 = d8r2 if inrange(d8r2,1,5)
label values educ5 d8r2

label define educ5_lbl ///
    1 "15 or younger" ///
    2 "16-19" ///
    3 "20+" ///
    4 "Still studying" ///
    5 "No full-time education"

label values educ5 educ5_lbl

*age drop, N=3 for below 15

gen age = d11 if d11 < 900
drop if age<15
drop if age>100

*sex
gen female = (d10 == 2) if inlist(d10,1,2,3)

*bill/income, N=142 missing 
gen bills3 = d60 if inrange(d60,1,3)
label values bills3 d60

*urban, N=8 missing 
gen urban3 = d25 if inrange(d25,1,3)
label values urban3 d25

encode iso2, gen(country_id)





*weight 

bysort country_id: egen sum_w1 = total(w1)
egen N = count(w1)
count
local n = r(N)

levelsof country_id, local(countries)
local K : word count `countries'


gen w1_norm = (w1 / sum_w1) * (`n' / `K')



* QB6 contact frequency -> reverse so higher = more contact ----
gen qb6_contact = .
replace qb6_contact = 5 if qb6 == 1   // daily
replace qb6_contact = 4 if qb6 == 2   // weekly+
replace qb6_contact = 3 if qb6 == 3   // monthly+
replace qb6_contact = 2 if qb6 == 4   // yearly+
replace qb6_contact = 1 if qb6 == 5   // less often/never
label define contact 1 "Less often/never" 2 "At least yearly" 3 "At least monthly" 4 "At least weekly" 5 "Daily"
label values qb6_contact contact
label var qb6_contact "Frequency of social interaction with immigrants (reverse-coded, higher=more contact)"
egen z_contact = std(qb6_contact)


gen log_gdp_pc = log(gdp_pc)
gen log_pop = log(pop)
egen z_pop = std(pop)
egen z_internet = std(internet_use)
egen z_country_sci = std(global_sci_2020_mean_log_cross_b)


*N=627 dropped for not having complete set of individual controls 
drop if missing(foreign_born, educ5, female, age, bills3, urban3, log_gdp_pc, log_pop, z_internet)


save eb_clean,replace


