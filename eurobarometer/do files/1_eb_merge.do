clear 

local eu27_efta_uk ///
 AT BE BG HR CY CZ DK EE FI FR DE EL GR HU IE IT LV LT LU MT NL PL PT RO SK SI ES SE ///
 CH NO IS LI GB UK




use "ZA7847_v1-0-0", clear

* Country code corrections
gen iso2 = isocntry
replace iso2 = "DE" if inlist(isocntry,"DE-W","DE-E")
replace iso2 = "GR" if isocntry == "EL"



*merge SCI 
merge m:1 iso2  using "global_sci_final_neil_2020.dta"

*drop countries not in eurobarometer 
drop if _merge==2 


save eb_sci,replace 

*merge WDI 

import delimited "wdi_indicators.csv", clear

keep if time=="2021"

rename populationtotalsppoptotl pop

gen gdp_pc = real(gdppercapitaconstantlcunygdppcap)
gen internet_use = real(individualsusingtheinternetofpop)



kountry countrycode, from(iso3c) to(iso2c)

list countryname countrycode  _ISO2C_ in 1/20


rename  _ISO2C_ iso2

keep iso2 gdp_pc internet_use pop

save "pop_gdp_2021.dta", replace


use eb_sci,clear
drop _merge
merge m:m iso2 using "pop_gdp_2021.dta"

*drop countries not in eurobarometer 
drop if _merge==2 

drop _merge


save eb_merged,replace

