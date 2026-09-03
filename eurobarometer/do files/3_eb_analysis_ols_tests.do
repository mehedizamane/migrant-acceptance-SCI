use eb_clean,clear

*country level specification and tests 
reg qb7_comfort_index ///
    c.z_country_sci##c.z_contact ///
	    c.log_gdp_pc c.log_pop c.z_internet ///
    c.age i.female i.educ5 i.bills3 i.urban3 i.foreign_born ///
    [pweight=w1_norm], vce(cluster country_id)

boottest c.z_contact#c.z_country_sci, ///
    cluster(country_id) reps(9999) seed(12345)
boottest c.z_contact, ///
    cluster(country_id) reps(9999) seed(12345)
	
boottest c.z_country_sci, ///
    cluster(country_id) reps(9999) seed(12345)