import excel "C:\Users\zhenf\Desktop\data & code\Stata\Duration.xlsx", sheet("Sheet1")
rename A Duration
save Duration, clear

import excel "C:\Users\zhenf\Desktop\data & code\Stata\Dur_datetime.xlsx", sheet("Sheet1") allstring
rename A date
rename B time
save datetime, clear

use datetime, clear
merge 1:1 _n using Duration
drop _merge
save data_in_use

gen hour = substr(time,1,2)
gen minute = substr(time,4,2)
gen second = substr(time,7,2)

destring hour, replace
destring minute, replace
destring second, replace

replace minute = minute + 1 if second != 0 /// Eliminate the effect of leap second
drop second

/// Generate dummy variables based on the time nodes - 9:30, 10:00, 11:00, 13:00, 14:00, 14:30, 15:00
gen I1 = (hour == 9)
gen I2 = (hour == 10)
gen I3 = (hour == 11)
gen I4 = (hour == 13)
gen I5 = ((hour == 14) & (minute < 30))
gen I6 = ((hour == 14) & (minute >= 30))

/// Calculate time period between durations and time nodes
gen minute_30 = minute - 30
replace minute_30 = minute if ((hour != 9) & (hour != 14))
replace minute_30 = minute_30 + 30 if minute_30 < 0
rename minute_30 t_k
gen t_k_2 = t_k * t_k
gen t_k_3 = t_k * t_k * t_k

save regression_data

/// Cubic spline interpolation
forvalues i = 1/6{
	use regression_data, clear
	
	keep if I`i' == 1
	regress Duration t_k t_k_2 t_k_3
	matrix A`i' = e(b)
	
	clear
	svmat A`i'
	rename A`i'1 d1`i'
	rename A`i'2 d2`i'
	rename A`i'3 d3`i'
	rename A`i'4 c`i'
	
	save coefficient_`i'
}

/// Calculate diurnally adjusted durations
forvalues i = 1/6{
	use regression_data, clear
	
	keep if I`i' == 1
	gen phi= A`i'[1,4] + A`i'[1,1] * t_k + A`i'[1,2] * t_k_2 + A`i'[1,3] * t_k_3
	gen ad_duration = Duration / phi
	
	save ad_duration_`i'
}

use ad_duration_1, clear

forvalues i = 2/6{
	append using ad_duration_`i'
}

save ad_duration

gen year = substr(date,7,4)
gen month = substr(date,4,2)
gen day = substr(date,1,2)

destring year, replace
destring month, replace
destring day, replace

gsort year month day hour minute

save data4model
