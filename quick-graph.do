capture drop n
capture drop beta

gen n = _n-$panelsize-1 if _n-$panelsize-1 <= $panelsize
gen beta = .
forvalues i = 1/$panelsize{
qui{
if `i'> 1 replace beta = _b[lagdiff`i'] if n ==-`i'
replace beta = _b[leaddiff`i'] if n == `i'
}
}
replace beta = _b[event] if n == 0

graph twoway connect beta n, xlabel(-$panelsize(12)$panelsize)
