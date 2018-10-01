clear all

set obs 21

set seed 83112



gen beta = rnormal(40,20) in 11
replace beta = beta[_n-1] + rnormal(10,10) in 12/21
replace beta = rnormal(0,20) in 1/10

gen corr_t = _n-11
gen rand = runiform()
sort rand
gen rand_t = _n-11 
/*
replace beta = 0 if start_t == -1

sort beta


gen corr_t = start_t
gen rand_t = start_t

forvalues i = 1/4{

	local time = start_t in -`i'

	if `i' == 1 {
		replace corr_t = `time' if start_t == 2
		replace corr_t = 2 in -`i'
		replace rand_t = `time' if rand_t == 5
		replace rand_t = 5 in -`i'
	}
	else if `i' == 2{
		replace corr_t = `time' if start_t == 3
		replace corr_t = 3 in -`i'
		replace rand_t = `time' if rand_t == 9
		replace rand_t = 9 in -`i'
	}
	else if `i' == 3{
		replace corr_t = `time' if start_t == 1
		replace corr_t = 1 in -`i'
		replace rand_t = `time' if rand_t == 2
		replace rand_t = 2 in -`i'
	}
	else{
		replace corr_t = `time' if start_t == 0
		replace corr_t = 0 in -`i'
		replace rand_t = `time' if rand_t == -3 
		replace rand_t = -3 in -`i'
	}
}
*/
summ beta
local max = r(max)

gen ub = beta+.9*`max'
gen lb = beta-.9*`max'


sort corr_t

#delimit ;
graph twoway line lb ub corr_t, lwidth(vthin vthin) lcolor(gs10 gs10)||line beta corr_t || , graphregion(color(white)) xline(0) xsize(7) legend(off)
name(correlated, replace) xtitle(Event Time) ytitle(Beta) title(This should (Perhaps) be meaningful)
;

sort rand_t;

graph twoway line lb ub rand_t, lwidth(vthin vthin) lcolor(gs10 gs10)||line beta rand_t|| , graphregion(color(white)) xline(0) xsize(7) legend(off)
name(random, replace) xtitle(Event Time) ytitle(Beta) title(This definitely should not)
;



graph combine correlated random, graphregion(color(white)) col(1) xsize(7) ysize(8)
note("Top graph: Pre-event dummies drawn from a random normal, mean 0 and " 
"SD 20.  Time zero drawn from a normal(40,20), and post-event times are "
"equal to the previous period coefficient plus an increase drawn from a "
"normal(10,10). Confidence Interval equal to +/- .9 of the largest draw." " "
"Bottom Graph: Same points, same confidence interval, sorted randomly.");

graph export graphs/inference_puzzle.pdf, replace;

#delimit cr
