clear all
macro drop _all




/* Multiple Event Study Monte Carlo: This is the basic program for
generating event study data and testing multiple methods for estimating
the treatment effect.

Our methods are:

1. Allow multiple event-time dummies to be turned on
2. Only consider the first event affecting each location
3. Duplicate observations so there is one observation per location per event

We allow for multiple data-generating processes for the outcome variable,
i.e. multiple different kinds of treatment effects and dynamics. These are:

o	No dynamics, all events same size
o	No dynamics, subsequent events have lesser effect
o	No dynamics, subsequent events have greater effect
o	No dynamics, effect of subsequent events depends on proximity to previous event



In this program, all specifications will have no maximum on the number of events, with
each event time drawn independently from a uniform distribution.

In this base version the event timing DGP is:
*/

global eventdesc "Any number of events, assigned by independent random draws"
global eventspec random_assignment_nomax


/* Set parameters*/

do parameters

/*Use this to set the gap between events*/

//local eventgap 36
local eventgap int(rand2)



/* Set the variance of the error term*/

local evar 30



tempfile temp



/*Now start the Monte Carlo Loop*/

forvalues iter = 1/$maxmc{
qui{ // Suppress all output
drop _all

set obs $obs

/*Generate a location identifier*/
gen location = ceil(_n/${maxtime})

/*And time*/
bysort location:  gen time = _n

/*Determine how many events*/

gen rand_events = runiform() if time == 1

by location: egen events = max(int(rand_events* ($maxtime-$panelsize*2)/6)+1)


summ events
global enum = r(max)

forvalues i = 1/$enum{
	gen rand`i' = runiform() if time == 1
	by location: egen event_r`i' = max($panelsize+int(rand`i'*(${maxtime}-$panelsize*2))) if events >=`i'
}

forvalues i = 1/$enum{
	by location: gen byte eventsize`i' = 1
}

/*Make sure we don't have any ties in event times*/

forvalues i = 1/$enum{
	local next = `i' + 1
	forvalues j = `next'/$enum{
		replace event_r`j' = event_r`j'+1 if event_r`j' == event_r`i'
	}
}


drop rand*

local missingtest 0
forvalues i = 1/$enum{
	local missingtest `missingtest' + (event_r`i'==240)
}

forvalues i = 1/$enum{
	egen event`i' = rowmin(event_r1-event_r$enum) if events >=`i'
	forvalues j = 1/$enum{
		replace event_r`j' = ${maxtime}*2 if event_r`j' == event`i' & event_r`j'!=. & `missingtest' <`i'
	}
}

drop event_r*	

	
/*Generate the error for the regression*/
gen epsilon = rnormal(0,`evar')

/*Generate Y*/

do y_dgp

/*Figure out how many y variables we have*/
global ynum
forvalues i = 1/500{
	capture confirm variable y`i'
	if _rc != 0 continue, break
	global ynum = `i'
}



do run_eventstudies

	

	/*Save the betas*/
	drop _all 
	local betas = $panelsize*2 + 1
	set obs `betas'

	gen t = _n-($panelsize + 1)
	forvalues y = 1/$ynum{
		foreach method in multiple ignore duplicate{	
			gen beta`method'`y' = 0
			forvalues i = 1/$panelsize{
				qui replace beta`method'`y' = ${lag`i'`method'`y'} if t == -`i'
				qui replace beta`method'`y' = ${lead`i'`method'`y'} if t== `i'
			}
			qui replace beta`method'`y' = ${event`method'`y'} if t == 0
		}
	}
	if `iter' > 1 append using `temp'
	save `temp', replace
	}
}


/*Plot the betas*/



collapse beta*, by(t)

local function1 y = 0 * (x<=0) + $rho*(x>0)
local function2 y = 0 * (x<=0) + ($rho/10)*x*(x>0&x<=10) + $rho*(x>10)
local function3 y = 0 * (x<=0) - ($rho/10)*x*(x>0&x<=10) - $rho*(x>10)
local function4 y = 0 * (x<=0) + ($rho/10)*x*(x>0&x<=10) + ($rho -($rho/15)*(x-10))*(x>10 & x<=20) + (1/3)*$rho*(x>20)
local function5 `function1'

local graphlist
forvalues y = 1/$ynum{

/*Estimate the slope of the average event dummies before an after
the event, and compare to what the slope should be*/
summ betamultiple`y' if t <0, meanonly
local offset = r(mean)

#delimit ;
graph twoway connect beta*`y' t , msymbol(S + T) lcolor(gs2 gs8 gs12) mcolor(gs2 gs8 gs12) ||
function `function`y'', range(-24 24) lpattern(dot) lcolor(gs0)||
function `function`y'' + `offset', range(-24 24) lpattern(dot) lcolor(gs0)||,
graphregion(color(white)) xsize(7) 
legend(label(1 Multiple Dummies On) label(2 Ignore Second Event) 
label(3 Duplicate Observations) col(1) order(1 2 3))
title("$eventdesc") 
subtitle("Outcome DGP: ${y`y'lab}")
;
#delimit cr

graph export graphs/y_compare_`y'.eps, replace
local graphlist `graphlist' graphs/y_compare_`y'.eps
}

instanttex using graphs/y_comparison.tex, fig(`graphlist') replace


