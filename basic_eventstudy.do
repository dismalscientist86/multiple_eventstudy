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



gen events = 1

summ events
global enum = r(max)

	gen rand = runiform() if time == 1
	by location: egen event = max($panelsize+int(rand*(${maxtime}-$panelsize*2))) 




drop rand*

gen timerandom =  rnormal(0,30) if location == 1
gen locationrandom = rnormal(0,30) if time == 1

bysort time: egen fe_t = max(timerandom)
bysort location: egen fe_x = max(locationrandom)

gen locationtrendrandom = rnormal(0,2) if time == 1
bysort location: egen trendslope = max(locationtrendrandom)

gen locationtrend = 0 //trendslope*time	

	
/*Generate the error for the regression*/
gen epsilon = rnormal(0,`evar')

/*Generate Y*/

gen y= (time>=event)*${rho} + fe_t + fe_x + locationtrend + epsilon

/*Figure out how many y variables we have*/

forvalues i = 1/$panelsize{
	gen lagdiff`i' = time -event == -`i'
	gen leaddiff`i' = time - event == `i'
}
replace leaddiff$panelsize = 1 if time-event>$panelsize
replace lagdiff$panelsize = 1 if time -event<-$panelsize

gen eventdummy = time == event

drop lagdiff1

regress y lagdiff* leaddiff* eventdummy i.time c.time#i.location i.location 


	/*Save the betas*/
	drop _all 
	local betas = $panelsize*2 + 1
	set obs `betas'

	gen t = _n-($panelsize + 1)
	gen beta = 0
	forvalues i = 1/$panelsize{
		if `i'>1 replace beta = _b[lagdiff`i'] if t==-`i'
		replace beta = _b[leaddiff`i'] if t==`i'
	}
	replace beta = _b[eventdummy] if t ==0

		gen iter = `iter'

	if `iter' > 1 append using `temp'
	save `temp', replace
	}
}

#delimit ;
twoway 
connect beta t if iter ==1 ||
connect beta t if iter ==2 ||
connect beta t if iter ==3 ||
connect beta t if iter ==4 ||
connect beta t if iter ==5 ||
;
/*Plot the betas*/


/*
collapse beta*, by(t)
twoway connect beta t

/*
local function1 y = 0 * (x<=0) + $rho*(x>0)

local graphlist
forvalues y = 1/$ynum{

/*Estimate the slope of the average event dummies before an after
the event, and compare to what the slope should be*/
summ betamultiple`y' if t <=0, meanonly
local offset = r(mean)

#delimit ;
graph twoway connect beta*`y' t , msymbol(S + T) ||
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


