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
o	Treatment effect increases over time
o	Treatment effect decreases over time
o	Treatment effect increases and then decreases


This is the basic version of the event study code, which draws event times independently
from a uniform distribution.

In this base version the event timing DGP is:
*/

global eventdesc "Any number of events, assigned by independent random draws"
global eventspec random_assignment_nomax


/* Set parameters*/

/*Number of locations*/

global locations = 2

/*Panel length*/

global maxtime = 20


/* The event time parameters*/

global rho = 100

/*The size of our balanced panel*/

global panelsize = 5

/*Set conditions on the regression*/

global inwindow on


set seed 11710


global obs = $maxtime*$locations



/* Set covariates for regression */
global covariates i.location i.time

/*Use this to set the gap between events*/

//local eventgap 36
local eventgap int(rand2)



/* Set the variance of the error term*/

local evar 30



tempfile temp



/*Now start the Monte Carlo Loop*/



drop _all

set obs $obs

/*Generate a location identifier*/
gen location = ceil(_n/${maxtime})

/*And time*/
bysort location:  gen time = _n

global eventdesc "Two events, assigned by independent random draws"


gen event1 = .
replace event1 = 7 if location == 1
replace event1 = 7 if location == 2
gen event2 = .
replace event2 = 15 if location == 1
replace event2 = 15 if location == 2
	
global enum = 2

/*Generate Y*/


/*NO DYNAMICS, ALL EFFECTS SAME SIZE*/
*Create list of events
local events
forvalues i=1/500{
	capture confirm variable event`i'
	if _rc !=0 continue, break
	local events `events' (time>=event`i')*$rho +
	}
gen y1 = (time>=event1)*$rho
global y`ynum'lab "No dynamics, all events same size"

/*Figure out how many y variables we have*/
global ynum = 1



/*Method 1: have multiple event dummies turned on*/



/*Generate event dummies*/
forvalues i = 1/$panelsize{
	gen byte leaddiff`i' = 0
	gen byte lagdiff`i' = 0
	forvalues e = 1/$enum {
		replace leaddiff`i' = 1 if (time - event`e'== `i')
		replace lagdiff`i' = 1 if (time - event`e' == -`i')
	}
}

/*With multiple event dummies turned on, we don't want an omitted category*/
gen byte event = 0
forvalues e = 1/$enum{
	replace event = 1 if time == event`e'
}

/*Set the end point dummies to stay on, and count up*/
sort location time

by location: replace leaddiff$panelsize = sum(leaddiff$panelsize)

gsort location - time

by location: replace lagdiff$panelsize = sum(lagdiff$panelsize)

sort location time

order location time event1 event2 lagdiff5 lagdiff4 lagdiff3 lagdiff2 lagdiff1 event leaddiff*


/*Run the regression*/

	regress y1 /*lagdiff5 lagdiff4*/ lagdiff* event leaddiff* if location == 1, nocons



/*Plot data to show why this is going wrong*/

gen eventtime = time - event1
gen eventtime2 = time - event2

local N =_N + $panelsize*2 + 1
set obs `N'

gen beta = .
replace eventtime = _N-_n - $panelsize if location == .

replace beta = 0 if location == . & eventtime == -1
replace beta = _b[event]  if location ==. & eventtime == 0
forvalues i = 1/$panelsize{
	if `i'> 1 replace beta = _b[lagdiff`i']  if location == . & eventtime == -`i'
	replace beta = _b[leaddiff`i'] if location == . & eventtime == `i'
}
regress
/*
#delimit ;
graph twoway 
connect y1 eventtime if location ==1 & abs(time - event1)<=$panelsize, jitter(1) ||
connect y1 eventtime2 if location ==1 & abs(time - event2)<=$panelsize, jitter(1) ||
connect y1 eventtime if location ==2 & abs(time - event1)<=$panelsize, jitter(1) ||
connect y1 eventtime2 if location ==2 & abs(time - event2)<=$panelsize, jitter(1) ||
connect beta eventtime if location == . ||
 if abs(eventtime)<=$panelsize,
legend(label(1 "Location 1, Event 1") label(2 "Location 1, Event 2")
label(3 "Location 2, Event 1") label(4 "Location 2, Event 2"))
graphregion(color(white)) name(eventtime, replace) xline(0)
title(Combined in Event Time) 
;

graph export graphs/duplicate_fail_eventtime.eps, replace;
graph twoway 
connect y1 time if location == 1, xline(5 7)
graphregion(color(white)) name(location1, replace)
title(Location 1)
;

graph twoway 
connect y1 time if location == 2, xline(10 14)
graphregion(color(white)) name(location2, replace)
title(Location 2)
;
#delimit cr


graph combine location1 location2, row(1) graphregion(color(white)) xsize(10)

graph export graphs/duplicate_fail_data.eps, replace

graph combine location1 location2 eventtime, col(1) graphregion(color(white)) ysize(12)

graph export graphs/duplicate_fail.eps, replace

instanttex using graphs/duplicate_failplot.tex, figure(graphs/duplicate_fail graphs/duplicate_fail_data graphs/duplicate_fail_eventtime) replace
