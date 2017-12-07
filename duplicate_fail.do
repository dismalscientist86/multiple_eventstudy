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
replace event1 = 5 if location == 1
replace event1 = 10 if location == 2
gen event2 = .
replace event2 = 7 if location == 1
replace event2 = 14 if location == 2
	
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
gen y1 = `events' 0
global y`ynum'lab "No dynamics, all events same size"

/*Figure out how many y variables we have*/
global ynum = 1



/*Method 3: Duplicate observations so we have one observation per event
per location per time*/

egen location_events = rownonmiss(event1-event$enum)

expand location_events, gen(event2flag)

sort location time event2flag
by location time: gen n = _n


forvalues e = 2/$enum{
	replace event1 = event`e' if event2flag & n == `e'
}

bysort location time: gen weight = 1/_N

/*Generate event dummies*/
forvalues i = 1/$panelsize{
	gen leaddiff`i' = (time - event1 == `i' )
	gen lagdiff`i' = (time - event1 == -`i' )
}
/*Set the end point dummies to stay on*/
replace leaddiff$panelsize = 1 if time - event1 > $panelsize & time - event1 !=.
replace lagdiff$panelsize = 1 if time - event1 < -$panelsize

/*With multiple event dummies turned on, we don't want an omitted category*/
gen event = (time ==event1)
drop lagdiff1

local conditions

if "$inwindow" == "on" local conditions if abs(time-event1)<=$panelsize

/*Run the regression*/

	regress y1 leaddiff* lagdiff* event  /*$covariates*/ [aw=weight] `conditions'
/*
	forvalues i=1/$panelsize{
		global lag`i'duplicate = _b[lagdiff`i']
		global lead`i'duplicate = _b[leaddiff`i']
	}
	global eventduplicate = _b[event]
	global lag1duplicate = 0
*/
/*Plot data to show why this is going wrong*/

gen eventtime = time - event1

local N =_N + $panelsize*2 + 1
set obs `N'

gen beta = .
replace eventtime = _N-_n - $panelsize if location == .

replace beta = _b[_cons] if location == . & eventtime == -1
replace beta = _b[event] + _b[_cons] if location ==. & eventtime == 0
forvalues i = 1/$panelsize{
	if `i'> 1 replace beta = _b[lagdiff`i'] + _b[_cons] if location == . & eventtime == -`i'
	replace beta = _b[leaddiff`i'] + _b[_cons] if location == . & eventtime == `i'
}

#delimit ;
graph twoway 
connect y1 eventtime if location ==1 & n == 1, jitter(1) ||
connect y1 eventtime if location ==1 & n == 2, jitter(1) ||
connect y1 eventtime if location ==2 & n == 1, jitter(1) ||
connect y1 eventtime if location ==2 & n == 2, jitter(1) ||
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
