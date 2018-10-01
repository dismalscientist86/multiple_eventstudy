clear all
macro drop _all




/* Multiple Event Study Monte Carlo: This is the basic program for
generating event study data and testing multiple methods for estimating
the treatment effect.

Our methods are:

1. Allow multiple event-time dummies to be turned on
2. Only consider the first event affecting each location
3. Duplicate observations so there is one observation per location per event

This program tests the robustness of the methods to failing to control for unobserved
trends, time effects, individual effects etc.


This is the basic version of the event study code, which draws event times independently
from a uniform distribution.

In this base version the event timing DGP is:
*/

global eventdesc "Any number of events, assigned by independent random draws"
global eventspec robust_unobservables


/* Set parameters*/

do parameters

/*Use this to set the gap between events*/

//local eventgap 36
local eventgap int(rand2)



/* Set the variance of the error term*/

local evar 30



tempfile temp


/*We will do this twice: once with fixed effects included, once without*/


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

by location: egen events = max(int(rand_events* ($maxtime-$panelsize*2)/6))

summ events
global enum = r(max)

forvalues i = 1/$enum{
	gen rand`i' = runiform() if time == 1
	by location: egen event_r`i' = max($panelsize+int(rand`i'*(${maxtime}-$panelsize*2))) if events >=`i'
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

*Base treatment effect is $rho
local ynum 0

/*Create fixed effects and trends*/

xtset location time

gen timerandom =  rnormal(sqrt(event1),30) if location == 1
gen locationrandom = rnormal(events^2,30) if time == 1

/*Create a complete polynomial trend*/
gen trend = (time/16-6.5)^2*(2*(time)^2-7*time-700)*(2*time+40)*.00003

bysort time: egen fe_t = max(timerandom)
bysort location: egen fe_x = max(locationrandom)


*Create list of events
local events
forvalues i=1/500{
	capture confirm variable event`i'
	if _rc !=0 continue, break
	local events `events' (time>event`i')*$rho +
	}

/*Generate Y variables with various combinations of unobservables*/
	
gen y`++ynum' = `events'  epsilon
global y`ynum'lab "No Unobservables"

gen y`++ynum' = `events' fe_x + epsilon
global y`ynum'lab "Locations with many events have higher outcomes"

gen y`++ynum' = `events' trend + epsilon
global y`ynum'lab "Complex Time Trend"

gen y`++ynum' = `events' fe_x + trend + epsilon
global y`ynum'lab "Trend and Individual Effects"



STOP 
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
	
		}
		qui replace betamultiple`y' = ${event`y'} if t == 0
	}
	if `iter' > 1 append using `temp'
	save `temp', replace
	}
}


/*Plot the betas*/


collapse beta*, by(t)


local function y = 0 * (x<=0) + $rho*(x>0)

local graphlist
forvalues y = 1/$ynum{

/*Estimate the slope of the average event dummies before an after
the event, and compare to what the slope should be*/
summ betamultiple`y' if t <=0, meanonly
local offset = r(mean)

#delimit ;
graph twoway connect beta*`y' t ||
function `function', range(-24 24) lpattern(dot) lcolor(gs0)||
function `function' + `offset', range(-24 24) lpattern(dot) lcolor(gs0)||,
graphregion(color(white)) xsize(7) 
legend(label(1 Multiple Dummies On) label(2 Ignore Second Event) 
label(3 Duplicate Observations) col(1) order(1 2 3))
title("$eventdesc") 
subtitle("Outcome DGP: ${y`y'lab}")
;
#delimit cr

graph export graphs/${eventspec}_`y'.eps, replace
local graphlist `graphlist' graphs/${eventspec}_`y'.eps
}



instanttex using graphs/${eventspec}.tex, fig(`graphlist') replace

