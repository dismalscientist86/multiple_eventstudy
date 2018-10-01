clear all
macro drop _all




/* Multiple Event Study Monte Carlo: This is the basic program for
generating event study data and testing multiple methods for estimating
the treatment effect.

Our methods are:

1. Allow multiple event-time dummies to be turned on
2. Only consider the first event affecting each location
3. Duplicate observations so there is one observation per location per event

This program focuses on one data generating process for the outcome variable,
with no dynamics post-event.  However, we generate 7 distinct event generating
processes.  These are:

o	2 events, drawn independently from a uniform distribution
o	Any number of events, each drawn indepentently from a uniform
o	One event with certainty, a second with some probability (times still draw
	from a uniform)
o	2 events, first is drawn from a uniform, time to second is drawn from a lognormal
o	2 events, first is drawn from a uniform, time to second is drawn from a lognormal
	plus the length of the event window (that is, no overlapping event windows)
o	2 events, second is 6 periods after the first
o	2 events, second is 5-7 periods after the first.





In this base version the event timing DGP is:
*/


global eventspec random_assignment_nomax


/* Set parameters*/

do parameters

/*Use this to set the gap between events*/

//local eventgap 36
local eventgap int(rand2)



/* Set the variance of the error term*/

local evar 30



tempfile temp

local graphlist

/*Now start the Monte Carlo Loop*/
foreach event_dgp in random_assignment_nomax  {
	global eventspec `event_dgp'

forvalues iter = 1/$maxmc{
qui{ // Suppress all output
drop _all

set obs $obs

/*Generate a location identifier*/
gen location = ceil(_n/${maxtime})

/*And time*/
bysort location:  gen time = _n

/*Generate events*/

do event_code/${eventspec}

forvalues i = 1/$enum{
	by location: gen byte eventsize`i' = 1
}

compress
	
/*Generate the error for the regression*/
gen epsilon = rnormal(0,`evar')

/*Generate Y*/

*Base treatment effect is $rho
local ynum 0

/*Create fixed effects*/

bysort time: egen e_at_t = total(time==event1)
bysort location: gen events_sofar = sum(e_at_t)

gen timerandom =  rnormal(0,30) if location == 1
gen locationrandom = rnormal(0,30) if time == 1

bysort time: egen fe_t = max(timerandom)
bysort location: egen fe_x = max(locationrandom)

/*NO DYNAMICS, ALL EFFECTS SAME SIZE*/
*Create list of events
local events
forvalues i=1/500{
	capture confirm variable event`i'
	if _rc !=0 continue, break
	local events `events' (time>=event`i')*$rho +
	}
gen y`++ynum' = `events' fe_x + fe_t + epsilon
global ylab "No dynamics, all events same size"

global ynum = 1

/*Change covariate list*/

global covariates i.location

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


gen iter = ceil(_n/49)
gen spec = "${eventspec}"
if "${eventspec}"=="random_assignment_nomax" save MC_data/event_dgps, replace
else{
	append using MC_data/event_dgps
	save MC_data/event_dgps, replace
}

collapse beta*, by(t)




/*Estimate the slope of the average event dummies before an after
the event, and compare to what the slope should be*/
summ betamultiple`y' if t <0, meanonly
local offset = r(mean)

#delimit ;
graph twoway connect beta*`y' t , msymbol(S + T) lcolor(gs2 gs8 gs12) mcolor(gs2 gs8 gs12) ||
function y = 0 * (x<=0) + $rho*(x>0), range(-24 24) lpattern(dot) lcolor(gs0)||
function y = 0 * (x<=0) + $rho*(x>0) + `offset', range(-24 24) lpattern(dot) lcolor(gs0)||,
graphregion(color(white)) xsize(7) 
legend(label(1 Multiple Dummies On) label(2 Ignore Second Event) 
label(3 Duplicate Observations) col(1) order(1 2 3))
title("$eventdesc") 
subtitle("Outcome DGP: ${ylab}")
;
#delimit cr

graph export graphs/no_fe.pdf, replace




}

