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



tempfile temp temp2

local graphlist

/*Now start the Monte Carlo Loop*/
forvalues maxevent = 2/2{

forvalues iter = 1/$maxmc{
qui{ // Suppress all output
drop _all

set obs $obs

/*Generate a location identifier*/
gen location = ceil(_n/${maxtime})

/*And time*/
bysort location:  gen time = _n

/*Generate events*/

/*Determine how many events*/

gen rand_events = runiform() if time == 1

//by location: egen events = max(int(rand_events* ($maxtime-$panelsize*2)/6)+1)

gen events = `maxevent'

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
gen y1 = (time>=event1)*$rho + fe_x + fe_t + epsilon
gen y2 = (time>=event1)*($rho/10)*min(time-event1+1,10) + fe_x + fe_t + epsilon
gen y3 = -(time>=event1)*($rho/10)*min(time-event1+1,10) + fe_x + fe_t + epsilon
gen y4 = (time>=event1)*($rho/10)*min(time-event1+1,10) -(time>event1 +10)*($rho/15)*min(time-event1-10+1,10) + fe_x + fe_t + epsilon
global ylab "No dynamics, all events same size"

global ynum = 4

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






/*Plot the betas*/

/*
gen iter = ceil(_n/49)
gen spec = "${eventspec}"
if "${eventspec}"=="random_assignment_nomax" save MC_data/event_dgps, replace
else{
	append using MC_data/event_dgps
	save MC_data/
	event_dgps, replace
}
*/
collapse beta*, by(t)

if `maxevent'>2 joinby t using `temp2'
qui save `temp2', replace
}


collapse beta*, by(t)

local function1 y = 0 * (x<0) + $rho*(x>=0)
local function2 y = 0 * (x<0) + ($rho/10)*(x+1)*(x>=0&x<=9) + $rho*(x>9)
local function3 y = 0 * (x<0) - ($rho/10)*(x+1)*(x>=0&x<=9) - $rho*(x>9)
local function4 y = 0 * (x<0) + ($rho/10)*(x+1)*(x>=0&x<=9) + $rho*(x>9 & x<=10) + ($rho -($rho/15)*(x-9))*(x>10 & x<=19) + (1/3)*$rho*(x>19)
local function5 `function1'

global y1lab "No dynamics, all events same size"

global y2lab "Treatment effect increases over time"

global y3lab "Treatment effect decreases over time"

global y4lab "Treatment happens gradually, then decays"

global eventdesc "Two events, assigned by independent random draws"

local graphlist
forvalues y = 1/$ynum{

/*Estimate the slope of the average event dummies before an after
the event, and compare to what the slope should be*/
summ betamultiple`y' if t <=0, meanonly
local offset = r(mean)

#delimit ;
graph twoway connect beta*`y' t , msymbol(S + T) lcolor(gs2 gs8 gs12) mcolor(gs2 gs8 gs12) ||
function `function`y'', range(-24 24) lpattern(dot) lcolor(gs0)||
/*function `function`y'' + `offset', range(-24 24) lpattern(dot) lcolor(gs0)||*/,
graphregion(color(white)) xsize(7) 
legend(label(1 Multiple Dummies On) label(2 Ignore Second Event) 
label(3 Duplicate Observations) col(1) order(1 2 3))
title("$eventdesc") 
subtitle("Outcome DGP: ${y`y'lab}")
;
#delimit cr

graph export graphs/no_effect_`y'.eps, replace
local graphlist `graphlist' graphs/no_effect_`y'.eps
}




instanttex using graphs/no_effect.tex, fig(`graphlist') replace
