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
local graphlist

forvalues censoring = -1/1{
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

by location: egen events = max(int(rand_events* ($maxtime+2*$panelsize*`censoring')/6)+1)



summ events
global enum = r(max)

forvalues i = 1/$enum{
	gen rand`i' = runiform() if time == 1
	by location: egen event_r`i' = max($panelsize*(-`censoring')+int(rand`i'*(${maxtime}+2*$panelsize*`censoring'))) if events >=`i'

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
	gen rand`i' = runiform() if time == 1
	by location: egen eventsize`i' = max(rand`i') if events>=`i'
}
drop rand*

	
/*Generate the error for the regression*/
gen epsilon = rnormal(0,`evar')

/*Generate Y*/

local ynum 0

/*Create fixed effects*/

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
	local events `events' (time>=event`i')*$rho*cond(missing(eventsize`i'),0,eventsize`i') +
	}
gen y`++ynum' = `events' fe_x + fe_t + epsilon
global y`ynum'lab "No dynamics, all events same size"

/*Figure out how many y variables we have*/
global ynum
forvalues i = 1/500{
	capture confirm variable y`i'
	if _rc != 0 continue, break
	global ynum = `i'
}


/*Method 1: have multiple event dummies turned on*/



/*Generate event dummies*/
forvalues i = 1/$panelsize{
	gen byte leaddiff`i' = 0
	gen byte lagdiff`i' = 0
	forvalues e = 1/$enum {
		replace leaddiff`i' = eventsize`e' if (time - event`e'== `i')
		replace lagdiff`i' = eventsize`e' if (time - event`e' == -`i')
	}
}

/*With multiple event dummies turned on, we don't want an omitted category*/
gen byte event = 0
forvalues e = 1/$enum{
	replace event = eventsize`e' if time == event`e'
}

/*Set the end point dummies to stay on, and count up*/
sort location time

by location: replace leaddiff$panelsize = sum(leaddiff$panelsize)

gsort location - time

by location: replace lagdiff$panelsize = sum(lagdiff$panelsize)

sort location time

/*
/*Set the end point dummies to stay on*/

replace leaddiff$panelsize = 1 if time - min(event1,event2) > $panelsize
replace lagdiff$panelsize = 1 if time - max(event2,event1) < -$panelsize

replace leaddiff$panelsize = 2 if time - max(event2,event1) > $panelsize
replace lagdiff$panelsize = 2 if time - min(event1,event2) < -$panelsize

*/


/*Construct a macro for the conditions on the regression*/

if "$inwindow" == "on" {
	local conditions if 1 == 0
	forvalues e = 1/$enum{
		local conditions `conditions' | abs(time-event`e')<=$panelsize
	}
}

/*Run the regression*/
forvalues y = 1/$ynum{
	regress y`y' leaddiff* lagdiff* event $covariates `conditions'


	/*Store the betas*/
	global eventmultiple`y' = _b[event]
	forvalues i=1/$panelsize{
		global lag`i'multiple`y' = _b[lagdiff`i']
		global lead`i'multiple`y' = _b[leaddiff`i']
	}
}

drop event lagdiff* leaddiff*

/*Method 2: Ignore succeeding events*/

/*Generate event dummies*/
forvalues i = 1/$panelsize{
	gen leaddiff`i' = (time - event1 == `i' )*eventsize1
	gen lagdiff`i' = (time - event1 == -`i' )*eventsize1
}
/*Set the end point dummies to stay on*/
replace leaddiff$panelsize = eventsize1 if time - event1 > $panelsize & time - event1 !=.
replace lagdiff$panelsize = eventsize1 if time - event1 < -$panelsize

gen event = time == event1
drop lagdiff1

local conditions

if "$inwindow" == "on" local conditions if abs(time-event1)<=$panelsize

/*Run the regression*/
forvalues y = 1/$ynum{
	regress y`y' leaddiff* lagdiff* event $covariates  `conditions'


	forvalues i=1/$panelsize{
		if `i' > 1 global lag`i'ignore`y' = _b[lagdiff`i']
		global lead`i'ignore`y' = _b[leaddiff`i']
	}
	global lag1ignore`y'= 0
	global eventignore`y' = _b[event]
}


/*Store the betas*/
//global event = _b[event]

drop lagdiff* leaddiff* event

/*Method 3: Duplicate observations so we have one observation per event
per location per time*/

egen location_events = rownonmiss(event1-event$enum)

expand location_events, gen(event2flag)

sort location time event2flag
by location time: gen n = _n


forvalues e = 2/$enum{
	replace event1 = event`e' if event2flag & n == `e'
	replace eventsize1 = eventsize`e' if event2flag & n == `e'
}

bysort location time: gen weight = 1/_N

/*Generate event dummies*/
forvalues i = 1/$panelsize{
	gen leaddiff`i' = (time - event1 == `i' )*eventsize1
	gen lagdiff`i' = (time - event1 == -`i' )*eventsize1
}
/*Set the end point dummies to stay on*/
replace leaddiff$panelsize = eventsize1 if time - event1 > $panelsize & time - event1 !=.
replace lagdiff$panelsize = eventsize1 if time - event1 < -$panelsize

/*With different sizes of event, we don't want an omitted category*/
gen event = (time ==event1)*eventsize1

local conditions

if "$inwindow" == "on" local conditions if abs(time-event1)<=$panelsize

/*Run the regression*/
forvalues y = 1/$ynum{
	regress y`y' leaddiff* lagdiff* event $covariates [aw=weight] `conditions'

	forvalues i=1/$panelsize{
		global lag`i'duplicate`y' = _b[lagdiff`i']
		global lead`i'duplicate`y' = _b[leaddiff`i']
	}
	global eventduplicate`y' = _b[event]
}

/*Store the betas*/



	

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
	if abs(${eventmultiple1}) >4000 replace betamultiple = betamultiple - ${eventmultiple1} 
	if `iter' > 1 append using `temp'
	save `temp', replace
	}
}


/*Plot the betas*/



collapse beta*, by(t)

local function1 y = 0 * (x<=0) + $rho*(x>0)




/*Estimate the slope of the average event dummies before an after
the event, and compare to what the slope should be*/
summ betamultiple`y' if t <0, meanonly
local offset = r(mean)

#delimit ;
graph twoway connect beta*`y' t , msymbol(S + T) lcolor(gs2 gs8 gs12) mcolor(gs2 gs8 gs12) ||
function `function1', range(-24 24) lpattern(dot) lcolor(gs0)||
function `function1' + `offset', range(-24 24) lpattern(dot) lcolor(gs0)||,
graphregion(color(white)) xsize(7) 
legend(label(1 Multiple Dummies On) label(2 Ignore Second Event) 
label(3 Duplicate Observations) col(1) order(1 2 3))
title("$eventdesc") 
subtitle("Outcome DGP: No Dynamics, Events Different Sizes")
;
#delimit cr

graph export graphs/censoring_`censoring'.eps, replace
local graphlist `graphlist' graphs/censoring_`censoring'.eps
}

instanttex using graphs/censoring.tex, fig(`graphlist') replace


