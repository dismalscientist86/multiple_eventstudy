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

gen events = 1

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
sort location time
by location: gen events_sofar = sum(e_at_t)

gen timerandom =  rnormal(0,30) if location == 1
gen locationrandom = rnormal(0,30) if time == 1

bysort time: egen fe_t = max(timerandom)
bysort location: egen fe_x = max(locationrandom)


gen y1 = /*(time>=event1)*$rho +*/ fe_x + fe_t + epsilon


global ynum = 1

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
//if `iter'== 148 STOP
/*Run the regression*/

	regress y leaddiff* lagdiff*   $covariates /*`conditions'*/, vce(cluster location)
	

	/*Store the betas*/
	global eventreject = 0 //(abs(_b[event]/_se[event])>=1.96)
	global eventbeta = 0
	forvalues i=1/$panelsize{
		global lag`i'reject = (abs(_b[lagdiff`i']/_se[lagdiff`i'])>=1.96)
		global lead`i'reject = (abs(_b[leaddiff`i']/_se[leaddiff`i'])>=1.96)
		global lag`i'beta = _b[lagdiff`i']
		global lead`i'beta = _b[leaddiff`i']
	}


	

	/*Save the betas*/
	drop _all 
	local betas = $panelsize*2 + 1
	set obs `betas'

	gen t = _n-($panelsize + 1)

		foreach var in beta reject {	
			gen `var' = 0
			forvalues i = 1/$panelsize{
				qui replace `var'= ${lag`i'`var'} if t == -`i'
				qui replace `var' = ${lead`i'`var'} if t== `i'
			}
			qui replace `var' = ${event`var'} if t == 0
		}
	
	egen total = total(reject)
	gen iter = `iter'
	if `iter' > 1 append using `temp'
	save `temp', replace
	}
}


/*Plot the betas*/

local plots
forvalues i = 1/$maxmc{
	local plots `plots' line beta t if iter == `i' ||
}


bysort t: egen rank = rank(beta)

sort iter t
by iter: egen lagrank = max(rank*(t==-24))

#delimit ;
graph twoway line beta t if lagrank == round(${maxmc}*.05)|| line beta t if lagrank ==round(${maxmc}*.95)|| line beta t if lagrank == round(${maxmc}*0.5)||, graphregion(color(white)) 
legend(label(1 "5th percentile") label(2 "95th percentile") label(3 "Median") row(1) title("Percentile of ${panelsize}-period lag dummes"))
;
#delimit cr


//graph twoway `plots', yline(0) graphregion(color(white))




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
collapse beta reject total, by(t)


graph twoway line reject t, yline(0.05) graphregion(color(white))

