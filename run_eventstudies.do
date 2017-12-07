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

gen event = (time == event1)
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
	global eventignore`y' = _b[event]
	global lag1ignore`y' = 0
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

/*With multiple event dummies turned on, we don't want an omitted category*/
gen event = (time ==event1)
drop lagdiff1

local conditions

if "$inwindow" == "on" local conditions if abs(time-event1)<=$panelsize

/*Run the regression*/
forvalues y = 1/$ynum{
	regress y`y' leaddiff* lagdiff* event $covariates [aw=weight] `conditions'

	forvalues i=1/$panelsize{
		if `i' > 1 global lag`i'duplicate`y' = _b[lagdiff`i']
		global lead`i'duplicate`y' = _b[leaddiff`i']
	}
	global eventduplicate`y' = _b[event]
	global lag1duplicate`y' = 0
}

/*Store the betas*/
//global event = _b[event]

