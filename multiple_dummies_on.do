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

drop lagdiff1

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



/*Run the regression*/

regress y leaddiff* lagdiff* event $covariates `conditions'


	/*Store the betas*/
	global eventmultiple`y' = _b[event]

	forvalues i=1/$panelsize{
		if `i'>1 global lag`i'multiple`y' = _b[lagdiff`i']
		global lead`i'multiple`y' = _b[leaddiff`i']
	}

	

