*This do-file creates the DGPs for y (event timing DGP coded elsewhere)
*Base treatment effect is $rho
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
	local events `events' (time>=event`i')*$rho +
	}
gen y`++ynum' = `events' fe_x + fe_t + epsilon
global y`ynum'lab "No dynamics, all events same size"

/*TREATMENT EFFECT INCREASES OVER TIME*/
*Create list of events
local events
forvalues i=1/500{
	capture confirm variable event`i'
	if _rc !=0 continue, break
	local j =`i'-1
	local events `events' /*(time>event`i')*$rho+*/cond(event`i',(time>=event`i')*($rho/10)*min(time-event`i'+1,10),0,0) +
	}
gen y`++ynum' = `events' fe_x + fe_t + epsilon
global y`ynum'lab "Treatment effect increases over time"

/*TREATMENT EFFECT DECREASES OVER TIME*/
*Create list of events
local events
forvalues i=1/500{
	capture confirm variable event`i'
	if _rc !=0 continue, break
	local j =`i'-1
	local events `events' /*(time>event`i')*$rho*/-cond(event`i',(time>=event`i')*($rho/10)*min(time-event`i'+1,10),0,0) +
	}
gen y`++ynum' = `events' fe_x + fe_t + epsilon
global y`ynum'lab "Treatment effect decreases over time"


/*TREATMENT EFFECT INCREASES, THEN DECREASES*/
*Create list of events
local events
forvalues i=1/500{
	capture confirm variable event`i'
	if _rc !=0 continue, break
	local j =`i'-1
	local events `events' cond(event`i',(time>=event`i')*($rho/10)*min(time-event`i'+1,10) -(time>event`i' +10)*($rho/15)*min(time-event`i'-10+1,10),0,0) +
	}
gen y`++ynum' = `events' fe_x + fe_t + epsilon
global y`ynum'lab "Treatment happens gradually, then decays"


/*NO DYNAMICS, SUBSEQUENT EVENTS HAVE NO EFFECT*/

gen y`++ynum' = (time>=event1)*$rho
global y`ynum'lab "No Dynamics, only first event has an effect"

/*
/*NO DYNAMICS, SUBSEQUENT EVENTS HAVE LESSER EFFECT*/
*Create list of events
local events
forvalues i=1/500{
	capture confirm variable event`i'
	if _rc !=0 continue, break
	local events `events' (time>event`i')*$rho/`i' +
	*Alternative code: local events `events' (time>event`i')*$rho/2 +
	}
gen y`++ynum' = `events' fe_x + fe_t +  epsilon
global y`ynum'lab "No dynamics, subsequent events have lesser effect"

/*NO DYNAMICS, SUBSEQUENT EVENTS HAVE GREATER EFFECT*/
*Create list of events
local events
forvalues i=1/500{
	capture confirm variable event`i'
	if _rc !=0 continue, break
	local events `events' (time>event`i')*$rho*`i' +
	*Alternative code: local events `events' (time>event`i')*$rho*2 +
	}
gen y`++ynum' = `events' fe_x + fe_t + epsilon
global y`ynum'lab "No dynamics, subsequent events have greater effect"

/*NO DYNAMICS, EFFECT OF SUBSEQUENT EVENTS DEPENDS ON PROXIMITY TO PREVIOUS EVENT*/
*Create list of events
local events (time>event1)*$rho +
forvalues i=2/500{
	capture confirm variable event`i'
	if _rc !=0 continue, break
	local j =`i'-1
	local eventdist=event`i'-event`j'
	local events `events' (time>event`i')*$rho*(`eventdist'/360) +
	}
gen y`++ynum' = `events' fe_x + fe_t + epsilon
global y`ynum'lab "No dynamics, effect of subsequent events depends on proximity to previous event"

*/
