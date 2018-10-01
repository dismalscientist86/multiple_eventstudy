/*Number of locations*/

global locations = 50

/*Panel length*/

global maxtime = 120


/* The event time parameters*/

global rho = 100

/*The size of our balanced panel*/

global panelsize = 24

/*The list of pre- and post- coefficients for a test command*/

forvalues i = 2/$panelsize{
	global lagtestlist ${lagtestlist} lagdiff`i'
}
forvalues i = 1/$panelsize{
	global leadtestlist ${leadtestlist} leaddiff`i'
}
global leadtestlist ${leadtestlist} event

/*Set conditions on the regression*/

global inwindow on

global conditions if abs(time-event1)<=$panelsize | abs(time-event2)<=$panelsize


set seed 11710


global obs = $maxtime*$locations


global maxmc = 500

/* Set covariates for regression */
global covariates i.location i.time
