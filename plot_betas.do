/*Plot the betas*/


collapse beta*, by(t)


local function1 y = 0 * (x<=0) + $rho*(x>0)
local function2 y = 0 * (x<=0) + ($rho/10)*x*(x>0&x<=10) + $rho*(x>10)
local function3 y = 0 * (x<=0) - ($rho/10)*x*(x>0&x<=10) - $rho*(x>10)
local function4 y = 0 * (x<=0) + ($rho/10)*x*(x>0&x<=10) + ($rho -($rho/15)*(x-10))*(x>10 & x<=20) + (1/3)*$rho*(x>20)

local graphlist
forvalues y = 1/$ynum{

/*Estimate the slope of the average event dummies before an after
the event, and compare to what the slope should be*/
summ betamultiple`y' if t <=0, meanonly
local offset = r(mean)

#delimit ;
graph twoway connect beta*`y' t , msymbol(S + T) ||
function `function`y'', range(-24 24) lpattern(dot) lcolor(gs0)||
function `function`y'' + `offset', range(-24 24) lpattern(dot) lcolor(gs0)||,
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


