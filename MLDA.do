/*This program sets up and creates an event study using Dobkin's MLDA data and the three different menthods
clear 
set more off
macro drop _all
set matsize 800

use extract_d_sandler.dta, clear
gen month=mofd(date)
gen monthofyear=month(date)
destring fips, replace
xtset fips month

* ******************************************************************** 
*Use the month of the law change, even when older individuals are grandfathered in
gen change_beer=D.MnAgeGfB
gen event=(D.MnAgeGfB!=0) if L.change_beer==0
replace event=0 if event==.
replace event=0 if month <180 

/*Count the number of events that have happened so far*/
bysort fips: gen eventcount=sum(event)

/*Create variables that capture the month of event*/
forvalues i=1/4{
	by fips : egen event`i' = min(month + 1000000 * (eventcount != `i'))
	replace event`i'=. if event`i'>30000
	}
forvalues i=1/4{
	by fips : egen eventsize`i' = total(change_beer* (eventcount == `i'))
	//by fips: egen eventsizedown`i' = min(change_beer + 100000*(eventcount!=`i'))
	//replace eventsize`i' = eventsizedown`i' if eventsizedown`i' <0 & eventsize`i' ==0
}
/*
by fips: egen bigevent = max(abs(change_beer))
by fips: egen biggestevent = min(month + 1000000*(abs(change_beer) == bigevent))
*/

egen bigeventsize = rowmax(eventsize*)
gen bigevent = .
forvalues i = 4(-1)1{
	replace bigevent = event`i' if eventsize`i' == bigeventsize
}

gen alt1 = .
gen altsize1 =.
replace alt1 = event1 if eventsize1!=bigeventsize
replace altsize1 = eventsize1 if eventsize1 !=bigeventsize
replace event1 = bigevent if event1 !=bigevent
forvalues i = 2/4{
	replace eventsize`i' = altsize1 if event`i' == event1
	replace event`i' = alt1 if event`i' == event1
}
	

save MLDA.dta, replace*/
 use MLDA, replace
 set more off
 gen y1=MVA_fars_evening_r_18_20
 global ynum=1
 global enum=4
 global panelsize=12
 gen time = month
 drop event
gen location = fips
 global eventspec MLDA
 global eventdesc "Events Are Changes in the Minimum Legal Drinking Age"
 global y1lab "Empirical"
 global covariates i.location i.time 
 global inwindow on
 
 xtset location time
 
 do run_eventstudies_agg.do
 
 /*Save the betas*/
	drop _all 
	local betas = $panelsize + 1
	set obs `betas'

	gen t = (_n-($panelsize/2 + 1))*2
	forvalues y = 1/$ynum{
		foreach method in multiple ignore duplicate{	
			gen beta`method'`y' = 0
			forvalues i = 2(2)$panelsize{
				qui replace beta`method'`y' = ${lag`i'`method'`y'} if t == -`i'
				qui replace beta`method'`y' = ${lead`i'`method'`y'} if t== `i'
			}
	
		}
		qui replace betamultiple`y' = ${eventmultiple`y'} if t == 0
		qui replace betaduplicate`y' = ${eventduplicate`y'} if t==0
	}

	summ betamultiple if t==0
	local adjust = r(mean)
	replace betamultiple = betamultiple-`adjust'
	//replace betaduplicate = betaduplicate - ${eventduplicate1}
	/*
#delimit ;
graph twoway connected /*betamultiple betad*/ betaig t, msymbol(S + T) ||,
graphregion(color(white)) xsize(7) 
legend(label(1 Multiple Dummies On (Adjusted)) label(3 Ignore Second Event) 
label(2 Duplicate Observations) col(1) order(1 2 3))
title("$eventdesc") ytitle("Deaths per 100,000 Person-Years")
subtitle("Outcome DGP: Empirical")
yline(0) xline(0)
note("Dependent Variable is the Motor Vehicle Fatality rate of 18-20 year olds and the unit of observation is a "
"state-month. An event should be interpreted as a 1 year increase in the state minimum drinking age. Multiple "
"Dummies coefficients are adjusted so that the coefficient for event time zero equals zero. Regression controls"
"for state and month fixed effects.")
;
#delimit cr
*/




#delimit ;
graph twoway connected /*betamultiple betad*/ betaig t, msymbol(S + T) ||,
graphregion(color(white)) xsize(7) 
legend(label(1 Multiple Dummies On (Adjusted)) label(3 Ignore  Event) 
label(2 Duplicate Observations) col(1) order(1 2 3))
ytitle("Deaths per 100,000 Person-Years") xtitle(Months since largest change)
 yscale(range(-4 4.5)) title(Largest Event)
yline(0) xline(0) name(ignore, replace) xlabel(-12(3)12)

;

graph twoway connected betamultiple /*betad betaig*/ t, msymbol(S + T) ||,
graphregion(color(white)) xsize(7) title(Multiple Dummies On)
legend(label(1 Multiple Dummies On (Adjusted)) label(3 Ignore Second Event) 
label(2 Duplicate Observations) col(1) order(1 2 3)) xtitle(Months since event)
ytitle("Deaths per 100,000 Person-Years") yscale(range(-4 4.5))
yline(0) xline(0) name(multiple, replace) xlabel(-12(3)12)
;

graph combine ignore multiple, xsize(10) title(${eventdesc})
note("Dependent Variable is the Motor Vehicle Fatality rate of 18-20 year olds and the unit of observation is a state-month. An event should be interpreted as a 1 year increase in the state minimum drinking age. Multiple"
"Dummies coefficients are adjusted so that the coefficient for event time zero equals zero. Regression controls for state and month fixed effects.")
;
graph export graphs/MLDA.eps, replace;
instanttex using graphs/MLDAgraph.tex, fig(graphs/MLDA) replace ;
 
