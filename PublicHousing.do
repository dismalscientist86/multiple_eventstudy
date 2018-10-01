/* This program sets up the public housing demolition data in order to run the 3 event study methods*/
/*capture log close
log using PublicHousing.log, replace
use ~/Public_Housing/OneObsMerge.dta, clear
set matsize 10000
macro drop _all
*Only looking within a quarter mile radius
keep if distance<=.25

*Create Event Dummies
gen demomonth=mofd(demo_start)
gen eventtime=month-demomonth
gen event=(eventtime==0)

*Collapse to one observation per block month	
collapse (sum) event units (mean) total,  by(month year monthofyear tract block)

xtset block month

/*Count the number of events*/
bysort block: gen eventcount=sum(event)

/*Create variables that capture the month of event*/
forvalues i=1/9{
	by block : egen event`i' = min(month + 1000000 * (eventcount != `i'))
	qui replace event`i'=. if event`i'>30000
	}
forvalues i=1/9{
	*by block : egen eventsize`i' = max(units* (eventcount == `i'))
	by block : egen eventsize`i' = max(1* (eventcount == `i'))
}
 
 
 save Public_Housing.dta, replace
 
 */
 set matsize 1000
 use Public_Housing, replace
 set more off
 gen y1=total
 global ynum=1
 global enum=9
 global panelsize=24
 gen time = month
 drop event
gen location = block
 global eventspec publichousing
 global eventdesc "Events Are Demolitions of Public Housing"
 global y1lab "Empirical"
 global covariates i.location i.time 
 global inwindow off
 
 xtset location time
 
 do run_eventstudies_agg.do
 
 /*Save the betas*/
	drop _all 
	local betas = $panelsize*2 + 1
	set obs `betas'

	gen t = _n-($panelsize + 1)
	forvalues y = 1/$ynum{
		foreach method in multiple ignore duplicate{	
			gen beta`method'`y' = .
			forvalues i = 2(2)$panelsize{
				qui replace beta`method'`y' = ${lag`i'`method'`y'} if t == -`i'
				qui replace beta`method'`y' = ${lead`i'`method'`y'} if t== `i'
			}
	
		}
		qui replace betamultiple`y' = ${eventmultiple`y'} if t == 0
		qui replace betaduplicate`y' = ${eventduplicate`y'} if t==0
		qui replace betaignore`y' = 0 if t==0
	}

	summ betamultiple if t==0
	local adjust = r(mean)
	replace betamultiple = betamultiple-`adjust'
	
	summ betaduplicate if t==0
	local adjust2 = r(mean)
	replace betaduplicate = betaduplicate-`adjust2'

*One graph 	
#delimit ;
graph twoway connect betamultiple betad  t, msymbol(S  T) mcolor(gs2 gs12) lcolor(gs2 gs12)
graphregion(color(white)) 
legend(label(1 Multiple Dummies On (Adjusted))  
label(2 Duplicate Observations (Adjusted)) col(1) order(1 2))
title("$eventdesc") xlabel(-$panelsize(12)$panelsize)
subtitle("Outcome DGP: Empirical")
xline(0, lwidth(medthick)) text(.10 0.25 "Demolition", place(ne) just(left) width(10)) 
xline(-6) text(.13 -5.75 "Eviction", place(ne) just(left) width(10))
yline(0)
note("Dependent Variable is the total crime committed within a 1/4 mile of the demolition, and the unit of "
"observation is a block-month. Coefficients are adjusted so that the coefficient for event time zero"  
"equals zero. Regression controls for state and year fixed effects.")
;
graph export graphs/publichousing.eps, replace;

*Three graphs ;
graph twoway connect betamultiple t, msymbol(S) mcolor(gs2) lcolor(gs2)
graphregion(color(white)) 
title("$eventdesc") xlabel(-$panelsize(12)$panelsize)
subtitle("Outcome DGP: Empirical")
xline(0, lwidth(medthick)) text(.10 0.25 "Demolition", place(ne) just(left) width(10)) 
xline(-6) text(.13 -5.75 "Eviction", place(ne) just(left) width(10))
yline(0)
;
graph save multiple, replace;
graph export graphs/ph_multiple.eps, replace;

graph twoway connect betad  t, msymbol(T) mcolor(gs12) lcolor(gs12) 
graphregion(color(white)) 
xlabel(-$panelsize(12)$panelsize)
xline(0, lwidth(medthick)) text(.10 0.25 "Demolition", place(ne) just(left) width(10)) 
xline(-6) text(.13 -5.75 "Eviction", place(ne) just(left) width(10))
yline(0)
;
graph save duplicate, replace;
graph export graphs/ph_duplicate.eps, replace;

gen yline=0 ;

graph twoway (connect betaignore  t, msymbol(+) mcolor(gs8) lcolor(gs8)) (line yline t, lcolor(black)) ,
graphregion(color(white)) title("$eventdesc") subtitle("Outcome DGP: Empirical")
xlabel(-$panelsize(12)$panelsize)
xline(0, lwidth(medthick)) text(.10 0.25 "Demolition", place(ne) just(left) width(10)) 
xline(-6) text(.13 -5.75 "Eviction", place(ne) just(left) width(10))
 yscale(r(-.1 .15)) ylabel(-.1(.05).15) ytitle("")
legend(label(1 Ignore Second Event) order(1))
note("Dependent Variable is the total crime committed within a 1/4 mile of the demolition, and the unit of "
"observation is a block-month. Regression controls for state and year fixed effects.")
;
graph save ignore, replace;
graph export graphs/ph_ignore.eps, replace;

graph combine multiple.gph duplicate.gph ignore.gph, ycommon xcommon cols(1) colfirst;

graph export graphs/publichousingcombo.eps, replace;
#delimit cr



instanttex using graphs/PublicHousinggraph.tex, fig(graphs/publichousing graphs/publichousingcombo graphs/ph_multiple graphs/ph_duplicate graphs/ph_ignore) replace
 
