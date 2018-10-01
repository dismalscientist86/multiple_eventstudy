/*set more off
quietly{
/*Read in the state data from the University of Kentucky Center for Poverty Research*/
insheet using UKCPR_National_Data_Set_12_14_11.csv, comma clear

/*Drop blanks*/
drop if state==.

/*Keep only variables of interest*/
keep stname state year employment population unemploymentrate federalminimumwage stateminimumwage povertyrate grossstateproduct percentlowincomeunisuredchildren afdctanfrecipients foodstampsnaprecipients

/*Rename long variables*/
rename unemploymentrate unemployment
rename federalminimumwage minwage_fed
rename stateminimumwage minwage_state
rename povertyrate poverty
rename percentlowincomeunisuredchildren uninsured_children
rename afdctanfrecipients AFDC_recipients
rename foodstampsnaprecipients foodstamp_recipients
rename grossstateproduct GSP


/*Remove commas*/
foreach var in AFDC_recipients foodstamp_recipients employment population GSP{
	replace `var' = subinstr(`var',",","",.)
}

qui destring, replace

gen GSP_PC = GSP/population

joinby state year using minwage_cps


/*Create a combined minimum wage variable*/
*gen minwage=max(minwage_fed, minwage_state)
gen minwage=minwage_state
/*Create a variable that indicates a change in the minimum wage*/
xtset state year
gen changewage=round(minwage-L.minwage,.01)
replace changewage=0 if year==1980 | (minwage_state<=minwage_fed & minwage_fed ==L1.minwage_fed)
gen event=changewage!=0
*Drop events that come from a COL or inflation adjustment (Information on COL adjustments from US Department of Labor)
replace event=0 if stname=="CO" & year>=2008
replace event=0 if stname=="FL" & year>=2007
replace event=0 if stname=="MO" & year>=2008
replace event=0 if stname=="MT" & year>=2011
replace event=0 if stname=="NV" & year>=2008
replace event=0 if stname=="OR" & year>=2004
replace event=0 if stname=="VT" & year>=2007
replace event=0 if stname=="WA" & year>=2001
gen percentchange=changewage/L.minwage

/*Count the number of events*/
bysort state: gen eventcount=sum(event)

/*Create variables that capture the year of event*/
forvalues i=1/19{
	by state : egen event`i' = min(year + 1000000 * (eventcount != `i'))
	replace event`i'=. if event`i'>3000
}
forvalues i = 1/19{
	by state : egen eventsize`i' = max(changewage* (eventcount == `i'))
}
 
}
 save minwage.dta, replace
 */
 use minwage, replace
 set more off
 gen y1=wage_male
 global ynum=1
 global enum=19
 global panelsize=4
 gen time = year
 drop event
gen location = state
 global eventspec minwage
 global eventdesc "Events Are Changes in the State Minimum Wage"
 global y1lab "Empirical"
 global covariates i.location i.time unemployment GSP_PC population
 global inwindow on
 
 xtset location time
 
 do run_eventstudies.do
 
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

	
	replace betamultiple = betamultiple-${lag1multiple1}
	
#delimit ;
graph twoway connect betamultiple betad t, msymbol(S + T) ||,
graphregion(color(white)) xsize(5) 
legend(label(1 Multiple Dummies On (Adjusted)) label(3 Ignore Second Event) 
label(2 Duplicate Observations) col(1) order(1 2))
yline(0) xline(0)
title("$eventdesc")
subtitle("Outcome DGP: Empirical")
note("Dependent Variable is the 10th percentile of wages earned by men 18-25, calculated from "
"the CPS, and the unit of observation is a state-year. An event should be interpreted as a 1 "
"dollar increase in the state minumum wage. Multiple Dummies coefficients are adjusted so "
"that the coefficient for event time zero equals zero. Regression controls for state and year "
"fixed effects, state population, GSP per capita, and the state unemployment rate.")
;
#delimit cr

graph export graphs/minwage.eps, replace

instanttex using graphs/minwagegraph.tex, fig(graphs/minwage) replace
 
