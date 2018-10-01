clear all


set obs 9


gen t = _n

gen eventtime = _n-5



gen fin = (eventtime == 0)*100 if eventtime <=-2 | eventtime==0

gen io = (eventtime ==0)*100 if eventtime <=0

gen pf = (eventtime>=0)*100

local graphs
foreach var in fin io pf{

#delimit ;
graph twoway connect `var' eventtime, ytitle(Outcome) cmissing(no) lcolor(gs4) mcolor(gs4)
xtitle(Periods after event) xline(0) graphregion(color(white))
;
#delimit cr
graph export graphs/event_example_`var'.eps, replace
local graphs `graphs' graphs/event_example_`var'
}

instanttex using event_example.tex, fig(`graphs') replace

