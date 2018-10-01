global eventdesc "Two events, time to second event follows lognormal"
local eventgap int(rand2)



/*We will randomly assign an event to each location such that we are guarranteed a balanced panel*/
gen rand1 = runiform() if time ==1

/*Generate a lognormal random variable to assign the gap between the second event and the first*/
gen rand2 = exp(rnormal(2,.5)) if time ==1

by location: egen event_r1 = max($panelsize+int(rand1*($maxtime-$panelsize*2))) 

/*The second event will be a fixed period after the first event*/
by location: egen event_r2 = max(event_r1 + `eventgap')
replace event_r2 = event_r2 + 1 if event_r2 == event_r1

global enum
forvalues i = 1/500{
	capture confirm variable event_r`i'
	if _rc != 0 continue, break
	global enum = `i'
}

forvalues i = 1/$enum{
	egen event`i' = rowmin(event_r1-event_r$enum)
	forvalues j = 1/$enum{
		replace event_r`j' = ${maxtime}*2 if event_r`j' == event`i'
	}
}
drop event_r*		
