global eventdesc "Max 2 events, second happens with probability .5"

/*We will randomly assign an event to each location such that we are guarranteed a balanced panel*/
gen rand1 = runiform() if time ==1

/*Generate a lognormal random variable to assign the gap between the second event and the first*/
gen rand2 = runiform() if time == 1 //exp(rnormal(2,.5)) if time ==1

/*Generate a variable to determine whether the second event happens*/

gen nextevent = runiform() if time == 1
by location: egen next = max(nextevent > .5)
gen events = 1 + next

by location: egen event_r1 = max($panelsize+int(rand1*(${maxtime}-$panelsize*2))) 

/*The second event will be a fixed period after the first event*/
by location: egen event_r2 = max($panelsize+int(rand2*(${maxtime}-$panelsize*2))) if events ==2 //max(event1 + `eventgap')



global enum
forvalues i = 1/500{
	capture confirm variable event_r`i'
	if _rc != 0 continue, break
	global enum = `i'
}

/*Now get the event variables ordered sequentially; this will make our lives easier later*/

/*Build a macro to test how many events have been ordered so far*/
local missingtest 0
forvalues i = 1/$enum{
	local missingtest `missingtest' + (event_r`i'==240)
}

/*Make event i the smallest event time that has not been ordered.  Then replace the event_r
variable to be very large.  But only do this once for each event*/
forvalues i = 1/$enum{
	egen event`i' = rowmin(event_r1-event_r$enum) if events >=`i'
	forvalues j = 1/$enum{
		replace event_r`j' = ${maxtime}*2 if event_r`j' == event`i' & event_r`j'!=. & `missingtest' <`i'
	}
}		
