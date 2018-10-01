global eventdesc "Two events, 5-7 period gap in between"
local eventgap int(rand2) +5



/*We will randomly assign an event to each location such that we are guarranteed a balanced panel*/
gen rand1 = runiform() if time ==1

/*Generate a lognormal random variable to assign the gap between the second event and the first*/
gen rand2 = runiform()*3 if time == 1 //exp(rnormal(2,.5)) if time ==1

by location: egen event1 = max($panelsize+int(rand1*(${maxtime}-$panelsize*2))) 

/*The second event will be a fixed period after the first event*/
by location: egen event2 = max(event1 + `eventgap')

global enum = 2
