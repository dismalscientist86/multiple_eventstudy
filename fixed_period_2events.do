global eventdesc "Two events, second is 6 periods after first"

/*We will randomly assign an event to each location such that we are guarranteed a balanced panel*/
gen rand1 = runiform() if time ==1

/*Generate a lognormal random variable to assign the gap between the second event and the first*/
gen rand2 = runiform() if time == 1 //exp(rnormal(2,.5)) if time ==1

by location: egen event1 = max($panelsize+int(rand1*(${maxtime}-$panelsize*2))) 
gen event2= event1 + 6


global enum = 2
