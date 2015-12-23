* upload the household level data 
use "~\q1q4iccctrl.dta", clear
* or...upload the individual level data
use "~\q2icc_nohh.dta", clear
* format the dataset
destring cluster, replace
destring dataid, replace
* get ICC for each exposure....sample code....
gllamm q11_tube, family(binomial) link(logit) i(cluster) eform adapt 
nlcom [clus1]_cons^2 / ( [clus1]_cons^2 +3.14*3.14/3)



