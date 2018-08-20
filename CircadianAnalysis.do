/***
Chronobiologic Rhythm of Death in Cancer Patients: a Population-Based Study
***/

/***
Authors:
Shing Fung Lee, Miguel Angel Luque-Fernandez (Senior author)
***/

/***
Copyright (c) 2018  

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NON INFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

***/

cd "your_path"
use "cancer_death_hongkong.dta"

// Time to string to extract time in minutes
gen str16  stime = DischargeDateyyyymmddHHMM
gen shours = substr(stime, 12, 2) //takes the first two digits
gen hours = real(shour) //reads first two digits as number
gen sminutes = substr(stime, 15, 2) //takes digits 4 and 5
gen minutes = real(sminutes)
*now add up
gen totminutes = 60*hours + minutes 

// Scatter plot of the number of deaths by time in hours and Prevalence Ratios
preserve 
rename dtime_group time24
collapse (sum) death, by(time24)
summ death, det

tw(scatter death time24)
//Testing for overdispersion: Score based t-test
glm death i.time24, family(poisson) link(log) eform nocons nolog baselevels
predict double mu, mu
//Score based t-test: H0 E(x|y) = Var(x|y)
gen double zrs =((d-mu)^2-d)/(mu)
//Negative binomial model Var(x|y) = E(x|y) + E(x|y)^2 x alpha. Alpha can be approximate using the score based t-test
reg zrs mu, nocons nohead //Beta is an approximation of alpha = .0349428
glm death i.time24, f(nb  .0349428) l(log) eform 
restore

// Scatter plot of the number of deaths by time in minutes (bimodal periodicity)
preserve 
collapse (sum) death, by(totminutes) //by(time24)
tw(scatter death totminutes)(mspline death totminutes)
restore

// Scatter plot of the number of deaths by time in hours (by 4 segments/day) and Prevalence Ratios
drop dtime_4groups
gen dtime_4groups = dtime_group
recode dtime_4groups 0=24
recode dtime_4groups 7/12=0 13/18=1 19/24=2 1/6=3
tab dtime_4groups

preserve 
rename dtime_4groups time24
collapse (sum) death, by(time24)
tw(scatter death time24)
glm death i.time24, family(poisson) link(log) eform nocons nolog baselevels
predict double mu, mu
//Score based t-test: H0 E(x|y) = Var(x|y)
gen double zrs =((d-mu)^2-d)/(mu)
//Negative binomial model Var(x|y) = E(x|y) + E(x|y)^2 x alpha. Alpha can be approximate using the score based t-test
reg zrs mu, nocons nohead //Beta is an approximation of alpha = .347243 
glm death ib2.time24, f(nb .347243) l(log) eform 
restore

//Circadian Models

*Trigonometric Model: One and two cycles per 24 hours of time
preserve 
rename dtime_group time24

collapse (sum) death, by(time24)
gen f=(time24 - 0.5)/24
gen w=(2*_pi)*f
gen sin=sin(w)
gen cos=cos(w)
gen w2=(4*_pi)*f
gen sin2=sin(w2)
gen cos2=cos(w2)

*cosinor two and one cycles in a period of 24 hours
xi: glm death sin cos sin2 cos2 time24, family(poisson) link(log) scale(x2)
xi: glm death sin cos time24, family(poisson) link(log) scale(x2) eform

*Amplitude
sum time24
return list
gen  N = r(N)
nlcom (sqrt(_b[sin]^2 + _b[cos]^2))
*Phase
nlcom  24*((atan(_b[sin]/_b[cos]))/(2*_pi))+1
*Circadian Test
display (sqrt((_b[sin]^2/N) + (_b[cos]^2/N)))
*p-value
display sqrt((sum(sin)/N)^2 + (sum(cos)^2/N)^2)
display N * (.04166667)^2
display exp(-.04166667)

restore
