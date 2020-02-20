clear all
cap log close
set more 1
set dp comma

use "df_probit_unmatched.dta", clear

global ylist mahog_area
global xlist ln_gdp_pc gdpfr_ag area_plant2 und5_m heart_m infecc_m neop_m suic_m traff_m pol_m  

* ----------
*  PROBIT
* ----------

probit $ylist $xlist
estimates store probit1

* Predicting propensity scores
predict pscore_mt4, pr 

* Plot
qui graph twoway (kdensity pscore_mt4 if mahog_area==1, msize(small)) ///
(kdensity pscore_mt4 if mahog_area==0, msize(small) lpattern(shortdash_dot)), ///
subtitle(, bfcolor(none)) ///
xtitle("p score (2000)", size(medlarge)) ///
xscale(titlegap(*7)) ytitle("Densidade", size(medlarge)) yscale(titlegap(*5)) ///
legend(pos(12) ring(0) col(1)) ///
legend( label(1 "Tratados") label(2 "Controles")) saving(BEFORE, replace)


* ----------
*  PSMATCH2
* ----------

set seed 1000

psmatch2 mahog_area, pscore(pscore_mt4) out (hom_tx) kernel com 
psgraph

qui graph twoway (kdensity pscore_mt4 [aweight=_weight] if mahog_area==1, msize(small)) ///
(kdensity pscore_mt4 [aweight=_weight] if mahog_area==0, msize(small) lpattern(shortdash_dot)), ///
subtitle(, bfcolor(none)) ///
xtitle(" p score (2010) ", size(medlarge)) ///
xscale(titlegap(*7)) ytitle("Densidade", size(medlarge)) yscale(titlegap(*5)) ///
legend(pos(12) ring(0) col(1)) ///
legend( label(1 "Tratados") label(2 "Controles")) saving(AFTER , replace)

graph combine BEFORE.gph AFTER.gph

save df_probit_matched.dta
