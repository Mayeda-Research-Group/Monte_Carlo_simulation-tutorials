/******************************************************************************/
/* GSA Workshop: Bias Analysis 												  */
/* November 13, 2019														  */
/* Elizabeth Rose Mayeda, UCLA Fielding School of Public Health				  */
/* Questions/Comments: ermayeda@ph.ucla.edu									  */
/*																			  */
/* Simulation example: Collider-stratification bias							  */
/*																			  */
/* Required files to run simulation:										  */
/* 1) simulation_workshop_colldier_bias_preamble.do:						  */ 
/*		-Sets parameter inputs for simulation 								  */
/* 2) simulation_workshop_colldier_bias_data_gen_analysis.do: 				  */
/*		-Generates and analyzes data and stores results						  */
/* 3) simulation_workshop_colldier_bias_run.do: 							  */ 
/*		-Runs simulation and stores results									  */
/*																			  */
/* Stata version 15															  */
/******************************************************************************/

/******************************************************************************/
/* The purpose of this file is to run multiple iterations of sample 		  */
/* and summarize and store results.											  */
/******************************************************************************/

/*******************************************************************************
The "simulate" command performs Monte Carlo-type simulations. 
The basic form of the code is: 
	"simulate <list of scalars from each iteration of sample generation>, /// 
	reps(B) seed(#): do <data_generation_and_analysis_do_file>" 
	 
This calls the data generation and analysis file B times to generate and 
analyze B data sets and stores the results from each iteration of sample 
generation as scalars. 
	-Note that the first lines of code call in and assign the value from each 
	 scalar to a variable with the same name. 
	-"reps(B)" specifies the number of iterations of sample generation. 
	-"seed(#)" specifies the random number seed. 

Simulate documentation: http://www.stata.com/manuals13/rsimulate.pdf 

The rest of the code summarizes results across the B iterations of sample 
generation stores the results from the B replications in a Stata data file  
and in an Excel file(one iteration of sample generation = one row in the file). 
*******************************************************************************/

set more off
clear all


/*Define working folder: update with your file path*/
cd "C:\Users\emayeda\Box\GSA2019_bias_workshop\Collider_bias_simulation_tutorial\collider_bias_simulation_Stata\Do_files"


/*Set working directories (for Windows--modify for Mac)*/
local logs		"..\Logs\"
local data		"..\Data\"
local plots 	"..\Plots\"


/*Save log file*/  
capture log close
log using "`logs'simulation_workshop_collider_bias_log", replace


/*Start timer*/
timer clear 2
timer on 2

*call preamble file
do 1_simulation_workshop_collider_bias_preamble.do

/*Specify desired number of iterations of sample generation*/
local B = $B

/*Create local variable for causal/true OR for A on Y specified in data
generation and analysis file*/
local causal_OR_AY = exp($b1)


/*Create a local macro variable for "scalar_X = variable_X".
We will this local macro in the simulate command. This will pull in scalars  
from the data generation and analysis file and store each scalar as a variable
with the same name.*/
local simlist ""
foreach x in OR_AY_S1 ub_OR_AY_S1 lb_OR_AY_S1 ///
	 OR_AY_all ub_OR_AY_all lb_OR_AY_all ///
	 mean_U ///
	 mean_U_A1_all mean_U_A0_all ///
	 mean_U_A1_S1 mean_U_A0_S1 ///
	 p_A p_Y p_S ///
	 p_S_A1 p_S_A0  ///
	 p_Y_A0 p_Y_A1 ///
	 p_Y_A0_S1 p_Y_A1_S1 ///
{ 
   local simlist "`simlist' `x'=`x'"
}


/***Run simulation***/
simulate `simlist', ///
reps(`B') seed(67208105): do 2_simulation_workshop_collider_bias_data_gen_analysis //replace with name of your data generation do file

/*Across all B replications, calculate and store mean value of each variable as a scalar*/ 
*Take the log of the ORs so that means are calculated correctly
*Round scalars to two decimal places
gen log_OR_AY_S1 = ln(OR_AY_S1)																						
gen log_OR_AY_all = ln(OR_AY_all)
foreach b in log_OR_AY_S1 ///
			 log_OR_AY_all ///
			 mean_U ///
			 mean_U_A1_all mean_U_A0_all ///
			 mean_U_A1_S1 mean_U_A0_S1 ///
			 p_A p_Y p_S ///
			 p_S_A1 p_S_A0  ///
			 p_Y_A0 p_Y_A1 ///
			 p_Y_A0_S1 p_Y_A1_S1 {
summarize `b', meanonly
scalar mean_`b' = round(r(mean),0.001)
local mean_`b' = round(r(mean),0.001)
}
scalar exp_mean_log_OR_AY_S1   = round(exp(mean_log_OR_AY_S1),0.01)
local exp_mean_log_OR_AY_S1    = round(exp(mean_log_OR_AY_S1),0.01)
scalar exp_mean_log_OR_AY_all  = round(exp(mean_log_OR_AY_all),0.01)
local exp_mean_log_OR_AY_all   = round(exp(mean_log_OR_AY_all),0.01)


/*For each sample, generate indicator variable for whether the 95% CI for 
the estimated OR_AY includes the causal/true OR_AY*/
gen covg_OR_AY_S1 = (lb_OR_AY_S1  < `causal_OR_AY' & ub_OR_AY_S1  > `causal_OR_AY')
gen covg_OR_AY_all = (lb_OR_AY_all  < `causal_OR_AY' & ub_OR_AY_all  > `causal_OR_AY')


/*Across all B replications, calculate and store the proportion of times the
95% CI includes the causal/true OR_AY (95% CI coverage)*/ 
foreach b in covg_OR_AY_S1 covg_OR_AY_all {
summarize `b', meanonly
scalar P_`b' = round(r(mean),0.001)
local P_`b' = round(r(mean),0.001)
}


/*Generate plots of estimated OR and 95% CI across B simulated samples*/
gen sample = _n

*Plot estimated OR S=1
twoway (rcap lb_OR_AY_S1 ub_OR_AY_S1 sample) ///
(scatter OR_AY_S1 sample, msymbol(circle) mcolor(navy) mfcolor(edkblue)) ///
(function y=1 , lwidth(thick) clcolor(red) lpat(solid) range(sample) ) ///
, ///
ytitle(estimated OR (95% CI))  ///
ylabel(0.6(.2)2.0) xscale(off) legend(off) ///
title(Memory complaints=1: Estimated OR and 95% CI from `B' simulated samples, size(medsmall) ) ///
subtitle(mean estimated OR = `exp_mean_log_OR_AY_S1'; 95% CI coverage = `P_covg_OR_AY_S1', size(medsmall))
	graph export `plots'plot_est_OR_CI_S1.png,replace

*Plot estimated OR for whole sample
local P_covg_OR_AY_all : di %4.3f `P_covg_OR_AY_all' /*fixing rounding problem in this local variable*/
twoway (rcap lb_OR_AY_all ub_OR_AY_all  sample) ///
(scatter OR_AY_all sample, msymbol(circle) mcolor(navy) mfcolor(edkblue)) ///
(function y=1 , lwidth(thick) clcolor(red) lpat(solid) range(sample))  ///
, ///
ytitle(estimated OR (95% CI))  ///
ylabel(0.6(.2)2.0) xscale(off) legend(off) ///
title(Whole population: Estimated OR and 95% CI from `B' simulated samples, size(medsmall) ) ///
subtitle(mean estimated OR = `exp_mean_log_OR_AY_all'; 95% CI coverage = `P_covg_OR_AY_all', size(medsmall))
	graph export `plots'plot_est_OR_CI_all.png,replace


/*Plot mean U across the B iterations of sample generation by S and overall*/
*Plot for S=1
twoway (histogram mean_U_A1_S1, start(-0.11) width(0.01) color(blue)) ///
       (histogram mean_U_A0_S1, start(-0.11) width(0.01) ///
	   fcolor(cranberry) lcolor(cranberry)), ytitle(, size(large)) ///
	   ylabel(, labsize(large)) xlabel(, labsize(large)) ///
	   xscale(range(-0.1 1.2)) xlabel(-0.2(.2)1.2) ///
	   yscale(range(0 25)) ylabel(0(5)25) ///
	   legend(order(1 "anxiety=1" 2 "anxiety=0" ) size(large)) ///
	   title(memory complaints=1: distribution of U,size(large)) graphregion(color(white)) ///
	   subtitle(mean U anxiety=1 = `mean_mean_U_A1_S1'; mean U anxiety=0 = `mean_mean_U_A0_S1')
	   
	   graph export `plots'histogram_S1.png,replace
	   
*Plot for whole sample	   
twoway (histogram mean_U_A1_all, start(-0.11) width(0.01) fcolor(blue) lcolor(blue)) ///
       (histogram mean_U_A0_all, start(-0.11) width(0.01) ///
	   fcolor(cranberry) lcolor(cranberry)), ytitle(, size(large)) ///
	   ylabel(, labsize(large)) xlabel(, labsize(large)) ///
	   xscale(range(-0.1 1.2)) xlabel(-0.2(.2)1.2) ///
	   yscale(range(0 25)) ylabel(0(5)25) ///
	   legend(order(1 "anxiety=1" 2 "anxiety=0" ) size(large)) ///
	   title(whole population: distribution of U,size(large)) graphregion(color(white)) ///
	   subtitle(mean U anxiety=1 = `mean_mean_U_A1_all'; mean U anxiety=0 = `mean_mean_U_A0_all')
	   
	   graph export `plots'histogram_all.png,replace
	

/***List results across the B iterations of sample generation***/

*Check proportions of people with: 
	*exposure (A) = 1
	*selection (S) = 1
	*outcome (Y) = 1
	*outcome (Y) = 1 by exposure (A)
scalar list mean_p_A mean_p_S mean_p_Y mean_p_Y_A1 mean_p_Y_A0 mean_mean_U 

*Check proportion of people with outcome (Y) = 1 by exposure (A) among S=1
scalar list mean_p_Y_A1_S1 mean_p_Y_A0_S1

*Check distribution of U by A among S=1
scalar list mean_mean_U_A1_S1 mean_mean_U_A0_S1 

*Check distribution of U by A among full sample
scalar list mean_mean_U_A1_all mean_mean_U_A0_all

*Check mean OR and 95% CI coverage in whole population (no bias anticipated)
scalar list exp_mean_log_OR_AY_all
scalar list P_covg_OR_AY_all

*Estimate of primary interest: Mean estimated ORs among S=1
scalar list exp_mean_log_OR_AY_S1

*Proportion of 95% CIs that include the causal/true OR 95% CI coverage)
scalar list P_covg_OR_AY_S1 					


/*Save results. One row = one iteration of sample generation*/	   
*Excel file
export excel using `data'simulation_workshop_collider_bias_results_each_replication, sheet("Results") sheetmodify firstrow(variables)

*Stata data file
save "`data'simulation_workshop_collider_bias_results_each_replication.dta", replace


/*End timer and display computational time (in seconds)*/
timer off 2
timer list 2
