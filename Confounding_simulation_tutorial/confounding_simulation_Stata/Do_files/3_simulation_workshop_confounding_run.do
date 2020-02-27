/******************************************************************************/
/* GSA Workshop: Bias Analysis 												  */
/* November 13, 2019														  */
/* Elizabeth Rose Mayeda, UCLA Fielding School of Public Health				  */
/* Questions/Comments: ermayeda@ph.ucla.edu									  */
/*																			  */
/* Simulation example: Bias from unmeasured confounding						  */
/*																			  */
/* Required files to run simulation:										  */
/* 1) simulation_workshop_confounding_preamble.do:							  */ 
/*		-Sets parameter inputs for simulation 								  */
/* 2) simulation_workshop_confounding_data_gen_analysis.do: 				  */
/*		-Generates and analyzes data and stores results						  */
/* 3) simulation_workshop_confounding_run.do: 								  */ 
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
cd "C:\Users\emayeda\Box\GSA2019_bias_workshop\Confounding_simulation_tutorial\confounding_simulation_Stata\Do_files"


/*Set working directories (for Windows--modify for Mac)*/
local logs		"..\Logs\"
local data		"..\Data\"
local plots 	"..\Plots\"


/*Save log file*/  
capture log close
log using "`logs'simulation_workshop_confounding_log", replace


/*Start timer*/
timer clear 1
timer on 1

*call preamble file
do 1_simulation_workshop_confounding_preamble.do

/*Specify desired number of iterations of sample generation*/
local B = $B //Number of iterations of sample generation

/*Create local variable for causal/true OR for effect of exposure on outcome*/
local true_OR_exposure_outcome = exp($b1)


/*Create a local macro variable for "scalar_X = variable_X".
We will this local macro in the simulate command. This will pull in scalars  
fromthe data generation and analysis file and store each scalar as a variable
with the same name.*/
local simlist ""
foreach x in mean_U p_exposure p_outcome ///
		OR_exposure_Uyes OR_exposure_Uno ///
		lb_OR_exposure_Uyes ub_OR_exposure_Uyes ///
		lb_OR_exposure_Uno ub_OR_exposure_Uno { 
   local simlist "`simlist' `x'=`x'"
}


/*Run simulation*/
simulate `simlist', ///
reps(`B') seed(67208113): do 2_simulation_workshop_confounding_data_gen_analysis //use name of your data generation do file


/*Across B replications, calculate and store mean value of each variable as a scalar and as a local variable*/ 
*Take the log of the ORs so that means are calculated correctly
*Round scalars to two decimal places
gen log_OR_exposure_Uyes = ln(OR_exposure_Uyes)																						
gen log_OR_exposure_Uno = ln(OR_exposure_Uno)
foreach b in mean_U ///
			 p_exposure ///
			 p_outcome ///
			 log_OR_exposure_Uyes ///
			 log_OR_exposure_Uno {
summarize `b', meanonly
scalar mean_`b' = round(r(mean),0.001)
local mean_`b' = round(r(mean),0.001)
}
scalar exp_mean_log_OR_exposure_Uyes = round(exp(mean_log_OR_exposure_Uyes),0.01)
local exp_mean_log_OR_exposure_Uyes  = round(exp(mean_log_OR_exposure_Uyes),0.01)
scalar exp_mean_log_OR_exposure_Uno  = round(exp(mean_log_OR_exposure_Uno),0.01)
local exp_mean_log_OR_exposure_Uno   = round(exp(mean_log_OR_exposure_Uno),0.01)


/*For each replication, generate indicator variable for whether the 95% CI for the weight status-death OR includes the causal/true OR*/
gen covg_OR_exposure_Uyes = (ub_OR_exposure_Uyes > `true_OR_exposure_outcome' & lb_OR_exposure_Uyes < `true_OR_exposure_outcome')
gen covg_OR_exposure_Uno = (ub_OR_exposure_Uno > `true_OR_exposure_outcome' & lb_OR_exposure_Uno < `true_OR_exposure_outcome')


/*Across B replications, calculate and store mean value of each variable as a scalar and as a local variable*/ 
foreach b in covg_OR_exposure_Uno covg_OR_exposure_Uyes {
summarize `b', meanonly
scalar P_`b' = round(r(mean),0.001)
local P_`b' = round(r(mean),0.001)
}


/*Generate plots of estimated OR and 95% CI across B simulated samples*/
gen sample = _n

*Plot estimated OR adjusted for U
local P_covg_OR_exposure_Uyes : di %4.3f `P_covg_OR_exposure_Uyes' /*fixing rounding problem in this local variable*/
twoway (rcap lb_OR_exposure_Uyes ub_OR_exposure_Uyes sample) ///
(scatter OR_exposure_Uyes sample, msymbol(circle) mcolor(navy) mfcolor(edkblue)) ///
(function y=1 , lwidth(thick) clcolor(red) lpat(solid) range(sample) ) ///
, ///
ytitle(estimated OR (95% CI)) ///
ylabel(0.6(.2)2.8) xscale(off) legend(off) ///
title(Estimated adjusted OR and 95% CI from `B' simulated samples, size(medium)) ///
subtitle(mean estimated OR = `exp_mean_log_OR_exposure_Uyes'; 95% CI coverage = `P_covg_OR_exposure_Uyes')
	graph export `plots'confounding_simulation_plot_est_OR_CI_adj.png,replace

*Plot estimated OR unadjusted for U
twoway (rcap lb_OR_exposure_Uno ub_OR_exposure_Uno sample) ///
(scatter OR_exposure_Uno sample, msymbol(circle) mcolor(navy) mfcolor(edkblue)) ///
(function y=1 , lwidth(thick) clcolor(red) lpat(solid) range(sample) ) ///
, ///
ytitle(estimated OR (95% CI)) ///
ylabel(0.6(.2)2.8) xscale(off) legend(off) ///
title(Estimated crude OR and 95% CI from `B' simulated samples, size(medium)) ////
subtitle(mean estimated OR = `exp_mean_log_OR_exposure_Uno'; 95% CI coverage = `P_covg_OR_exposure_Uno')
	graph export `plots'confounding_simulation_plot_est_OR_CI_crude.png,replace


/*List results  across the B iterations of sample generation*/
*Across B replications, average mean value of U, proportion of people with exposure=1, proportion of people with outcome=1
scalar list mean_mean_U mean_p_exposure mean_p_outcome

*Across B replications, average estimated OR for exposure on outcome, adjusting for U
scalar list exp_mean_log_OR_exposure_Uyes

*Across B replications, average estimated OR for exposure on outcome, not adjusting for U
scalar list exp_mean_log_OR_exposure_Uno

*Proportion of 95% CIs for OR for exposure on outcome that include the causal/true OR 
scalar list P_covg_OR_exposure_Uyes P_covg_OR_exposure_Uno 


/***export data to Excel***/
export excel using `data'simulation_workshop_confounding_results_each_replication, sheet("Results") sheetmodify firstrow(variables)

/***save data***/
save "`data'simulation_workshop_confounding_results_each_replication.dta", replace

/*End timer and display computational time (in seconds)*/
timer off 1
timer list 1
