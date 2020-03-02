/******************************************************************************/
/* Tutorial: Using Monte Carlo Simulations for Quantitative Bias Analysis     */
/* Last updated: March 2, 2020				                      */
/* Elizabeth Rose Mayeda, UCLA Fielding School of Public Health		      */
/* Questions/Comments: ermayeda@ph.ucla.edu				      */
/*						  			      */
/* Simulation example: Bias from unmeasured confounding			      */
/*									      */
/* Required files to run simulation:					      */
/* 1) simulation_workshop_confounding_preamble.do:			      */ 
/*		-Sets parameter inputs for simulation 			      */
/* 2) simulation_workshop_confounding_data_gen_analysis.do: 		      */
/*		-Generates and analyzes data and stores results		      */
/* 3) simulation_workshop_confounding_run.do: 				      */ 
/*		-Runs simulation and stores results			      */
/*									      */
/* Stata version 15							      */
/******************************************************************************/

/******************************************************************************/
/* The purpose of this file is to generate and analyze data and store the     */
/* parameter estimates for one iteration of sample generation		      */
/******************************************************************************/

*Start with a blank data set (no observations)
set more off
clear
*set seed 67208113 	//For multiple iterations of sample generation, the seed 
			//is set in the run simulation file
				
				
/*Define working folder: update with your file path*/
cd "C:\Users\emayeda\Box\GSA2019_bias_workshop\Confounding_simulation_tutorial\confounding_simulation_Stata\Do_files"
			
			
*call preamble file
do 1_simulation_workshop_confounding_preamble.do


/************************************************************/
/********	Step A. Create blank data set	     ********/
/************************************************************/

set obs $num_obs //creates blank dataset with desired # observations
gen id = _n


/************************************************************/
/********	Step B. Set parameters		     ********/
/************************************************************/

*Parameters for odds of exposure
local g0 = $g0 //log odds [p/(1-p)] of exposure for ref group (U=0)
local g1 = $g1 //effect of one-unit higher U on log odds of exposure 

*Parameters for odds of outcome
local b0 = $b0 //log odds [p/(1-p)] of outcome for ref group (U=0,exp=0)
local b1 = $b1 //effect of exposure on log odds of outcome
local b2 = $b2 //effect of one-unit higher U on log odds of outcome 


/************************************************************/
/********	Step C. Generate data		     ********/
/************************************************************/

*Generate U, where U~N(0,1)
gen U = rnormal()

*Generate exposure
*First, generate P(exposure)
*Next, assign each person exposure status based on P(exposure)
gen Pexposure = exp(`g0' + `g1'*U )/(1 + exp(`g0' + `g1'*U))			
gen exposure= runiform()<Pexposure

*Generate outcome
*First, generate P(outcome)
*Next, assign each person outcome status based on P(exposure)
gen Poutcome = exp(`b0' + `b1'*exposure + `b2'*U)/(1 + exp(`b0' + `b1'*exposure + `b2'*U))
gen outcome= runiform()<Poutcome


/************************************************************/
/******** 	End data generation steps	     ********/
/************************************************************/


/*************************************************************/
/********	Step D. Analyze data & store results  ********/
/*************************************************************/

*Store results as scalar variables	 
*Scalar variables store single numbers or stings, not variables in the dataset								

/**** Check distributions and store results ***/

*Check mean value of U and proportion of people with exposure=1 and outcome=1
foreach x in U exposure outcome {
summarize `x'
scalar p_`x' = round(r(mean),0.001)
}
scalar mean_U = p_U
scalar drop p_U 


/*** Look at associations and store results ***/

*Estimated OR & 95% CI for exposure-outcome association
*Adjusting for U (confounder)
logistic outcome i.exposure c.U 
	matrix list r(table)
	matrix matrix1 = r(table)
	scalar OR_exposure_Uyes = matrix1[1,2]
	scalar lb_OR_exposure_Uyes = matrix1[5,2]
	scalar ub_OR_exposure_Uyes = matrix1[6,2]
	
*Estimated OR & 95% CI for exposure-outcome association
*Without adjusting for U (confounder)
logistic outcome i.exposure 
	matrix list r(table)
	matrix matrix2 = r(table)
	scalar OR_exposure_Uno = matrix2[1,2]
	scalar lb_OR_exposure_Uno = matrix2[5,2]
	scalar ub_OR_exposure_Uno = matrix2[6,2]

	
*Look at stored results
scalar list
