/******************************************************************************/
/* Tutorial: Using Monte Carlo Simulations for Quantitative Bias Analysis     */
/* Last updated: March 2, 2020				 	              */
/* Elizabeth Rose Mayeda, UCLA Fielding School of Public Health		      */
/* Questions/Comments: ermayeda@ph.ucla.edu				      */
/*							   		      */
/* Simulation example: Collider-stratification bias			      */
/*								  	      */
/* Required files to run simulation:					      */
/* 1) simulation_workshop_colldier_bias_preamble.do:			      */ 
/*		-Sets parameter inputs for simulation 			      */
/* 2) simulation_workshop_colldier_bias_data_gen_analysis.do: 		      */
/*		-Generates and analyzes data and stores results		      */
/* 3) simulation_workshop_colldier_bias_run.do: 			      */ 
/*		-Runs simulation and stores results			      */
/*																			  */
/* Stata version 15															  */
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
cd "C:\Users\emayeda\Box\GSA2019_bias_workshop\Collider_bias_simulation_tutorial\collider_bias_simulation_Stata\Do_files"
			
			
*call preamble file
do 1_simulation_workshop_collider_bias_preamble.do					
	
	
/*******************************************************/
/********	Step A. Create blank data set	********/
/*******************************************************/

set obs $num_obs //creates blank dataset with desired # observations
gen id = _n


/*******************************************************/
/********	Step B. Set parameters		********/
/*******************************************************/

*Specify prevalence of A (exposure)
local P_A = $P_A 

*Parameters for odds of S (selection)
local g0 = $g0 	//log odds of S for ref group (A=0 and U=0) 
local g1 = $g1 	//log OR for effect of A on log odds of selection (OR=5.0)	
local g2 = $g2 	//log OR for effect of U on log odds of selection (OR=5.0)
local g3 = $g3 	//log OR for interaction between A and U on  (OR=1.0)


*Parameters for odds of Y (outcome)
local b0 = $b0 	//log odds of Y for ref group (A=0, U=0)
local b1 = $b1 	//log OR for effect of A on log odds of Y (OR=1.0) 
local b2 = $b2 	//log OR for effect of U on log odds of Y (OR=5.0)


/*******************************************************/
/********	Step C. Generate data		********/
/*******************************************************/

*Generate A, a binary variable drawn from a Bernoulli distribution with P(A=1) = "P_A"
gen A = 0
replace A = runiform()<`P_A'


*Generate U, a continuous variable where U~N(0,1)
gen U = rnormal() 

*Generate S
*First, generate P(S)
*Next, assign each person S status based on P(S)
gen P_S = exp(`g0' + `g1'*A + `g2'*U + `g3'*U*A)/(1 + exp(`g0' + `g1'*A + `g2'*U + `g3'*U*A))			
gen S = runiform()<P_S

*Generate Y 
*First, generate P(Y)
*Next, assign each person Y status based on P(Y)
gen P_Y = exp(`b0' + `b1'*A + `b2'*U)/(1 + exp(`b0' + `b1'*A + `b2'*U))		
gen Y = runiform()<P_Y

/*******************************************************/
/******** 	End data generation steps	********/
/*******************************************************/


/*******************************************************/
/********	Step D. Look at data 		********/
/*******************************************************/

/********Check distributions of variables and store results********/

*Check proportions of people with: exposure (A) = 1, selection (S) = 1, outcome (Y) = 1
*Check mean U
foreach x in A U S Y {
summarize `x'
	scalar p_`x' = round(r(mean),0.001)
}
scalar mean_U = p_U
	scalar drop p_U

*Check proportion of people with selection (S) = 1 by exposure (A)
summarize S if A==0, meanonly
	scalar p_S_A0 = round(r(mean),0.001)
summarize S if A==1, meanonly
	scalar p_S_A1 = round(r(mean),0.001)

*Check proportion of people with outcome (Y) = 1 by exposure (A)
summarize Y if A==1, meanonly
	scalar p_Y_A1 = round(r(mean),0.001)
summarize Y if A==0, meanonly
	scalar p_Y_A0 = round(r(mean),0.001)
	
*Check proportion of people with outcome (Y) = 1 by exposure (A) among S=1
summarize Y if (A==1 & S==1), meanonly
	scalar p_Y_A1_S1 = round(r(mean),0.001)
summarize Y if (A==0 & S==1), meanonly
	scalar p_Y_A0_S1 = round(r(mean),0.001)

*Look at mean U by A
sum U if (A==0), meanonly
	scalar mean_U_A0_all = r(mean)
sum U if (A==1), meanonly
	scalar mean_U_A1_all = r(mean)
	
*Look at mean U by A among S=1
sum U if (A==0 & S==1), meanonly
	scalar mean_U_A0_S1 = r(mean)
sum U if (A==1 & S==1), meanonly
	scalar mean_U_A1_S1 = r(mean)


/********Look at associations and store results********/

*Check ORs for A and Y (whole population). No bias anticipated
logistic Y A
	matrix list r(table)
	matrix matrix2 = r(table)
	scalar OR_AY_all = matrix2[1,1]
	scalar lb_OR_AY_all = matrix2[5,1]
	scalar ub_OR_AY_all = matrix2[6,1]
	
*Estimates of primary interest: Estimated ORs for A and Y among S=1 (store ORs and 95% CI limits)
logistic Y A if (S==1)
	matrix list r(table)
	matrix matrix3 = r(table)
	scalar OR_AY_S1 = matrix3[1,1]
	scalar lb_OR_AY_S1 = matrix3[5,1]
	scalar ub_OR_AY_S1 = matrix3[6,1]


*Look at stored results
scalar list
