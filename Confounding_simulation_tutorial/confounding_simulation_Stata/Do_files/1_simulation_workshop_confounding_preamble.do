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
/* The purpose of this file is to set the parameter values (stored as 		  */
/* global macro variables) to be used in the data-generation and analysis 	  */ 
/* file. This way, you can have separate files for different parameter sets   */ 
/* and use one file for the data-generation and analysis code. 		 		  */
/* If a specific path is not included for a given causal structure of 		  */ 
/* interest, set the parameter to 0.										  */
/******************************************************************************/


/*******************************************************/
/******** Set parameters						********/
/*******************************************************/
global B = 1000         		//number of iterations of sample generation
global num_obs = 5000   		//number of observations in each sample

*Parameters for odds of exposure
global g0 = ln(0.20/(1-0.20))	//log odds [p/(1-p)] of exposure for ref group (U=0)
global g1 = ln(2.0) 			//effect of one-unit higher U on log odds of exposure 

*Parameters for odds of outcome
global b0 = ln(0.10/(1-0.10)) 	//log odds [p/(1-p)] of outcome for ref group (U=0,exp=0)
global b1 = ln(1.0) 			//effect of exposure on log odds of outcome
global b2 = ln(2.0) 			//effect of one-unit higher U on log odds of outcome 


