#******************************************************************************#
# SER 2019: Using Using Monte Carlo Simulations for Quantitative Bias Analysis #
# February ___, 2019                                                           #
# Crystal Shaw                                                                 #
# Questions/Comments: c.shaw@ucla.edu                                          #
#                                                                              #
# Simulation example: Bias from unmeasured confounding                         #
#                                                                              #
# Required files to run simulation:                                            #
#   (1) confounding_sim_par.R:                                                 #
#       -Sets parameter inputs for the simulation                              #
#   (2) confounding_sim_script.R:                                              #
#       -Generates and analyzes data and stores results                        #
#   (3) confounding_sim_run.R:                                                 #
#       -Runs simulation and stores results                                    #                                                              
#******************************************************************************#

#******************************************************************************#
# This file generates and analyzes data and stores the parameter estimates for #
# one iteration of sample generation.                                          #
#******************************************************************************#

#---- Package loading & options ----
set.seed(67208113)
  *Start with a blank data set (no observations)
set more off
clear
*set seed 67208113 	//For multiple iterations of sample generation, the seed 
//is set in the run simulation file


/************************************************************/
  /********	Step A. Create blank data set	  		 ********/
  /************************************************************/
  
  set obs 5000 //creates blank dataset with desired # observations
gen id = _n


/************************************************************/
  /********	Step B. Set parameters					 ********/
  /************************************************************/
  *call preamble file
do 1_SER_simulation_workshop_2018_confounding_preamble.do

*Parameters for odds of exposure
local g0 = $g0 //log odds [p/(1-p)] of exposure for ref group (U=0)
local g1 = $g1 //effect of one-unit higher U on log odds of exposure 

*Parameters for odds of outcome
local b0 = $b0 //log odds [p/(1-p)] of outcome for ref group (U=0,exp=0)
local b1 = $b1 //effect of exposure on log odds of outcome
local b2 = $b2 //effect of one-unit higher U on log odds of outcome 


/************************************************************/
  /********	Step C. Generate data					 ********/
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
                    /******** 	End data generation steps		   		 ********/
                    /************************************************************/
                    
                    
                    /*************************************************************/
                    /********	Step D. Analyze data & store results     ********/
                    /*************************************************************/
                    
                    *Store results as scalar variables	 
                    *Scalar variables store single numbers or  stings, not variables in the dataset								
                    
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
matrix logit1 = r(table)
scalar OR_exposure_Uyes = logit1[1,2]
scalar lb_OR_exposure_Uyes = logit1[5,2]
scalar ub_OR_exposure_Uyes = logit1[6,2]

*Estimated OR & 95% CI for exposure-outcome association
*Without adjusting for U (confounder)
logistic outcome i.exposure 
matrix list r(table)
matrix logit2 = r(table)
scalar OR_exposure_Uno = logit2[1,2]
scalar lb_OR_exposure_Uno = logit2[5,2]
scalar ub_OR_exposure_Uno = logit2[6,2]


*Look at stored results
scalar list
