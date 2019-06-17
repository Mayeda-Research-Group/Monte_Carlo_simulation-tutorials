#******************************************************************************#
# Using Using Monte Carlo Simulations for Quantitative Bias Analysis           #
# Date                                                                         #
# Crystal Shaw, PhD Student UCLA department of Biostatistics                   #
# Questions/Comments: c.shaw@ucla.edu                                          #
#                                                                              #
# Simulation example: Bias from unmeasured confounding                         #
#                                                                              #
# Required files to run simulation:                                            #
#   (1) confounding_sim_par.R:                                                 #
#       -Sets parameter inputs for the simulation                              #
#   (2) confounding_sim_script.R:                                              #
#       -Generates and analyzes data and returns results                       #
#   (3) confounding_sim_run.R:                                                 #
#       -Runs simulation and stores results                                    #                                                              
#******************************************************************************#

#******************************************************************************#
# This file sets the parameter values (stored as global variables) to be used  #
# in the data generation and analysis file.  This allows for different files   #
# for each parameter set but a single file for data-generation and analysis    #
# code. If a specific path is not included for a given causal structure of     #
# interest, set the parameter to 0.                                            #
#******************************************************************************#

#******************************************************************************#
#                             SET PARAMETERS                                   #
#******************************************************************************#
B = 1000         #number of iterations of sample generation
num_obs = 5000   #number of observations in each sample

#---- Parameters for odds of exposure ----
g0 = log(0.20/(1-0.20))	#log odds [p/(1-p)] of exposure for ref group (U = 0)
g1 = log(2.0) 			    #effect of one-unit higher U on log odds of exposure 

#---- Parameters for odds of outcome ----
b0 = log(0.10/(1-0.10)) #log odds [p/(1-p)] of outcome for ref group 
                        #(U = 0, exp = 0)
b1 = log(1.0) 			    #effect of exposure on log odds of outcome
b2 = log(2.0) 		      #effect of one-unit higher U on log odds of outcome 


