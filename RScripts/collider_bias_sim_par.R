#******************************************************************************#
# Using Using Monte Carlo Simulations for Quantitative Bias Analysis           #
# Date                                                                         #
# Crystal Shaw                                                                 #
# Questions/Comments: c.shaw@ucla.edu                                          #
#                                                                              #
# Simulation example: Collider-stratification bias                             #
#                                                                              #
# Required files to run simulation:                                            #
#   (1) collider_bias_sim_par.R:                                               #
#       -Sets parameter inputs for the simulation                              #
#   (2) collider_bias_sim_script.R:                                            #
#       -Generates and analyzes data and returns results                       #
#   (3) colider_bias_sim_run.R:                                                #
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
#---- Specify the model parameters ----
samp_size <- 5000

#Specify the prevalence of A (exposure)
p_A <- 0.2

#Parameters for the odds of S (selection)
g0 <- log(0.10/(1 - 0.10)) #log odds of S for ref group (A = 0 and U = 0)
g1 <- log(5.0) #log OR for effect of A on log odds of selection (OR = 5.0)
g2 <- log(5.0) #log OR for effect of U on log odds of selection (OR = 5.0)
g3 <- log(1.0)  #log OR for interaction between A and U on    (OR = 5.0)

#Parameters for the odds of Y (outcome)
b0 <- log(0.05/(1 - 0.05)) #log odds of Y for ref group (A = 0, U = 0, S = 0)
b1 <- log(1.0) #log OR for effect of A on log odds of Y (OR = 1.0)
b2 <- log(5.0) #log OR for effect of U on log odds of Y (OR = 5.0)

