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
# This file generates and analyzes data and returns the parameter estimates for#
# one iteration of sample generation.                                          #
#******************************************************************************#

#******************************************************************************#
# Package quick reference and resources                                        #
#                                                                              #
# pacman: clean way to import many R packages at one time                      #
#         (https://www.rdocumentation.org/packages/pacman/versions/0.5.0)      #
#                                                                              #
# here: prevents the issues associated with specifying file paths              #
#       (https://www.rdocumentation.org/packages/here/versions/0.1)            #
#                                                                              #
# tidyverse: a suite of data wrangling, analysis, and visualization packages   #
#            (https://www.rdocumentation.org/packages/tidyverse/versions/1.2.1)#
#                                                                              #
#******************************************************************************#

#******************************************************************************#
# A note on custom functions...                                                #
#   It is best practice to save each user-defined function in a separate       #
#   RScript for debugging and testing purposes. For simplicity in this tutorial# 
#   and because the functions are so small, however, we chose to include them  #
#   in the "#---- Custom functions ----" portion of this RScript.              #
#******************************************************************************#

#---- Package Loading and Options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("here", "tidyverse")

#Using standard notation (as opposed to scientific), rounded to three 
#decimal places
options(scipen = 999)
options(digits = 3)

#---- Custom Functions ----
p_S <- function(A, U){
  exp(g0 + g1*A + g2*U + g3*U*A)/(1 + exp(g0 + g1*A + g2*U + g3*U*A))
}

p_Y <- function(A, U){
  exp(b0 + b1*A + b2*U)/(1 + exp(b0 + b1*A + b2*U))
}

#---- The simulation function ----
collider_sim <- function(){
  #Generating the data
  #Comment about tibbles
  #Comment on each line # about what it does
  #Creating IDs and exogenous variables
  #Generating S and Y
  obs <- tibble("id" = seq(from = 1, to = samp_size, by = 1), 
                "A" = rbinom(n = samp_size, size = 1, p = p_A), 
                "U" = rnorm(n = samp_size, mean = 0, sd = 1), 
                "probS" = p_S(A, U), 
                "probY" = p_Y(A, U), 
                "S" = rbinom(n = samp_size, size = 1, prob = probS), 
                "Y" = rbinom(n = samp_size, size = 1, prob = probY))
  
  #Estimates of primary interest:  Estimated ORs for A and Y among S = 1
  #Store ORs and 95% CI limits
  selected <- obs %>% filter(S == 1)
  sel_model <- glm(Y ~ A, family = binomial(), data = selected)
  sel_OR <- exp(sel_model$coefficients["A"])
  ci_95 <- exp(confint(sel_model, "A", level = 0.95))
  return(c(sel_OR, ci_95))
}

