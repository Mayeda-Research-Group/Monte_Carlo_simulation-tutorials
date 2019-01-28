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

p_load("tidyverse")

#Using standard notation (as opposed to scientific)
options(scipen = 999)

#---- Custom Functions ----
#Probability of selection
p_S <- function(A, U){
  exp(g0 + g1*A + g2*U + g3*U*A)/(1 + exp(g0 + g1*A + g2*U + g3*U*A))
}

#Probability of exposure
p_Y <- function(A, U){
  exp(b0 + b1*A + b2*U)/(1 + exp(b0 + b1*A + b2*U))
}

#---- The simulation function ----
collider_sim <- function(){
  #---- Creating the dataset ----
  # We generate a tibble to store all the data. Tibbles are just specific types 
  # of dataframes... they only print the first 10 rows in the console when 
  # called (unless more rows are specified) and list the datatype under each 
  # column name. I find them nifty for debugging issues with datatypes
  
  # What's happening in each column...
  #   "id": Numbering the observations
  #   "A": Flip a coin with success probability p_A for each observation
  #   "U": Creating num_obs number of standard normal random variables
  #   "prob_S": Applying the p_S function to the column of A values and column
  #             of U values to determine probability of selection
  #   "prob_Y": Applying the p_Y function to the column of A values and column
  #             of U values to determine probability of exposure
  #   "S": Flip a coin with success probability prob_S for each observation
  #   "Y": Flip a coin with success probability prob_y for each observation
  obs <- tibble("id" = seq(from = 1, to = num_obs, by = 1), 
                "A" = rbinom(n = num_obs, size = 1, p = p_A), 
                "U" = rnorm(n = num_obs, mean = 0, sd = 1), 
                "probS" = p_S(A, U), 
                "probY" = p_Y(A, U), 
                "S" = rbinom(n = num_obs, size = 1, prob = probS), 
                "Y" = rbinom(n = num_obs, size = 1, prob = probY))
  
  #---- Analyze data ----
  #Check the mean value of U (should be very close, if not exactly, 0)
  #Calculate the proportion of people with exposure (A) = 1
  #Calculate the proportion of people with selection (S) = 1
  #Calculate the proportion of people with outcome (Y) = 1
  
  mean_U <- mean(obs$U)
  p_A <- mean(obs$A)
  p_S <- mean(obs$S)
  p_Y <- mean(obs$Y)
  
  #Check proportion of people with selection (S) = 1 by exposure (A)
  #Save as a scalar (rather than a 1 x 1 dataframe)
  p_S_A0 <- obs %>% filter(A == 0) %>% summarise_at("S", mean) %>% as.numeric()
  p_S_A1 <- obs %>% filter(A == 1) %>% summarise_at("S", mean) %>% as.numeric()
  
  #Check proportion of people with outcome (Y) = 1 by exposure (A)
  #Save as a scalar (rather than a 1 x 1 dataframe)
  p_Y_A0 <- obs %>% filter(A == 0) %>% summarise_at("Y", mean) %>% as.numeric()
  p_Y_A1 <- obs %>% filter(A == 1) %>% summarise_at("Y", mean) %>% as.numeric()
  
  #Check proportion of people with outcome (Y) = 1 by exposure (A) among S = 1
  #Save as a scalar (rather than a 1 x 1 dataframe)
  p_Y_A0_S1 <- obs %>% filter(A == 0 & S == 1) %>% summarise_at("Y", mean) %>% 
    as.numeric()
  p_Y_A1_S1 <- obs %>% filter(A == 1 & S == 1) %>% summarise_at("Y", mean) %>% 
    as.numeric()
  
  #Look at mean U by A
  #Save as a scalar (rather than a 1 x 1 dataframe)
  mean_U_A0_all <- obs %>% filter(A == 0) %>% summarise_at("U", mean) %>% 
    as.numeric()
  mean_U_A1_all <- obs %>% filter(A == 1) %>% summarise_at("U", mean) %>% 
    as.numeric()
  
  #Look at mean U by A among S = 1
  #Save as a scalar (rather than a 1 x 1 dataframe)
  mean_U_A0_S1 <- obs %>% filter(A == 0 & S == 1) %>% 
    summarise_at("U", mean) %>% as.numeric()
  mean_U_A1_S1 <- obs %>% filter(A == 1 & S == 1) %>% 
    summarise_at("U", mean) %>% as.numeric()
  
  #---- Look at associations and return results ----
  #Check ORs for A and Y (whole population). No bias anticipated
  no_bias_model <- glm(Y ~ A, family = binomial(), data = obs)
  
  OR_AY_all <- exp(no_bias_model$coefficients[["A"]])
  ci_95_no_bias <- exp(confint(no_bias_model, "A", level = 0.95))
  ub_OR_AY_all <- ci_95_no_bias[["97.5 %"]]
  lb_OR_AY_all <- ci_95_no_bias[["2.5 %"]]
  
  #Estimates of primary interest:  Estimated ORs for A and Y among S = 1
  #Return ORs and 95% CI limits
  selected <- obs %>% filter(S == 1)
  sel_model <- glm(Y ~ A, family = binomial(), data = selected)
  
  OR_AY_S1 <- exp(sel_model$coefficients[["A"]])
  ci_95_selected <- exp(confint(sel_model, "A", level = 0.95))
  ub_OR_AY_S1 <- ci_95_selected[["97.5 %"]]
  lb_OR_AY_S1 <- ci_95_selected[["2.5 %"]]
  
  #Values to return
  return(c("OR_AY_S1" = OR_AY_S1, "OR_AY_all" = OR_AY_all, "mean_U" = mean_U, 
           "mean_U_A1_all" = mean_U_A1_all, "mean_U_A0_all" = mean_U_A0_all, 
           "mean_U_A1_S1" = mean_U_A1_S1, "mean_U_A0_S1" = mean_U_A0_S1, 
           "p_A" = p_A, "p_Y" = p_Y, "p_S" = p_S, "p_S_A1" = p_S_A1, 
           "p_S_A0" = p_S_A0, "p_Y_A0" = p_Y_A0, "p_Y_A1" = p_Y_A1, 
           "p_Y_A0_S1" = p_Y_A0_S1, "p_Y_A1_S1" = p_Y_A1_S1, 
           "ub_OR_AY_all" = ub_OR_AY_all, "lb_OR_AY_all" = lb_OR_AY_all,  
           "ub_OR_AY_S1" = ub_OR_AY_S1, "lb_OR_AY_S1" = lb_OR_AY_S1))
}

