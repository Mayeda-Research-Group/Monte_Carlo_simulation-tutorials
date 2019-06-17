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
#       -Generates and analyzes data and stores results                        #
#   (3) confounding_sim_run.R:                                                 #
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

#---- Package loading & options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse")

#Using standard notation (as opposed to scientific)
options(scipen = 999)

#For multiple iterations of sample generation, the seed is set in the 
#collider_bias_sim_run.R simulation file
#set.seed(4418)

#---- Custom functions ----
#Probability of exposure
prob_exp <- function(U){
  exp(g0 + g1*U)/(1 + exp(g0 + g1*U))
}

#Probability of outcome
prob_out <- function(U, exposure){
  exp(b0 + b1*exposure + b2*U)/(1 + exp(b0 + b1*exposure + b2*U))
}

#---- The simulation function ----
confounding_sim <- function(){
  
  #---- Creating the dataset ----
  # We generate a tibble to store all the data. Tibbles are just specific types 
  # of dataframes... they only print the first 10 rows in the console when 
  # called (unless more rows are specified) and list the datatype under each 
  # column name. I find them nifty for debugging issues with datatypes
  
  # What's happening in each column...
  #   "ID": Numbering the observations
  #   "U": Creating num_obs number of standard normal random variables
  #   "Pexposure": Applying the prob_exp function to the column of U values to 
  #                get a probability of exposure for each observation
  #   "exposure": Flip a coin with success probability Pexposure for each 
  #               observation
  #   "Poutcome": Applying the prob_out function to the column of U values and 
  #               column of exposure values to get a probability of outcome for 
  #               each observation
  #   "outcome": Flip a coin with success probability Poutcome for each 
  #              observation
  
  dataset <- tibble("ID" = seq(from = 1, to = num_obs, by = 1), 
                    "U" = rnorm(n = num_obs, mean = 0, sd = 1), 
                    "Pexposure" = prob_exp(U),
                    "exposure" = rbinom(n = num_obs, size = 1, 
                                        prob = Pexposure),  
                    "Poutcome" = prob_out(U, exposure), 
                    "outcome" = rbinom(n = num_obs, size = 1, 
                                       prob = Poutcome))
  
  #---- Analyze data ----
  #Check the mean value of U (should be very close, if not exactly, 0)
  #Calculate the proportion of people with exposure = 1
  #Calculate the proportion of people with outcome = 1
  
  mean_U <- mean(dataset$U)
  p_exposure <- mean(dataset$exposure) 
  p_outcome <- mean(dataset$outcome)        
  
  #---- Look at associations and return results ----
  
  #Estimated OR & 95% CI for exposure-outcome association
  #Adjusting for U (confounder)
  model_U_yes <- glm(outcome ~ exposure + U, family = binomial(), 
                     data = dataset)
  
  OR_exposure_Uyes <- exp(model_U_yes$coefficients[["exposure"]])
  ci_95_Uyes <- exp(confint(model_U_yes, "exposure", level = 0.95))
  ub_OR_exposure_Uyes <- ci_95_Uyes[["97.5 %"]]
  lb_OR_exposure_Uyes <- ci_95_Uyes[["2.5 %"]]
  
  #Estimated OR & 95% CI for exposure-outcome association
  #Without adjusting for U (confounder)
  model_U_no <- glm(outcome ~ exposure, family = binomial(), 
                    data = dataset)
  
  OR_exposure_Uno <- exp(model_U_no$coefficients[["exposure"]])
  ci_95_Uno <- exp(confint(model_U_no, "exposure", level = 0.95))
  ub_OR_exposure_Uno <- ci_95_Uno[["97.5 %"]]
  lb_OR_exposure_Uno <- ci_95_Uno[["2.5 %"]]
  
  #Values to return
  return(c("mean_U" = mean_U, "p_exposure" = p_exposure, 
           "p_outcome" = p_outcome, "OR_exposure_Uyes" = OR_exposure_Uyes, 
           "ub_OR_exposure_Uyes" = ub_OR_exposure_Uyes, 
           "lb_OR_exposure_Uyes" = lb_OR_exposure_Uyes,  
           "OR_exposure_Uno" = OR_exposure_Uno, 
           "ub_OR_exposure_Uno" = ub_OR_exposure_Uno, 
           "lb_OR_exposure_Uno" = lb_OR_exposure_Uno))
}


