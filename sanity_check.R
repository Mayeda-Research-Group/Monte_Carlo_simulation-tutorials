#---- Package Loading and Options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse", "MASS", "ggplot2", "latex2exp")

#---- Specifying options ----
#Using standard notation (as opposed to scientific), rounded to three 
#decimal places
options(scipen = 999)
options(digits = 3)

#---- Specifying the source file ----
source("stroke_sim_par.R")
source("stroke_sim_script.R")

#---- Generating the data ----
#Creating IDs and exogenous variables
obs <- tibble("id" = seq(from = 1, to = samp_size, by = 1), 
              "A" = rbinom(n = samp_size, size = 1, p = p_A), 
              "U" = rnorm(n = samp_size, mean = 0, sd = 1), 
              "probS" = p_S(A, U), 
              "probY" = p_Y(A, U), 
              "S" = rbinom(n = samp_size, size = 1, prob = probS), 
              "Y" = rbinom(n = samp_size, size = 1, prob = probY))

#---- Look at the data ----
#Check the proportion of people with 
#exposure (A) = 1
#selection (S) = 1
#outcome (Y) = 1

grand_means <- obs %>% summarise_at(vars(A, U, S, Y), mean)

#Check the proportion of people with 
#selection (S) = 1
#outcome (Y) = 1
#and check the mean of U
#by exposure (A)

exposure_means <- obs %>% group_by(A) %>% summarise_at(vars(S, Y, U), mean)

#Check the proportion of people with outcome (Y) = 1 and mean of U by 
#exposure (A) among S = 1

selection_means <- obs %>% filter(S == 1) %>% group_by(A) %>% 
  summarize_at(vars(Y, U), mean)

#Check ORs for A and Y (whole population). No bias anticipated
pop_model <- glm(Y ~ A, family = binomial(), data = obs)
pop_ORs <- exp(pop_model$coefficients)
