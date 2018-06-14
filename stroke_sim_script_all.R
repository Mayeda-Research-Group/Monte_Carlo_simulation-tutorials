#---- Package Loading and Options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse", "MASS")

#Using standard notation (as opposed to scientific), rounded to three 
#decimal places
options(scipen = 999)
options(digits = 3)

set.seed(4418)

#---- Step 1: Setting the parameters ----
#Specify the prevalence of A (exposure)
p_A <- 0.2

#Parameters for the odds of S (selection)
g0 <- log(0.10/(1 - 0.10)) #log odds of S for ref group (A = 0 and U = 0)
g1 <- log(5.0) #log OR for effect of A on log odds of selection (OR = 5.0)
g2 <- log(5.0) #log OR for effect of U on log odds of selection (OR = 5.0)
g3 <- log(1.0)  #log OR for interaction between A and U on    (OR = 5.0)

#Parameters for the odds of Y (outcome)
#Is S = 0 for reference group correct?
b0 <- log(0.05/(1 - 0.05)) #log odds of Y for ref group (A = 0, U = 0, and S = 0)
b1 <- log(1.0) #log OR for effect of A on log odds of Y (OR = 1.0)
b2 <- log(5.0) #log OR for effect of U on log odds of Y (OR = 5.0)

#---- Step 2: Defining the probability mechanisms for S and Y ----
p_S <- function(A, U){
  exp(g0 + g1*A + g2*U + g3*U*A)/(1 + exp(g0 + g1*A + g2*U + g3*U*A))
}

p_Y <- function(A, U){
  exp(b0 + b1*A + b2*U)/(1 + exp(b0 + b1*A + b2*U))
}

#---- Step 3: Generating the data ----
#Creating IDs and exogenous variables
obs <- tibble("id" = seq(from = 1, to = 5000, by = 1), 
              "A" = rbinom(n = 5000, size = 1, p = p_A), 
              "U" = rnorm(n = 5000, mean = 0, sd = 1)) %>%
#Generating S and Y
  mutate("probS" = p_S(A, U), "probY" = p_Y(A, U)) %>%
  mutate("S" = rbinom(n = 5000, size = 1, prob = probS), 
         "Y" = rbinom(n = 5000, size = 1, prob = probY))

#---- Step 4: Look at the data ----
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

#---- Look at the associations and store the results ----
#Check ORs for A and Y (whole population). No bias anticipated
pop_model <- glm(Y ~ A, family = binomial(), data = obs)
pop_ORs <- exp(pop_model$coefficients)

#Estimates of primary interest:  Estimated ORs for A and Y among S = 1
#Store ORs and 95% CI limits
selected <- obs %>% filter(S == 1)
sel_model <- glm(Y ~ A, family = binomial(), data = selected)
sel_ORs <- exp(sel_model$coefficients)
ci_95 <- exp(confint(sel_model, "A", level = 0.95))




