#---- Package Loading and Options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse", "MASS")

#Using standard notation (as opposed to scientific), rounded to three 
#decimal places
options(scipen = 999)
options(digits = 3)

#---- Defining the probability mechanisms for S and Y ----
p_S <- function(A, U){
  exp(g0 + g1*A + g2*U + g3*U*A)/(1 + exp(g0 + g1*A + g2*U + g3*U*A))
}

p_Y <- function(A, U){
  exp(b0 + b1*A + b2*U)/(1 + exp(b0 + b1*A + b2*U))
}

#---- The simulation function ----
stroke_sim <- function(){
  #Generating the data
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

