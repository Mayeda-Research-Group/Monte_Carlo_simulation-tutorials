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
#Is S = 0 for reference group correct?
b0 <- log(0.05/(1 - 0.05)) #log odds of Y for ref group (A = 0, U = 0, and S = 0)
b1 <- log(1.0) #log OR for effect of A on log odds of Y (OR = 1.0)
b2 <- log(5.0) #log OR for effect of U on log odds of Y (OR = 5.0)

