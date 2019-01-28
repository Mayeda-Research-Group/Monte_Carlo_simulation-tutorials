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
# The purpose of this file is to run multiple iterations of sample generation  # 
# and summarize and store the results.                                         #
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
# latex2exp: converts latex to expressions that can be included in plot text   #
#           (https://www.rdocumentation.org/packages/latex2exp/versions/0.4.0) #
#******************************************************************************#

#---- Package Loading and Options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("here", "tidyverse", "latex2exp")

#Using standard notation (as opposed to scientific)
options(scipen = 999)

#---- Set the seed ----
set.seed(4418)

#---- Source Files ----
source(here("RScripts", "collider_bias_sim_par.R"))     #parameter file
source(here("RScripts", "collider_bias_sim_script.R"))  #simulation script

#---- Run the simulation ----  
#Start timer
start_time <- Sys.time()

#Create variable for causal/true ORs for A on Y specified in data generation 
#and analysis file 
causal_OR_AY = 1.0

#Run simulation 
#   supressMessages hides messages from fitting the logistic regression
#   %>% is a "pipe"... it feeds the output of one command to the next
#   t() transposes the results so that we have one simulation result per row of
#   the matrix
#   convert this to a dataframe so we can use the tidyverse functions
sim_data <- suppressMessages(replicate(B, collider_sim())) %>% t() %>%
  as.data.frame()

#---- Numerical summaries ----
#Across B replications, calculate and store mean value of selected variables 
#Round to three decimal places

mean_results <- sim_data %>% 
  select("OR_AY_S1", "OR_AY_all", "mean_U", "mean_U_A1_all", "mean_U_A0_all", 
         "mean_U_A1_S1", "mean_U_A0_S1", "p_A", "p_Y", "p_S", "p_S_A1", 
         "p_S_A0", "p_Y_A0", "p_Y_A1", "p_Y_A0_S1", "p_Y_A1_S1") %>%
  colMeans() %>% round(3)

#For each replication, generate indicator variable for whether the 95% CI for 
#the OR includes the causal/true OR
#Appends this to the simulation results dataset
sim_data %<>% 
  mutate("covg_OR_AY_S1" = 
           if_else(ub_OR_AY_S1 > causal_OR_AY & 
                     lb_OR_AY_S1 < causal_OR_AY, 1, 0), 
         "covg_OR_AY_all" = 
           if_else(ub_OR_AY_all > causal_OR_AY & 
                     lb_OR_AY_all < causal_OR_AY, 1, 0))

#Average coverage probabilities across B replications
coverage_prob <- sim_data %>% 
  select("covg_OR_AY_S1", "covg_OR_AY_all") %>% 
  colMeans(., na.rm = TRUE) %>% round(3)

#---- Produce and save plots ----
plot_title <- paste("Average OR for Selected = ", 
                    round(avg_OR, 2), "; ",   
                    round(coverage*100, 0), 
                    "% of 95% CI include the true OR = 1", 
                    sep = "")
CI_plot <- 
  ggplot(data = conf_ints, aes(x = seq(from = 1, to = nrow(conf_ints), by = 1), 
                               y = A)) + 
  geom_errorbar(width = 10, alpha = 0.5, color = "#A9A9A9",
                ymin = conf_ints$L, ymax = conf_ints$U) +
  geom_point(size = 2, alpha = 0.75) +
  geom_hline(aes(yintercept = 1), size = 1.5) + 
  geom_hline(aes(yintercept = avg_OR), size = 1.5, lty = 2, alpha = 0.75, 
             color = "#4ABDAC") +
  labs(x = " ", 
       y = TeX("$\\widehat{OR}_{AY|S = 1}$")) +
  theme_minimal() + 
  ggtitle(plot_title) + 
  theme(plot.title = element_text(size = 15)) + 
  theme(axis.title.y = element_text(size = 12)) + 
  scale_y_continuous(limits = c(0.2, 1.5), 
                     breaks = c(seq(from = 0.2, to = 1.5, by = 0.2))) + 
  scale_x_continuous(breaks = NULL)
  

ggsave(filename = here("Plots", "CI_95_plot.jpeg"), 
                       plot = CI_plot, width = 8, height = 6, dpi = 300, 
                       units = "in", device = 'jpeg')

#End timer
stop_time <- Sys.time()

#---- Display run time ----
stop_time - start_time



