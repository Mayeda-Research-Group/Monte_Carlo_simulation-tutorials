#******************************************************************************#
# Using Monte Carlo Simulations for Quantitative Bias Analysis                 #
# Last updated: March 2, 2020                                                  #
# Tutorial designed by Elizabeth Rose Mayeda, Assistant Professor UCLA         #
# R code by Crystal Shaw, PhD Student UCLA Department of Biostatistics         #
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
# magrittr: contains the %<>% (double pipe) which allows efficient             #
#           reassignment/replacement of newly created objects to old ones      #
#           (https://www.rdocumentation.org/packages/magrittr/versions/1.5)    #
#******************************************************************************#

#---- Package Loading and Options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("here", "tidyverse", "magrittr")

set.seed(20190130)

#Using standard notation (as opposed to scientific), 
options(scipen = 999)

#---- Source Files ----
source(here("Confounding_simulation_tutorial", "confounding_simulation_R", 
            "RScripts", "confounding_sim_par.R"))     #call parameter file
source(here("Confounding_simulation_tutorial", "confounding_simulation_R", 
            "RScripts", "confounding_sim_script.R"))  #call simulation script

#---- Run the simulation ----  
#Start timer
start_time <- Sys.time()
  
#Create variable for causal/true ORs for effect of exposure on outcome 
#specified in data generation and analysis file 
true_OR_exposure_outcome = 1.0

#Run simulation 
#   supressMessages hides messages from fitting the logistic regression
#   %>% is a "pipe"... it feeds the output of one command to the next
#   t() transposes the results so that we have one simulation result per row of
#   the matrix
#   convert this to a dataframe so we can use the tidyverse functions
sim_data <- suppressMessages(replicate(B, confounding_sim())) %>% t() %>%
  as.data.frame()

#---- Numerical summaries ----
#Take the log of the ORs so that means are calculated correctly
sim_data %<>% mutate("log_OR_exposure_Uyes" = log(OR_exposure_Uyes), 
                     "log_OR_exposure_Uno" = log(OR_exposure_Uno))

#Across B replications, calculate and store mean value of selected variables 
#Round to two decimal places

mean_results <- sim_data %>% select("mean_U", "p_exposure", "p_outcome", 
                                    "log_OR_exposure_Uyes", 
                                    "log_OR_exposure_Uno") %>%
  colMeans() %>% round(2)

#For each replication, generate indicator variable for whether the 95% CI for 
#the OR includes the causal/true OR
#Appends this to the simulation results dataset
sim_data %<>% 
  mutate("covg_OR_exposure_Uyes" = 
           if_else(ub_OR_exposure_Uyes > true_OR_exposure_outcome & 
                     lb_OR_exposure_Uyes < true_OR_exposure_outcome, 1, 0), 
         "covg_OR_exposure_Uno" = 
           if_else(ub_OR_exposure_Uno > true_OR_exposure_outcome & 
                     lb_OR_exposure_Uno < true_OR_exposure_outcome, 1, 0))
                  
#Average coverage probabilities across B replications
coverage_prob <- sim_data %>% 
  select("covg_OR_exposure_Uyes", "covg_OR_exposure_Uno") %>%
  colMeans(., na.rm = TRUE) %>% round(3)

#---- Visualizations ----
#Generate plots of estimated OR and 95% CI across B simulated samples
#Plot estimated OR adjusted for U
plot_title_adj <- 
  paste("Estimated adjusted OR and 95% CI from", B, "simulated samples", 
        sep = " ")

plot_subtitle_adj <- 
  paste("mean estimated OR = ", 
        round(exp(mean_results[["log_OR_exposure_Uyes"]]), 2), 
        "; 95% CI coverage probability = ", 
        coverage_prob[["covg_OR_exposure_Uyes"]],
        sep = "")

#To create horizontal lines in the plot with legend
h_lines <- tibble(x = c(-Inf, Inf), "True OR" = true_OR_exposure_outcome, 
                  "Estimated OR" = 
                    round(exp(mean_results[["log_OR_exposure_Uyes"]])), 2) %>%
  gather(key = "line_type", value = "value", 
         c("True OR", "Estimated OR")) %>% 
  mutate_at("line_type", as.factor)

CI_adjusted_plot <- 
  ggplot(data = sim_data, aes(x = seq(from = 1, to = nrow(sim_data), by = 1), 
                              y = OR_exposure_Uyes)) + 
  geom_errorbar(width = 10, alpha = 0.5, color = "#A9A9A9",
                ymin = sim_data$lb_OR_exposure_Uyes, 
                ymax = sim_data$ub_OR_exposure_Uyes) +
  geom_point(size = 2, alpha = 0.75) +
  geom_line(data = h_lines, aes(x, value, linetype = line_type, 
                                color = line_type), size = 1.5) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_color_manual(values = c("#4ABDAC", "black")) + 
  theme_minimal() + 
  theme(legend.background = element_rect("light gray"), 
        legend.position = c(0.90, 0.875), legend.title = element_blank()) + 
  labs(title = plot_title_adj, 
       subtitle = plot_subtitle_adj) + 
  ylab("estimated OR (95% CI)") + xlab("") +
  theme(plot.title = element_text(size = 15)) + 
  theme(axis.title.y = element_text(size = 12)) + 
  coord_cartesian(ylim = c(0.6, 2.8)) +
  scale_y_continuous(breaks = c(seq(from = 0.6, to = 3, by = 0.2))) + 
  scale_x_continuous(breaks = NULL)

ggsave(filename = here("Confounding_simulation_tutorial", 
                       "confounding_simulation_R", 
                       "Plots", "plot_est_OR_CI_adj.jpeg"), 
       plot = CI_adjusted_plot, width = 8, height = 6, dpi = 300, units = "in", 
       device = 'jpeg')                            

#Generate plots of estimated OR and 95% CI across B simulated samples
#Plot estimated OR unadjusted for U
plot_title_unadj <- 
  paste("Estimated crude OR and 95% CI from", B, "simulated samples", 
        sep = " ")

plot_subtitle_unadj <- paste("mean estimated OR = ", 
                           round(exp(mean_results[["log_OR_exposure_Uno"]]), 3), 
                           "; 95% CI coverage probability = ", 
                           coverage_prob[["covg_OR_exposure_Uno"]],
                           sep = "")

h_lines <- tibble(x = c(-Inf, Inf), "True OR" = true_OR_exposure_outcome, 
                  "Estimated OR" = 
                    round(exp(mean_results[["log_OR_exposure_Uno"]]), 2)) %>%
  gather(key = "line_type", value = "value", 
         c("True OR", "Estimated OR")) %>% 
  mutate_at("line_type", as.factor)

CI_unadjusted_plot <- 
  ggplot(data = sim_data, aes(x = seq(from = 1, to = nrow(sim_data), by = 1), 
                              y = OR_exposure_Uno)) + 
  geom_errorbar(width = 10, alpha = 0.5, color = "#A9A9A9",
                ymin = sim_data$lb_OR_exposure_Uno, 
                ymax = sim_data$ub_OR_exposure_Uno) +
  geom_point(size = 2, alpha = 0.75) +
  geom_line(data = h_lines, aes(x, value, linetype = line_type, 
                                color = line_type), size = 1.5) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_color_manual(values = c("#4ABDAC", "black")) + 
  theme_minimal() + 
  theme(legend.background = element_rect("light gray"), 
        legend.position = c(0.90, 0.875), legend.title = element_blank()) + 
  labs(title = plot_title_unadj, 
       subtitle = plot_subtitle_unadj) + 
  ylab("estimated OR (95% CI)") + xlab("") +
  theme(plot.title = element_text(size = 15)) + 
  theme(axis.title.y = element_text(size = 12)) + 
  coord_cartesian(ylim = c(0.6, 2.8)) +
  scale_y_continuous(breaks = c(seq(from = 0.6, to = 2.8, by = 0.2))) + 
  scale_x_continuous(breaks = NULL)

ggsave(filename = here("Confounding_simulation_tutorial", 
                       "confounding_simulation_R", 
                       "Plots", "plot_est_OR_CI_crude.jpeg"), 
       plot = CI_unadjusted_plot, width = 8, height = 6, dpi = 300, 
       units = "in", device = 'jpeg')     

#---- Save results as .csv ----
write_csv(sim_data, here("Confounding_simulation_tutorial", 
                         "confounding_simulation_R", 
                         "Data", "confounding_results_each_replication.csv"))

#---- Display numerical results ----
#Across B replications, average mean value of U, 
#proportion of people with exposure = 1, proportion of people with outcome = 1
mean_results[c("mean_U", "p_exposure", "p_outcome")]

#Across B replications, average estimated OR for exposure on outcome, 
#adjusting for U
mean_OR_exposure_Uyes = exp(mean_results["log_OR_exposure_Uyes"]) #Save value
names(mean_OR_exposure_Uyes) = c("mean_OR_exposure_Uyes")         #Name value
mean_OR_exposure_Uyes                                             #Display value

#Across B replications, average estimated OR for exposure on outcome, 
#not adjusting for U
mean_OR_exposure_Uno = exp(mean_results["log_OR_exposure_Uno"]) #Save value
names(mean_OR_exposure_Uno) = c("mean_OR_exposure_Uno")         #Name value
mean_OR_exposure_Uno                                            #Display value

#Proportion of 95% CIs for OR for exposure on outcome that include the 
#causal/true OR 

coverage_prob

#End timer
stop_time <- Sys.time()

#---- Display run time ----
stop_time - start_time
           
