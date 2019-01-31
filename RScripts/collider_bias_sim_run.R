#******************************************************************************#
# Using Using Monte Carlo Simulations for Quantitative Bias Analysis           #
# Date                                                                         #
# Crystal Shaw, PhD Student UCLA department of Biostatistics                   #
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

#---- Visualizations ----
#Generate plots of estimated OR and 95% CI across B simulated samples
#Plot estimated OR for selected
plot_title_S1 <- paste("Memory complaints = 1: Estimated OR and 95% CI from", 
                       B, "simulated samples", 
                       sep = " ")

plot_subtitle_S1 <- paste("mean estimated OR = ", 
                          mean_results[["OR_AY_S1"]], 
                          "; 95% CI coverage probability = ", 
                          coverage_prob[["covg_OR_AY_S1"]],
                          sep = "")

CI_S1_plot <- 
  ggplot(data = sim_data, aes(x = seq(from = 1, to = nrow(sim_data), by = 1), 
                              y = OR_AY_S1)) + 
  geom_errorbar(width = 10, alpha = 0.5, color = "#A9A9A9",
                ymin = sim_data$lb_OR_AY_S1, 
                ymax = sim_data$ub_OR_AY_S1) +
  geom_point(size = 2, alpha = 0.75) +
  geom_hline(aes(yintercept = 1), size = 1.5) + 
  geom_hline(aes(yintercept = mean_results[["OR_AY_S1"]]), size = 1.5, 
             lty = 2, alpha = 0.75, color = "#4ABDAC") +
  theme_minimal() + 
  labs(title = plot_title_S1, 
       subtitle = plot_subtitle_S1) + 
  ylab("estimated OR (95% CI)") + xlab("") +
  theme(plot.title = element_text(size = 15)) + 
  theme(axis.title.y = element_text(size = 12)) + 
  scale_y_continuous(limits = c(0.2, 2), 
                     breaks = c(seq(from = 0.2, to = 2, 
                                    by = 0.2))) + 
  scale_x_continuous(breaks = NULL)

ggsave(filename = here("Plots", "plot_est_OR_CI_S1.jpeg"), 
       plot = CI_S1_plot, width = 8.25, height = 6, dpi = 300, units = "in", 
       device = 'jpeg')                            

#Generate plots of estimated OR and 95% CI across B simulated samples
#Plot estimated OR for total population
plot_title_all <- paste("Whole population: Estimated OR and 95% CI from", 
                        B, "simulated samples", 
                        sep = " ")

plot_subtitle_all <- paste("mean estimated OR = ", 
                           mean_results[["OR_AY_all"]], 
                           "; 95% CI coverage probability = ", 
                           coverage_prob[["covg_OR_AY_all"]], sep = "")

CI_all_plot <- 
  ggplot(data = sim_data, aes(x = seq(from = 1, to = nrow(sim_data), by = 1), 
                              y = OR_AY_all)) + 
  geom_errorbar(width = 10, alpha = 0.5, color = "#A9A9A9",
                ymin = sim_data$lb_OR_AY_all, 
                ymax = sim_data$ub_OR_AY_all) +
  geom_point(size = 2, alpha = 0.75) +
  geom_hline(aes(yintercept = 1), size = 1.5) + 
  geom_hline(aes(yintercept = mean_results[["OR_AY_all"]]), size = 1.5, 
             lty = 2, alpha = 0.75, color = "#4ABDAC") +
  theme_minimal() + 
  labs(title = plot_title_all, 
       subtitle = plot_subtitle_all) + 
  ylab("estimated OR (95% CI)") + xlab("") +
  theme(plot.title = element_text(size = 15)) + 
  theme(axis.title.y = element_text(size = 12)) + 
  scale_y_continuous(limits = c(0.2, 2), 
                     breaks = c(seq(from = 0.2, to = 2, 
                                    by = 0.2))) + 
  scale_x_continuous(breaks = NULL)

ggsave(filename = here("Plots", "plot_est_OR_CI_all.jpeg"), 
       plot = CI_all_plot, width = 8, height = 6, dpi = 300, units = "in", 
       device = 'jpeg')

#Make a long dataset for histogram plotting
hist_data <- sim_data %>% 
  select(mean_U_A0_all, mean_U_A1_all, mean_U_A0_S1, mean_U_A1_S1) %>% 
  gather()

#Plot mean U across the B iterations of sample generation by S and overall
#Plot for S = 1
U_hist_S1 <- 
  ggplot(data = hist_data %>% 
           filter(key == "mean_U_A0_S1" | key == "mean_U_A1_S1"), 
         aes(x = value, fill = key)) + geom_histogram() + theme_minimal() +
  scale_fill_manual(name = "", values = c("#A9A9A9", "#4ABDAC"), 
                    labels = c("A = 0","A = 1")) + xlab("mean U") + 
  ggtitle("Distributions of mean U for Selected (S = 1)")

ggsave(filename = here("Plots", "histogram_S1.jpeg"), 
       plot = U_hist_S1, width = 8, height = 6, dpi = 300, units = "in", 
       device = 'jpeg')

#Plot for whole population
U_hist_all <- 
  ggplot(data = hist_data %>% 
           filter(key == "mean_U_A0_all" | key == "mean_U_A1_all"), 
         aes(x = value, fill = key)) + geom_histogram() + theme_minimal() +
  scale_fill_manual(name = "", values = c("#A9A9A9", "#4ABDAC"), 
                    labels = c("A = 0","A = 1")) + xlab("mean U") + 
  ggtitle("Distributions of mean U for population")

ggsave(filename = here("Plots", "histogram_all.jpeg"), 
       plot = U_hist_all, width = 8, height = 6, dpi = 300, units = "in", 
       device = 'jpeg')

#---- Save results as .csv ----
write_csv(sim_data, here("Data", "collider_bias_results_each_replication.csv"))

#---- Display numerical results ----
#Check proportions of people with: 
# exposure (A) = 1
# selection (S) = 1
# outcome (Y) = 1
# outcome (Y) = 1 by exposure (A)
mean_results[c("p_A", "p_S", "p_Y", "p_Y_A1", "p_Y_A0", "mean_U")]

#Check proportion of people with outcome (Y) = 1 by exposure (A) among S = 1
mean_results[c("p_Y_A1_S1", "p_Y_A0_S1")]

#Check distribution of U by A among S = 1
mean_results[c("mean_U_A1_S1", "mean_U_A0_S1")]

#Check distribution of U by A among full sample
mean_results[c("mean_U_A1_all", "mean_U_A0_all")]

#Check ORs and 95% CI coverage in whole population (no bias anticipated)
mean_results["OR_AY_all"]
coverage_prob["covg_OR_AY_all"]

#Estimate of primary interest: Estimated ORs among S = 1
mean_results["OR_AY_S1"]

#Proportion of 95% CIs that include the causal/true OR 95% CI coverage
coverage_prob["covg_OR_AY_S1"]

#---- Display run time ----
#End timer
stop_time <- Sys.time()

stop_time - start_time



