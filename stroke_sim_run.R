#---- Package Loading and Options ----
if (!require("pacman")) 
  install.packages("pacman", repos='http://cran.us.r-project.org')

p_load("tidyverse", "MASS", "ggplot2", "latex2exp")

#Using standard notation (as opposed to scientific), rounded to three 
#decimal places
options(scipen = 999)
options(digits = 3)

set.seed(4418)

#---- Specify the parameter file ----
source("stroke_sim_par.R")
source("stroke_sim_script.R")

#---- Performing the simulation ----
#suppressMessages command silences the output from fitting the GLM
suppressMessages(conf_ints <- as_data_frame(t(replicate(1000, stroke_sim()))))
colnames(conf_ints) <- c("A", "L", "U")

#Finding the coverage probability and average estimated OR
avg_OR <- conf_ints %>% summarize_at("A", mean)
coverage <- conf_ints %>% mutate(capture = L < 1 & U > 1) %>%  
  summarize_at("capture", mean)

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
  

ggsave(filename = "CI_95_plot.jpeg", plot = CI_plot, width = 8, height = 6, 
       dpi = 300, units = "in", device = 'jpeg')



