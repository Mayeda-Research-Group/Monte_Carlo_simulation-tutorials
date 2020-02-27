This repository provides Stata and R code for quantifying bias in estimation of the effect of a binary exposure (A) on a binary outcome (Y) arising from conditioning on a collider (S), where S is a common effect of the exposure and an unmeasured determinant of the outcome (U). 
Input parameters can be modified by the user.  

**Steps for running the simulation code in Stata**
1) For each causal structure of interest, save two Stata do-files: a. the data generation and analysis file and b. the run simulation file.
-The data generation and analysis file generates one data set, analyzes the generated data, and stores the parameter estimates from analyses as scalars. 
-The run simulation file calls the data generation and analysis file B times to generate and analyze B data sets, stores the results from each iteration of data set generation, summarizes the results across the B data sets, and exports the results for each iteration of sample generation in a Statafile (one sample = one row in the file).

2) Update the parameter values in "Step 1: Set parameters" of the data generation and analysis file to reflect the desired input paramters for investigation. 
-For example, to modify the  simulation so there is a causal effect of A on Y, modify local parameter "b1" in the data generation and analysis file. 

3) Run B iterations of data generation and summarize results from statistical models estimating the effect of A on Y across the  B iterations of sample generation by running the run simulation Stata do-file.
   
The data-generating process is outlined graphically in the file Collider_bias_simulation_data_generation.pdf.

**Steps for running the simulation code in R**
1) For each causal structure of interest, save two RScripts: a. the data generation and analysis file and b. the run simulation file.
-The data generation and analysis file generates one data set, analyzes the generated data, and stores the parameter estimates from analyses as scalars. 
-The run simulation file calls the data generation and analysis file B times to generate and analyze B data sets, stores the results from each iteration of data set generation, summarizes the results across the B data sets, and exports the results for each iteration of sample generation in a .csv (one sample = one row in the file).

2) Update the parameter values in "Step 1: Set parameters" of the data generation and analysis file to reflect the desired input paramters for investigation. 
-For example, to modify the  simulation so there is a causal effect of A on Y, modify local parameter "b1" in the data generation and analysis file. 

3) Run B iterations of data generation and summarize results from statistical models estimating the effect of A on Y across the  B iterations of sample generation by running the run simulation Stata do-file.
   
The data-generating process is outlined graphically in the file Collider_bias_simulation_data_generation.pdf.
