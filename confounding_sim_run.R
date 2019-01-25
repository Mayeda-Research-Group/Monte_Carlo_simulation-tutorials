/****************************************************************************/
  /* SER 2018: Using Monte Carlo Simulations for Quantitative Bias Analysis	*/
  /* June 19, 2018															*/
  /* Hailey Banack & Elizabeth Rose Mayeda									*/
  /* Questions/Comments: ermayeda@ph.ucla.edu									*/
  /*																			*/
  /* Simulation example: Bias from unmeasured confounding						*/
  /*																			*/
  /* Required files to run simulation:										*/
  /* 1) SER_simulation_workshop_2018_confounding_preamble.do:					*/ 
  /*		-Sets parameter inputs for simulation 								*/
  /* 2) SER_simulation_workshop_2018_confounding_data_gen_analysis.do: 		*/
  /*		-Generates and analyzes data and stores results						*/
  /* 3) SER_simulation_workshop_2018_confounding_run.do: 						*/ 
  /*		-Runs simulation and stores results									*/
  /****************************************************************************/
  
  /****************************************************************************/
  /* The purpose of this file is to run multiple iterations of sample 		*/
  /* and summarize and store results.											*/
  /****************************************************************************/
  
  /****************************************************************************
  The "simulate" command performs Monte Carlo-type simulations. 
The basic form of the code is: 
  "simulate <list of scalars from each iteration of sample generation>, /// 
reps(B) seed(#): do <data_generation_and_analysis_do_file>" 

This calls the data generation and analysis file B times to generate and 
analyze B data sets and stores the results from each iteration of sample 
generation as scalars. 
-Note that the first lines of code call in and assign the value from each 
scalar to a variable with the same name. 
-"reps(B)" specifies the number of iterations of sample generation. 
-"seed(#)" specifies the random number seed. 

Simulate documentation: http://www.stata.com/manuals13/rsimulate.pdf 

The rest of the code summarizes results across the B iterations of sample 
generation stores the results from the B replications in a Stata data file  
and in an Excel file(one iteration of sample generation = one row in the file). 
******************************************************************************/
  
  
  set more off
clear all

/*Save log file*/  
  capture log close
log using "SER_simulation_workshop_2018_confounding_log", replace


/*Start timer*/
  timer clear 1
timer on 1

/*Specify desired number of iterations of sample generation*/
  local B = 1000 //Number of iterations of sample generation

/*Create local variable for causal/true ORs for effect of exposure on outcome
specified in data generation and analysis file*/
  local true_OR_exposure_outcome = 1.0


/*Create a local macro variable for "scalar_X = variable_X".
We will this local macro in the simulate command. This will pull in scalars  
fromthe data generation and analysis file and store each scalar as a variable
with the same name.*/
  local simlist ""
foreach x in mean_U p_exposure p_outcome ///
  OR_exposure_Uyes OR_exposure_Uno ///
  lb_OR_exposure_Uyes ub_OR_exposure_Uyes ///
  lb_OR_exposure_Uno ub_OR_exposure_Uno { 
    local simlist "`simlist' `x'=`x'"
  }


/*Run simulation*/
  simulate `simlist', ///
reps(`B') seed(67208113): do 2_SER_simulation_workshop_2018_confounding_data_gen_analysis //use name of your data generation do file


/*Across B replications, calculate and store mean value of each variable as a scalar*/ 
  *Round scalars to two decimal places
foreach b in mean_U p_exposure p_outcome OR_exposure_Uyes OR_exposure_Uno {
  summarize `b', meanonly
  scalar mean_`b' = round(r(mean),0.01)
}


/*For each replication, generate indicator variable for whether the 95% CI for the weight status-death OR includes the causal/true OR*/
  gen covg_OR_exposure_Uyes = (ub_OR_exposure_Uyes > `true_OR_exposure_outcome' & lb_OR_exposure_Uyes < `true_OR_exposure_outcome')
                               gen covg_OR_exposure_Uno = (ub_OR_exposure_Uno > `true_OR_exposure_outcome' & lb_OR_exposure_Uno < `true_OR_exposure_outcome')
                               
                               
                               /*Across B replications, calculate and store mean value of each variable as a scalar*/ 
                               foreach b in covg_OR_exposure_Uno covg_OR_exposure_Uyes {
                               summarize `b', meanonly
                               scalar P_`b' = round(r(mean),0.001)
                               }
                               
                               
                               /*Generate plots of estimated OR and 95% CI across B simulated samples*/
                               gen sample = _n
                               
                               *Plot estimated OR adjusted for U
                               twoway (rcap lb_OR_exposure_Uyes ub_OR_exposure_Uyes sample) ///
                               (scatter OR_exposure_Uyes sample, msymbol(circle) mcolor(navy) mfcolor(edkblue)) , ///
                               ytitle(estimated OR (95% CI)) yline(1, lwidth(thick) lcolor(red)) ///
                               ylabel(0.6(.2)2.8) xscale(off) legend(off) ///
                               subtitle(Estimated adjusted OR and 95% CI from `B' simulated samples)
graph save plot_est_OR_CI_adj, replace
graph export plot_est_OR_CI_adj.png,replace

*Plot estimated OR unadjusted for U
twoway (rcap lb_OR_exposure_Uno ub_OR_exposure_Uno sample) ///
  (scatter OR_exposure_Uno sample, msymbol(circle) mcolor(navy) mfcolor(edkblue)) , ///
  ytitle(estimated OR (95% CI)) yline(1, lwidth(thick) lcolor(red)) ///
  ylabel(0.6(.2)2.8) xscale(off) legend(off) ///
  subtitle(Estimated crude OR and 95% CI from `B' simulated samples)
           graph save plot_est_OR_CI_crude, replace
           graph export plot_est_OR_CI_crude.png,replace
           
           
           /*List results */
           *Across B replications, average mean value of U, proportion of people with exposure=1, proportion of people with outcome=1
           scalar list mean_mean_U mean_p_exposure mean_p_outcome
           
           *Across B replications, average estimated OR for exposure on outcome, adjusting for U
           scalar list mean_OR_exposure_Uyes
           
           *Across B replications, average estimated OR for exposure on outcome, not adjusting for U
           scalar list mean_OR_exposure_Uno
           
           *Proportion of 95% CIs for OR for exposure on outcome that include the causal/true OR 
           scalar list P_covg_OR_exposure_Uyes P_covg_OR_exposure_Uno 
           
           /***export data to Excel***/
           export excel using SER_simulation_workshop_2018_confounding_results_each_replication, sheet("Results") sheetmodify firstrow(variables)
           
           /***save data***/
           save "SER_simulation_workshop_2018_confounding_results_each_replication.dta", replace
           
           /*End timer and display computational time (in seconds)*/
           timer off 1
           timer list 1
           