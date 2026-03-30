# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
#library(tarchetypes)

# Set target options:
tar_option_set(
  packages = c("glmnet","tidyverse","dummy") # Packages that your targets need for their tasks.
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source('./Helper_functions.R')
# tar_source("other_functions.R") # Source other scripts as needed.


# --- CONFIGURATION ---
filtering <- TRUE   # Change to FALSE if you want to skip association filtering
alpha <- c(0.3, 0.4, 0.5, 0.7, 0.9, 1.0)
output<-'Example/Output'
# After checking drpoc_csv:
run_RN    <- TRUE   # Normalization switch 
run_IMP   <- TRUE   # Imputation switch


# --- Replace the target list below with your own ---
list(
  
  # 1. Load the raw data and phenotype file
  tar_target(name = data_file,  "Example/Example_Data/data.csv",  format = "file"),
  tar_target(name = phe_file, "Example/Example_Data/pheno.csv", format = "file"),
  tar_target(name = data,  read.csv(data_file)),
  tar_target(name = phefile, read.csv(phe_file)),
  
  # 2. Pre-processing steps:
  ## Check for normality and missingness
  tar_target(name = dproc, command = check_NM(data, id='id')),
  
  ## Save check results as csv
  tar_target(dproc_csv, command = {
    
    # Define the file path relative to the project root
    file_path <- paste0(output,"/check_results.csv")
    
    # Create directory if it doesn't exist 
    if (!dir.exists(output)) dir.create(output, recursive = TRUE)  
    
    write.csv(dproc, file_path, row.names = FALSE)
    file_path },
    format = "file"),
  
  ## Rank-normalize and half minimum imputation (optional)
  tar_target(name = data_proc, 
             command = RN_IMP(data, id = 'id', miss_prop = 0.25,
                              dproc, RN = run_RN, IMP = run_IMP), 
             format = 'rds'),
  
  # 3. Define analysis parameters
  tar_target(name = outcomes, command = c("htn_incident")), 
  tar_target(name = cat_covs, command = c("sex", "race")),
  tar_target(name = num_covs, command = c("age", "bmi")),
  
  ## Filtering non-associated features (optional)
  tar_target(filtered_res,  
             command = if (filtering) {
               filter(data_proc, phefile = phefile, id='id', 
                      outcome = outcomes, categorical_covariates = cat_covs, 
                      numeric_covariates = num_covs) } else {
                 NULL # If filtering is FALSE, return NULL so ENselect knows not to look for results
                        }),
  
  # 4. Run elastic net
  tar_target(name = selected_features, 
             command = ENselect(data_proc = data_proc, phefile = phefile,
                                filtered_res = filtered_res, alpha = alpha,
                                id = "id", outcome = outcomes,
                                categorical_covariates = cat_covs,
                                numeric_covariates = num_covs)),
  
  # Save results to csv
  tar_target(save_selected_reatues,  command = {
    
    ENres_path <- paste0(output,"/Selected_features.csv")
    
    write.csv(selected_features, ENres_path, row.names = FALSE)
    ENres_path }, 
    format = "file")
  
)


