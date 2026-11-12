#%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Generate Conditions ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Make sure environment has single threaded applications
Sys.setenv(BLIS_NUM_THREADS = "1")
Sys.setenv(AOCL_NUM_THREADS = "1")
Sys.setenv(OMP_NUM_THREADS = "1")
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)

# Create directory
dir.create("./results", showWarnings = FALSE)

# Load packages
library(psychNetwoRk); library(parallel); library(tidyverse)

# Load source functions
source("./source/methods.R")
source("./source/results.R")

# Set conditions
conditions <- expand.grid(
  BLOCKS = c(2, 4, 8),
  NODES = c(4, 8, 12),
  WITHIN = c(0.85, 0.95),
  BETWEEN = c(0.20, 0.40),
  DIFFUSION = c(0.30, 0.50),
  N = c(500, 1000, 2500, 5000)
)

# Number of conditions
n_conditions <- nrow(conditions)
n_digits <- latentFactoR:::digits(n_conditions) - 1

# Initialize results
results <- vector("list", n_conditions)

# Loop over results
for(i in seq_len(n_conditions)){

  # Obtain condition
  condition <- conditions[i,]

  # Set up condition data frame
  condition_df <- do.call(
    rbind.data.frame, lapply(seq_len(10), function(x){condition})
  )

  # Set condition number
  condition_number <- latentFactoR:::format_integer(i, n_digits)

  # Load condition
  load(file = paste0("./conditions/condition_", condition_number, ".RData"))

  # Generate replicate results
  replicate_results <- mclapply(
    replicate, function(x){

      # Obtain result for no error
      no_error_results <- obtain_results(x$no_error, condition_df)
      no_error_results$ERROR <- FALSE

      # Obtain result for error
      error_results <- obtain_results(x$error, condition_df)
      error_results$ERROR <- TRUE

      # Combine results
      return(rbind.data.frame(no_error_results, error_results))

    }, mc.cores = 50, mc.preschedule = FALSE
  )

  # Save condition
  save(replicate_results, file = paste0("./results/condition_", condition_number, ".RData"))

  # Print update
  cat(paste0("\r Run ", condition_number, " of ", n_conditions, " completed."))

}