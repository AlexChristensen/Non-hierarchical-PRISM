#%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Generate Conditions ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Provide password once, then keep sudo alive for as long as this R process lives
system("sudo -v")
system(sprintf(
  "( while kill -0 %d 2>/dev/null; do sudo -n -v; sleep 60; done ) >/dev/null 2>&1 &",
  Sys.getpid()))
system("sudo sh -c 'echo 3 > /proc/sys/vm/drop_caches'")


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
  BLOCKS = c(1, 3, 5, 7),
  NODES = c(4, 8, 12),
  WITHIN = c(0.85, 0.95),
  BETWEEN = c(0.10, 0.20),
  DIFFUSION = c(0.10, 0.20)
)

# Set for unidimensional conditions
## Remove WITHIN
conditions <- conditions[
  !((conditions$BLOCKS == 1) & (conditions$WITHIN == 0.85)),
]
## Remove BETWEEN
conditions <- conditions[
  !((conditions$BLOCKS == 1) & (conditions$BETWEEN == 0.20)),
]
## Remove DIFFUSION
conditions <- conditions[
  !((conditions$BLOCKS == 1) & (conditions$DIFFUSION == 0.20)),
]
## Set BETWEEN and DIFFUSION to zero
conditions$BETWEEN[conditions$BLOCKS == 1] <- 0
conditions$DIFFUSION[conditions$BLOCKS == 1] <- 0

# Number of conditions
n_conditions <- nrow(conditions)
n_digits <- latentFactoR:::digits(n_conditions) - 1

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

  # Set up membership
  membership <- rep(seq_len(condition$BLOCKS), each = condition$NODES)

  # Generate replicate results
  replicate_results <- mclapply(
    replicate, function(x){

      # Obtain base results
      base_results <- do.call(
        rbind.data.frame, lapply(x, function(y){
          obtain_results(y$base, membership, condition_df)
        })
      )

      # Add error
      base_results$ERROR <- FALSE

      # Obtain error results
      error_results <- do.call(
        rbind.data.frame, lapply(x, function(y){
          obtain_results(y$error, membership, condition_df)
        })
      )

      # Add error
      error_results$ERROR <- TRUE

      # Combine results
      return(rbind.data.frame(base_results, error_results))

    }, mc.cores = 50, mc.preschedule = FALSE
  )

  # Check for failed replicates (try-error or NULL from killed workers)
  failed <- !vapply(replicate_results, is.data.frame, logical(1))
  if(any(failed)){
    warning(
      "Condition ", condition_number, ": failed replicates ",
      paste(which(failed), collapse = ", "), call. = FALSE
    )
  }

  # Save condition
  save(replicate_results, file = paste0("./results/condition_", condition_number, ".RData"))

  # Remove replicate
  rm(replicate_results)
  gc(verbose = FALSE)
  if((i %% 12) == 0){
    system("sudo sh -c 'echo 3 > /proc/sys/vm/drop_caches'")
  }

  # Print update
  cat(paste0("\r Run ", condition_number, " of ", n_conditions, " completed."))

}

# Load results
FILES <- list.files("./results")
sbm_results <- do.call(rbind.data.frame, lapply(FILES, function(file){
  load(paste0("./results/", file))
  return(do.call(rbind.data.frame, replicate_results))
}))

# Save results
save(sbm_results, file = "sbm_results.RData")

sbm_results$NDIM <- as.numeric(sbm_results$NDIM)
sbm_results$BLOCKS <- as.numeric(sbm_results$BLOCKS)

# Summary of results
sbm_results %>%
  filter(!grepl("Error", METHOD)) %>%
  group_by(METHOD, BLOCKS, NODES, ERROR) %>%
  summarize(
    PC = mean(NDIM == BLOCKS, na.rm = TRUE),
    MAE = mean(abs(NDIM - BLOCKS), na.rm = TRUE),
    MBE = mean(NDIM - BLOCKS, na.rm = TRUE)
  ) %>% as.data.frame() %>% print(digits = 2)
