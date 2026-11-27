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
dir.create("./conditions", showWarnings = FALSE)

# Load packages
library(psychNetwoRk); library(parallel); library(tidyverse)

# Load source functions
source("./source/generate_condition.R")

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

# Remove row names
row.names(conditions) <- NULL

# Number of conditions
n_conditions <- nrow(conditions)
n_digits <- latentFactoR:::digits(n_conditions) - 1

# Initialize results
results <- vector("list", n_conditions)

# Loop over results
for(i in seq_len(n_conditions)){

  # Obtain condition
  condition <- conditions[i,]

  # Set condition number
  condition_number <- latentFactoR:::format_integer(i, n_digits)

  # Generate replicate results
  replicate <- mclapply(
    seq_len(100), function(j){
      generate_condition(condition)
    }, mc.cores = 50, mc.preschedule = FALSE
  )

  # Save condition
  save(replicate, file = paste0("./conditions/condition_", condition_number, ".RData"))

  # Remove replicate
  rm(replicate)
  gc(verbose = FALSE)
  if((i %% 12) == 0){
    system("sudo sh -c 'echo 3 > /proc/sys/vm/drop_caches'")
  }

  # Print update
  cat(paste0("\r Run ", condition_number, " of ", n_conditions, " completed."))

}
