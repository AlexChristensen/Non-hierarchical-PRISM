#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Generate Error for Conditions ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

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
dir.create("./error", showWarnings = FALSE)

# Load packages
library(latentFactoR); library(parallel); library(tidyverse)

# Load source functions
source("./source/generate_error.R")
source("./source/results.R")

# Set conditions
conditions <- expand.grid(
  NFAC = c(3, 5, 7),
  NVAR = c(4, 8, 12),
  LOAD = c(0.50, 0.60, 0.70),
  CORF = c(0.25, 0.50, 0.75)
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

  # Set condition number
  condition_number <- latentFactoR:::format_integer(i, n_digits)

  # Create condition data frame
  condition_df <- do.call(
    rbind.data.frame, lapply(seq_len(40), function(x){condition})
  )

  # Generate replicate results
  replicate <- mclapply(
    seq_len(100), function(j){
      return(cbind.data.frame(condition_df, generate_condition(condition)))
    }, mc.cores = 50, mc.preschedule = FALSE
  )

  # Check for failed replicates (try-error or NULL from killed workers)
  failed <- !vapply(replicate, is.data.frame, logical(1))
  if(any(failed)){
    warning(
      "Condition ", condition_number, ": failed replicates ",
      paste(which(failed), collapse = ", "), call. = FALSE
    )
    replicate <- replicate[!failed]
  }

  # Summarize replicate
  replicate_df <- do.call(rbind.data.frame, replicate)
  rm(replicate)

  # Save result
  save(replicate_df, file = paste0("./error/condition_", condition_number, ".RData"))

  # Remove replicate
  rm(replicate_df)
  gc(verbose = FALSE)
  if((i %% 12) == 0){
    system("sudo sh -c 'echo 3 > /proc/sys/vm/drop_caches'")
  }

  # Print update
  cat(paste0("\r Run ", condition_number, " of ", n_conditions, " completed."))

}

# Load results
FILES <- list.files("./error")
results <- do.call(rbind.data.frame, lapply(FILES, function(file){
  load(paste0("./error/", file))
  return(replicate_df)
}))

# Ensure numeric
results[,c(7:12)] <- lapply(results[,c(7:12)], as.numeric)

# Ensure logical
results[,13:14] <- lapply(results[,13:14], as.logical)

# Set up name
error_results <- results

# Save results
save(error_results, file = "error_results.RData")

# Print results
error_results %>%
  group_by(method, metric, value) %>%
  summarize(
    mean_cfi = mean(CFI, na.rm = TRUE),
    # sd_cfi = sd(CFI, na.rm = TRUE),
    min_cfi = min(CFI, na.rm = TRUE),
    max_cfi = max(CFI, na.rm = TRUE),
    mean_tli = mean(TLI, na.rm = TRUE),
    sd_tli = sd(TLI, na.rm = TRUE),
    min_tli = min(TLI, na.rm = TRUE),
    max_tli = max(TLI, na.rm = TRUE),
    mean_rmsea = mean(RMSEA, na.rm = TRUE),
    # sd_rmsea = sd(RMSEA, na.rm = TRUE),
    min_rmsea = min(RMSEA, na.rm = TRUE),
    max_rmsea = max(RMSEA, na.rm = TRUE),
    mean_srmr = mean(SRMR, na.rm = TRUE),
    sd_srmr = sd(SRMR, na.rm = TRUE),
    min_srmr = min(SRMR, na.rm = TRUE),
    max_srmr = max(SRMR, na.rm = TRUE),
    mean_res = mean(max_residual, na.rm = TRUE),
    sd_res = sd(max_residual, na.rm = TRUE),
    min_res = min(max_residual, na.rm = TRUE),
    max_res = max(max_residual, na.rm = TRUE),
    recoverable = mean(recoverable, na.rm = TRUE)
  ) %>% as.data.frame() %>% print(digits = 3)
