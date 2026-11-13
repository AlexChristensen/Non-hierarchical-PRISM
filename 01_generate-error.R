#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Generate Error for Conditions ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Make sure environment has single threaded applications
Sys.setenv(BLIS_NUM_THREADS = "1")
Sys.setenv(AOCL_NUM_THREADS = "1")
Sys.setenv(OMP_NUM_THREADS = "1")
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)

# Create directory
dir.create("./error", showWarnings = FALSE)

# Load packages
library(psychNetwoRk); library(parallel); library(tidyverse)

# Load source functions
source("./source/generate_error.R")

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

# Path to the saved file for a condition
condition_path <- function(i){
  paste0(
    "./error/condition_",
    latentFactoR:::format_integer(i, n_digits), ".RData"
  )
}

# Check whether a condition was already completed
# (file must exist AND load cleanly with a non-empty data frame)
condition_completed <- function(i){

  # Get file
  file <- condition_path(i)

  # No file means not completed
  if(!file.exists(file)){
    return(FALSE)
  }

  # Attempt to load (corrupt or partial files return FALSE)
  tryCatch(
    {
      load_env <- new.env()
      load(file, envir = load_env)
      is.data.frame(load_env$replicate_df) && nrow(load_env$replicate_df) > 0
    },
    error = function(e){FALSE}
  )

}

# Determine which conditions remain (allows restarting where it left off)
completed <- vapply(seq_len(n_conditions), condition_completed, logical(1))
cat(paste0(
  " Found ", sum(completed), " of ", n_conditions, " conditions completed; ",
  "resuming with ", sum(!completed), " remaining.\n"
))

# Loop over results
for(i in which(!completed)){

  # Obtain condition
  condition <- conditions[i,]

  # Set condition number
  condition_number <- latentFactoR:::format_integer(i, n_digits)

  # Create condition data frame
  condition_df <- do.call(
    rbind.data.frame, lapply(seq_len(20), function(x){condition})
  )

  # Generate replicate results
  replicate <- mclapply(
    seq_len(100), function(j){

      # Regenerate with error
      data_condition <- list()
      class(data_condition) <- "try-error"

      # Loop until good condition
      while(is(data_condition, "try-error")){
        data_condition <- try(generate_condition(condition), silent = TRUE)
      }

      # Combine results
      return(cbind.data.frame(condition_df, data_condition))

    }, mc.cores = 50, mc.preschedule = FALSE
  )

  # Summarize replicate
  replicate_df <- do.call(rbind.data.frame, replicate)
  rm(replicate)

  # Save result (write to hidden temporary file, then rename so that
  # stopping mid-save never leaves a partial condition file behind)
  temporary_file <- paste0("./error/.condition_", condition_number, ".tmp")
  save(replicate_df, file = temporary_file)
  file.rename(temporary_file, condition_path(i))

  # Remove replicate
  rm(replicate_df)
  gc(verbose = FALSE)

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
results[,10:14] <- lapply(results[,10:14], as.numeric)

# Ensure logical
results[,15:17] <- lapply(results[,15:17], as.logical)

# Set up name
error_results <- results

# Save results
save(error_results, file = "error_results.RData")


# Print results
results %>%
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
    mean_walktrap = mean(walktrap, na.rm = TRUE),
    mean_louvain = mean(louvain, na.rm = TRUE),
    mean_prism = mean(prism, na.rm = TRUE)
  ) %>% as.data.frame() %>% print(digits = 3)
