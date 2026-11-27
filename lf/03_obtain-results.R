#%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Obtain Results ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%#

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
library(latentFactoR); library(parallel); library(tidyverse)

# Load source functions
source("./source/methods.R")
source("./source/results.R")

# Set conditions
conditions <- expand.grid(
  NFAC = c(1, 3, 5, 7),
  NVAR = c(4, 8, 12),
  LOAD = c(0.50, 0.60, 0.70),
  CORF = c(0.30, 0.60)
)

# Set for unidimensional conditions
conditions <- conditions[
  !((conditions$NFAC == 1) & (conditions$CORF %in% c(0.60))),
]

# Number of conditions
n_conditions <- nrow(conditions)
n_digits <- latentFactoR:::digits(n_conditions) - 1

# Number of methods (rows returned per dataset by `obtain_results`)
n_methods <- 10

# Loop over results
for(i in 21){# seq_len(n_conditions)){

  # Obtain condition
  condition <- conditions[i,]

  # Set up condition data frame
  condition_df <- do.call(
    rbind.data.frame, lapply(seq_len(n_methods), function(x){condition})
  )

  # Set condition number
  condition_number <- latentFactoR:::format_integer(i, n_digits)

  # Load condition
  load(file = paste0("./conditions/condition_", condition_number, ".RData"))

  # Generate replicate results
  replicate_results <- mclapply(
    seq_along(replicate), function(j){

      # Obtain result
      result <- do.call(
        rbind.data.frame, lapply(replicate[[j]], obtain_results, condition_df)
      )

      # Remove row names
      row.names(result) <- NULL

      # Return result
      return(result)

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
  rm(replicate, replicate_results)
  gc(verbose = FALSE)
  if((i %% 12) == 0){
    system("sudo sh -c 'echo 3 > /proc/sys/vm/drop_caches'")
  }

  # Print update
  cat(paste0("\r Run ", condition_number, " of ", n_conditions, " completed."))

}

# Load results
FILES <- list.files("./results")
lf_results <- do.call(rbind.data.frame, lapply(FILES, function(file){
  load(paste0("./results/", file))
  return(do.call(rbind.data.frame, replicate_results))
}))

# Save results
save(lf_results, file = "lf_results.RData")

# Summary of results
lf_results %>%
  group_by(METHOD, NFAC, ERROR, OVERLAP) %>%
  summarize(
    PC = mean(NDIM == NFAC, na.rm = TRUE),
    MAE = mean(abs(NDIM - NFAC), na.rm = TRUE),
    MBE = mean(NDIM - NFAC, na.rm = TRUE),
    ARI = mean(ARI, na.rm = TRUE)
  ) %>% as.data.frame() %>% print(digits = 2)
