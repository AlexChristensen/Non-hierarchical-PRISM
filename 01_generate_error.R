#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Generate Error for Conditions ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Provide password
system("sudo sh -c 'echo 3 > /proc/sys/vm/drop_caches'")

# Make sure environment has single threaded applications
Sys.setenv(BLIS_NUM_THREADS = "1")
Sys.setenv(AOCL_NUM_THREADS = "1")
Sys.setenv(OMP_NUM_THREADS = "1")
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)

# Create directory
dir.create("./error")

# Load packages
library(psychNetwoRk); library(parallel)

# Load source functions
source("./source/generate_error.R")

# Set conditions
conditions <- expand.grid(
  BLOCKS = c(2, 4, 8),
  NODES = c(4, 8, 12),
  WITHIN = c(0.80, 0.90),
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
#
# # Load results
# FILES <- list.files("./error")
# results <- do.call(rbind.data.frame, lapply(FILES, function(file){
#   load(paste0("./error/", file))
#   return(replicate_df)
# }))
#
# # Print results
# results %>%
#   group_by(overlap, target_cfi) %>%
#   summarize(
#     mean_recoverable = mean(recoverable, na.rm = TRUE),
#     sd_recoverable = sd(recoverable, na.rm = TRUE),
#     min_recoverable = min(recoverable, na.rm = TRUE),
#     max_recoverable = max(recoverable, na.rm = TRUE),
#     mean_cfi = mean(cfi, na.rm = TRUE),
#     sd_cfi = sd(cfi, na.rm = TRUE),
#     min_cfi = min(cfi, na.rm = TRUE),
#     max_cfi = max(cfi, na.rm = TRUE),
#     mean_res = mean(max_abs_res, na.rm = TRUE),
#     sd_res = sd(max_abs_res, na.rm = TRUE),
#     min_res = min(max_abs_res, na.rm = TRUE),
#     max_res = max(max_abs_res, na.rm = TRUE)
#   ) %>% as.data.frame() %>% print(digits = 3)
