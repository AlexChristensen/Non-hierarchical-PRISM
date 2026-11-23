#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Partial Correlations Distribution ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Load packages
library(latentFactoR); library(parallel); library(tidyverse)

# Load source functions
source("./source/distributions.R")

# Set conditions
conditions <- expand.grid(
  NFAC = c(1, 2, 5, 8),
  NVAR = c(4, 8, 12),
  LOAD = c(0.50, 0.70),
  CROSS = c(0.05, 0.075),
  CORF = c(0.00, 0.30, 0.60)
)

# Remove extra unidimensional conditions
conditions <- conditions[
  !((conditions$NFAC == 1) & (conditions$CORF %in% c(0.30, 0.60))),
]

# Set distributions
distributions <- c("beta", "exp", "gamma", "lnorm", "norm", "weibull")

# Set number of conditions
n_conditions <- nrow(conditions)

# Set replicates list
replicates <- vector("list", n_conditions)

# Loop over conditions
for(i in seq_len(n_conditions)){

  # Set condition
  condition <- conditions[i,]

  # Set up condition data frame
  condition_df <- do.call(rbind.data.frame, lapply(seq_len(6), function(x){condition}))

  # Loop over replicates
  replicates[[i]] <- do.call(
    rbind.data.frame, mclapply(
      seq_len(100), function(j){

        # Simulate structure
        structure <- simulate_factors(
          factors = condition$NFAC, variables = condition$NVAR,
          loadings = condition$LOAD, cross_loadings = condition$CROSS,
          correlations = condition$CORF, sample_size = 100
        )

        # Obtain partial correlation matrix
        P <- EGAnet:::cor2pcor(structure$population_correlation)

        # Set lower triangle
        P_lower <- abs(P[lower.tri(P)])

        # Loop over distributions
        output <- lapply(
          distributions, function(distr){
            get_parameters(x = P_lower, distr = distr)
          })

        # Collapse into data frame
        output_df <- do.call(rbind.data.frame, output); row.names(output_df) <- NULL

        # Add replicate
        combined_df <- cbind.data.frame(condition_df, output_df)
        combined_df$REPLICATE <- j

        # Return data frame
        return(combined_df)

      }, mc.cores = 10, mc.preschedule = FALSE
    )
  )

  # Set condition number
  conditon_number <- EGAnet:::format_integer(i, 1)

  # Print progress
  cat(
    paste0(
      "\r Condition ", conditon_number,
      " of ", n_conditions, " completed."
    )
  )

}

# Compile lists
replicate_df <- do.call(rbind.data.frame, replicates)

# Save data frame
save(replicate_df, file = "replicate_df.RData")

# Checkpoint ----

# Load data frame
load(file = "replicate_df.RData")

# Set distributions
distributions <- c("beta", "exp", "gamma", "lnorm", "norm", "weibull")
n_distributions <- length(distributions)

# Check for maximum log-likelihood
output <- replicate_df %>%
  group_by(NFAC, NVAR, LOAD, CROSS, CORF, REPLICATE) %>%
  summarize(logLik = which.max(logLik))

# Check table
frequency <- tabulate(output$logLik)

# Print proportions
structure(
  round(frequency / sum(frequency), 3),
  names = distributions
)

# Compute minimal differences between likelihoods
differences <- matrix(
  0, nrow = n_distributions, ncol = n_distributions,
  dimnames = list(distributions, distributions)
)

# Loop over distributions
for(i in seq_len(n_distributions)){

  # Store value
  first <- replicate_df$logLik[replicate_df$distr == distributions[i]]

  # Compute differences with second value
  for(j in seq_len(n_distributions)){
    differences[i,j] <- round(mean((first - replicate_df$logLik[replicate_df$distr == distributions[j]]) > 2, na.rm = TRUE), 3)
  }

}

# Rename distributions
colnames(differences) <- row.names(differences) <- c(
  "Beta", "Exponential", "Gamma", "Log-Normal", "Normal", "Weibull"
)

# Construct heatmap
ll_heat <- EGAnet:::ggheatmap(differences, type = "full", color = "lightgrey") +
  scale_fill_gradient2(low = "white", high = "#5DC09A") +
  geom_text(aes(label = EGAnet:::format_decimal(Values, 3))) +
  labs(
    x = "Worse Fitting Distribution",
    y = expression("Better Fitting Distribution (" * Delta * "LL" > 2 * ")"),
    title = "Pairwise Log-Likelihood Comparison",
    fill = "Proportion"
  ) +
  theme(
    panel.background = element_blank(), axis.ticks = element_blank(),
    legend.title = element_text(hjust = 0.5),
    plot.title = element_text(size = 12, hjust = 0.5),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
  ) +
  guides(fill = guide_colorbar(ticks.colour = "darkgrey"))

# Print plot
ll_heat

# Check margins
margin_df <- data.frame(
  better_than = rowSums(differences),
  worse_than = colSums(differences),
  ratio = rowSums(differences) / colSums(differences)
)

# Print
print(margin_df, digits = 2)

# Summary of parameters
weibull_df <- replicate_df[replicate_df$distr == "weibull",]

# Remove extremes
weibull_df <- weibull_df[weibull_df$param1 < 10,]

# Shape
hist(weibull_df$param1, breaks = "FD")
weibull_df %>%
  group_by(NFAC, LOAD) %>%
  reframe(tibble::as_tibble_row(as.list(summary(param1))))

# Scale
hist(weibull_df$param2, breaks = "FD")
weibull_df %>%
  group_by(NFAC, LOAD) %>%
  reframe(tibble::as_tibble_row(as.list(summary(param2))))
