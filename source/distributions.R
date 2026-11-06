# Load packages
library(fitdistrplus); library(compiler)

# Function get distribution parameters
get_parameters <- function(x, distr)
{

  # Attempt
  attempt <- try(
    fitdist(data = x, distr = distr, method = "mle"), silent = TRUE
  )

  # Return estimates
  if(is(attempt, "try-error")){
    return(data.frame(distr = distr, param1 = NA, param2 = NA, logLik = NA))
  }else(
    return(
      data.frame(
        distr = distr, param1 = attempt$estimate[1],
        param2 = ifelse(
          length(attempt$estimate) > 1, attempt$estimate[2], NA
        ), logLik = logLik(attempt)
      )
    )
  )


}

# Compile function
get_parameters <- cmpfun(get_parameters)
