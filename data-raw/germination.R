# ===================================================================
# DATA GENERATION SCRIPT
# This script uses a detailed simulation based on parameters from
# Haj Sghaier et al. (2022) to create the 'germination' dataset.
# ===================================================================

# Set a global seed for overall reproducibility of the script
set.seed(2025)

# --- Part 1: Original Parameters and Core Function ---

# Define All Target Parameters from the Manuscript
original_params <- data.frame(
  temperature = c(5, 10, 15, 20, 25, 30, 35),
  start_day = c(13, 7, 4, 4, 3, 4, 9),
  end_day = c(25, 20, 16, 16, 15, 16, 20),
  slope = c(0.9846, 3.6848, 4.9765, 3.8693, 6.5766, 1.1384, 0.1560),
  intercept = c(-13.627, -34.076, -28.116, -17.012, -35.621, -1.2833, -1.0571),
  r_squared = c(0.9585, 0.9501, 0.9307, 0.9524, 0.9564, 0.9237, 0.9574),
  germination_rate = c(0.88, 0.90, 0.92, 0.98, 0.94, 0.91, 0.55)
)

# Core Data Generation Function
generate_hybrid_data <- function(params) {
  set.seed(as.integer(params$temperature * 10))
  n_replicates <- 4
  n_seeds_per_dish <- 20
  n_samples_total <- n_replicates * n_seeds_per_dish

  time_data <- round(runif(n_samples_total, min = params$start_day, max = params$end_day))

  perfect_growth <- params$intercept + params$slope * time_data
  perfect_growth_variance <- var(perfect_growth)
  if (params$r_squared < 1 && params$r_squared > 0) {
    residual_variance <- perfect_growth_variance * (1 - params$r_squared) / params$r_squared
  } else {
    residual_variance <- 0
  }
  residual_std_dev <- sqrt(residual_variance)
  residuals <- rnorm(n = n_samples_total, mean = 0, sd = residual_std_dev)
  actual_growth <- perfect_growth + residuals
  actual_growth[actual_growth < 0] <- 0

  total_germinated <- round(n_samples_total * params$germination_rate)
  status_vector <- c(rep(1, total_germinated), rep(0, n_samples_total - total_germinated))
  shuffled_status <- sample(status_vector)

  replicate_vector <- rep(1:n_replicates, each = n_seeds_per_dish)
  actual_growth[shuffled_status == 0] <- 0

  data.frame(
    temperature = params$temperature,
    replicate = replicate_vector,
    time = time_data,
    growth = actual_growth,
    germinated = shuffled_status
  )
}

# Generate data for the original 7 temperature points
original_data_list <- lapply(1:nrow(original_params), function(i) generate_hybrid_data(original_params[i, ]))
original_sim_data <- do.call(rbind, original_data_list)


# --- Part 2: Simulation of Interpolated Data ---

new_temps <- c(7.5, 12.5, 17.5, 22.5, 27.5, 32.5)
interp_slope <- approx(original_params$temperature, original_params$slope, xout = new_temps)$y
interp_intercept <- approx(original_params$temperature, original_params$intercept, xout = new_temps)$y
interp_r_squared <- approx(original_params$temperature, original_params$r_squared, xout = new_temps)$y
interp_germ_rate <- approx(original_params$temperature, original_params$germination_rate, xout = new_temps)$y
interp_start_day <- approx(original_params$temperature, original_params$start_day, xout = new_temps)$y
interp_end_day <- approx(original_params$temperature, original_params$end_day, xout = new_temps)$y

interpolated_params <- data.frame(
  temperature = new_temps,
  slope = interp_slope,
  intercept = interp_intercept,
  r_squared = interp_r_squared,
  germination_rate = interp_germ_rate,
  start_day = interp_start_day,
  end_day = interp_end_day
)

interpolated_data_list <- lapply(1:nrow(interpolated_params), function(i) generate_hybrid_data(interpolated_params[i, ]))
interpolated_sim_data <- do.call(rbind, interpolated_data_list)


# --- Part 3: Combine and Save ---

# Combine both datasets for the final analysis
germination <- rbind(original_sim_data, interpolated_sim_data)

# Save the final 'germination' object to 'data/germination.rda'
usethis::use_data(germination, overwrite = TRUE)
