#' Title
#'
#' @param model_parameters list of model parameters
#' @param simulation_steps number of steps to run the simulation for (in model time) (if dt = 1/12, this should be 12 * years to run simulation for)
#'
#' @importFrom odin.dust odin_dust
#' @return simulated data over time as data frame
#' @export
#'
simulate_model <- function(model_parameters, simulation_steps){
  WF_PPxSero <- odin.dust::odin_dust("model.R") # imports model

  WFmodel_ppxSero <- WF_PPxSero$new(pars = model_parameters,
                                    time = 0,
                                    n_particles = 10L,
                                    n_threads = 4L,
                                    seed = 1L, deterministic = TRUE) # can also run this in stochastic mode with deterministic = FALSE

  simulated_data_long <- list()
  simulated_sum_data_long <- list()
  list_of_particle_list <- list()
  for (i in 0:simulation_steps) {
    simMean_all <- rowMeans(WFmodel_ppxSero$run(i)[-(1:(model_parameters$GPSC_no+1)),])
    simMean_matrix <- matrix(simMean_all, ncol = model_parameters$sero_no, nrow = model_parameters$GPSC_no, byrow = FALSE)
    simulated_data_long[[i+1]] <- simMean_matrix
    simSum_all <- rowSums(WFmodel_ppxSero$run(i)[-(1:(model_parameters$GPSC_no+1)),])
    simSum_matrix <- matrix(simSum_all, ncol = model_parameters$sero_no, nrow = model_parameters$GPSC_no, byrow = FALSE)
    simulated_sum_data_long[[i+1]] <- simSum_matrix
    list_of_particle_list[[i+1]] <-  WFmodel_ppxSero$run(i)[(2:(model_parameters$GPSC_no+1)),]
  }
  simulated_data_over_time_clust <- sapply(simulated_data_long, rowSums)
  simulated_data_over_time_clust_df <- as.data.frame(simulated_data_over_time_clust)
  simulated_data_over_time_clust_df
}

#' Title
#'
#' @param state values of model
#' @param observed values of data
#' @param pars parameters
#'
#' @return likelihood
#'
likelihood <- function(state, observed, pars = NULL) {
  result <- 0
  data_size <- sum(unlist(observed[as.character(1:(length(unlist(observed))-4))]))
  model_size = sum(unlist(state[-1, , drop = TRUE]))
  exp_noise <- 1e6
  data_vals <- unlist(observed[as.character(1:(length(unlist(observed))-4))])
  model_vals <- rep(0, length(unlist(observed))-4)
  data_missing <- FALSE
  for (i in 1:(length(unlist(observed))-4)){
    state_name <- paste("sum_clust", i, sep = "")
    model_vals[i] <- state[state_name, , drop = TRUE]
    if (is.na(observed[[as.character(i)]])) {
      data_missing <- TRUE
    }
  }
  models_vals_err <- model_vals + stats::rexp(n = length(model_vals), rate = exp_noise)
  if(data_missing){
    ll_obs <- 0
  }
  else{
    ll_obs <- stats::dmultinom(x = (data_vals), prob = models_vals_err/model_size, log = TRUE)
  }
  result <- ll_obs
  result
}

#' Title
#'
#' @param data data to fit model to (GPSC counts)
#' @param fixed_parameters input parameters from data, such as start population
#' @param steps_mcmc1 number of mcmc steps for first run
#' @param steps_mcmc2 number of mcmc steps for second run (should be a lot higher for a good fit ~10,000 but then needs to be run on an hpc and will take ~24 hours)
#'
#' @import mcstate
#' @importFrom odin.dust odin_dust
#' @return fitted mcmc object
#' @export
#'
fit_model_to_data <- function(data, fixed_parameters, steps_mcmc1 = 10, steps_mcmc2 = 50){
# #' @importFrom mcstate particle_filter_data pmcmc_parameters pmcmc pmcmc_thin particle_deterministic pmcmc_control
  # convert input data to mcstate object
  fitting_data <- data.frame("year" = 1:ncol(data), t(data))
  names(fitting_data) <- c("year", as.character(1:(fixed_parameters$GPSC_no)))
  fitting_data <- mcstate::particle_filter_data(data = fitting_data,
                                                       time = "year",
                                                       rate = 1 / fixed_parameters$dt,
                                                       initial_time = 0)

  WF <- odin.dust::odin_dust("model.R")
  det_filter <- mcstate::particle_deterministic$new(data = fitting_data,
                                           model = WF,
                                           compare = likelihood)
  make_transform <- function(m) {
    function(theta) {
      as_double_mtx <- function(x){
        sapply(x,as.double)
      }
      c(lapply(m, as_double_mtx), as.list(theta))
    }
  }
  index <- function(info) {
    list(run = c(sum_clust = info$index$Pop_tot),
         state = c(Pop = info$index$Pop))
  }
  proposal_matrix <- diag(c(exp(1), 0.1, exp(1), 0.1))
  mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f", runif(n=1, min=-10, max=0), min = -1000, max = 0, prior = function(a) 1/a), mcstate::pmcmc_parameter("prop_f", runif(n=1, min=0, max=1), min = 0, max = 1, prior = function(a) a), mcstate::pmcmc_parameter("m", runif(n=1, min=-10, max=0), min = -1000, max = 0, prior = function(a) 1/a), mcstate::pmcmc_parameter("v", runif(n=1, min=0, max=1), min = 0, max = 1, prior = function(a) a)), proposal_matrix, make_transform(fixed_parameters))
  mcmc_pars$initial()
  det_filter <- mcstate::particle_deterministic$new(data = fitting_data,
                                           model = WF,
                                           index = index,
                                           compare = likelihood)

  n_steps <- 5 # this is a work-around (do a non-parallel run for five steps, which sets up the model and ensures that the parallel mcmc runs)
  n_burnin <- 0

  control <- mcstate::pmcmc_control(
    n_steps,
    save_state = TRUE,
    save_trajectories = TRUE,
    progress = TRUE,
    adaptive_proposal = TRUE,
    n_chains = 1)
  det_pmcmc_run <- mcstate::pmcmc(mcmc_pars, det_filter, control = control)

  n_steps <- steps_mcmc1
  n_burnin <- 0

  control <- mcstate::pmcmc_control(
    n_steps,
    save_state = TRUE,
    save_trajectories = TRUE,
    progress = TRUE,
    adaptive_proposal = TRUE,
    n_chains =4, n_workers = 4,
    n_threads_total = 4)

  det_pmcmc_run <- mcstate::pmcmc(mcmc_pars, det_filter, control = control)
  processed_chains <- mcstate::pmcmc_thin(det_pmcmc_run, burnin = 0.2 * steps_mcmc1, thin = 1)
  parameter_mean_hpd <- apply(processed_chains$pars, 2, mean)

  det_mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f",  parameter_mean_hpd[1], min = -1000, max = 0, prior = function(a) 1/a), mcstate::pmcmc_parameter("prop_f",  parameter_mean_hpd[2], min = 0, max = 1, prior = function(a) a), mcstate::pmcmc_parameter("m",  parameter_mean_hpd[3], min = -1000, max = 0, prior = function(a) 1/a), mcstate::pmcmc_parameter("v",  parameter_mean_hpd[4], min = 0, max = 1, prior = function(a) a)), proposal_matrix, make_transform(fixed_parameters))


  det_filter <- mcstate::particle_deterministic$new(data = fitting_data,
                                           model = WF,
                                           index = index,
                                           compare = likelihood)
  n_steps <- 5
  n_burnin <- 0

  control <- mcstate::pmcmc_control(
    n_steps,
    save_state = TRUE,
    save_trajectories = TRUE,
    progress = TRUE,
    adaptive_proposal = TRUE,
    n_chains = 1)
  det_pmcmc_run <- mcstate::pmcmc(det_mcmc_pars, det_filter, control = control)

  n_steps <- steps_mcmc2
  n_burnin <- 0

  control <- mcstate::pmcmc_control(
    n_steps,
    save_state = TRUE,
    save_trajectories = TRUE,
    progress = TRUE,
    adaptive_proposal = TRUE,
    n_chains = 4, n_workers = 4, n_threads_total = 4)
  det_pmcmc_run2 <- mcstate::pmcmc(det_mcmc_pars, det_filter, control = control)

  return(det_pmcmc_run2)
}

#' Title
#'
#' @param sim_parameters list of parameter values for v, sigma_f, prop_f, m
#' @param time_steps number of time steps (in months) to run the model for
#'
#' @importFrom mcstate particle_filter_data pmcmc_parameters pmcmc pmcmc_thin
#' @importFrom odin.dust odin_dust
#' @return fitted mcmc object
#' @export
#'
simulate_example <- function(sim_parameters = list("v" = 0.081, "sigma_f" = log(0.035), "prop_f" = 0.30, "m" = log(0.013)), time_steps = 36){
  # Set parameters for model
  vacc_time <- 4 # time of vaccination, in years after start of dataset
  PPsero_startpop <- readRDS("data/PPsero_startpop.rds") # start population for model (GPSC x serotypes)
  no_GPSC <- nrow(PPsero_startpop) # number of GPSCs
  #freq_sero
  no_sero <- ncol(PPsero_startpop) # number of serotypes
  Genotypes_matrix <- readRDS("data/Genotypes_matrix.rds")
  delta_test <- readRDS("data/gene_delta_ranking.rds") # delta statistic (computed as in Corander et al.)
  # calculates the changes in gene frequencies of the first time point (which are assumed to be at equilibrium) to the last time point
  # Genes that change the least in frequency have a low delta statistic value, genes that change more have a higher value.
  # then apply rank() function to determine order of genes
  # (NFDS model will determine a cut-off for genes to be under NFDS or not based on this order)
  gene_no_test <- length(delta_test) # number of intermediate-frequency genes
  migMatr_test_mtx <- readRDS("data/PPsero_mig.rds") # immigration matrix (rows are GPSCs, columns are serotypes)
  # compute immigration matrix by assigning all GPSC-serotype combinations that appear at least once in the data set equal probabilities of immigration (i.e. one over the number of GPSC-serotype combinations)
  vaccTypes_test <- readRDS(file = "data/SeroVT.rds") # information on which types are affected by vaccine (1 = serotype affected by vaccine, 0 = serotypes not affected)
  dt_test <- 1/12 # determines how many generations there are between two data points (here: 12 generations per year)

  example_params <- list(dt = dt_test, GPSC_no = no_GPSC, sero_no = no_sero, gene_no = gene_no_test, Pop_start = PPsero_startpop, Pop_eq = rowSums(PPsero_startpop), capacity = sum(PPsero_startpop), Genotypes = Genotypes_matrix, delta = (delta_test), Pop_mig_dist = migMatr_test_mtx, vaccTypes = vaccTypes_test, vacc_time = 4, v = sim_parameters$v, sigma_f = sim_parameters$sigma_f, prop_f = sim_parameters$prop_f, m = sim_parameters$m)

  simulate_model(model_parameters = example_params, simulation_steps = time_steps)
}

#' Title
#'
#' @param sim_parameters list of parameter values for v, sigma_f, prop_f, m
#' @param steps_mcmc1 number of mcmc steps for first run
#' @param steps_mcmc2 number of mcmc steps for second run (should be a lot higher for a good fit ~10,000 but then needs to be run on an hpc and will take ~24 hours)
#' @param time_steps_for_sim number of time steps (in months) for simulating data
#'
#' @importFrom mcstate particle_filter_data pmcmc_parameters pmcmc pmcmc_thin
#' @importFrom odin.dust odin_dust
#' @return fitted mcmc object
#' @export
#'
fit_example_to_sim_data <- function(sim_parameters = list("v" = 0.081, "sigma_f" = log(0.035), "prop_f" = 0.30, "m" = log(0.013)), steps_mcmc1 = 10, steps_mcmc2 = 50, time_steps_for_sim = 72){
  # simulate example data
  sim_data <- simulate_example(time_steps = time_steps_for_sim)

  # Set parameters for model
  vacc_time <- 4 # time of vaccination, in years after start of dataset
  PPsero_startpop <- data.frame(readRDS("data/PPsero_startpop.rds")) # start population for model (GPSC x serotypes)
  no_GPSC <- nrow(PPsero_startpop) # number of GPSCs
  #freq_sero
  no_sero <- ncol(PPsero_startpop) # number of serotypes
  Genotypes_matrix <- as.data.frame(readRDS("data/Genotypes_matrix.rds"))
  delta_test <- readRDS("data/gene_delta_ranking.rds") # delta statistic (computed as in Corander et al.)
  # calculates the changes in gene frequencies of the first time point (which are assumed to be at equilibrium) to the last time point
  # Genes that change the least in frequency have a low delta statistic value, genes that change more have a higher value.
  # then apply rank() function to determine order of genes
  # (NFDS model will determine a cut-off for genes to be under NFDS or not based on this order)
  gene_no_test <- length(delta_test) # number of intermediate-frequency genes
  migMatr_test_mtx <- data.frame(readRDS("data/PPsero_mig.rds")) # immigration matrix (rows are GPSCs, columns are serotypes)
  # compute immigration matrix by assigning all GPSC-serotype combinations that appear at least once in the data set equal probabilities of immigration (i.e. one over the number of GPSC-serotype combinations)
  vaccTypes_test <- readRDS(file = "data/SeroVT.rds") # information on which types are affected by vaccine (1 = serotype affected by vaccine, 0 = serotypes not affected)
  dt_test <- 1/12 # determines how many generations there are between two data points (here: 12 generations per year)

  example_params <- list(dt = dt_test,
                         GPSC_no = no_GPSC,
                         sero_no = no_sero,
                         gene_no = gene_no_test,
                         Pop_start = PPsero_startpop,
                         Pop_eq = rowSums(PPsero_startpop),
                         capacity = sum(PPsero_startpop),
                         Genotypes = Genotypes_matrix,
                         delta = delta_test,
                         Pop_mig_dist = migMatr_test_mtx,
                         vaccTypes = vaccTypes_test,
                         vacc_time = vacc_time)

  Stubentiger::fit_model_to_data(data = sim_data, fixed_parameters = example_params, steps_mcmc1 = steps_mcmc1, steps_mcmc2 = steps_mcmc2)
}
