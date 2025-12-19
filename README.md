<img src="STUBENTIGER_logo.png" alt="" width="200"/>

# Stubentiger
(STUdying Balancing Evolution (Nfds) To Investigate GEnome Replacement)

This is an R package to simulate and fit a population genetic model with parameters for vaccination effectiveness, negative frequency-dependent selection (NFDS) and immigration.

## Installation

Install Stubentiger's dependencies, `odin.dust`, `dust`, and `mcstate`:
```R
install.packages(
  "odin.dust", "dust", "mcstate",
  repos = c("https://mrc-ide.r-universe.dev", "https://cloud.r-project.org"))
```

Install `devtools`, and use it to install Stubentiger from its GitHub repository:
```R
install.packages("devtools")
devtools::install_github("bacpop/Stubentiger")
```

## Usage
```R
library(Stubentiger)
simulate_example() # simulate model with example inputs and default parameters
fit_example_to_sim_data() # simulates data and fits model to simulated data, using example inputs and default parameters
```

To fit to your own data, first follow the steps described in the ProcessDataForModel vignette. Then to simulate run
```R
simulate_model(model_parameters, simulation_steps)
```
and for fitting/inference run
```R
fit_model_to_data(data, fixed_parameters, steps_mcmc1 = 10, steps_mcmc2 = 50)
```
