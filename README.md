<img src="https://github.com/bacpop/Stubentiger/blob/add_logo/STUBENTIGER_logo.png" alt="" width="200"/>

# Stubentiger 
(STUdying Balancing Evolution (Nfds) To Investigate GEnome Replacement)

This is an R package to simulate and fit a population genetic model with parameters for vaccination effectiveness, negative frequency-dependent selection (NFDS) and immigration.

## Installation
```R
# install devtools which allows users to install SBMLtoOdin
install.packages("devtools")
library(devtools)

# install odin, odin.dust, dust, and mcstate

# install Stubentiger from our github repository
install_github("bacpop/Stubentiger")
library(Stubentiger)
```

## Usage
```R
simulate_example() # simulate model with example inputs and default parameters
fit_example_to_sim_data() # simulates data and fits model to simulated data, using example inputs and default parameters

# to fit to your own data, follow the steps described in the ProcessData.Rmd file on github.com/bacpop/NFDS_Model/tree/development/ModelExample 
# then run (for simulation)
simulate_model(model_parameters, simulation_steps)
# or run (for fitting / inference)
fit_model_to_data(data, fixed_parameters, steps_mcmc1 = 10, steps_mcmc2 = 50)
```
