##############################################################
######################## Simulation ##########################
##############################################################

########################## Note ##############################
## Choose the working directory as the directory where 
## this code file locates.
##############################################################

########################## Note ##############################
## Please first install required R and Python packages.
## R: System_preparation.R
## Python: pip install -r requirements.txt
##############################################################



rm(list = ls())
source('../../TreeTC_class/mvnorm.R', echo=F)

# Generate simulation data -----------------------------------------------------
num_sources = 10
num_obs_each = 10

set.seed(22)

dims = 2
DataOrigin = vector(mode = "list", length = num_sources)

for (s in 1:num_sources) {
  if (s %% 2 == 1) {
    tmpdat = rbind(rmvnorm(num_obs_each, mean = rep(0, dims), sigma = diag(2.0^2, dims, dims)),
                   rmvnorm(num_obs_each, mean = c(-5, 5), sigma = diag(1.0^2, dims, dims)),
                   rmvnorm(num_obs_each, mean = c(-5, 7.5), sigma = diag(0.5^2, dims, dims)),
                   rmvnorm(num_obs_each, mean = c(-7.5, 5), sigma = diag(0.5^2, dims, dims)),
                   rmvnorm(num_obs_each, mean = c(-5, -5), sigma = diag(1.0^2, dims, dims)),
                   rmvnorm(num_obs_each, mean = c(5, 5), sigma = diag(1.0^2, dims, dims)))
  } else {
    tmpdat = rbind(rmvnorm(num_obs_each, mean = rep(0, dims), sigma = diag(2.0^2, dims, dims)),
                   rmvnorm(num_obs_each, mean = c(5, -5), sigma = diag(1.0^2, dims, dims)),
                   rmvnorm(num_obs_each, mean = c(5, -7.5), sigma = diag(0.5^2, dims, dims)),
                   rmvnorm(num_obs_each, mean = c(7.5, -5), sigma = diag(0.5^2, dims, dims)),
                   rmvnorm(num_obs_each, mean = c(-5, -5), sigma = diag(1.0^2, dims, dims)),
                   rmvnorm(num_obs_each, mean = c(5, 5), sigma = diag(1.0^2, dims, dims)))
  }
  DataOrigin[[s]] = tmpdat
}

testData <- Reduce(rbind, DataOrigin)


save(list = ls(), file = "Sim_data.RData")
