library(iterators)
library(parallel)
library(foreach)
library(doParallel)

library(ape)
library(TiPS)

# Create Model ------------------------------------------------------------

reactions <- c(
  "S [beta*S*(I+i)/N] -> E",
  "E [sigma*E] -> I",
  "I [gamma*I] -> R",
  "R [s*R] -> 0",
  "R [(1-s)*R] -> U"
)
model_simu <- build_simulator(reactions)


# Faster Simulation -------------------------------------------------------
# Time of simulation is 5 years = 500 time units
TIME <- seq(from=0, to=500, length.out=700)
mTIME <- (TIME[1:699]+TIME[2:700])/2
dT = 0.005

beta_t <- function(B, b, m, T=100, p=0, MIN_TRANS=5e-3){
  # Return the transmission signal beta_t evaluated in m
  beta_t <- B*(1+b*sin(2*pi*(m/T +p)))
  beta_t[beta_t<MIN_TRANS] <- MIN_TRANS
  
  return(beta_t)
}

rbnorm <- function(n, mean=0, std=1, min= NULL, max=NULL){
  # Return a vector of size n which stores realizations of 
  # a normal distribution truncated between min and max
  
  random_vect <- rnorm(n, mean=mean, sd=std)
  
  if (is.null(min) ||  is.null(max)){
    return(random_vect)
  }
  
  for(i in seq_len(n)){
    while (random_vect[i]< min || random_vect[i]>max)
      random_vect[i] <- rnorm(1, mean=mean, sd=std)
  }
  return(random_vect)
}

simulate_dynamic <- function(params, N=100000){
  # Return a simulated tree from params which stores a value for each parameters
  
  params <- data.frame(params)
  beta_t <- beta_t(params["mean_signal",], params["amplitude",], mTIME)
  
  incubation_period <- rbnorm(1, mean=5.9, std=0.5, min=3, max=15)*100/365
  incubation_rate <- 1/incubation_period
  infectious_period <- rbnorm(1, mean=7, std=0.5, min=3, max=10)*100/365
  recovery_rate <- 1/infectious_period
  
  theta <- list(N=N, beta=beta_t, sigma=incubation_rate, gamma=recovery_rate, 
                i=params["import_param",], s=params["sampling_proba",])

  
  # Simulation of the trajectory --------------------------------------------
  # capture.output is used to prevent messages to be printing
  capture.output(traj_i <- model_simu(paramValues = theta,
                                      initialStates = c(S=N-1, E=0, I=1, R=0, U=0),
                                      times = TIME,
                                      method = "approximate",
                                      tau=dT,
                                      seed = 282875),
                 file= nullfile(), type="message")
  
  if (length(traj_i)==0){
    # if the simulation of the trajectory has failed
    return(FALSE)
  }
  
  # Filter seasonal time series ---------------------------------------------
  X = !duplicated(round(traj_i$traj$Time,0)) # register indices where time point is not duplicated
  P = spectrum(traj_i$traj$I[X], log='no', plot=FALSE)
  ind = which.max(P$spec)
  
  if (ind != 5){
    return(FALSE)
  }
  
  # capture.output is used to prevent messages to be printing
  capture.output(tree_i <- simulate_tree(simuResults = traj_i,
                          deme = c("I", "R"),
                          root = "I",
                          isFullTrajectory = TRUE, # dead do generate leaves
                          nTrials = 20, # integer, number of unsuccessful simulation trials before giving up
                          addInfos = FALSE),
                 file=nullfile(), type="message")
  
  nb_tips <- length(tree_i$tip.label)
  if (nb_tips<60 || nb_tips>600){
    # tree_i == list(), in case of this error : 
    # Cannot sample compartment R, the number of individuals is not sufficient.
    return(FALSE)
  }

  tree_i$node.label <- NULL
  
  return(tree_i)
}
  

# Set random parameters and generate data ---------------------------------

path = "" 

numCores <- detectCores()
# print(paste("Number of cores:", numCores))

generate_faster <- function(nb_data, first){
  # Simulate nb_data random parameters and generate corresponding trees
  # first is a boolean, if first is false, new trees are added to the existing file
  
  params <- data.frame(
    mean_signal = runif(nb_data, 0.1, 0.7),
    amplitude = runif(nb_data, 0.1, 1),
    import_param = runif(nb_data, 0, 15)
  )
  params$sampling_proba <- rep(0.02, nb_data)

  registerDoParallel(numCores)
  D = nb_data/numCores
  trees <- foreach(i=1:numCores, .combine=cbind) %dopar% {
    apply(params[(1+(i-1)*D):(i*D), ], 1, simulate_dynamic)
    }

  ind <- sapply(trees, isFALSE)

  write.tree(phy=trees[!ind], file=paste(path, "Trees.newick", sep="/"), append=!first)
  print(paste(length(trees[!ind]), "valid trees"))
  write.table(params[!ind, ], paste(path,"params.csv", sep="/"), append=!first, row.names=FALSE, col.names=first)
  print(paste(nrow(params[!ind,]), "valid params!"))

}

generate_faster(nb_data=80000, first=FALSE)
