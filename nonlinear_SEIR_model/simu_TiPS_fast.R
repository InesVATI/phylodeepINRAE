library(iterators)
library(parallel)
library(foreach)
library(doParallel)

library(ape)
library(TiPS)
path = "/home/ivati/Documents/data_TiPS1" 

# Create Model ------------------------------------------------------------

# reactions <- c(
#   "S [beta*S*(I+i)/N] -> E",
#   "0 [mu*(S+E+I+R)] -> S",
#   "S [mu*S] -> 0",
#   "E [sigma*E] -> I",
#   "E [mu*E] -> 0",
#   "I [mu*I] -> 0",
#   "I [gamma*I] -> R",
#   "R [mu*R] -> 0"
# )
# model_simu <- build_simulator(reactions)

reactions1 <- c(
  "S [beta*S*(I+i)/N] -> E",
  "E [sigma*E] -> I",
  "I [gamma*I] -> R",
  "R [s*R] -> 0",
  "R [(1-s)*R] -> U"
)
model1_simu <- build_simulator(reactions1)



# Faster Simulation -------------------------------------------------------
# Time of simulation is 5 years = 500 time units
TIME <- seq(from=0, to=500, length.out=700)
mTIME <- (TIME[1:699]+TIME[2:700])/2
dT = 0.009

beta_t <- function(B, b, m, T=100, p=0, MIN_TRANS=5e-4, plot=FALSE){
  beta_t <- B*(1+b*sin(2*pi*(m/T +p)))
  beta_t[beta_t<MIN_TRANS] <- MIN_TRANS
  
  if(plot){
    plot(m, beta_t, col="dark red", main=paste("mean=", B, "b=",b, "p=",p))
  }
  
  return(beta_t)
}

simulate_dynamic <- function(params, N=8000){

  params <- data.frame(params)
  beta_t <- beta_t(params["mean_signal",], params["amplitude",], mTIME)
  
  theta <- list(N=N, beta=beta_t, sigma=params["incubation_rate",], 
               gamma=params["recovery_rate",], i=params["import_param",], s=params["sampling_proba",])

  
  # Simulation of the trajectory --------------------------------------------

  capture.output(traj_i <- model1_simu(paramValues = theta,
                                      initialStates = c(S=N-1, E=0, I=1, R=0, U=0),
                                      times = TIME,
                                      method = "approximate",
                                      tau=dT,
                                      seed = 282875),
                 file= nullfile(), type="message")
  
  
  if (length(traj_i)==0){
    return(FALSE)
  }
  if (max(traj_i$traj$Time)<450){
    return(FALSE)
  }

  capture.output(tree_i <- simulate_tree(simuResults = traj_i,
                          deme = c("I", "R"),
                          root = "I",
                          isFullTrajectory = TRUE, # dead do generate leaves
                          nTrials = 10, # integer, number of unsuccessful simulation trials before giving up
                          addInfos = FALSE),
                 file=nullfile(), type="message")
  
  if (length(tree_i$tip.label)<90 || length(tree_i$tip.label)>900){
    # tree_i == list(), in case of this error : 
    #Cannot sample compartment R, the number of individuals is not sufficient.
    return(FALSE)
  }
  
  tree_i$node.label <- NULL
  
  return(tree_i)
}
  

# Set random parameters and generate data ---------------------------------

numCores <- detectCores()
# print(paste("Number of cores:", numCores))

generate_faster <- function(nb_data, first){

  params <- data.frame(
    mean_signal = runif(nb_data, 0.01, 0.35),
    amplitude = runif(nb_data, 0.5, 8),
    incubation_rate = runif(nb_data, 0.05, 1),
    recovery_rate = runif(nb_data, 0.02, 0.5)
  )
  # params$import_param = rep(2, nb_data) #fixed import parameters
  params$import_param = runif(nb_data, 0, 10) # varying import parameters
  params$sampling_proba <- rep(0.2, nb_data)

  registerDoParallel(numCores)
  trees <- foreach(i=1:4, .combine=cbind) %dopar% {
    apply(params[(1+(i-1)*nb_data/4):(i*nb_data/4), ], 1, simulate_dynamic)
    }

  ind <- sapply(trees, isFALSE)

  write.table(params[!ind, ], paste(path,"params_TiPS.csv", sep="/"), append=!first, row.names=FALSE, col.names=first)
  write.tree(phy=trees[!ind], file=paste(path, "Trees_TiPS.newick", sep="/"), append=!first)

}

start = Sys.time()
print(start)
generate_faster(nb_data=130000, first=FALSE)
print(Sys.time()-start)
