library(ape)
library(TiPS)

# Define the simulator SEIR ------------------------------------------------

reactions1 <- c(
  "S [beta*S*(I+i)/N] -> E",
  "E [sigma*E] -> I",
  "I [gamma*I] -> R",
  "R [s*R] -> 0",
  "R [(1-s)*R] -> U"
)
model1_simu <- build_simulator(reactions1)


# Simulate dynamics and trees ---------------------------------------------

MIN_TRANS <- 5e-3
# TIME <- seq(from=0, to=800, length.out=1200)
# mTIME <- (TIME[1:1199]+TIME[2:1200])/2
TIME <- seq(from=0, to=500, length.out=700)
mTIME <- (TIME[1:699]+TIME[2:700])/2

rbnorm <- function(n, mean=0, std=1, min= NULL, max=NULL){
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
# rbnorm(3, min= 0, max=2)

beta <- function(B, b, m, T=100, p=0, plot=FALSE){
  
  beta_t <- B*(1+b*sin(2*pi*(m/T +p)))
  beta_t[beta_t<=MIN_TRANS] <- MIN_TRANS
  
  if (plot){
    plot(m, beta_t, col="dark red")
  }
  
  return(beta_t)
}
N=100000
s=0.02

incubation_period <- rbnorm(1, mean=5.9, std=0.5, min=3, max=15)*100/365
sigma <- 1/incubation_period 
sigma=0.56
infectious_period <- rbnorm(1, mean=7, std=0.5, min=3, max=10)*100/365
gamma <- 1/infectious_period 
gamma=0.5

# mean_trans <- runif(1, 0.3, 0.6) # 0.4
mean_trans = 0.58
# b <- runif(1, 0.15, 1) #0.5
b =0.3
p <- 0 #runif(1, -0.5, 0.5)
# i <- runif(1, 0, 15)
i=2

# Simulate traj -----------------------------------------------------------

method = "approximate"
tau=0.005
traj1 <- model1_simu(paramValues = list(N=N, beta=beta(B=mean_trans, b=b, m=mTIME, p=p, plot=TRUE), sigma=sigma, gamma=gamma,i=i, s=s),
                     initialStates = c(S=N-1, E=0, I=1, R=0, U=0),
                     times=TIME,
                     method=method,
                     tau=tau,
                     seed=282875)
plot(traj1$traj$Time, traj1$traj$I, type='l', main=paste("Infected individuals method", method, "dT", tau))
# plot(traj1$traj$Time, traj1$traj$E, type='l', main=paste("E individuals method", method, "dT", tau))
# plot(traj1)
# title(main=paste("B", round(mean_trans, 2), "b", round(b,2), "i",round(i,2), "sig", round(sigma,2), "gam",round(gamma,2)))
# without given sampling dates (this can be done isFullTrajectory option)
tree <- simulate_tree(simuResults = traj1,
                      deme = c("I","R"),
                      root = "I",
                      isFullTrajectory = TRUE, # dead do generate leaves
                      nTrials = 20, # integer, number of unsuccessful simulation trials before giving up
                      addInfos = TRUE)
print(length(tree$tip.label))

inode_cols <- ifelse(grepl(x=tree$node.label, pattern="R"), "blue", "red")
tips_cols <- ifelse(grepl(x=tree$tip.label, pattern="R"), "blue", "red")
plot.phylo(tree, root.edge = T, no.margin = F, show.tip.label = F)
tiplabels(pch=20, col=tips_cols)
nodelabels(pch=20, col=inode_cols)



# Save some trajectories --------------------------------------------------
n=5
write.csv(traj1$traj[,-c(2,3,4, 8)], paste('traj_',n, '.csv', sep=""))
write.tree(tree, paste('tree_',n, '.newick', sep=""))

# tree=read.tree('tree_0.newick')


# Filter right trajectories -----------------------------------------------
# traj = read.csv('traj_2.csv')
traj = traj1$traj

# plot(traj$Time, traj$I, type='l')
# plot(traj$I, type='l')

X = !duplicated(round(traj$Time,0))
plot(traj$I[X], type='l', main="0")

# library(TSA)
# f = 5/length(traj$I[X])
# f
P = spectrum(traj$I, log='no', plot=FALSE)
P = spectrum(traj$I[X], log='no', plot=FALSE)
# P$spec
# P = spectrum(traj$I[X], span=4,  log='no')
# spx <- P$freq
# spy <- P$spec
# plot (spy~spx, subset=spx<=0.005, type='l')
ind = which.max(P$spec)
ind
# P$spec[5]
# P$freq[ind]
# round(P$freq[ind], 4) == round(f, 4)

ord <- order(P$spec)
t = as.numeric(length(ord))
5 %in% ord[(t-2):t]

library(zoo)
sI = zooreg(traj$I[X])
# Y = c(traj$I[1])
# init=traj$Time[1]
# last=length(traj$Time)-1
# for(i in 1:last){
#   print(init + 0.001*i)
#   Y = c(Y, traj$I[traj$Time == init + 0.009*i])
# }
# plot(Y)

library(xts)
ITS = xts(traj$I, order.by=as.POSIXct('1970-01-01')+traj$Time)
seir_jacobian = jacobian.net(data=ITS, m=3:3, timelapse="VARIABLE")
summary(seir_jacobian)
exponent <- lyapunov.max(data=seir_jacobian, B=100, doplot = TRUE)
summary(exponent)