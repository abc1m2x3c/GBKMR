# R programs for the article "Generalized Bayesian Kernel Machine Regression"
This is an example to implement Generalized Bayesian Kernel Machine Regression (GBKMR) for variable selection and effect estimation. The manuscript has been submitted to Biostatistics for review. Before use, please download the following R files in this repository  to the working directory: gbkmr.R, setparams.R, simgendata.R, src.R. 

Install missing packages or load the package if already installed    
```
packages <- c('invgamma', 'mvtnorm', 'truncnorm')
package.check <- sapply(
  packages,
  FUN = function(x) {
    if (!requireNamespace(x, quietly = TRUE)) {
      install.packages(x, dependencies = TRUE)
    }
    library(x, character.only = TRUE)
  }
)
```
        
Now we are ready to fit GBKMR model following steps below.

## Clean memory and set seed to guarantee the repetition
```
rm(list=ls())
set.seed(12345)
```
## Generate simulated data in the simulation of the paper 
```
source('simgendata.R')
n=200   # sample size
p=12    # number of candidate Z's
p0=1    # number of X variables
beta.true <- 0.1
family='gaussian' # Data type of the variable
pattern='nonlinear' # How X\beta+h(Z) is associated with g(\mu)
                    # In the manuscript, we provide three patterns
                    #   'linear' g(\mu)=X\beta+Z1+1.5*Z2+2*Z3
                    #   'quadratic' g(\mu)=X\beta+0.75*(Z1+Z2)^2+Z3
                    #   'nonlinear' g(\mu)=X\beta+2*sin(Z1*Z2*pi)+Z3
data=sim.gen.data(n,p,p0,beta.true,family,pattern) # The output of the function are a list ensembling:
                                                   # X an n by p0 matrix
                                                   # Y a length-n vector 
                                                   # Z an n by p matrix, where the first three columns Z1, Z2, and Z3
                                                   #   are truly associated with outcome Y through the pattern equation
```
For fitting GBKMR model using your own data, please ensemble the data in a list via the following code
```
data=list(X=X,y=y,Z=Z)
```

## Set start, prior, proposal distribution parameters
```
source('setparams.R')
start.params=set.start.param(data,family,delta.start=rep(FALSE,p)) # The starting status of 12 candidates Z1 to Z12 are 
                                                                   # set to be all negative 
prior.params=set.prior.param(data)
prop.params=set.prop.param()
```

## Variable selection using GBKMR
```
source('gbkmr.R')
res=gbkmr.varselect(data=data,family=family,start.params=start.params,prior.params=prior.params,prop.params=prop.params,T=1000,multijump=TRUE,r.prior.dist='invuniform')
# Function gbkmr.varselect
# Input:
# T: number of maximum iterations in the chain
# multijump: an option to inncrease the possibility of jumping to other statuses of delta 
#            designed to decrease the chance of trapping in local maxima trap
# r.prior.dist: prior distribution of parameter r, can be 'uniform' or 'invuniform'
#               The current version separate the prior option of r from the rest of prior.params
# Output:
# A list ensembling a record of parameter estimates in each iteration
# and some diagnostic parameters such as acceptance rate of Metroplis-Hasting Algorithm.
# The most related outcome is pip, which is a length-p vector indicating the posterior inclusion probability of each
# candidate variable
```

The output of the function will be shown as below.
```
[1] "Iteration 10 out of 1000 completed"
...
[1] "Iteration 980 out of 1000 completed"
[1] "Iteration 990 out of 1000 completed"
[1] "Iteration 1000 out of 1000 completed"
> res$pip
 [1] 1.000 1.000 1.000 0.000 0.000 0.000 0.000 0.000 0.272 0.000 0.000
[12] 0.000
```
The correct variables Z1 to Z3 are selected since they have 100% PIP, distinctly higher than the rest variables.

## Effect estimation using GBKMR
We can further use function gbkmr.pred to estimate the relationship between Z's and the outcome y. Below is a code example to estimate the relationship between z1 and y when z2=0.5 and z3=0.
```
# First, set the delta.start to be the status selected by gbkmr.varselect
start.params=set.start.param(data,family,delta.start=c(TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE))

# Second, construct the design matrix. In this example code, we estimate the relationship between Z1 (in a grid ranging
# from -1 to 1) and y when Z2=0.5 and Z3=0. Since we generated the data following pattern='nonlinear', we are actually
# estimating 2*sin(z1*pi/2)
z1_grid=seq(-1,1,length=50) #a grid of z1 used to predict f(z1)
if (pattern=='nonlinear'){
  dsgn.mat <- matrix(0, nrow = length(z1_grid), ncol = ncol(Z) - 2, byrow = TRUE)
  Z.new <- cbind(z1 = z1_grid, 0.5,dsgn.mat)  
}

# Third, fit the model and plot figures where dotted red is truth and solid black is estimation
res=gbkmr.pred(data,family,Z.new,start.params,prior.params,prop.params,T=1000,r.prior.dist='invuniform')
plot(z1_grid,res$h.est,type='l',xlab='z1',ylab='Estimate')
lines(z1_grid,2*sin(z1_grid*pi/2),col='red',lty=2)

```

<img src="https://github.com/abc1m2x3c/GBKMR/blob/main/example.png" width="400" height="400">

