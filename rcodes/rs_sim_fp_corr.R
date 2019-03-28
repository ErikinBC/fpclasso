###################################################################################################
# -------- SCRIPT TO SIMULATE FALSE POSITIVE CONTROL FOR DIFFERENT CORRELATION LEVELS ----------- #
###################################################################################################

##########################################
###### ---- (0) PRELIMINARIES ---- #######

rm(list=ls())

# Libraries for optparse
library(stringr,logical.return=F,warn.conflicts=F,quietly=T,verbose=F)
library(optparse,logical.return=F,warn.conflicts=F,quietly=T,verbose=F)

# Optional linux argument
option_list <- list(
  make_option('--dir_codes',type='character',default=getwd(),help='Path to code directory [default WD]'),
  make_option('--dir_output',type='character',default=getwd(),help='Path to output directory [default WD]')
);

# Parse and assign
opt <- parse_args(OptionParser(option_list=option_list))
dir.codes <- opt$dir_codes
dir.output <- opt$dir_output

# dir.codes='/home/erik/Documents/projects/FPC_Lasso/current/codes'
# dir.output='/home/erik/Documents/projects/FPC_Lasso/current/output'

# Remaining packages
pckgs <- c('data.table','cowplot','forcats','glmnet','survival')
for (pp in pckgs) { library(pp,character.only = T) }

# Load in the data generating functions
source(file.path(dir.codes,'rs_dgp.R'))
# Load in the FPC-Lasso algorithm
source(file.path(dir.codes,'rs_fpc_wrapper.R'))

############################################
###### ---- (1) RUN SIMULATIONS ---- #######

# Number of non-zero elements
nsupp <- function(x) { length(x[x!=0]) }

# Number of simulations to run for each experiment type
nsim <- 250
n <- 100
p <- 1000
k <- 5
nfp <- c(1,3,5,10)

yfam.seq <- c('gaussian','binomial','cox')
xfam.seq <- c('gaussian','binomial','exponential')
corr.seq <- c(0,0.25,0.5)

lst.store <- list()
jj <- 0
for (yf in yfam.seq) {
  if (yf == 'gaussian') {
    s2 <- sqrt(k^2/2)
  } else {
    s2 <- 0
  }
  for (xf in xfam.seq) {
    for (corr in corr.seq) {
      jj <- jj + 1
      print(sprintf('iteration: %i ---- y-dist: %s, x-dist: %s, correlation: %0.2f',jj,yf,xf,corr))
      # get lambda sequence
      lam.seq <- sapply(nfp,function(fp) qnorm(1-fp/(2*p)))
      # Temporary storage
      storage <- vector('list',length(nsim))
      for (ii in seq(nsim)) {
        if (ii %% 5 == 0) print(ii)
        set.seed(ii)
        # generate data design
        X1 <- dgp.X(n=n,p=p-k,xfam=xf, corr = corr, standardize = T, adj_skew = T,plant=ii)
        X2 <- dgp.X(n=n,p=k,xfam=xf, corr = 0.0, standardize = T, adj_skew = T,plant=ii+1)
        X <- cbind(X2, X1)
        yX <- dgp.X_2_y(X,yfam=yf,b0=1,k=k,s2=s2,a=0.5,plant=ii)
        if (yf == 'cox') {
          y <- Surv(time=yX[,1],event=yX[,2])
          X <- yX[,-(1:2)]
        } else {
          y <- yX[,1]
          X <- yX[,-1]
        }
        # Fit the models
        mdl.seq <- lapply(lam.seq,function(lam) fpc.lasso(x=X,y=y,family=yf,lambda=lam))
        # Store the lambda values
        tmp.store <- rbindlist(lapply(mdl.seq,function(ll) data.frame(lam_fpc=ll$lam.fpc,lam_lasso=ll$lam.glmnet)))
        stopifnot(nrow(tmp.store) == length(lam.seq))
        tmp.store <- data.table(lam_target=lam.seq, tmp.store)
        if (yf == 'cox') {
          lst.bhat <- lapply(mdl.seq, function(ll) as.vector(coef(object=ll$mdl,s=ll$lam.glmnet)))  
        } else {
          lst.bhat <- lapply(mdl.seq, function(ll) coef(object=ll$mdl,s=ll$lam.glmnet)[-1] )  
        }
        # Calculate number of FPs
        fp.seq <- unlist(lapply(lst.bhat,function(ll) nsupp(ll[-(1:k)]) ))
        # Calculate the number of true positives
        tp.seq <- unlist(lapply(lst.bhat,function(ll) nsupp(ll[1:k]) ))
        # Append
        tmp.store <- data.table(tmp.store, expected=nfp, fp=fp.seq, tp=tp.seq)
        # store
        storage[[ii]] <- tmp.store
      }
      storage <- melt(rbindlist(storage),measure.vars=c('fp','tp'),variable.name='type')
      storage[,`:=` (corr=corr,xfam=xf,yfam=yf) ]
      lst.store[[jj]] <- storage
    }
  }
}

storage.all <- rbindlist(lst.store)
storage.all[,list(mu=mean(value)),by=list(xfam,yfam,corr,expected,type)]
fwrite(x=storage.all,file=file.path(dir.output,'storage_fp_corr.csv'))

