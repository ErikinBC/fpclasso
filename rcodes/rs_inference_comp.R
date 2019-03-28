#################################################################################
# -------- SCRIPT TO COMPARE INFERENCE PROPERTIES OF KNOCKOFFS + SI ----------- #
#################################################################################

rm(list=ls())

# Libraries for optparse
library(stringr,logical.return=F,warn.conflicts=F,quietly=T,verbose=F)
library(optparse,logical.return=F,warn.conflicts=F,quietly=T,verbose=F)

# Optional linux argument
option_list <- list(
  make_option('--dir_codes',type='character',default=getwd(),help='Path to R code directory[default WD]'),
  make_option('--dir_output',type='character',default=getwd(),help='Path to output directory [default WD]'),
  make_option('--dir_data',type='character',default=getwd(),help='Path to data directory to load survival/class/reg datasets [default WD]')
);

# Parse and assign
opt <- parse_args(OptionParser(option_list=option_list))
dir.codes <- opt$dir_codes
dir.output <- opt$dir_output
dir.data <- opt$dir_data
stopifnot(dir.exists(dir.codes),dir.exists(dir.output),dir.exists(dir.data))

# dir.codes="/home/erik/Documents/projects/FPC_Lasso/current/codes"
# dir.output="/home/erik/Documents/projects/FPC_Lasso/current/output"
# dir.data="/home/erik/Documents/projects/FPC_Lasso/current/output/data"

# Remaining packages
pckgs <- c('data.table','survival','glmnet','selectiveInference','knockoff','Rfast')
for (pp in pckgs) { library(pp,character.only = T) }

# Load in the wrapper
source(file.path(file.path(dir.codes),'rs_fpc_wrapper.R'))
source(file.path(file.path(dir.codes),'rs_dgp.R'))

nsupp <- function(x) { length(x[x!=0]) }

###################################################################
############## ----- (1) LOAD IN THE DATA ------ ##################

# Survival
load(file.path(dir.data,"surv_datasets.RData"))
ds.surv <- names(lst.surv)
df.surv.size <- rbindlist(lapply(lst.surv,function(ll) data.frame(n=nrow(ll$X),p=ncol(ll$X))))
df.surv.size[, `:=` (td=unlist(lapply(lst.surv,function(ll) ncol(ll$So)==3)), dataset=names(lst.surv)) ]
setcolorder(df.surv.size,'dataset')
# glmnet only accept non-time dependent choices
df.surv.size <- df.surv.size[td == F]
# Set E(FP)=1
df.surv.size <- df.surv.size[p >= 10 & n <= 1000]

# Simulation parameters
nsim <- 250
n <- 100
p <- 100
max.steps <- min(n,p)-1
s <- 5
fdr_alpha <- 0.1
frac <- 0.5
b0 <- 0.5


################################################################
############## ----- (2) FPC VS G'SELL ------ ##################

cn_gsell <- c('gsell_tp','gsell_fp','fpc_tp','fpc_fp','efp')
simstore_gsell <- data.frame(matrix(NA,nrow=nsim,ncol=length(cn_gsell),dimnames=list(NULL,cn_gsell)))
for (ii in seq(nsim)) {
  if (ii %% 5 == 0) print(ii)
  # Generate the data
  tmp.yX <- dgp.Xy(n=n,p=p,yfam='gaussian',xfam='gaussian',b0=b0,k=s,s2=1, plant=ii)
  tmp.y <- tmp.yX[,1]
  tmp.X <- scale(tmp.yX[,-1])
  
  # Run G-sell
  tmp.fsfit <- selectiveInference::fs(x=tmp.X,y=tmp.y,maxsteps = max.steps)
  tmp.out <- fsInf(tmp.fsfit,sigma = 1)
  tmp.steps <- forwardStop(tmp.out$pv, alpha=fdr_alpha)
  stopifnot(tmp.steps < max.steps) # Ensure full path not explored
  tmp.sel <- tmp.out$vars[seq(tmp.steps)]
  tmp.tp.gsell <- sum(tmp.sel %in% seq(s))
  tmp.fp.gsell <- length(tmp.sel) - tmp.tp.gsell
  tmp.efp.gsell <- fdr_alpha * length(tmp.sel)
  
  if (length(tmp.sel)==0) {
    print('No variables selected; jumping loop iteration')
    next
  }
  
  # Run FPC Lasso
  tmp.lam.fpc <- qnorm(1-tmp.efp.gsell/(2*p))
  tmp.fpc <- fpc.lasso(x=tmp.X,y=tmp.y,family='gaussian',lambda = tmp.lam.fpc)
  tmp.fpc.bhat <- coef(tmp.fpc$mdl,s=tmp.fpc$lam.glmnet)[-1]
  tmp.tp.fpc <- nsupp(tmp.fpc.bhat[seq(s)])
  tmp.fp.fpc <- nsupp(tmp.fpc.bhat[-seq(s)])
  
  # Store
  simstore_gsell[ii,] <- c(tmp.tp.gsell, tmp.fp.gsell, tmp.tp.fpc, tmp.fp.fpc, tmp.efp.gsell)
}

simstore_gsell <- data.table(simstore_gsell)
simstore_gsell[, iter := seq(.N)]
simstore_gsell <- melt(simstore_gsell,id.vars = c('iter','efp'))
simstore_gsell[, `:=` (approach=str_split_fixed(variable,'\\_',2)[,1],
                       measure=str_split_fixed(variable,'\\_',2)[,2], variable=NULL)]

simstore_gsell[,list(nsel=mean(value),efp=mean(efp)),by=list(approach,measure)][order(measure)]

#############################################################################
############## ----- (3) FPC VS SELECTIVE INFERENCE ------ ##################

fun.sim.si <- function(family,n,p,b0,s,frac,nsim,fdr_alpha) {
  cn_si <- c('si_tp','si_fp','fpc_tp','fpc_fp','efp')
  simstore_si <- data.frame(matrix(NA,nrow=nsim,ncol=length(cn_si),dimnames=list(NULL,cn_si)))
  
  ii <- 0
  tic <- 0
  while(tic < nsim) {
    ii <- ii + 1
    if (ii %% 5 == 0) print(ii)
    # Generate the data
    if (family=='gaussian') {
      tmp.yX <- dgp.Xy(n=n,p=p,yfam=family,xfam='gaussian',b0=b0,k=s, plant=ii, s2=1)  
    } else {
      tmp.yX <- dgp.Xy(n=n,p=p,yfam=family,xfam='gaussian',b0=b0,k=s, plant=ii)
    }
    
    if (family=='cox') {
      tmp.y <- Surv(time=tmp.yX[,1],event=tmp.yX[,2])
      tmp.X <- scale(tmp.yX[,-(1:2)])
    } else {
      tmp.y <- tmp.yX[,1]
      tmp.X <- scale(tmp.yX[,-1])
    }
    
    # Lambda the lambda max and choose some fixed threshold
    tmp.lammax <- glmnet(x=tmp.X,y=tmp.y,family=family,nlambda = 2,lambda.min.ratio = 0.999,standardize = F)$lambda[2]
    tmp.lamtest <- frac * tmp.lammax
    tmp.mdl <- glmnet(x=tmp.X,y=tmp.y,family=family,lambda = tmp.lamtest,standardize = F)
    tmp.beta <- as.vector(coef(tmp.mdl))
    if (family == 'gaussian') {
      tmp.beta <- tmp.beta[-1]
    }
    if (family=='cox') {
      tmp.si <- fixedLassoInf(x=tmp.X,y=tmp.y[,1],beta=tmp.beta, lambda=n*tmp.lamtest,family=family,status=tmp.y[,2])
    } else if (family=='gaussian') {
      tmp.si <- fixedLassoInf(x=tmp.X,y=tmp.y,beta=tmp.beta, lambda=n*tmp.lamtest,family=family,
                              sigma=estimateSigma(x=tmp.X,y=tmp.y)$sigmahat)
    } else {
      tmp.si <- fixedLassoInf(x=tmp.X,y=tmp.y,beta=tmp.beta, lambda=n*tmp.lamtest,family=family)
    }
    tmp.pv <- tmp.si$pv
    # Number of significant variables
    tmp.sel <- as.vector(tmp.si$vars)[tmp.pv < fdr_alpha]
    # Expected numer of false positives (number of rejected Null's times the type-1 error rate)
    tmp.efp.si <- length(tmp.sel)*fdr_alpha
    if (length(tmp.sel) == 0) {
      print('No variables selected, skipping')
      next
    } else {
      tic <- tic + 1
    }
    tmp.tp.si <- sum(tmp.sel %in% seq(s))
    tmp.fp.si <- length(tmp.sel) - tmp.tp.si
    
    # Run FPC Lasso
    tmp.lam.fpc <- qnorm(1-tmp.efp.si/(2*p))
    tmp.fpc <- fpc.lasso(x=tmp.X,y=tmp.y,family=family,lambda = tmp.lam.fpc)
    tmp.fpc.bhat <- coef(tmp.fpc$mdl,s=tmp.fpc$lam.glmnet)[-1]
    tmp.tp.fpc <- nsupp(tmp.fpc.bhat[seq(s)])
    tmp.fp.fpc <- nsupp(tmp.fpc.bhat[-seq(s)])
    
    # Store
    simstore_si[tic,] <- c(tmp.tp.si, tmp.fp.si, tmp.tp.fpc, tmp.fp.fpc, tmp.efp.si)
    
  }
  simstore_si <- data.table(simstore_si)[!is.na(si_tp)]
  simstore_si[, iter := seq(.N)]
  simstore_si <- melt(simstore_si,id.vars = c('iter','efp'))
  simstore_si[, `:=` (approach=str_split_fixed(variable,'\\_',2)[,1],
                      measure=str_split_fixed(variable,'\\_',2)[,2], variable=NULL)]
  simstore_si[, family := family ]
}

# Run FPC vs SI for three different family types
simstore_si <- fun.sim.si(family='gaussian',n=n,p=p,b0=b0,s=s,frac=frac,nsim=nsim,fdr_alpha=fdr_alpha)
# simstore_si <- rbindlist(simstore_si)
simstore_si[,list(nsel=mean(value),efp=mean(efp)),by=list(approach,measure)][order(measure,approach)]

# # Note, again, some of the EFP is not valid because SI doesn't have uniform p-values
# with(dcast(simstore_si[measure=='tp'],'iter~approach',value.var='value'),wilcox.test(fpc,si))
# with(dcast(simstore_si[measure=='tp'],'iter~approach',value.var='value'),t.test(fpc,si))

###################################################################
############## ----- (4) FPC VS KNOCKOFFS ------ ##################

# KNOCKOFFS DOES NOT SELECT ANYTHING!?!?

cn_knock <- c('knock_tp','knock_fp','fpc_tp','fpc_fp','efp')
simstore_knock <- data.frame(matrix(NA,nrow=nsim,ncol=length(cn_knock),dimnames=list(NULL,cn_knock)))

ii <- 0
tic <- 0
while(tic < nsim) {
  ii <- ii + 1
  if (ii %% 1 == 0) print(ii)
  # Generate the data
  tmp.yX <- dgp.Xy(n=n,p=floor(n/2),yfam='gaussian',xfam='gaussian',b0=b0,k=s,s2=1, plant=ii)
  tmp.y <- tmp.yX[,1]
  tmp.X <- scale(tmp.yX[,-1])

  # Run Knockoffs
  tmp.fixed <- create.fixed(X=tmp.X, method='sdp')
  tmp.W <- stat.glmnet_coefdiff(X=tmp.fixed$X, X_k = tmp.fixed$Xk, y=tmp.y, nfolds=10, family="gaussian")
  tmp.thresh <- knockoff.threshold(tmp.W, fdr=0.2, offset=0)
  tmp.sel <- which(tmp.W >= tmp.thresh)
  if (length(tmp.sel) == 0) {
    print('No variables selected, skipping')
    next
  } else {
    tic <- tic + 1
  }
  tmp.tp.knock <- sum(tmp.sel %in% seq(s))
  tmp.fp.knock <- length(tmp.sel) - tmp.tp.knock
  tmp.efp.knock <- fdr_alpha * length(tmp.sel)
  
  # Run FPC Lasso
  tmp.lam.fpc <- qnorm(1-tmp.efp.knock/(2*p))
  tmp.fpc <- fpc.lasso(x=tmp.X,y=tmp.y,family='gaussian',lambda = tmp.lam.fpc)
  tmp.fpc.bhat <- coef(tmp.fpc$mdl,s=tmp.fpc$lam.glmnet)[-1]
  tmp.tp.fpc <- nsupp(tmp.fpc.bhat[seq(s)])
  tmp.fp.fpc <- nsupp(tmp.fpc.bhat[-seq(s)])
  # Store
  simstore_knock[tic,] <- c(tmp.tp.knock, tmp.fp.knock, tmp.tp.fpc, tmp.fp.fpc, tmp.efp.knock)
}
simstore_knock <- data.table(simstore_knock)
simstore_knock[, iter := seq(.N)]
simstore_knock <- melt(simstore_knock,id.vars = c('iter','efp'))
simstore_knock[, `:=` (approach=str_split_fixed(variable,'\\_',2)[,1],
                       measure=str_split_fixed(variable,'\\_',2)[,2], variable=NULL)]

simstore_knock[,list(nsel=mean(value),efp=mean(efp)),by=list(approach,measure)][order(measure)]

################################################################
############## ----- (5) REAL DATASETS ------ ##################

# cn.real <- c('fpc_nsel','fpc_efp','fpc_etp','fpc_fdr')
# simstore_real <- data.frame(matrix(NA,nrow=nrow(df.surv.size),ncol=length(cn.real),dimnames=list(NULL,cn.real)))

simstore_real <- vector('list',nrow(df.surv.size))

ii <- 0
for (ds in df.surv.size$dataset) {
  ii <- ii + 1
  tmp.X <- lst.surv[[ds]]$X
  tmp.So <- lst.surv[[ds]]$So
  # Create orthogonal design matrix with X
  tmp.Xpca <- prcomp(tmp.X)$x
  # Create a skew-adjusted version
  tmp.Xpca <- apply(tmp.Xpca,2,skew.adj)
  # Scale
  tmp.Xpca <- standardise(tmp.Xpca)
  tmp.p <- ncol(tmp.Xpca)
  tmp.n <- nrow(tmp.Xpca)
  print(sprintf('Iteration %i --- Survival dataset: %s, n = %i, p = %i',ii,ds,tmp.n, tmp.p))
  
  # (1) Fit FPC-lasso with a single FP target
  efp.seq <- seq(1,12)
  tmp.store <- vector('list',length(efp.seq))
  for (jj in seq_along(efp.seq)) {
    tmp.efp <- efp.seq[jj]
    lam.star <- qnorm(1-tmp.efp/(2*tmp.p))
    mdl.fpc <- fpc.lasso(x=tmp.Xpca, y=tmp.So, family='cox',lambda=lam.star,nl=50)
    if (abs(mdl.fpc$lam.fpc - lam.star) > 0.01) { print('Breaking loop no lambda found'); next }
    bhat.fpc <- as.vector(coef(mdl.fpc$mdl, s=mdl.fpc$lam.glmnet))
    nsel.fpc <- nsupp(bhat.fpc)
    tmp.etp <- nsel.fpc - tmp.efp
    tmp.fdr <- tmp.efp / nsel.fpc
    tmp.store[[jj]] <- data.frame(efp=tmp.efp, etp=tmp.etp, nsel=nsel.fpc, fdr=tmp.fdr)
  }
  tmp.store <- rbindlist(tmp.store)[fdr == min(fdr)]
  tmp.store <- tmp.store[etp == max(etp)]
  tmp.store[, dataset := ds]
  simstore_real[[ii]] <- tmp.store
}

simstore_real <- rbindlist(simstore_real)  

###################################################################
############## ----- (6) SAVE THE RESULTS ------ ##################

inference.store <- list(gsell=simstore_gsell, si=simstore_si, knock=simstore_knock,real=simstore_real)
save(inference.store, file=file.path(dir.output,'inference_store.RData'))
