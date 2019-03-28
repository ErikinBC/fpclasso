#################################################################################################
############# ----------- EMPIRICAL PROOF OF REGULARITY CONDITIONS --------------- ##############
#################################################################################################

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

# Load in the other packages
pckgs <- c('data.table','glmnet','survival')
for (pp in pckgs) { library(pp,character.only = T) }

# Load in the support functions
source(file.path(dir.codes,'rs_fpc_wrapper.R'))

###################################################################
############## ----- (1) LOAD IN THE DATA ------ ##################

# Survival
load(file.path(dir.data,"surv_datasets.RData"))
ds.surv <- names(lst.surv)

# Classification
load(file.path(dir.data,"class_datasets.RData"))
ds.class <- names(lst.class)

# Regression
load(file.path(dir.data,"reg_datasets.RData"))
ds.reg <- names(lst.reg)

#####################################################################
############## ----- (2) SURVIVAL DATASETS  ------ ##################

store.surv <- vector('list',length(ds.surv))
for (ii in seq_along(ds.surv)) {
  tmp.ds <- ds.surv[ii]
  tmp.td <- ncol(lst.surv[[tmp.ds]]$So)==2
  tmp.n <- nrow(lst.surv[[tmp.ds]]$So)
  tmp.p <- ncol(lst.surv[[tmp.ds]]$X)
  if (tmp.td & (tmp.n < 2500)) {
    print(sprintf('Dataset %i of %i: %s, n: %i, p: %i',ii,length(ds.surv),tmp.ds,tmp.n, tmp.p))
    tmp.So <- lst.surv[[tmp.ds]]$So
    tmp.So[,1] <- tmp.So[,1]+1
    tmp.X <- lst.surv[[tmp.ds]]$X
    if (class(tmp.X) != 'matrix') {
      tmp.X <- as.matrix(tmp.X)
    }
    tmp.y <- risksets(tmp.So)
    # Fit model
    tmp.mdl <- glmnet(x=tmp.X,y=tmp.So,family='cox', nlambda=50, lambda.min.ratio = 0.01)
    # Calculate residuals
    tmp.eta <- predict(tmp.mdl,newx = tmp.X)
    tmp.s2 <- apply(tmp.eta,2,function(ee) sqrt(sum(resfun.lst$cox(eta=ee,y=tmp.y)^2)) )
    tmp.lam.lasso <- tmp.mdl$lambda
    tmp.lam.fpc <- (tmp.n * tmp.lam.lasso) / tmp.s2
    tmp.df <- data.table(dataset=tmp.ds, n=tmp.n, p=tmp.p , lam_fpc=tmp.lam.fpc, lam_lasso=tmp.lam.lasso)
    store.surv[[ii]] <- tmp.df
  } else {
    print('Time dependent survival or n>2500 ~~ skipping')
  }
}
df.surv <- data.table(type='survival',rbindlist(store.surv))

###########################################################################
############## ----- (3) CLASSIFICATION DATASETS  ------ ##################

store.class <- vector('list',length(ds.class))
for (ii in seq_along(ds.class)) {
  tmp.ds <- ds.class[ii]
  tmp.n <- nrow(lst.class[[tmp.ds]]$X)
  tmp.p <- ncol(lst.class[[tmp.ds]]$X)
  print(sprintf('Dataset %i of %i: %s, n: %i, p: %i',ii,length(ds.class),tmp.ds,tmp.n, tmp.p))
  tmp.y <- lst.class[[tmp.ds]]$y
  tmp.X <- lst.class[[tmp.ds]]$X
  # Fit model
  tmp.mdl <- glmnet(x=tmp.X,y=tmp.y,family='binomial', nlambda=50, lambda.min.ratio = 0.01)
  # Calculate residuals
  tmp.eta <- predict(tmp.mdl,newx = tmp.X)
  tmp.s2 <- apply(tmp.eta,2,function(ee) sqrt(sum(resfun.lst$binomial(eta=ee,y=tmp.y)^2)) )
  tmp.lam.lasso <- tmp.mdl$lambda
  tmp.lam.fpc <- (tmp.n * tmp.lam.lasso) / tmp.s2
  tmp.df <- data.table(dataset=tmp.ds, n=tmp.n, p=tmp.p , lam_fpc=tmp.lam.fpc, lam_lasso=tmp.lam.lasso)
  store.class[[ii]] <- tmp.df
}
df.class <- data.table(type='class',rbindlist(store.class))
df.class[,list(ma=max(lam_fpc),mi=min(lam_fpc)),by=dataset]

#######################################################################
############## ----- (4) REGRESSION DATASETS  ------ ##################

store.reg <- vector('list',length(ds.reg))
for (ii in seq_along(ds.reg)) {
  tmp.ds <- ds.reg[ii]
  tmp.n <- nrow(lst.reg[[tmp.ds]]$X)
  tmp.p <- ncol(lst.reg[[tmp.ds]]$X)
  print(sprintf('Dataset %i of %i: %s, n: %i, p: %i',ii,length(ds.reg),tmp.ds,tmp.n, tmp.p))
  tmp.y <- lst.reg[[tmp.ds]]$y
  tmp.X <- lst.reg[[tmp.ds]]$X
  # Fit model
  tmp.mdl <- glmnet(x=tmp.X,y=tmp.y,family='gaussian', nlambda=50, lambda.min.ratio = 0.01)
  # Calculate residuals
  tmp.eta <- predict(tmp.mdl,newx = tmp.X)
  tmp.s2 <- apply(tmp.eta,2,function(ee) sqrt(sum(resfun.lst$gaussian(eta=ee,y=tmp.y)^2)) )
  tmp.lam.lasso <- tmp.mdl$lambda
  tmp.lam.fpc <- (tmp.n * tmp.lam.lasso) / tmp.s2
  tmp.df <- data.table(dataset=tmp.ds, n=tmp.n, p=tmp.p , lam_fpc=tmp.lam.fpc, lam_lasso=tmp.lam.lasso)
  store.reg[[ii]] <- tmp.df
}
df.reg <- data.table(type='reg',rbindlist(store.reg))
df.reg[,list(ma=max(lam_fpc),mi=min(lam_fpc)),by=dataset]

#####################################################################
############## ----- (5) COMBINE AND CHECK  ------ ##################

# Combine all
df.merge <- rbind(df.surv, df.class, df.reg)
# scale lam_lasso from 0,1
df.merge[,lam_scale := (lam_lasso-min(lam_lasso))/(max(lam_lasso) - min(lam_lasso)),by=dataset]
# Calculate the step
df.merge[,diff_lam_fpc := as.numeric(lam_fpc - shift(lam_fpc,1)),by=list(dataset)]
# Flag any violations
df.merge[, violation := as.numeric(df.merge$diff_lam_fpc) > 0]

# Save the data for later
fwrite(x=df.merge,file=file.path(dir.output,'df_lamseq.csv'))



