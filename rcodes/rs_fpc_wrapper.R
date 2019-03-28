###########################################################################################
############# ----------- WRAPPER FOR FPC LASSO USING GLMNET --------------- ##############
###########################################################################################

###############################################################################
#  ------------------ (1) Survival/Support functions ------------------------ #

# Create the Ymatrix (Note! Uses Breslow approximation)
risksets <- function(So) {
  n <- nrow(So)
  Y <- matrix(0,nrow=n, ncol=n)
  if (ncol(So) == 2) {
    endtime <- So[,1]
    event <- So[,2]
    for (i in seq(n)) {
      Y[i,] <- endtime[i] >= endtime
    }
  } else {
    starttime <- So[,1]
    endtime <- So[,2]
    event <- So[,3]
    for (i in seq(n)) {
      Y[i,] <- (endtime[i] >= endtime) & (starttime[i] < endtime)
    }
  }
  diag(Y) <- event
  return(Y)
}

# Partial likelihood P matrix
Pfun <- function(Y,tY,eta) {
  diag(Y) <- diag(tY) <- 1
  rsk <- exp(eta)
  haz <- as.vector( tY %*% rsk )
  Pmat <- outer(rsk,haz,'/') * Y
  return(Pmat)
}

# Sigmoid function
sigmoid <- function(x) { 1 / (1+exp(-x)) }

# Soft-thresholding
softfun <- function(x,t) { sign(x)*pmax(abs(x)-t,0) }

#######################################################################
#  ------------------ (2) Residual Functions ------------------------ #

resfun.gaussian <- function(eta,y) { as.vector(y - eta) }
resfun.binomial <- function(eta,y) { as.vector(y - sigmoid(eta)) }
resfun.cox <- function(eta,y) {
  dd <- diag(y)
  as.vector(dd - (Pfun(y,t(y),eta) %*% dd))
}

resfun.lst <- list(gaussian=resfun.gaussian, binomial=resfun.binomial, cox=resfun.cox)

rm(list=c('resfun.gaussian','resfun.binomial','resfun.cox'))

###################################################################
#  ------------------ (3) Loss Functions ------------------------ #

lfun.lasso.gaussian <- function(eta,y,bhat,lambda) {
  0.5 * mean( (y - eta)^2 ) + lambda * sum( abs(bhat[-1]) )
}
lfun.lasso.binomial <- function(eta,y,bhat,lambda) {
  py <- sigmoid(eta)
  nll <- -1*mean( y*log(py) + (1-y)*log(1-py) ) + lambda * sum( abs(bhat[-1]) )
  return(nll)
}
lfun.lasso.cox <- function(eta,y,bhat,lambda) {
  dd <- diag(y)
  pp <- diag(Pfun(y,t(y), eta))
  pp <- pp[dd == 1]
  -mean(log(pp)) + lambda * sum( abs(bhat) )
}

lfun.lasso.lst <- list(gaussian=lfun.lasso.gaussian, binomial=lfun.lasso.binomial, cox=lfun.lasso.cox)

rm(list=c('lfun.lasso.gaussian','lfun.lasso.binomial','lfun.lasso.cox'))

############################################################
#  ------------------ (4) WRAPPER ------------------------ #

# x=tmp.Xpca;y=tmp.So
# family='cox';lambda=lam.star;nl=20;nlr=0.01;lamseq=NULL;retmdl=FALSE
fpc.lasso <- function(x,y,family,lambda,nl=20,nlr=0.01) {
  if (class(x) != 'matrix') { x <- as.matrix(x) }
  # Data size
  n <- nrow(x)
  p <- ncol(x)
  
  # Sanity checks
  stopifnot(family %in% c('gaussian','binomial','cox'))
  if (family == 'cox') {
    stopifnot(class(y) == 'Surv')
    y2 <- risksets(y)
  } else {
    y2 <- y
  }
  
  # Assign the family functions
  resfun <- resfun.lst[[family]]
  lfun <- lfun.lasso.lst[[family]]
  
  # glmnet wrapper
  fun.glmnet <- function(nl, nlr,lamseq=NULL,retmdl=FALSE) {
    if ( is.null(lamseq) ) {
      mdl.glmnet <- glmnet(x=x,y=y,family=family,alpha=1,nlambda=nl,lambda.min.ratio=nlr)
    } else {
      mdl.glmnet <- glmnet(x=x,y=y,family=family,alpha=1,lambda = lamseq)
    }
    lam.classical <- mdl.glmnet$lambda
    n.lam <- length(lam.classical)
    res.l2 <- apply(predict(mdl.glmnet,newx=x),2,function(eta) sqrt(sum(resfun(eta, y2)^2)) )
    lam.fpc <- (n * lam.classical) / res.l2
    
    idx.lam <- which.min((lam.fpc - lambda)^2)
    df.tab <- data.table(classical=lam.classical, fpc=lam.fpc, nearest=ifelse(lam.fpc==lam.fpc[idx.lam],T,F))
    if (retmdl) {
      return(list(df=df.tab, mdl=mdl.glmnet))
    } else {
      return(df.tab)
    }
  }
  
  # (1) recursive wrapper: keep dividing ratio by 10 until we find model
  recur <- T
  lam.seq <- NULL
  jj <- 0
  tab.store <- vector('list',4)
  while (recur) {
    jj <- jj + 1
    tab.mdl <- fun.glmnet(nl=nl,nlr=nlr,lamseq=lam.seq)
    if (jj == 1) lam.max <- min(tab.mdl$classical) / nlr
    recur <- which(tab.mdl$nearest) == nrow(tab.mdl)
    lam.seq <- exp(seq(log(lam.max*nlr),log(lam.max*(nlr/10)),length.out = nl))
    nlr <- nlr / 10
    tab.store[[jj]] <- tab.mdl
    if (jj > 4) { print('Error!!! Resursion applied more than 4 times!'); break  }
  }
  tab.mdl <- rbindlist(tab.store)
  if ( max(tab.mdl$fpc) < lambda ) { # Check if lambda.max < fpc.max
    ret.mdl <- fun.glmnet(nl=3,nlr=0.99, retmdl=TRUE)
    ret.list <- list(lam.fpc=tab.mdl[1]$fpc, lam.glmnet=tab.mdl[1]$classical, ratio=lambda /tab.mdl[1]$fpc,
                     mdl=ret.mdl$mdl)
    return(ret.list)
  }
  if (all(tab.mdl$fpc > lambda)) {
    ret.mdl <- glmnet(x,y,family,lambda=min(tab.mdl$classical))
    ret.list <- list(lam.fpc=tail(tab.mdl,1)$fpc, lam.glmnet=tail(tab.mdl,1)$classical,
                     ratio=tail(tab.mdl,1)$classical/lam.max, mdl=ret.mdl)
    return(ret.list)
  }
  
  # (2) Refine model
  lam.seq <- exp(seq(log(tab.mdl[fpc < lambda][which.max(fpc)]$classical),
                     log(tab.mdl[fpc > lambda][which.min(fpc)]$classical),length.out=nl))
  ret.mdl <- fun.glmnet(lamseq = lam.seq, retmdl=TRUE)
  
  df.star <- ret.mdl$df[which.min(abs(fpc - lambda))]
  lam.classical.star <- df.star$classical
  lam.fpc.star <- df.star$fpc
  nlr.star <- lam.classical.star / lam.max
  # Return
  ret.list <- list(lam.fpc=lam.fpc.star, lam.glmnet=lam.classical.star, ratio=lam.classical.star/lam.max,
                   mdl=ret.mdl$mdl)
  return(ret.list)
}


# # Preconditioning
# if (precond) {
#   mF <- with(svd(x), u %*% diag(1/d) %*% t(u))
#   xF <- mF %*% x
#   yF <- as.vector(mF %*% y)
#   # standardize
#   yF <- yF - mean(yF)
#   xF <- scale(xF)
#   yF2 <- sqrt(sum(yF^2))
#   bhat.lam <- as.vector(t(xF) %*% yF) / yF2
#   # bhat.lam <- softfun(bhat.lam, lambda)
#   return(bhat.lam)
# }

# ##################################################################
# #  ------------------ (5) Data Examples ------------------------ #
# 
# # Continuous
# X.gaussian <- model.matrix(~.,data=ElemStatLearn::ozone[,-1])[,-1]
# y.gaussian <- log(ElemStatLearn::ozone[,1])
# 
# # Binary
# X.binomial <- model.matrix(~., data=ElemStatLearn::SAheart[,-10])[,-1]
# y.binomial <- ElemStatLearn::SAheart$chd
# 
# # Right-censored
# X.cox <- model.matrix(~factor(trt)+celltype+karno+diagtime+age+factor(prior),data=veteran)[,-1]
# y.cox <- with(veteran, Surv(time=time,event=status))
# 
# # Test
# fpc.lasso(X.gaussian, y.gaussian, 'gaussian',lambda=1,nl=10,nlr=0.2)
# fpc.lasso(X.binomial, y.binomial, 'binomial',lambda=2,nl=10,nlr=0.2)
# fpc.lasso(X.cox, y.cox, 'cox',lambda=2,nl=10,nlr=0.2)

