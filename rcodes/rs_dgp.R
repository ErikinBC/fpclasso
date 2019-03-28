#############################################################
# -------- FUNCTIONS FOR DATA GENEATING PROCESS ----------- #
#############################################################

sigmoid <- function(x) { 1/(1+exp(-x)) }
activation <- function(x,t=0.5) { ifelse(x > t, 1, 0) }

######################################################
# -------- FUNCTIONS FOR GENERATING SNSs ----------- #
######################################################

skewness <- function (x, na.rm = FALSE, type = 3) {
  if (any(ina <- is.na(x))) {
    if (na.rm) 
      x <- x[!ina]
    else return(NA)
  }
  if (!(type %in% (1:3))) 
    stop("Invalid 'type' argument.")
  n <- length(x)
  x <- x - mean(x)
  y <- sqrt(n) * sum(x^3)/(sum(x^2)^(3/2))
  if (type == 2) {
    if (n < 3) 
      stop("Need at least 3 complete observations.")
    y <- y * sqrt(n * (n - 1))/(n - 2)
  }
  else if (type == 3) 
    y <- y * ((1 - 1/n))^(3/2)
  y
}

skew.score <- function(c, x) (skewness(log(x + c)))^2

skew.adj <- function(x) {
  x <- (x-mean(x))/sd(x)
  xmin <- min(x)
  xmax <- max(x)
  cstar <- optimise(f=skew.score,lower=-xmin+1e-6,upper=max(-xmin,xmax)+10,x=x)$minimum
  xstar <- log(x + cstar)
  if (abs(skewness(xstar)) > abs(skewness(x))) {
    return(x)
  } else {
    return(xstar)  
  }
}

# Data generation function
dgp.dist <- function(dist,n,a=NULL,b=NULL) {
  if (dist == 'gaussian') {
    x <- rnorm(n,mean=a,sd=b)
  } else if (dist == 'exponential') {
    x <- rexp(n,a)
  } else if (dist == 'gamma') {
    x <- rgamma(n,a,b)
  } else if (dist == 'beta') {
    x <- rbeta(n,a,b)
  } else if (grepl('binomial',dist)) {
    x <- rbinom(n,1,a)
  } else {
    print('Error! Please pick a valid distribution')
  }
  return(x)
}

# model: x'y / || y ||_2

# SNS function (scale X only)
sns.fun <- function(d.x,d.y,n,a.x,b.x,a.y,b.y) {
  x <- dgp.dist(dist=d.x,n=n,a=a.x,b=b.x)
  x <- (x - mean(x))/sd(x)
  y <- dgp.dist(dist=d.y,n=n,a=a.y,b=b.y)
  y <- y - mean(y)
  sns <- sum(x*y)/sqrt(sum(y*y))
  return(sns)
}

# SNS function with skew adjustment
sns.skew.fun <- function(d.x,d.y,n,a.x,b.x,a.y,b.y) {
  x <- dgp.dist(dist=d.x,n=n,a=a.x,b=b.x)
  cstar <- optimise(f=skew.score,interval=c(abs(min(x)),10),x=x)$minimum
  x <- log(x+cstar)
  x <- (x-mean(x))/sd(x)
  y <- dgp.dist(dist=d.y,n=n,a=a.y,b=b.y)
  y <- y - mean(y)
  sns <- sum(x*y)/sqrt(sum(y*y))
  return(sns)
}

###################################################################
# -------- FUNCTIONS FOR GENERATING REGRESSION MODELS ----------- #
###################################################################

# n=100;p=1000;xfam='gaussian';corr=0.25;standardize=T;adj_skew=T;plant=3;inflation=0.1
dgp.X <- function(n,p,xfam,corr=0,standardize=T,adj_skew=T,plant=1,inflation=0.1) {
  stopifnot(xfam %in% c('gaussian','binomial','exponential'))
  stopifnot( corr >= 0)
  # ---- generate X ---- #
  set.seed(plant)
  if (xfam == 'gaussian') {
    if (corr != 0) {
      sig <- sqrt(1/corr - 1)
      X <- matrix(rep(rnorm(n),p),ncol=p)
      U <- matrix(rnorm(n*p,sd = sig),ncol=p)
      V <- X + U
    } else {
      V <- matrix(rnorm(n*p),ncol=p)
    }
  } else if (xfam == 'binomial') {
    if (corr != 0) {
      X <- matrix(rep(rnorm(n),p),ncol=p)
      check <- T
      z <- log(corr/(1-corr))
      jj <- 0
      while (check) {
        jj <- jj + 1
        sig <- sqrt(1/sigmoid(z) - 1)
        U <- matrix(rnorm(n*p,sd = sig),ncol=p)
        V <- X + U
        V <- activation(V)
        chat <- cor(V[,sample(p,floor(p*0.1))])
        cc <- mean(chat[lower.tri(chat)])
        diff <- cc - corr
        if (abs(diff) > 0.01) {
          if (sign(diff) == -1) {
            z <- z + inflation
          } else {
            z <- z - inflation
          }
        } else {
          check <- F
        }
      }
    } else {
      V <- matrix(rbinom(n*p,1,0.5),ncol=p)
    }
    
  } else if (xfam == 'exponential') {
    if (corr != 0) {
      lam <- 1 / sqrt(1/corr - 1)
      X <- matrix(rep(rexp(n,rate=1),p),ncol=p)
      U <- matrix(rexp(n*p,rate=lam),ncol=p)
      V <- X + U
    } else {
      V <- matrix(rexp(n*p),ncol=p)
    }
  } else {
    print('Pick a valid distribution please!'); break
  }
  # skew adjust
  if (adj_skew) {
    V <- apply(V,2,skew.adj)
  }
  # scale
  if (standardize) {
    V <- scale(V)
  }
  return(V)
}

# ---- FUNCTION THAT MAPS AN OUTPUTTED X TO A RESPONSE VECTOR

# x=X;yfam=yf;b0=1;k=k;s2=s2
dgp.X_2_y <- function(x,yfam,b0,k,s2=0,a=NULL,plant=1) {
  set.seed(plant*1234)
  n <- nrow(x)
  p <- ncol(x)
  # --- generate y --- #
  if (yfam == 'gaussian') {
    ydgp <- function(eta) { eta }
  } else if (yfam == 'binomial') {
    ydgp <- function(eta) { rbinom(n,1,1/(1+exp(-eta))) }
  } else if (yfam == 'cox') {
    ydgp <- function(eta) {
      # generate survival time
      tt <- rexp(n,rate=exp(eta))
      ratestar <- log(2)/quantile(tt,0.75)
      dd <- rexp(n,ratestar)
      is.event <- ifelse(tt < dd,1,0)
      is.time <- ifelse(tt < dd, tt, dd)
      y <- cbind(time=is.time,event=is.event)
      return(y)
    }
  } else {
    print('Pick a valid distribution please!'); break
  }
  # Generate
  bvec <- c(rep(b0,k),rep(0,p-k))
  eta <- as.vector(x %*% bvec)
  noise <- rnorm(n,0,sqrt(s2))
  y <- ydgp(eta + noise)
  
  # Return
  return(cbind(y,x))
  
}

# ---- FUNCTION TO GENERATE GLM AND exponential MATRIX ----- #
# n=100;p=5;yfam='cox';xfam='exponential';a=0.25;standardize=T;adj_skew=T;b0=1;k=5;s2=0;plant=1
dgp.Xy <- function(n,p,yfam,xfam,b0,k,s2=0,a=NULL,plant=1,standardize=T,adj_skew=T) {
  set.seed(plant)
  # ---- generate X ---- #
  if (xfam == 'gaussian') {
    xdgp <- function(n,a) rnorm(n)
  } else if (xfam == 'binomial') {
    xdgp <- function(n,a) rbinom(n,1,a)
  } else if (xfam == 'exponential') {
    xdgp <- function(n,a)  rexp(n,a)
  } else {
    print('Pick a valid distribution please!'); break
  }
  x <- matrix(xdgp(n*p,a),ncol=p)
  if (adj_skew) {
    x <- apply(x,2,skew.adj)
  }
  if (standardize) {
    x <- scale(x)
  }
  
  # --- generate y --- #
  if (yfam == 'gaussian') {
    ydgp <- function(eta) { eta }
  } else if (yfam == 'binomial') {
    ydgp <- function(eta) { rbinom(n,1,1/(1+exp(-eta))) }
  } else if (yfam == 'cox') {
    ydgp <- function(eta) {
      # generate survival time
      tt <- rexp(n,rate=exp(eta))
      ratestar <- log(2)/quantile(tt,0.75)
      dd <- rexp(n,ratestar)
      is.event <- ifelse(tt < dd,1,0)
      is.time <- ifelse(tt < dd, tt, dd)
      y <- cbind(time=is.time,event=is.event)
      return(y)
    }
  } else {
    print('Pick a valid distribution please!'); break
  }
  # Generate
  bvec <- c(rep(b0,k),rep(0,p-k))
  eta <- as.vector(x %*% bvec)
  noise <- rnorm(n,0,sqrt(s2))
  eta <- eta + noise  
  y <- ydgp(eta)
  
  # Return
  return(cbind(y,x))
  
}