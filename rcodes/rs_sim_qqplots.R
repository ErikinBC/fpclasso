#####################################################################################
# -------- SCRIPT TO SIMULATE QQ PLOTS FOR VARIETY OF A/B DISTRIBUTIONS ----------- #
#####################################################################################

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

# Remaining packages
pckgs <- c('data.table')
for (pp in pckgs) { library(pp,character.only = T) }

print(dir.output)

# Load in the data generating functions
source(file.path(dir.codes,'rs_dgp.R'))

#######################################################################################

start_time <- Sys.time()

# x,y are of length n, and we generate nrep SNSs
n <- 100
nrep <- 1000

# Loop through each
dist.seq <- c('gaussian','exponential','gamma','beta','binomial1','binomial2')
param.seq <- list(gaussian=list(a=0,b=1), exponential=list(a=1,b=NULL), gamma=list(a=5,b=1/2),
                  beta=list(a=2,b=1),binomial1=list(a=0.5,b=NULL),binomial2=list(a=0.2,b=NULL))

lst.store <- list()
set.seed(1234)
i <- 0
for (d1 in dist.seq) {
  for (d2 in dist.seq) {
    i <- i + 1
    print(sprintf('Dist x: %s, dist y: %s, iteration: %i of %i',d1,d2,i,length(dist.seq)^2))
    params.d1 <- param.seq[[d1]]
    params.d2 <- param.seq[[d2]]
    sns.hat <- rep(NA,nrep)
    sns.skew.hat <- rep(NA,nrep)
    for (k in seq(nrep)) {
      if (k %% 5000 == 0) print(k)
      sns.hat[k] <- sns.fun(d.x=d1,d.y=d2,n=n,a.x=params.d1$a,b.x=params.d1$b,a.y=params.d2$a,b.y=params.d2$b)
      sns.skew.hat[k] <- sns.skew.fun(d.x=d1,d.y=d2,n=n,a.x=params.d1$a,b.x=params.d1$b,a.y=params.d2$a,b.y=params.d2$b)
    }
    lst.store[[i]] <- data.table(scale=sns.hat,skew=sns.skew.hat,d_x=d1,d_y=d2)
    # names(lst.store)[i] <- str_c(d1,d2,sep='_')
  }
}

all.store <- melt(rbindlist(lst.store),id.vars=c('d_x','d_y'),variable.name='adj')
all.store <- all.store[order(d_x,d_y,adj,value)]
all.store[,pseq := seq(1e-6,1-1e-6,length.out=.N) , by=list(d_x,d_y,adj)]
all.store[,qseq := qnorm(pseq) , by=list(d_x,d_y,adj)]

# Store the data
fwrite(x=all.store,file=file.path(dir.output,'sim_qq_rwsns.csv'))

end_time <- Sys.time()

print(sprintf('Run time was: %0.2f minutes',as.numeric(difftime(end_time,start_time,units='mins'))))

