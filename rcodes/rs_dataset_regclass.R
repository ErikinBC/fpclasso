##########################################################################
# -------- SCRIPT TO COMBINE PYTHON DATASETS IN .RDATA FILES ----------- #
##########################################################################

rm(list=ls())

# Libraries for optparse
library(stringr,logical.return=F,warn.conflicts=F,quietly=T,verbose=F)
library(optparse,logical.return=F,warn.conflicts=F,quietly=T,verbose=F)

# Optional linux argument
option_list <- list(
  make_option('--dir_output',type='character',default=getwd(),help='Path to output directory [default WD]')
);

# Parse and assign
opt <- parse_args(OptionParser(option_list=option_list))
dir.output <- opt$dir_output
stopifnot(dir.exists(dir.output))
# dir.output='/home/erik/Documents/projects/FPC_Lasso/current/output/data'
# ii=1

# Remaining packages
pckgs <- c('data.table','forcats')
for (pp in pckgs) { library(pp,character.only = T) }

#####################################################
###### -------- SUPPORT FUNCTIONS ---------- ########

# Proportion of unique values
fun.prop.unique <- function(X) { apply(X, 2, function(cc) length(unique(cc))) / nrow(X) }

# Return the frequency of the second most frequent column
fun.cat.bal <- function(X) {
  tmp.p <- apply(X,2,function(cc) sort(prop.table(table(cc)),decreasing = T)[2] )
  tmp.p <- ifelse(is.na(tmp.p), 0, tmp.p)
  return(tmp.p)  
}

# Return the ratio of the first to second most frequent value
fun.cat.ratio <- function(X) {
  tmp.p <- apply(X,2,function(cc) sort(prop.table(table(cc)),decreasing = T))
  if (class(tmp.p)=='matrix') {
    tmp.r <- apply(tmp.p,2,function(cc) cc[1]/cc[2] )
  } else if (class(tmp.p)=='list') {
    tmp.r <- unlist(lapply(tmp.p, function(ll) ll[1]/ll[2]))  
  } else {
    stopifnot(F)
  }
  return(tmp.r)
}


##################################################
###### -------- CLASSIFICATION ---------- ########

# Get the list of the files
fn.class <- str_subset(str_subset(list.files(dir.output),'^class'),'\\.csv$')
stopifnot(length(fn.class)>10)

# Loop through and load the classification datasets
lst.class <- vector('list',length(fn.class))
for (ii in seq_along(lst.class)) {
  ff <- fn.class[ii]
  print(sprintf('Iteration %i: %s',ii,ff))
  tmp.df <- fread(file.path(dir.output,ff))
  if(colnames(tmp.df)[1]=='V1'){ tmp.df[, V1 := NULL]   }
  colnames(tmp.df) <- str_replace_all(tolower(colnames(tmp.df)),'\\s','_')
  tmp.y <- tmp.df$y
  tmp.df[, y:= NULL]
  tmp.n <- nrow(tmp.df)
  # Continuous: anything with >5% unique values
  tmp.prop <- fun.prop.unique(tmp.df)
  idx.cat <- which(tmp.prop < 0.05)
  cn.cat <- names(idx.cat)
  if (length(cn.cat) > 0) {
    # Convert to factors
    tmp.df[, (cn.cat) := lapply(.SD, function(ll) as.factor(as.character(ll))),.SDcols=cn.cat]
    # Aggregate if any factor is less than 5%
    bal.cat <- fun.cat.bal(X=tmp.df[,cn.cat,with=F])
    drop.cat <- names(which(bal.cat < 0.05))
    if (length(drop.cat) > 0) { tmp.df[, (drop.cat) := NULL] }
    keep.cat <- setdiff(cn.cat, drop.cat)
    if (length(keep.cat) >0 ) {
      tmp.df[, (keep.cat) := lapply(.SD, function(ll) fct_lump(ll,prop=0.05)),.SDcols=keep.cat]
    }
  }
  tmp.X <- model.matrix(~.,data=tmp.df)[,-1]
  # Get ratio of first to second most common value
  tmp.ratio <- fun.cat.ratio(tmp.X)
  # Drop anything that has a balance worse than 5% or < 100 observations
  idx.keep <- (tmp.ratio < 20) | (floor((1/tmp.ratio) * tmp.n) > 100)
  tmp.X <- tmp.X[,idx.keep]
  # Save in the list
  lst.class[[ii]] <- list(y=tmp.y, X=tmp.X)
}
# Assign nice names
names(lst.class) <- str_remove_all(fn.class,'class\\_|\\.csv')

##############################################
###### -------- REGRESSION ---------- ########

# Get the list of the files
fn.reg <- str_subset(str_subset(list.files(dir.output),'^reg'),'\\.csv$')
stopifnot(length(fn.reg)>10)

# Loop through and load the regification datasets
lst.reg <- vector('list',length(fn.reg))
for (ii in seq_along(lst.reg)) {
  ff <- fn.reg[ii]
  print(sprintf('Iteration %i: %s',ii,ff))
  tmp.df <- fread(file.path(dir.output,ff))
  if(colnames(tmp.df)[1]=='V1'){ tmp.df[, V1 := NULL]   }
  colnames(tmp.df) <- str_replace_all(tolower(colnames(tmp.df)),'\\s','_')
  tmp.y <- tmp.df$y
  tmp.df[, y:= NULL]
  tmp.n <- nrow(tmp.df)
  # Continuous: anything with >5% unique values
  tmp.prop <- fun.prop.unique(tmp.df)
  idx.cat <- which(tmp.prop < 0.05)
  cn.cat <- names(idx.cat)
  if (length(cn.cat) > 0) {
    # Convert to factors
    tmp.df[, (cn.cat) := lapply(.SD, function(ll) as.factor(as.character(ll))),.SDcols=cn.cat]
    # Aggregate if any factor is less than 5%
    bal.cat <- fun.cat.bal(X=tmp.df[,cn.cat,with=F])
    drop.cat <- names(which(bal.cat < 0.05))
    if (length(drop.cat) > 0) { tmp.df[, (drop.cat) := NULL] }
    keep.cat <- setdiff(cn.cat, drop.cat)
    if (length(keep.cat) >0 ) {
      tmp.df[, (keep.cat) := lapply(.SD, function(ll) fct_lump(ll,prop=0.05)),.SDcols=keep.cat]
    }
  }
  tmp.X <- model.matrix(~.,data=tmp.df)[,-1]
  # Get ratio of first to second most common value
  tmp.ratio <- fun.cat.ratio(tmp.X)
  # Drop anything that has a balance worse than 5% or < 100 observations
  idx.keep <- (tmp.ratio < 20) | (floor((1/tmp.ratio) * tmp.n) > 100)
  tmp.X <- tmp.X[,idx.keep]
  # Save in the list
  lst.reg[[ii]] <- list(y=tmp.y, X=tmp.X)
}
# Assign nice names
names(lst.reg) <- str_remove_all(fn.reg,'reg\\_|\\.csv')

#######################################################
###### -------- SAVE DATA AND CLEAN ---------- ########

n.class <- unlist(lapply(lst.class,function(ll) nrow(ll$X)))
p.class <- unlist(lapply(lst.class,function(ll) ncol(ll$X)))
n.reg <- unlist(lapply(lst.reg,function(ll) nrow(ll$X)))
p.reg <- unlist(lapply(lst.reg,function(ll) ncol(ll$X)))

tmp.dim <- rbind(data.table(dataset=names(lst.reg),n=n.reg,p=p.reg,type='regression'),
      data.table(dataset=names(lst.class),n=n.class,p=p.class,type='classifcation'))
print(tmp.dim[order(type,-n)])

save(lst.class, file=file.path(dir.output,'class_datasets.RData'))
save(lst.reg, file=file.path(dir.output,'reg_datasets.RData'))
print(sprintf('Class list: %s',format(object.size(lst.class),'MB')))
print(sprintf('Class list: %s',format(object.size(lst.reg),'MB')))

# Remove the CSVs
sapply(file.path(dir.output,fn.class),function(ss) file.remove(ss) )
sapply(file.path(dir.output,fn.reg),function(ss) file.remove(ss) )

###################################################
###### -------- POST PROCESSING ---------- ########

# Load the class again
load(file=file.path(dir.output,'class_datasets.RData'))

# Remove the parity datasets
lst.class <- lst.class[-str_which(names(lst.class),'parity')]
# Encode y as 0/1
for (ss in names(lst.class)) {
  print(ss)
  lst.class[[ss]]$y <- as.numeric(as.factor(lst.class[[ss]]$y))-1
}

# Re-save
save(lst.class, file=file.path(dir.output,'class_datasets.RData'))




