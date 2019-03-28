#############################################################################################################
############# ----------- SCRIPT TO GENERATE FIGURES FOR DIFFERENT SIMULATIONS --------------- ##############
#############################################################################################################

rm(list=ls())

# Libraries for optparse
library(stringr,logical.return=F,warn.conflicts=F,quietly=T,verbose=F)
library(optparse,logical.return=F,warn.conflicts=F,quietly=T,verbose=F)

# Optional linux argument
option_list <- list(
  make_option('--dir_output',type='character',default=getwd(),help='Path to output directory [default WD]'),
  make_option('--dir_figs',type='character',default=getwd(),help='Path to store figures [default WD]')
);

# Parse and assign
opt <- parse_args(OptionParser(option_list=option_list))
dir.output <- opt$dir_output
dir.figs <- opt$dir_figs
stopifnot(dir.exists(dir.output),dir.exists(dir.figs))

# dir.output="/home/erik/Documents/projects/FPC_Lasso/current/output"
# dir.figs="/home/erik/Documents/projects/FPC_Lasso/current/output/figures"

# Load in the other packages
pckgs <- c('data.table','stringr','forcats','cowplot')
for (pp in pckgs) { library(pp,character.only = T) }

# default ggplot color palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

####################################################
###### ---- (1) rsSNS QQ-PLOTS ------------ ########

# Load
all.store <- fread(file=file.path(dir.output,'sim_qq_rwsns.csv'))

d_map <- c(beta='Beta(2,1)',binomial1='Bin(0.5)',binomial2='Bin(0.2)',
           exponential='Exp(1)',gamma='Gamma(5,0.5)',gaussian='N(0,1)')

# Now make the ggplots
gg.sns.qq <-
  ggplot(all.store[adj == 'skew'],aes(x=qseq,y=value)) + 
  geom_point(size=1,alpha=0.5) + 
  facet_grid(d_y~d_x,labeller = labeller(d_y=d_map, d_x = d_map)) + 
  geom_abline(slope=1,color='red',alpha=0.5) + 
  scale_x_continuous(limits=c(-5,5),breaks = seq(-5,5,2.5)) + 
  scale_y_continuous(limits=c(-5,5),breaks = seq(-5,5,2.5)) + 
  labs(x='Theoretical quantiles',y='Sample quantiles')

# Add the titles
gg.sns.qq <- plot_grid(plot_grid(ggplot() + draw_label(label=expression('Distribution of:  ' * X[p]),fontface='bold',size=12),
                                 gg.sns.qq,ncol=1,rel_heights = c(1,20)),
                       ggplot() + draw_label(label=expression('Distribution of:  ' * hat(epsilon)[-p]),fontface='bold',size=12,angle=270),
                       rel_widths = c(20,1))

# save
save_plot(filename=file.path(dir.figs,'gg_qq_rwsns.png'),plot=gg.sns.qq,base_height=12,base_width=12)

#####################################################################
###### ---- (4) EMPIRICAL REGULARITY CONDITIONS ------------ ########

# Load in the regularity condition data
df.merge <- fread(file.path(dir.output,'df_lamseq.csv'))

bh <- bw <- 15

# Survival plot
gg.survival <- ggplot(df.merge[type=='survival'],aes(x=lam_scale,y=log(lam_fpc))) + 
  geom_point(color=gg_color_hue(3)[1]) + 
  facet_wrap(~dataset,scales = 'free') + 
  theme(legend.position = 'none') + 
  background_grid(major='xy',minor='none') + 
  labs(y=expression(log(lambda^FPC)), x=expression(frac(lambda-lambda[min],lambda[max]-lambda[min])))
save_plot(filename = file.path(dir.figs,'gg_lamseq_surv.png'),plot=gg.survival,base_height=bh,base_width=bw)

# Classification plot
gg.class <- ggplot(df.merge[type=='class'],aes(x=lam_scale,y=log(lam_fpc))) + 
  geom_point(color=gg_color_hue(3)[2]) + 
  facet_wrap(~dataset,scales = 'free') + 
  theme(legend.position = 'none') + 
  background_grid(major='xy',minor='none') + 
  labs(y=expression(log(lambda^FPC)), x=expression(frac(lambda-lambda[min],lambda[max]-lambda[min])))
save_plot(filename = file.path(dir.figs,'gg_lamseq_class.png'),plot=gg.class,base_height=bh,base_width=bw)

# Regression plot
gg.reg <- ggplot(df.merge[type=='reg'],aes(x=lam_scale,y=log(lam_fpc))) + 
  geom_point(color=gg_color_hue(3)[3]) + 
  facet_wrap(~dataset,scales = 'free') + 
  theme(legend.position = 'none') + 
  background_grid(major='xy',minor='none') + 
  labs(y=expression(log(lambda^FPC)), x=expression(frac(lambda-lambda[min],lambda[max]-lambda[min])))
save_plot(filename = file.path(dir.figs,'gg_lamseq_reg.png'),plot=gg.reg,base_height=bh,base_width=bw)

####################################################################
###### ---- (5A) FALSE POSITIVE CONTROL: INDEP ------------ ########

# Plotting parameters
pd <- position_dodge(0.8)
d1 <- c(gaussian='Gaussian',binomial='Binomial',cox='Cox-PH')
d2 <- c(gaussian='Gaussian',binomial='Binomial',exponential='Exponential')

# Load data
storage.indep <- fread(file.path(dir.output,'storage_fp_indep.csv'))

# Aggregate to the average
sim.indep <- storage.indep[,list(actual=mean(value),q25=quantile(value,1/4),q75=quantile(value,3/4)),by=list(expected,p,xfam,yfam,type)]
sim.indep[,`:=` (xfam = lvls_reorder(xfam,c(3,1,2)), yfam = lvls_reorder(yfam,c(3,1,2)))]

# Create the false positive plot
gg.fp.indep <-
  ggplot(sim.indep[type=='fp'],aes(x=factor(expected),y=actual,color=factor(p))) + 
  geom_point(position = pd,size=3,shape=17) + 
  geom_linerange(aes(ymin=q25,ymax=q75),position=pd) + 
  facet_grid(xfam~yfam,labeller = labeller(xfam=d2,yfam=d1)) +
  scale_y_continuous(breaks = c(1,3,5,10)) +
  background_grid(major='xy',minor='none') +
  labs(x='Expected # FPs',y='Actual (average)',subtitle='Lines shows interquartile range') +
  scale_color_manual(name='Dimensionality: ',values=c('lightblue','blue','darkblue')) + 
  theme(legend.position = 'bottom',legend.justification = 'center',plot.subtitle = element_text(size=10))
# Add labels
gg.fp.indep <- plot_grid(plot_grid(ggplot() + draw_label(label='Response class dist.',size=12,fontface = 'bold'),
                      gg.fp.indep,ncol=1,rel_heights = c(1,20)),
            ggplot() + draw_label(label='Design matrix dist.',size=12,angle = 270,fontface = 'bold'),
            rel_widths = c(20,1))
# save
save_plot(filename = file.path(dir.figs,'gg_fp_indep.png'),plot=gg.fp.indep,base_height = 10,base_width = 10)

# Create the power plot
gg.tp.indep <-
  ggplot(sim.indep[type=='tp'],aes(x=factor(expected),y=actual,color=factor(p))) + 
  geom_point(position = pd,size=3,shape=17) + 
  geom_linerange(aes(ymin=q25,ymax=q75),position=pd) + 
  facet_grid(xfam~yfam,labeller = labeller(xfam=d2,yfam=d1)) +
  # scale_y_continuous(breaks = c(1,3,5,10)) +
  background_grid(major='xy',minor='none') +
  labs(x='Expected # FPs',y='Actual TPs (average)',subtitle='Lines shows interquartile range') +
  scale_color_manual(name='Dimensionality: ',values=c('lightgreen','chartreuse3','darkgreen')) + 
  theme(legend.position = 'bottom',legend.justification = 'center',plot.subtitle = element_text(size=10))
gg.tp.indep <- plot_grid(plot_grid(ggplot() + draw_label(label='Response class dist.',size=12,fontface = 'bold'),
                                   gg.tp.indep,ncol=1,rel_heights = c(1,20)),
                         ggplot() + draw_label(label='Design matrix dist.',size=12,angle = 270,fontface = 'bold'),
                         rel_widths = c(20,1))
# save
save_plot(filename = file.path(dir.figs,'gg_tp_indep.png'),plot=gg.tp.indep,base_height = 10,base_width = 10)

###################################################################
###### ---- (5B) FALSE POSITIVE CONTROL: CORR ------------ ########

# Load data
storage.corr <- fread(file.path(dir.output,'storage_fp_corr.csv'))

# Aggregate to the average
sim.corr <- storage.corr[,list(actual=mean(value),q25=quantile(value,1/4),q75=quantile(value,3/4)),
                         by=list(expected,corr,xfam,yfam,type)]
sim.corr[,`:=` (xfam = lvls_reorder(xfam,c(3,1,2)), yfam = lvls_reorder(yfam,c(3,1,2)))]

# Create the false positive plot
gg.fp.corr <-
  ggplot(sim.corr[type=='fp'],aes(x=factor(expected),y=actual,color=factor(corr))) + 
  geom_point(position = pd,size=3,shape=17) + 
  geom_linerange(aes(ymin=q25,ymax=q75),position=pd) + 
  facet_grid(xfam~yfam,labeller = labeller(xfam=d2,yfam=d1)) +
  scale_y_continuous(breaks = c(1,3,5,10)) +
  background_grid(major='xy',minor='none') +
  labs(x='Expected # FPs',y='Actual (average)',subtitle='Lines shows interquartile range') +
  scale_color_manual(name='Columnwise correlation: ',values=c('pink','red','darkred')) + 
  theme(legend.position = 'bottom',legend.justification = 'center',plot.subtitle = element_text(size=10))
gg.fp.corr <- plot_grid(plot_grid(ggplot() + draw_label(label='Response class dist.',size=12,fontface = 'bold'),
                      gg.fp.corr,ncol=1,rel_heights = c(1,20)),
            ggplot() + draw_label(label='Design matrix dist.',size=12,angle = 270,fontface = 'bold'),
            rel_widths = c(20,1))

# save
save_plot(filename = file.path(dir.figs,'gg_fp_corr.png'),plot=gg.fp.corr,base_height = 10,base_width = 10)

# NOTE!!!! NO CHANGE IN POWER [THIS IS A WEAKNESS AS THERE IS INFORMATION LEFT ON THE TABLE....]



