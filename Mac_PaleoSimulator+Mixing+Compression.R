# Simulate a paleo record with mixing and compression of core layers using
#  the pseudocode in 'Pseudocode_Mixing+Compression_2016-08-25.pdf'
# Copyright Stephen R. Carpenter 2017

rm(list = ls())
graphics.off()

library('zoo')
library('tsDyn')

# Functions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Solve the model given parameters ----------------------------------------------------------

# 2D rate funtion for simulation; deterministic parts omitting the load
dW.noload = function(W) {
  Rden=mq + W^q
  # deterministic part of dW/dt w/o load 
  rate = -(s+h)*W + (r*M.frozen*W^q/Rden)
  return(rate)
}

Trajectory = function(lamda,nt) {  # lamda is a vector of load over time; unlisted parameters are implicit
  # Preliminaries
  env.noise = rnorm(nt)
  Wt=rep(0,nt)
  Wt[1]=W0 
  
  for(i in 2:nt)  {
    rate = dW.noload(Wt[i-1])
    Wnext = Wt[i-1] + (lamda[i] + rate)*dt + dtnoise*CV*Wt[i-1]*env.noise[i]
    Wt[i] = max(Wnext,0.1)  # set a floor for water P
  }
  return(Wt)
}

# Useful statistics ----------------------------------------------------------------------------------

# Compute autocorrelation time
ACtime = function(x) {
  zzz=acf(x,lag=1,plot=FALSE)
  zz=zzz$acf[2]
  ACt=-1/log(zz)
  return(ACt)
}

# Compute autocorrelation
AClag1 = function(x) {
  zzz=acf(x,lag=1,plot=FALSE)
  zz=zzz$acf[2]
  return(zz)
}

# End functions >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Main program

# Deterministic Parameters 
b = 0.002 # this rate is increased from previous papers; 
#           burial is 0.001 in Carpenter & Brock Ecology Letters 2006
h = 1-exp(-0.29) # export coef. from WSC model (Motew et al. Ecosystems 2017)
m = 4 # half-saturation for recycle; wide range 1-15 in Carpenter & Lathrop Ecosystems 2008
r = 0.019 # max. recycling rate; ~0.002 in Carpenter & Lathrop 2008; 0.019 in Carpenter & Brock 2006
q = 4 # 4 is near posterior mode in Carpenter & Lathrop Ecosystems 2008
mq = m^q  
s = 1-h # sedimentation from WSC model (Motew et al. Ecosystems 2017)

# Freeze sediment P at a constant value that is close to the bifurcation point
# This converts the 3-dimensional model of Carpenter & Brock 2006 Ecology Letters 
#  to a one-dimensional model that is easier to use for our purposes.
M.frozen = 400

# Noise process
# CV is multiplied by current mean to obtain a standard deviation
CV = 0.3 # Values of 0.3-0.4 are sedate; flickering starts around 0.45 and gets larger with larger CV

# Simulation control
nt.all = 10000
dt = 0.1 # time step in years
dtnoise = sqrt(dt) # stochastic part time step

# Initial conditions
W0 = 1

# Set up load (lamda) sequence
lamda.vec = seq(0.8,2,length.out=nt.all)

# Calculate a time series --------------------------------------------------------------------------
Wsim.all = Trajectory(lamda.vec,nt.all)

# Plot
Tvec.all = (1:nt.all)
quartz()
par(mfrow=c(1,1),mar=c(5, 4.3, 2, 2) + 0.1, cex.axis=1.6,cex.lab=1.6)
plot(Tvec.all,Wsim.all,type='l',lwd=2,col='blue',xlab='time step',ylab='Water P',
     main='Time Series at step dt')

# Extract annual sample  ----------------------------------------------------------------------------
n.year = round(1/dt)
ny = nt.all/n.year
Wsim.all.mat = matrix(Wsim.all,nr=ny,nc=n.year,byrow=T)
Wsim.ann = apply(Wsim.all.mat,1,mean)
T.ann = seq(10,nt.all,by=n.year)

# Overplot annual sample on original data
quartz()
par(mfrow=c(1,1),mar=c(5, 4.3, 2, 2) + 0.1, cex.axis=1.6,cex.lab=1.6)
plot(Tvec.all,Wsim.all,type='l',lwd=2,col='blue',xlab='time step',ylab='Water P',
     main='Annual sample red, full sample blue')
points(T.ann,Wsim.ann,type='l',lwd=1,col='red')

# Extract a subset of 400 years, around the switchpoint  ------------------------------------------------

# Use SETAR to find a point before or near the switchpoint.
fsetar = setar(Wsim.ann,nthresh=1,mL=2,mH=2)
print('Summary of setar fit',quote=F)
print(summary(fsetar))
quartz()
plot(fsetar,ask=F) # use ask=F to turn off plots

# setar switchpoint
par.setar=fsetar$coefficients
thresh = par.setar[7]
print('',quote=F)
print(c('SETAR switchpoint for Water P = ',thresh),quote=F)

# find first point between years 300 and 900 where WaterP crosses threshold going up
W300 = Wsim.ann[301:(ny-101)]
W300pre = Wsim.ann[299:(ny-103)] # it is more robust to check slope over two intervals
dW300 = W300-W300pre
upcross = rep(0,length(dW300))
upcross = ifelse(W300>thresh & W300pre<thresh & dW300>0,1,0)
iswitch=300+which.max(upcross)

print(c('setar estimate of switchpoint',iswitch),quote=F)

quartz()
par(mfrow=c(1,1),mar=c(5, 4.3, 2, 2) + 0.1, cex.axis=1.6,cex.lab=1.6)
plot(T.ann/10,Wsim.ann,type='l',lwd=2,col='blue',xlab='time step',ylab='Water P',
     main='setar estimate of switchpoint')
abline(v=iswitch,lwd=3,col='red')
abline(h=thresh,lwd=2,lty=2,col='magenta')

lo.sub = iswitch-300
hi.sub = iswitch+99
Wsim = Wsim.ann[lo.sub:hi.sub] 
nt=length(Wsim)
Tvec = (1:nt)

quartz()
par(mfrow=c(1,1),mar=c(5, 4.3, 2, 2) + 0.1, cex.axis=1.6,cex.lab=1.6)
plot(Tvec,Wsim,type='l',lwd=2,col='blue',xlab='time step',ylab='Water P',
     main='Annual subsample')

# Check time series statistics on the annual sample ************************************************
# Compute rolling window variance ------------------------------------------------------------
winlen = 30
var.WP = rollapply(Wsim,winlen,var,align='right')
Tvec.var = (winlen:nt) 
# Compute rolling window autocorrelation ----------------------------------------------------
winlen = 30
AC.WP = rollapply(Wsim,winlen,ACtime,align='right')
# Add to plot
Tvec.AC = (winlen:nt)

# Plot time series and typical resilience indicators
quartz(width=6,height=6)
par(mfrow=c(3,1),mar=c(4, 4.3, 2, 2) + 0.1, cex.axis=1.5,cex.lab=1.5)
plot(Tvec,Wsim,type='l',lwd=2,col='blue',xlab='time step',ylab='Water P',log='y',
     main='Statistics for Series Before Mixing & Compression')
abline(v=300)
plot(Tvec.var,var.WP,type='l',lwd=2,col='red',xlab='time step',ylab='Variance',log='y')
abline(v=300)
plot(Tvec.AC,AC.WP,type='l',lwd=2,col='magenta',xlab='time step',ylab='AC time')
abline(v=300)

# Obtain "core" sample and analyze it ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Mix the time series ####

sigma.mix = 2.5  # assume a sigma for mixing  

# Build the weight vector by integrating over intervals that match annual time steps
# This is done 'by hand' and could be improved
LL = c(-0.5,0.5,1.5,2.5,3.5,4.5) # Lower limits of intervals for weights
UL = c(0.5,1.5,2.5,3.5,4.5,5) # Upper limits of intervals for weights
# Compute the weights on the right hand side of the normal curve, counting the center interval
Weight.right = rep(0,length(LL))
# Function for normal density given mixing s.d.
Ndens = function(x) {dnorm(x,mean=0,sd=sigma.mix)}
# Integrate for intervals
for(i in 1:length(LL)) {
  W.int = integrate(Ndens,lower=LL[i],upper=UL[i])
  Weight.right[i] = W.int$value
}
# Make the left side
Weight.left=rev(Weight.right[2:length(LL)])
# Make the full weight vector for mixing, standardized to sum to one
Weight.mix = c(Weight.left,Weight.right) / sum(c(Weight.left,Weight.right))
NW = length(Weight.mix)

print('',quote=F)
print('Mixing Weights',quote=F)
print(Weight.mix)

# Mix the time series Wsim corresponding to nt time steps in Tvec
Wmix = rep(0,nt) # Vector to hold the mixed values

# Mix the time steps that will use the full weight vector
Nhalf = length(LL)-1 # number of time steps that will use only part of the weight vector
for(i in (Nhalf+1):(nt-Nhalf)) {
  lo = i-5
  hi = i+5
  simvec = Wsim[lo:hi]
  Wmix[i] = simvec%*%Weight.mix
}

# Mix the steps that use only part of the weight vector

# At the start of the time series
for(i in 1:Nhalf) {
  hi = i+Nhalf
  simvec = Wsim[1:hi]
  lo = (Nhalf+1)-(i-1)
  Wvec = Weight.mix[lo:NW]
  Wmix[i] = (simvec%*%Wvec)/sum(Wvec)
}

# At the end of the time series
count = 0
for(i in (nt-Nhalf):nt) {
  count = count+1
  lo = i-Nhalf
  simvec = Wsim[lo:nt]
  hi = NW-(count-1)
  Wvec = Weight.mix[1:hi]
  Wmix[i] = (simvec%*%Wvec)/sum(Wvec)
}

# Plot unmixed and mixed time series
quartz()
par(mfrow=c(1,1),mar=c(5, 4.3, 2, 2) + 0.1, cex.axis=1.6,cex.lab=1.6)
plot(Tvec,Wsim,type='l',lwd=2,col='blue',xlab='time step',ylab='Water P',
     main='True Series blue, Mixed Series red')
points(Tvec,Wmix,type='l',lwd=2,col='red')
legend('topleft',legend=c('mixed','not mixed'),
       lwd=c(1,2),col=c('red','blue'),
       cex=1.5,bty='n')

# Compress the time series into core slices ####
# Assume compression parameters; see 'Pseudocode_Mixing+Compression_2016-08-25.doc'
# Analysis of ~100 cores indicated provided the following summary statistics:
# a0: 1st Qu. = 1.2770, Median = 2.2400, 3rd Qu. = 4.3510  
# a1: 1st Qu. = 0.01501, Median =  0.03409, 3rd Qu. = 0.08815 

# Example shows a case of moderate compression
a0 = 2 # assume 2 years / cm at the top of the core
a1 = 0.025  # less compression than Lake Champlain (a1 ~ 0.1)

# Compute core length in cm
L.core = (1/a1)*( log(a1*nt + a0) - log(a0))

print('',quote=F)
print(c(nt,'annual time steps compressed to core of',round(L.core,2),'cm'),quote=F)

# Discard non-integer fraction and end of core
L.core = floor(L.core)

# Average tracer concentrations for each 1 cm slice of core
Wmix.rev = rev(Wmix) # reverse the time series to start at top of core
# Wmix.rev = rev(Wsim) # for no mixing case
core = rep(0,L.core)
t0 = 0 # start at the top of the core
t1 = 0
# save matrix of t0 and t1 values to check results
t01mat = matrix(0,nr=L.core,nc=2)

# Do the first slice separately to accommodate starting at exactly 0
# Compute new t1 from forumula in 'Pseudocode_Mixing+Compression_2016-08-25.doc'
t1 = (1/a1)*( (a1*t0 + a0)*exp(a1) - a0)
# Save t0 and t1
t01mat[1,]=c(t0,t1)
# average the tracer concentrations between t0 and t1
y.slice = ceiling(t1)-floor(t0)
wts = rep(1,y.slice)
x = Wmix.rev[1:ceiling(t1)]
wts[1] = ceiling(t0)-t0 # fraction of the first year in the slice
wts[y.slice] = t1 - floor(t1) # fraction of the final year in the slice
core[1] = x%*%wts/sum(wts)

# Do the remaining slices
for(i in 2:L.core) {  # loop over core slices
  t0 = t1 # previous t1 is now t0
  # Compute new t1 from forumula in 'Pseudocode_Mixing+Compression_2016-08-25.doc'
  t1 = (1/a1)*( (a1*t0 + a0)*exp(a1) - a0)
  # Save t0 and t1
  t01mat[i,]=c(t0,t1)
  # average the tracer concentrations between t0 and t1
  y.slice = ceiling(t1)-floor(t0)+1
  wts = rep(1,y.slice)
  x = Wmix.rev[floor(t0):ceiling(t1)]
  wts[1] = ceiling(t0)-t0 # fraction of the first year in the slice
  wts[y.slice] = t1 - floor(t1) # fraction of the final year in the slice
  core[i] = x%*%wts/sum(wts)
}

# Compute times at the center of each core slice
t01mat[1,1]=1.e-3 # replace the zero to prevent problem with log
T.core = exp( 0.5*(log(t01mat[,1])+log(t01mat[,2])) )
# reverse the sequence of times to correspond with the simulated series
T.core = nt - T.core

# Overplot core samples with mixing and compression and true time series
quartz()
par(mfrow=c(1,1),mar=c(5, 4.3, 2, 2) + 0.1, cex.axis=1.6,cex.lab=1.6)
plot(Tvec,Wsim,type='l',lwd=2,col='blue',xlab='time step',ylab='Water P',
     main='True Series blue, Mixed & Compressed Series red')
points(T.core,core,type='p',pch=19,col='red',cex=1.3)
abline(v=300)

# Save core sample and true time series for other analyses
# Tvec and Wsim are time and true biomass; 
# T.core and core are the mixed & compressed core sample
# 400 years, switch occurs year 300
# REMOVE # TO ACTIVATE SAVE
#save(Tvec,Wsim,T.core,core,file='Core_simulated.Rdata')

# Compute resilience statistics ####
# for the part of the time series before the jump
core.pre = subset(core,subset=(T.core<300))
n.pre = length(core.pre)
T.core.pre = T.core[(L.core-n.pre+1):L.core]

winlen = length(core.pre)/1.95
var.core = rollapply(core.pre,winlen,var,align='right')
AC.core = rollapply(core.pre,winlen,AClag1,align='right')

quartz()
par(mfrow=c(2,1),mar=c(5, 4.3, 2, 2) + 0.1, cex.axis=1.6,cex.lab=1.6)
plot(T.core.pre[winlen:n.pre],var.core,type='l',lwd=2,col='red',xlab='time before switch',ylab='Variance',
     main='Time Steps Before Switch are Plotted')
plot(T.core.pre[winlen:n.pre],AC.core,type='l',lwd=2,col='magenta',xlab='time before switch',ylab='AC')

