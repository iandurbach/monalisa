# This is code to create the various density plots from a Poisson point process 
# with linear trend, which are used in the Introduction.
# ---------------------------------------------------------------- dlb Oct 2018

library(secr)
library(scrmlebook)
library(sp)

usage.sigma = 5
error.sigma = 2.5
dmax = 60
errmult = 6 # ratio ofmax error to error.sigma
max.error.sigma = errmult*error.sigma
buffer = 30


traps = make.grid()
mask = make.mask(traps,buffer=buffer)
plot(mask)
plot(traps,add=TRUE)
summary(mask)
boundpoly = Polygons(list(Polygon(attributes(mask)$boundingbox)),"boundary")
boundary = SpatialPolygons(list(boundpoly))

tbbox = data.frame(x=c(min(traps$x),max(traps$x),max(traps$x),min(traps$x)),
                   y=c(min(traps$y),min(traps$y),max(traps$y),max(traps$y)))
trapboundpoly  = Polygons(list(Polygon(tbbox)),"trapboundary")
trapboundary = SpatialPolygons(list(trapboundpoly))
traplim = bbox(trapboundary)

D.linear.x = function(mask,Dleft,Dright) {
  xrange = range(mask$x)
  dx = diff(xrange)
  slope = (Dright-Dleft)/dx
  D = Dleft + slope*(mask$x-xrange[1])
}

distance = function(to,from) return(sqrt((to[1]-from[1])^2 + (to[2]-from[2])^2))

addnormal = function(ac,mask,sigma,lambda=1) {
  dists = unlist(apply(mask,1,distance,from=ac))
  dens = lambda*exp(-dists^2/(2*sigma^2))
  dens = dens/sum(dens)
  return(dens)
}


# Creat liner trend ac density surface and add to mask
covariates(mask)$D.0 = D.linear.x(mask,0,100)
# Make density surface a SpatialPixelsDataFrame, and create boundary polygon
spdf.D = SpatialPixelsDataFrame(as.matrix(mask),data=data.frame(D=covariates(mask)$D.0))
Dbbox = bbox(spdf.D)
Dbdf = data.frame(x=c(Dbbox[1,1],Dbbox[1,2],Dbbox[1,2],Dbbox[1,1]),
                  y=c(Dbbox[2,1],Dbbox[2,1],Dbbox[2,2],Dbbox[2,2]))
Dboundpoly  = Polygons(list(Polygon(Dbdf)),"bufferboundary")
Dboundary = SpatialPolygons(list(Dboundpoly))
# Make SpatialPixelsDataFrame for surface only inside trapboundary
inspdf.D = spdf.D[trapboundary,]

# Generate realisation of the point process
set.seed(1)
pop = sim.popn(covariates(mask)$D.0,mask,model2D="IHP",Ndist="fixed")
N = dim(pop)[1]
# make this population a SpatialPointsDataFrame
spop = SpatialPointsDataFrame(as.matrix(pop),data=data.frame(size=rep(1,N)))
# make the members inside trapboudary a SpatialPointsDataFrame
inpop=spop[trapboundary,]

# Create usage density surface (V slow!) and add do mask
ausage = rep(0,dim(mask)[1])
for(i in 1:dim(pop)[1]) ausage = ausage + addnormal(pop[i,],mask,sigma=usage.sigma) # This is slow!
covariates(mask)$usage = ausage
# Make usage surface a SpatialPixelsDataFrame
spdf.usage = SpatialPixelsDataFrame(as.matrix(mask),data=data.frame(D=covariates(mask)$usage))
inspdf.usage = spdf.usage[trapboundary,]


# ERROR-ONLY MAPS
# ==================
# Create unifrom error density surface (V slow!) and add do mask
acerr = rep(0,dim(mask)[1])
system.time(for(i in 1:dim(pop)[1]) acerr = acerr + addnormal(pop[i,],mask,sigma=error.sigma)) # This is slow!
covariates(mask)$acerr = acerr
# Make usage surface a SpatialPixelsDataFrame
spdf.acerr = SpatialPixelsDataFrame(as.matrix(mask),data=data.frame(D=covariates(mask)$acerr))
inspdf.acerr = spdf.acerr[trapboundary,]
plot(trapboundary)
plot(inspdf.acerr,col=terrain.colors(40),what="image",add=TRUE)
plot(inpop,pch=19,cex=0.25,col="darkgray",add=TRUE)

# Create unifrom error density surface with max.error.sigma (V slow!) and add do mask
maxacerr = rep(0,dim(mask)[1])
system.time(for(i in 1:dim(pop)[1]) maxacerr = maxacerr + addnormal(pop[i,],mask,sigma=max.error.sigma)) # This is slow!
covariates(mask)$maxacerr = maxacerr
# Make usage surface a SpatialPixelsDataFrame
spdf.maxacerr = SpatialPixelsDataFrame(as.matrix(mask),data=data.frame(D=covariates(mask)$maxacerr))
inspdf.maxacerr = spdf.maxacerr[trapboundary,]
plot(trapboundary)
plot(inspdf.maxacerr,col=terrain.colors(40),what="image",add=TRUE)
plot(inpop,pch=19,cex=0.25,col="darkgray",add=TRUE)

# Create increasing error density surface, up to max.error.sigma (V slow!) and add to mask
acerrtrend = rep(0,dim(mask)[1])
d = sigmad = rep(0,dim(pop)[1])
for(i in 1:dim(pop)[1]) d[i] = sqrt(sum((pop[i,]-errpts[1,])^2))
# make sigma increase up to sill of max.error.sigma at dmax:
for(i in 1:dim(pop)[1]) sigmad[i] = min(max.error.sigma,error.sigma + d[i]/dmax*(max.error.sigma-error.sigma))
system.time(for(i in 1:dim(pop)[1]) acerrtrend <- acerrtrend + addnormal(pop[i,],mask,sigma=sigmad[i])) # This is slow!
covariates(mask)$acerrtrend = acerrtrend
# Make usage surface a SpatialPixelsDataFrame
spdf.acerrtrend = SpatialPixelsDataFrame(as.matrix(mask),data=data.frame(D=covariates(mask)$acerrtrend))
inspdf.acerrtrend = spdf.acerrtrend[trapboundary,]
plot(trapboundary)
plot(inspdf.acerrtrend,col=terrain.colors(40),what="image",add=TRUE)
plot(inpop,pch=19,cex=0.25,col="gray",add=TRUE)


# Look at relationship between estimate of realised activity centre density, and x
dat = data.frame(x=inspdf.D@coords[,"x"], D=inspdf.D@data$D)
plot(dat)
lmfit = lm(D~x,data=dat)
summary(lmfit)

dat = data.frame(x=inspdf.acerr@coords[,"x"], D=inspdf.acerr@data$D)
plot(dat)
lmfit = lm(D~x,data=dat)
summary(lmfit)

dat = data.frame(x=inspdf.maxacerr@coords[,"x"], D=inspdf.maxacerr@data$D)
plot(dat)
lmfit = lm(D~x,data=dat)
summary(lmfit)

dat = data.frame(x=inspdf.acerrtrend@coords[,"x"], D=inspdf.acerrtrend@data$D)
plot(dat)
lmfit = lm(D~x,data=dat)
summary(lmfit)

# ERROR + USAGE MAPS
# ==================
# Create unifrom error + usage density surface (V slow!) and add do mask
acuseerr = rep(0,dim(mask)[1])
system.time(for(i in 1:dim(pop)[1]) acuseerr = acuseerr + addnormal(pop[i,],mask,sigma=error.sigma+usage.sigma)) # This is slow!
covariates(mask)$acuseerr = acuseerr
# Make usage surface a SpatialPixelsDataFrame
spdf.acuseerr = SpatialPixelsDataFrame(as.matrix(mask),data=data.frame(D=covariates(mask)$acuseerr))
inspdf.acuseerr = spdf.acuseerr[trapboundary,]
plot(trapboundary)
plot(inspdf.acuseerr,col=terrain.colors(40),what="image",add=TRUE)
plot(inpop,pch=19,cex=0.25,col="darkgray",add=TRUE)

# Create unifrom error + usage density surface with max.error.sigma (V slow!) and add do mask
maxacuseerr = rep(0,dim(mask)[1])
system.time(for(i in 1:dim(pop)[1]) maxacuseerr = maxacuseerr + addnormal(pop[i,],mask,sigma=max.error.sigma+usage.sigma)) # This is slow!
covariates(mask)$maxacuseerr = maxacuseerr
# Make usage surface a SpatialPixelsDataFrame
spdf.maxacuseerr = SpatialPixelsDataFrame(as.matrix(mask),data=data.frame(D=covariates(mask)$maxacuseerr))
inspdf.maxacuseerr = spdf.maxacuseerr[trapboundary,]
plot(trapboundary)
plot(inspdf.maxacuseerr,col=terrain.colors(40),what="image",add=TRUE)
plot(inpop,pch=19,cex=0.25,col="darkgray",add=TRUE)

# Create increasing error + usage density surface, up to max.error.sigma (V slow!) and add do mask
acuseerrtrend = rep(0,dim(mask)[1])
d = sigmad = rep(0,dim(pop)[1])
for(i in 1:dim(pop)[1]) d[i] = sqrt(sum((pop[i,]-errpts[1,])^2))
# make sigma increas up to sill of max.error.sigma at dmax:
for(i in 1:dim(pop)[1]) sigmad[i] = min(max.error.sigma+usage.sigma, error.sigma + d[i]/dmax*(max.error.sigma-error.sigma)+usage.sigma)
system.time(for(i in 1:dim(pop)[1]) acuseerrtrend <- acuseerrtrend + addnormal(pop[i,],mask,sigma=sigmad[i])) # This is slow!
covariates(mask)$acuseerrtrend = acuseerrtrend
# Make usage surface a SpatialPixelsDataFrame
spdf.acuseerrtrend = SpatialPixelsDataFrame(as.matrix(mask),data=data.frame(D=covariates(mask)$acuseerrtrend))
inspdf.acuseerrtrend = spdf.acuseerrtrend[trapboundary,]
plot(trapboundary)
plot(inspdf.acuseerrtrend,col=terrain.colors(40),what="image",add=TRUE)
plot(inpop,pch=19,cex=0.25,col="gray",add=TRUE)


# Example points to illustrate estimation error
# (Does not work well - don't think will use it)
# ---------------------------------------------
traplim.x = bbox(trapboundary)[1,]
traplim.y = bbox(trapboundary)[2,]
errpts = data.frame(x=c(mean(traplim.x),mean(traplim.x),mean(traplim.x),traplim.x[1]   ,traplim.x[2]   ,
                        traplim.x[1],traplim.x[1],traplim.x[2],traplim.x[2]),
                    y=c(mean(traplim.y),traplim.y[1]   ,traplim.y[2]   ,mean(traplim.y),mean(traplim.y),
                        traplim.y[1],traplim.y[2],traplim.y[1],traplim.y[2]))

# Uniform estimation error
errsunif = rep(0,dim(mask)[1])
for(i in 1:dim(errpts)[1]) errsunif = errsunif + addnormal(errpts[i,],mask,sigma=error.sigma) # This is slow!
covariates(mask)$errsunif = errsunif
spdf.errsunif = SpatialPixelsDataFrame(as.matrix(mask),data=data.frame(D=covariates(mask)$errsunif))
inspdf.errsunif = spdf.errsunif[trapboundary,]
# Plotting
plot(trapboundary)
plot(inspdf.errsunif,col=terrain.colors(40),what="image",add=TRUE)
points(errpts,pch="+")

# Linearly increasing estimation error
errtrend = rep(0,dim(mask)[1])
d = sigmad = rep(0,dim(errpts)[1])
for(i in 1:dim(errpts)[1]) d[i] = sqrt(sum((errpts[i,]-errpts[1,])^2))
for(i in 1:dim(errpts)[1]) sigmad[i] = error.sigma + d[i]/dmax*(max.error.sigma-error.sigma)
for(i in 1:dim(errpts)[1]) {
  adderr = addnormal(errpts[i,],mask,sigma=sigmad[i]) # This is slow!
  adderr = adderr/max(adderr)  # scale so all have same height
  errtrend = errtrend + adderr
}
# Plotting
covariates(mask)$errtrend = errtrend
spdf.errtrend = SpatialPixelsDataFrame(as.matrix(mask),data=data.frame(D=covariates(mask)$errtrend))
inspdf.errtrend = spdf.errtrend[trapboundary,]
plot(trapboundary)
plot(inspdf.errtrend,col=terrain.colors(40),what="image",add=TRUE)
#contour(inspdf.errtrend,add=TRUE)
points(errpts,pch="+")


# Plot various kinds of densities
# ===============================

#pdf(file="output/sigmas.pdf",h=2,w=6)
#par(mfrow=c(1,3),mar=c(1,4,1,1))
#plot(c(x0,xmax),rep(error.sigma,2),type="l",xaxt="n",ylim=c(0,max.error.sigma),ylab=expression(sigma))
#segments(c(x0+buffer,xmax-buffer),rep(0,2),c(x0+buffer,xmax-buffer),rep(max.error.sigma,2),lty=2)
#par(mar=c(1,1,1,1))
#plot(c(x0,xmax),rep(max.error.sigma,2),type="l",xaxt="n",ylim=c(0,max.error.sigma))
#segments(c(x0+buffer,xmax-buffer),rep(0,2),c(x0+buffer,xmax-buffer),rep(max.error.sigma,2),lty=2)
#plot(x,sigma,type="l",xaxt="n",ylim=c(0,max.error.sigma))
#segments(c(x0+buffer,xmax-buffer),rep(0,2),c(x0+buffer,xmax-buffer),rep(max.error.sigma,2),lty=2)




pdf(file="output/densities.pdf",h=2,w=6)
par(mfrow=c(1,3),mar=c(0.25,0.25,1.5,0.25))
# Plot the ac density surface
plot(trapboundary,main="(a)")
plot(inspdf.D,col=heat.colors(40),what="image",add=TRUE)
# Plot the realised acs
par(mar=c(0.25,0.25,1.5,0.25))
plot(trapboundary,main="(b)")
plot(inpop,pch=19,cex=0.25,add=TRUE)
# Plot the usage density surface
par(mar=c(0.25,0.25,1.5,0.25))
plot(trapboundary,main="(c)")
plot(inspdf.usage,col=heat.colors(40),what="image",add=TRUE)
plot(inpop,pch=19,cex=0.25,col="gray",add=TRUE)
dev.off()

# Plot sigmas and activity centre density maps with observation error
x0 = bbox(boundary)[1,1]
xmax = bbox(boundary)[1,2]
xmid = errpts[1,1]
x = c(x0,(xmid-dmax),xmid,(xmid+dmax),xmax)
sigma = c(max.error.sigma,max.error.sigma,error.sigma,max.error.sigma,max.error.sigma)

pdf(file="output/acesterr.pdf",h=2,w=6)
# Plot with small error
par(mfrow=c(1,3),mar=c(0.25,0.25,1.5,0.25))
plot(trapboundary,main="(a)")
plot(inspdf.acerr,col=heat.colors(40),what="image",add=TRUE)
plot(inpop,pch=19,cex=0.25,col="gray",add=TRUE)
# Plot with large error
par(mar=c(0.25,0.25,1.5,0.25))
plot(trapboundary,main="(b)")
plot(inspdf.maxacerr,col=heat.colors(40),what="image",add=TRUE)
plot(inpop,pch=19,cex=0.25,col="gray",add=TRUE)
# Plot with small-to-large error
par(mar=c(0.25,0.25,1.5,0.25))
plot(trapboundary,main="(c)")
plot(inspdf.acerrtrend,col=heat.colors(40),what="image",add=TRUE)
plot(inpop,pch=19,cex=0.25,col="gray",add=TRUE)
dev.off()

pdf(file="output/sigmas.pdf",h=1,w=6)
par(mfrow=c(1,3),mar=c(1,2,1.5,2))
#plot(c(x0,xmax),rep(error.sigma,2),type="l",xaxt="n",ylim=c(0,max.error.sigma),xlab="",ylab=expression(sigma),main="(d)")
plot(c(x0,xmax),rep(error.sigma,2),type="l",xaxt="n",ylim=c(0,max.error.sigma),xlab="",ylab="",main="(d)")
segments(c(x0+buffer,xmax-buffer),rep(0,2),c(x0+buffer,xmax-buffer),rep(max.error.sigma,2),lty=2)
#plot(c(x0,xmax),rep(max.error.sigma,2),type="l",xaxt="n",xlab="",ylab=expression(sigma),ylim=c(0,max.error.sigma),main="(e)")
plot(c(x0,xmax),rep(max.error.sigma,2),type="l",xaxt="n",xlab="",ylab="",ylim=c(0,max.error.sigma),main="(e)")
segments(c(x0+buffer,xmax-buffer),rep(0,2),c(x0+buffer,xmax-buffer),rep(max.error.sigma,2),lty=2)
#plot(x,sigma,type="l",xaxt="n",xlab="",ylab=expression(sigma),ylim=c(0,max.error.sigma),main="(f)")
plot(x,sigma,type="l",xaxt="n",xlab="",ylab="",ylim=c(0,max.error.sigma),main="(f)")
segments(c(x0+buffer,xmax-buffer),rep(0,2),c(x0+buffer,xmax-buffer),rep(max.error.sigma,2),lty=2)
dev.off()

# Now redo with labels that allow all three plots to be put in one, and as png instead of pdf:

jpeg(file="output/densities2.jpg",h=2,w=6,units="in",res=720)
par(mfrow=c(1,3),mar=c(0.25,0.25,1.5,0.25))
# Plot the ac density surface
plot(trapboundary,main="(a)")
plot(inspdf.D,col=heat.colors(40),what="image",add=TRUE)
# Plot the realised acs
par(mar=c(0.25,0.25,1.5,0.25))
plot(trapboundary,main="(b)")
plot(inpop,pch=19,cex=0.25,add=TRUE)
# Plot the usage density surface
par(mar=c(0.25,0.25,1.5,0.25))
plot(trapboundary,main="(c)")
plot(inspdf.usage,col=heat.colors(40),what="image",add=TRUE)
plot(inpop,pch=19,cex=0.25,col="gray",add=TRUE)
dev.off()


jpeg(file="output/acesterr2.jpg",h=2,w=6,units="in",res=720)
# Plot with small error
par(mfrow=c(1,3),mar=c(0.25,0.25,1.5,0.25))
plot(trapboundary,main="(d)")
plot(inspdf.acerr,col=heat.colors(40),what="image",add=TRUE)
plot(inpop,pch=19,cex=0.25,col="gray",add=TRUE)
# Plot with large error
par(mar=c(0.25,0.25,1.5,0.25))
plot(trapboundary,main="(e)")
plot(inspdf.maxacerr,col=heat.colors(40),what="image",add=TRUE)
plot(inpop,pch=19,cex=0.25,col="gray",add=TRUE)
# Plot with small-to-large error
par(mar=c(0.25,0.25,1.5,0.25))
plot(trapboundary,main="(f)")
plot(inspdf.acerrtrend,col=heat.colors(40),what="image",add=TRUE)
plot(inpop,pch=19,cex=0.25,col="gray",add=TRUE)
dev.off()


jpeg(file="output/sigmas2.jpg",h=1,w=6,units="in",res=720)
par(mfrow=c(1,3),mar=c(1,2,1.5,2))
#plot(c(x0,xmax),rep(error.sigma,2),type="l",xaxt="n",ylim=c(0,max.error.sigma),xlab="",ylab=expression(sigma),main="(d)")
plot(c(x0,xmax),rep(error.sigma,2),type="l",xaxt="n",ylim=c(0,max.error.sigma),xlab="",ylab="",main="(g)")
segments(c(x0+buffer,xmax-buffer),rep(0,2),c(x0+buffer,xmax-buffer),rep(max.error.sigma,2),lty=2)
#plot(c(x0,xmax),rep(max.error.sigma,2),type="l",xaxt="n",xlab="",ylab=expression(sigma),ylim=c(0,max.error.sigma),main="(e)")
plot(c(x0,xmax),rep(max.error.sigma,2),type="l",xaxt="n",xlab="",ylab="",ylim=c(0,max.error.sigma),main="(h)")
segments(c(x0+buffer,xmax-buffer),rep(0,2),c(x0+buffer,xmax-buffer),rep(max.error.sigma,2),lty=2)
#plot(x,sigma,type="l",xaxt="n",xlab="",ylab=expression(sigma),ylim=c(0,max.error.sigma),main="(f)")
plot(x,sigma,type="l",xaxt="n",xlab="",ylab="",ylim=c(0,max.error.sigma),main="(i)")
segments(c(x0+buffer,xmax-buffer),rep(0,2),c(x0+buffer,xmax-buffer),rep(max.error.sigma,2),lty=2)
dev.off()


pdf(file="output/acuseesterr.pdf",h=2,w=6)
par(mfrow=c(1,3),mar=c(0.25,0.25,1.5,0.25))
# Plot with small error
plot(trapboundary,main="(a)")
plot(inspdf.acuseerr,col=heat.colors(40),what="image",add=TRUE)
plot(inpop,pch=19,cex=0.25,col="gray",add=TRUE)
# Plot with large error
par(mar=c(0.25,0.25,1.5,0.25))
plot(trapboundary,main="(b)")
plot(inspdf.maxacuseerr,col=heat.colors(40),what="image",add=TRUE)
plot(inpop,pch=19,cex=0.25,col="gray",add=TRUE)
# Plot with small-to-large error
par(mar=c(0.25,0.25,1.5,0.25))
plot(trapboundary,main="(c)")
plot(inspdf.acuseerrtrend,col=heat.colors(40),what="image",add=TRUE)
plot(inpop,pch=19,cex=0.25,col="gray",add=TRUE)
dev.off()


# Do an SCR survey to illustrate prediction error size change
# ===========================================================
scrtraps = make.grid(nx=4,ny=4,spacex=10,detector="count",origin=c(35,35))
plot(boundary)
plot(scrtraps,add=TRUE)
plot(trapboundary,add=TRUE)
set.seed(2)
simch = sim.capthist(scrtraps,pop,detectpar=list(g0=1,sigma=10),noccasions=1)
set.seed(2)
simch_notrenum = sim.capthist(scrtraps,pop,detectpar=list(g0=1,sigma=10),noccasions=1, renumber = FALSE)
set.seed(93)
simch2 = sim.capthist(scrtraps,pop,detectpar=list(g0=1,sigma=10),noccasions=1)
set.seed(93)
simch2_notrenum = sim.capthist(scrtraps,pop,detectpar=list(g0=1,sigma=10),noccasions=1, renumber = FALSE)

founddet <- FALSE
thisseed <- 1
while(!founddet){
  thisseed <- thisseed + 1
  set.seed(thisseed)
  thisch = sim.capthist(scrtraps,pop,detectpar=list(g0=1,sigma=10),noccasions=1, renumber = FALSE)
  if(c('95') %in% row.names(thisch)){
    founddet = setequal(c(13), which(thisch["95",1,] > 0)) 
    }
}
# 32,93 for det13
# 5 for det14
# 2 for det15
# 20 for det16


xx2 <- simch2_notrenum[,1,]
row.names(simch_notrenum)
row.names(simch2_notrenum)
which(row.names(simch_notrenum) == "95")
which(row.names(simch2_notrenum) == "95")
which(row.names(simch_notrenum) == "61")
which(row.names(simch2_notrenum) == "61")

summary(simch)
plot(simch,tracks=TRUE,border=0)
simfit = secr.fit(simch,mask=mask)
fxtot = fx.total(simfit,mask=mask)
simfit2 = secr.fit(simch2,mask=mask)
fxtot2 = fx.total(simfit2,mask=mask)
#plotcovariate(fxtot,covariate="D.sum",what="image")
#plot(scrtraps,add=TRUE)

# prob of ch 1 / ch 2
g0 <- predict(simfit)[2,2]
sigma <- predict(simfit)[3,2]
d2 <- sqrt(sum((pop[95,] - scrtraps[13,])^2))
d1 <- sqrt(sum((pop[95,] - scrtraps[15,])^2))
p2 <- g0 * exp(-d2^2 / (2*sigma^2))
p1 <- g0 * exp(-d1^2 / (2*sigma^2))
p2
p1
p1/p2

pdf(file="output/screrr.pdf",h=4,w=8)
par(mar=c(1,1,1,1), mfrow = c(1,2))

#dev.off()
#pdf(file="output/screrr.pdf",h=4,w=4)
plot(trapboundary)
fxi.contour(simfit2,i=c(20,5),nx=200,add=TRUE,drawlabels=FALSE)
sf1 <- fxi.contour(simfit2,i=c(20,5),nx=200,add=TRUE,drawlabels=FALSE, plt = FALSE, output = "sf")
plot(scrtraps,add=TRUE)
ch1 = simfit2$capthist[20,1,]
detected1 = which(ch1>0)
points(scrtraps$x[detected1],scrtraps$y[detected1],pch=15)
ch2 = simfit2$capthist[5,1,]
detected2 = which(ch2>0)
points(scrtraps$x[detected2],scrtraps$y[detected2],pch=17)
text(pop[95,], "A")
text(pop[61,], "B")

plot(trapboundary)
fxi.contour(simfit,i=c(20,6),nx=200,add=TRUE,drawlabels=FALSE)
plot(scrtraps,add=TRUE)
ch1 = simfit$capthist[20,1,]
detected1 = which(ch1>0)
points(scrtraps$x[detected1],scrtraps$y[detected1],pch=15)
ch2 = simfit$capthist[6,1,]
detected2 = which(ch2>0)
points(scrtraps$x[detected2],scrtraps$y[detected2],pch=17)
text(pop[95,], "A")
text(pop[61,], "B")
dev.off()
