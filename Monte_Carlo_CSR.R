#Monte Carlo Simulations for testing Complete Spatial Randomness

#runs monte carlo simulations of randomly placed trees, 
#Calculates Ripley's k function for each simulation.
#Develops an envelope of maximum and minimum kest at different distances
#if K of true tree distribution is less than monte carlo extreme, they're evenly spaced, csr=1
#if K is more than monte carlo extremes, then they're clustered, csr=2
#if K is at some distances above and below the extremes, they're likely even at some spatial scale, and clustered at another csr=3
#if K never left the monte carlo envelope, they're random, csr=0

library(spatstat)
library(maptools)
library(raster)

#area of extent (often in square m)
area=400
#or a SpatialPolygon of the area
#area=readShapeSpatial("mytrees.shp)

#Planar Point Pattern of your trees
mydata=readShapeSpatial("myarea.shp")
trees=as.ppp(mydata)

#Number of trees (see below for multiple runs)
ntree=nrow(trees)

enveloper=function(area, ntree){  
  if (is.numeric(area)){
    sample_area=extent(0,sqrt(area),0,sqrt(area))
    sample_area=as(sample_area, 'SpatialPolygons')
  }else{
    sample_area=area
  }
  n=200
  Km = matrix(nrow=n, ncol=513)
  for(i in 1:n){
    Pmc = runifpoint(ntree,win=sample_area)    # Generate random points
    Kmc = Kest(Pmc, nsim=n, border="border")   # Compute K
    Km[i,] = Kmc$iso             
  }
  list(maxes=apply((Km),2,max), mins=apply((Km),2,min))
}
envelopes=enveloper(area,ntree)

kest=Kest(trees, correction="border")
kest$border[is.nan(kest$border)]<-0

plot(envelopes$maxes, pch=20)
points(envelopes$mins, pch=20)
points(kest$border, col='green')

if (sum((kest$border-env$mins)<0) > 0 && sum((env$maxes-kest$border)<0) < 1){
  #evenly spaced
  csr=1
}
#some above the envelope max, none below
if (sum((kest$border-env$mins)<0) < 1 && sum((env$maxes-kest$border)<0) > 0){
  #clustered
  csr=2
}
#some above and below the envelope
if (sum((kest$border-env$mins)<0) > 0 && sum((env$maxes-kest$border)<0) > 0){
  #clustered and spaced at different scales
  csr=3
}
#if observed were all above lower envelope, and all less than the upper envelope
if (sum((kest$border-env$mins)<0) == 0 && sum((env$maxes-kest$border)<0) == 0){
  #randomly distributed
  csr=0
}

print(csr)
###############################################################################################


#With multiple runs, and different numbers of trees per plot
#Generates a bunch of envelopes at once for different numbers of trees

library(doParallel)
#Number of processor cores to use for parallel processesing
#Default is the number of cores you have in your machine minus 1.
n= detectCores(all.tests = FALSE, logical = TRUE)-1

area=400

#run 200 monte carlo simulations of tree locations within an area. Calculate 
enveloper=function(area){
  cl=makeCluster(n)
  registerDoParallel(cl)
  envelopes=foreach(j=2:50, .combine='rbind', .inorder=TRUE) %dopar% {
    library("spatstat")
    library("maptools")
    library("raster")
    sample_area=extent(0,sqrt(area),0,sqrt(area))
    sample_area=as(sample_area, 'SpatialPolygons')
    n=100
    Km = matrix(nrow=n, ncol=513)
    for(i in 1:n){
      Pmc = runifpoint(j,win=sample_area)    # Generate random points
      Kmc = Kest(Pmc, nsim=n, border="border")   # Compute K
      Km[i,] = Kmc$iso             
    }
    list(maxes=apply((Km),2,max), mins=apply((Km),2,min))
  }
  stopCluster(cl)
  return(envelopes)
}
envelopes=enveloper(400)

kest=Kest(trees, correction="border")
kest$border[is.nan(kest$border)]<-0

env=envelopes[trees$n-1,]
#some below envelope min, none above
if (sum((kest$border-env$mins)<0) > 0 && sum((env$maxes-kest$border)<0) < 1){
  #evenly spaced
  csr=1
}
#some above the envelope max, none below
if (sum((kest$border-env$mins)<0) < 1 && sum((env$maxes-kest$border)<0) > 0){
  #clustered
  csr=2
}
#some above and below the envelope
if (sum((kest$border-env$mins)<0) > 0 && sum((env$maxes-kest$border)<0) > 0){
  #clustered and spaced at different scales
  csr=3
}
#if observed were all above lower envelope, and all less than the upper envelope
if (sum((kest$border-env$mins)<0) == 0 && sum((env$maxes-kest$border)<0) == 0){
  #randomly distributed
  csr=0
}

