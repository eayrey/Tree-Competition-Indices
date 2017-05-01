#dat input into each of these functions must have SP, DBH, X and Y
#If there are multiple plots, run the functions iteratively in a for loop

#Simpson Species Index
SimpsonS=function(dat){
  frequencies=table(dat$SPP)
  SIM_1=frequencies*(frequencies-1)
  SIM=1 - (sum(SIM_1)/ (nrow(Book1)*(nrow(Book1)-1)))
  SIM
}

#Simpson Diameter Index with Xcm diam classes
x=3
SimpsonD=function(dat){
  D_Classes=round(dat$DBH/x)*x
  frequencies=table(D_Classes)
  SIM_1=frequencies*(frequencies-1)
  SIM=1 - (sum(SIM_1)/ (nrow(dat)*(nrow(dat)-1)))
  SIM
}

#GINI based on diameters
GINI=function(dat){
  n=nrow(dat)
  mu=mean(dat$DBH)
  N=n * (n - 1)
  ox <- dat$DBH[order(dat$DBH)]
  dsum <- drop(crossprod(2 * 1:n - n - 1,  ox))
  GINI = dsum / (mu * N)
  GINI
}

Distances=as.matrix(dist(dat$X, dat$Y, method='e', upper=TRUE))

#TD Diameter Differentation Index, NO SPATIAL COMPONENT
TD=function(dat){
  #divide each tree's diameter by all the others in its plot
  biggers=lapply(dat$DBH,function(x){x/dat$DBH})
  biggers=do.call(rbind,biggers)
  #remove self comparisons
  diag(biggers)=NA
  #We only want the smaller trees divided by bigger trees
  biggers[biggers>1]=NA
  #1 minus the ratios for some reason
  TDI=1-biggers
  #summed then divided by total tree number
  TD=sum(TDI,na.rm=TRUE)/nrow(dat)
  TD
}

#Mingling Index (MI), using the 4 nearest neighbors of each tree
MI=function(dat){
  diag(Distances)=NA
  neighbors=data.frame()
  if (nrow(dat)>4){
    for (r in 1:nrow(Distances)){
      one_tree=Distances[r,]
      neighbors=rbind(neighbors,order(one_tree)[1:4])
    }
  }else{
    x=nrow(dat-1)
    for (r in 1:nrow(Distances)){
      one_tree=Distances[r,]
      neighbors=rbind(neighbors,order(one_tree)[1:x])
    }
  }
  neighborsI=data.frame()
  for (t in 1:nrow(neighbors)){
    the_sp=dat$SPP[t]
    same_sp=t(apply(neighbors[t,],1,function(x){dat$SPP[x]==the_sp}))
    neighbor_diversity=(1*!same_sp)/length(same_sp)
    neighborsI=rbind(neighborsI,sum(neighbor_diversity))
  }
  M=mean(neighborsI[,1])
  M
}

#Clark-Evans Index
area=400
CE=function(dat){
  mean_nn=mean(nndist(dat$X,dat$Y))
  R_EXP=1/(2*sqrt(nrow(dat)/area))
  CE = mean_nn / R_EXP
  CE
}


