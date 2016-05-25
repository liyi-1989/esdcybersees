library(ncdf4)
library(MASS)

# Load the GCM data
get_gcm_train=function(mylon,mylat){
  cat("In position: ",mylon,mylat,"\n")
  files=list.files("./gcm_train")
  pr=NULL
  for(i in files){
    cat("Read in file ",i,"\n")
    ncin = nc_open(paste0("./gcm_train/",i))
    lon = ncvar_get(ncin,"lon")
    lat = ncvar_get(ncin,"lat")
    Pr=ncvar_get(ncin,"pr")
    pr=c(pr,Pr[((mylon-0.1)<lon)&(lon<(mylon+0.1)),((mylat-0.1)<lat)&(lat<(mylat+0.1)),])
  }
  prm=matrix(pr,365,length(pr)/365)
  return(prm)
}

get_gcm_test=function(mylon,mylat){
  cat("In position: ",mylon,mylat,"\n")
  files=list.files("./gcm_test")
  pr=NULL
  for(i in files){
    cat("Read in file ",i,"\n")
    ncin = nc_open(paste0("./gcm_test/",i))
    lon = ncvar_get(ncin,"lon")
    lat = ncvar_get(ncin,"lat")
    Pr=ncvar_get(ncin,"pr")
    pr=c(pr,Pr[((mylon-0.1)<lon)&(lon<(mylon+0.1)),((mylat-0.1)<lat)&(lat<(mylat+0.1)),])
  }
  prm=matrix(pr,365,length(pr)/365)
  return(prm)
}
# Load the observation data, unusful, use matlab do this.
get_obs=function(mylon,mylat){
  cat("In position: [",mylon,"]*[",mylat,"]\n")
  pr=NULL
  for(i in 1979:2003){
    cat("Read in file ",i,"\n")
    ncin = nc_open(paste0("/media/liyi/新加卷/downscaling/pr_",i,".nc"))
    lon = ncvar_get(ncin,"lon")
    lat = ncvar_get(ncin,"lat")
    Pr=ncvar_get(ncin,"precipitation_amount") #[day,lon,lat] 
    tmp=Pr[,((mylon[1]-0.01)<lon)&(lon<(mylon[2]+0.01)),((mylat[1]-0.01)<lat)&(lat<(mylat[2]+0.01))]
    tmp=apply(tmp,1,mean)
    if(length(tmp)==366){
      tmp[31+28]=mean(tmp[(31+28):(31+29)])
      tmp=tmp[-(31+29)]
    }
    pr=c(pr,tmp)
    rm(Pr)
  }
  prm=matrix(pr,365,length(pr)/365)
  return(prm)
}


# Y=get_obs(mylon,mylat)

## Density estimation

# plot(rank(Xim)/25,rank(Yim)/25)
# plot(density(Xim))
# x0=2
# density(Xim, from=x0, to=x0, n=1)$x
# density(Xim, from=x0, to=x0, n=1)$y
# 
# f2=kde2d(rank(Xim)/25,rank(Yim)/25)
# image(f2)
# persp(f2, phi = 30, theta = 20, d = 5)
# contour(f2)
# 
# kde2d(rank(Xim)/25,rank(Yim)/25,lims=c(0.5,0.5,0.3,0.3),n=1)



kde1d_point=function(x0,x_given){
  density(x_given, from=x0, to=x0, n=1)$y
}

kde2d_point=function(x0,y0,x_given,y_given){
  kde2d(x_given,y_given,lims=c(x0,x0,y0,y0),n=1)$z
}

max_cond=function(x_fix,x_given,y_given){
  ygrid=seq(from=min(y_given),to=max(y_given),length=100)
  yt=ygrid
  u=rank(x_given)/25
  v=rank(y_given)/25
  for(i in 1:length(ygrid)){
    yt[i]=kde2d_point(rank(c(x_fix,x_given))[1]/26,rank(c(ygrid[i],y_given))[1]/26,u,v)*kde1d_point(x_fix,x_given)
  }
  #plot(ygrid,yt,type="o")
  #ygrid[yt==max(yt)]
  ygrid[which.max(yt)]
}
