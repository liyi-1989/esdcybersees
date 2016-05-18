setwd("~/Documents/esd")
library(ncdf4)
library(MASS)
# 1. Read in the GCM data
MON=1:365
names(MON)=c(rep("JAN",31),rep("FEB",28),rep("MAR",31),rep("APR",30),rep("MAY",31),rep("JUN",30),rep("JUL",31),rep("AUG",31),rep("SEP",30),rep("OCT",31),rep("NOV",30),rep("DEC",31))

get_gcm=function(mylon,mylat){
  cat("In position: ",mylon,mylat,"\n")
  files=list.files("./gcm")
  pr=NULL
  for(i in files){
    cat("Read in file ",i,"\n")
    ncin = nc_open(paste0("./gcm/",i))
    lon = ncvar_get(ncin,"lon")
    lat = ncvar_get(ncin,"lat")
    Pr=ncvar_get(ncin,"pr")
    pr=c(pr,Pr[((mylon-0.1)<lon)&(lon<(mylon+0.1)),((mylat-0.1)<lat)&(lat<(mylat+0.1)),])
  }
  prm=matrix(pr,365,length(pr)/365)
  return(prm)
}

X1=get_gcm(-113.75+180,43)
X2=get_gcm(-113.75+180,45)
X3=get_gcm(-116.25+180,45)
X4=get_gcm(-116.25+180,43)
X=(X1+X2+X3+X4)/4 # kg/m^2/s = mm/hr *3600
X=X4
#rm(X1,X2,X3,X4)
colnames(X)=1975:2005
rownames(X)=1:365
X=X[,5:29] * 3600 * 24 # mm # 365 days * 25 years
#http://weather.unisys.com/wxp/wxp5/Users/units.lup
Xi=X[MON[names(MON)=="JAN"],] # 31 days * 25 years
Xim=apply(Xi,2,sum) # 1 month * 25 years 

plot(ecdf(Xi))
plot(ecdf(Xim))
plot(1979:2003,Xim,type="o")

# 2. Read in the observation data
mylon=c(-116.25,-113.75)
mylat=c(43,45)

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

#Y=get_obs(mylon,mylat)
Y=read.csv("obs_avg.csv",header = F)
colnames(Y)=1979:2003
rownames(Y)=1:365
Yi=Y[MON[names(MON)=="JAN"],] # 31 days * 25 years
Yim=apply(Yi,2,sum) # 1 month * 25 years 

plot(ecdf(Yim))
plot(1979:2003,Yim,type="o")

# 3. Quantile Mapping
Z=quantile(Yim,ecdf(Xim)(Xim))
plot(ecdf(Z))
plot(1979:2003,Z,type="o")

# Quantile Mapping
plot(sort(Xim),rank(sort(Xim))/length(Xim),type="o",xlim=c(0,max(Yim)))
lines(sort(Yim),rank(sort(Yim))/length(Yim),type="o",col="red")
# Bias Correction Result
plot(1979:2003,Xim,type="o",ylim=c(0,max(Yim)))
lines(1979:2003,Yim,type="o",col="red")
lines(1979:2003,Z,type="o",col="blue")


plot(rank(Xim)/25,rank(Yim)/25)
plot(density(Xim))

x0=2
density(Xim, from=x0, to=x0, n=1)$x
density(Xim, from=x0, to=x0, n=1)$y

f2=kde2d(rank(Xim)/25,rank(Yim)/25)
image(f2)
persp(f2, phi = 30, theta = 20, d = 5)
contour(f2)

kde2d(rank(Xim)/25,rank(Yim)/25,lims=c(0.5,0.5,0.3,0.3),n=1)

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

Xim2=Xim
for(i in 1:length(Xim)){
  Xim2[i]=max_cond(Xim[i],Xim,Yim)
}

plot(1979:2003,Xim,type="o",ylim=c(0,max(Yim)))
lines(1979:2003,Yim,type="o",col="red")
lines(1979:2003,Z,type="o",col="blue")
lines(1979:2003,Xim2,type="o",col="green")

sqrt(mean((Yim-Z)^2))
sqrt(mean((Yim-Xim2)^2))





