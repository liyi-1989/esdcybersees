#setwd("~/Documents/esdcybersees")
source('functions.R')
MON=1:365
names(MON)=c(rep("JAN",31),rep("FEB",28),rep("MAR",31),rep("APR",30),rep("MAY",31),rep("JUN",30),
             rep("JUL",31),rep("AUG",31),rep("SEP",30),rep("OCT",31),rep("NOV",30),rep("DEC",31))
Mon=data.frame(mon=as.character(unique(names(MON))),num=1:12,days=c(31,28,31,30,31,30,31,31,30,31,30,31))
################## 1. Read in the GCM data(train) ##################
mylon=c(-116.25,-113.75)
mylat=c(43,45)
mon=1 # Month number: Jan=1
# 1.1 training
# http://nomads.gfdl.noaa.gov:8080/DataPortal/cmip5.jsp RCP6 & historical 
# X1=get_gcm(-113.75+180,43); X2=get_gcm(-113.75+180,45); X3=get_gcm(-116.25+180,45); 
X=get_gcm_train(mylon[1]+180,mylat[1]) # X=(X1+X2+X3+X4)/4 # kg/m^2/s = mm/hr *3600
colnames(X)=1975:2005; rownames(X)=1:365 # rm(X1,X2,X3,X4)
X=X[,5:29] * 3600 * 24 # mm # 365 days * 25 years # http://weather.unisys.com/wxp/wxp5/Users/units.lup

# 1.2 testing
X_test=get_gcm_test(mylon[1]+180,mylat[1])
colnames(X_test)=2006:2015; rownames(X_test)=1:365
X_test=X_test * 3600 * 24
################## 2. Read in the observation data(train) ##################
ncin = nc_open(paste0("./obs_in_grid/","obs_in_grid_",mylon[1],"-",mylat[1],".nc"))
lon = ncvar_get(ncin,"longitude")
lat = ncvar_get(ncin,"latitude")
t1=ncvar_get(ncin,"time")
Pr=ncvar_get(ncin,"precipitation") # 60*48*[365*37]  
# 2.1 testing
t1_test=t1[-(1:(365*27))]
Pr_test=Pr[,,-(1:(365*27))]
# 2.2 training
t1=t1[1:(365*25)]
Pr=Pr[,,1:(365*25)]

############# Month-specific Data Manipulation ############# 
# GCM
Xi=X[MON[names(MON)==Mon[mon,1]],] # 31 days * 25 years
Xim=apply(Xi,2,sum) # 1 month * 25 years 
# OBS
Pr1=Pr[,,is.element(t1,MON[names(MON)==Mon[mon,1]])]
Yi=matrix(apply(Pr1,3,mean),Mon[mon,3],25); colnames(Yi)=1979:2003
Yim=apply(Yi,2,sum) 


Xi_test=X_test[MON[names(MON)==Mon[mon,1]],]
Xim_test=apply(Xi_test,2,sum)

Pr1_test=Pr_test[,,is.element(t1_test,MON[names(MON)==Mon[mon,1]])]
Yi_test=matrix(apply(Pr1_test,3,mean),Mon[mon,3],10); colnames(Yi_test)=2006:2015
Yim_test=apply(Yi_test,2,sum) 


# Y=read.csv("obs_avg.csv",header = F); Y=t(Y); colnames(Y)=1979:2003; rownames(Y)=1:365
# Yi=Y[MON[names(MON)=="JAN"],] # 31 days * 25 years
# Yim=apply(Yi,2,sum) # 1 month * 25 years 
################## 3-1 Quantile Mapping ##################
# training
Qi=quantile(Yim,ecdf(Xim)(Xim))
# testing
# Xim_test=Xim
Qi_test=quantile(Yim,ecdf(Xim)(Xim_test))
################## 3-2 Copula ##################
Ci=Xim
for(i in 1:length(Xim)){
  Ci[i]=max_cond(Xim[i],Xim,Yim)
}
# testing
Ci_test=Xim_test
for(i in 1:length(Xim_test)){
  Ci_test[i]=max_cond(Xim_test[i],Xim,Yim)
}
################## 4. Spatial disaggregation ##################
rmse1=rmse2=rmse3=rep(0,10)

for(k in 1:10){
  T1=apply(Pr1_test[,,(1:Mon[mon,3])+(k-1)*Mon[mon,3]],c(1,2),sum) # Test truth (not available in real case) (for a given year)
  
  M=apply(Pr1,c(1,2),sum)/25 # Local Historical Monthly (accumulation) mean
  D1=Ci_test[k]*M/mean(Ci) #Copula; Downscale 1979 Jan, actually apply to training data (which need to try in testing data)
  D2=Qi_test[k]*M/mean(Qi) #QM Downscale 1979 Jan, actually apply to training data (which need to try in testing data)
  D3=Xim_test[k]*M/mean(Xim) #GCM
  
  rmse1[k]=sqrt(mean((D1-T1)^2)) 
  rmse2[k]=sqrt(mean((D2-T1)^2)) 
  rmse3[k]=sqrt(mean((D3-T1)^2)) 
}

plot(2006:2015,rmse1,col="black",type="o",ylim = c(0,max(c(rmse1,rmse2,rmse3))),ylab="RMSE",
     main=paste("Downscaling on Testing Set in",Mon[mon,1]))
lines(2006:2015,rmse2,col="blue",type="o")
lines(2006:2015,rmse3,col="red",type="o")
legend("top",c("Nonparametric Copula","Quantile Mapping","GCM"),lty=1,pch=1,col=c("black","blue","red"),cex=0.5)
mean(rmse1<rmse2)



################## 5. Make Plots ##################
# 5.1 Quantile Mapping
# plot(ecdf(Z)); plot(ecdf(Xim)); plot(ecdf(Yim))
plot(sort(Xim),rank(sort(Xim))/length(Xim),type="o",xlim=c(0,max(Yim)),xlab="Preicipitation",ylab="ECDF",
     main="Quantile Mapping (Monthly accumulation, 1979 - 2003)",col="blue")
lines(sort(Yim),rank(sort(Yim))/length(Yim),type="o",col="red")
legend("bottomright",c("GCM","Obs"),lty=1,pch=1,col=c("blue","red"))
# 5.2 Bias Correction
plot(1979:2003,Xim,type="o",ylim=c(0,max(Yim)),col="blue",xlab="Years",ylab="Precipitation(mm)",main="Bias Correction")
lines(1979:2003,Yim,type="o",col="red")
lines(1979:2003,Qi,type="o",col="green")
lines(1979:2003,Ci,type="o",col="black")
legend("top",c("GCM","Obs","Quantile Mapping","Nonparametric Copula"),lty=1,pch=1,col=c("blue","red","green","black"))
################## 6. Calculate Error ##################
# 6.1 BS step
# sqrt(mean((Yim-Xim)^2))
# sqrt(mean((Yim-Qi)^2))
# sqrt(mean((Yim-Ci)^2))

sqrt(mean((Yim_test-Xim_test)^2))
sqrt(mean((Yim_test-Qi_test)^2))
sqrt(mean((Yim_test-Ci_test)^2))

# 6.2 SD step
sqrt(mean((D1-T1)^2)) #[1] 93.93788
sqrt(mean((D2-T1)^2)) #[1] 50.06246
sqrt(mean((D3-T1)^2)) #[1] 60.90534


mean(abs(D1-T1)) #[1] 63.46105
mean(abs(D2-T1)) #[1] 33.90981
mean(abs(D3-T1)) #[1] 41.21309
