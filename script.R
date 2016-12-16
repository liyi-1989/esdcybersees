#setwd("~/Documents/esdcybersees")
source('functions.R')
MON=1:365
names(MON)=c(rep("JAN",31),rep("FEB",28),rep("MAR",31),rep("APR",30),rep("MAY",31),rep("JUN",30),
             rep("JUL",31),rep("AUG",31),rep("SEP",30),rep("OCT",31),rep("NOV",30),rep("DEC",31))
Mon=data.frame(mon=as.character(unique(names(MON))),num=1:12,days=c(31,28,31,30,31,30,31,31,30,31,30,31))
################## 1. Read in the GCM data(train) ##################
mylon=c(-116.25,-113.75)
#mylon=c(-71.25,-68.75)
mylat=c(43,45)
mylon=c(-73.75,-71.25)
mylat=c(41,43)

mylon=c(-118.75,-116.25)
mylat=c(35,37)

# 1.1 training
# http://nomads.gfdl.noaa.gov:8080/DataPortal/cmip5.jsp RCP6 & historical 
# X1=get_gcm(-113.75+180,43); X2=get_gcm(-113.75+180,45); X3=get_gcm(-116.25+180,45); 
X=get_gcm_train(mylon[2]+360,mylat[2]) # X=(X1+X2+X3+X4)/4 # kg/m^2/s = mm/hr *3600
colnames(X)=1975:2005; rownames(X)=1:365 # rm(X1,X2,X3,X4)
X=X[,5:29] * 3600 * 24 # mm # 365 days * 25 years # http://weather.unisys.com/wxp/wxp5/Users/units.lup

# 1.2 testing
X_test=get_gcm_test(mylon[2]+360,mylat[2])
colnames(X_test)=2006:2015; rownames(X_test)=1:365
X_test=X_test * 3600 * 24
################## 2. Read in the observation data(train) ##################
# mylon=mylon-1.25
# mylat=mylat-1
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

mon=2 # Month number: Jan=1
# GCM
Xi=X[MON[names(MON)==Mon[mon,1]],] # 31 days * 25 years
Xim=apply(Xi,2,sum) # 1 month * 25 years 

Xi_test=X_test[MON[names(MON)==Mon[mon,1]],]
Xim_test=apply(Xi_test,2,sum)

# OBS
Pr1=Pr[,,is.element(t1,MON[names(MON)==Mon[mon,1]])]
Yi=matrix(apply(Pr1,3,mean, na.rm = TRUE),Mon[mon,3],25); colnames(Yi)=1979:2003
Yim=apply(Yi,2,sum) 

Pr1_test=Pr_test[,,is.element(t1_test,MON[names(MON)==Mon[mon,1]])]
Yi_test=matrix(apply(Pr1_test,3,mean, na.rm = TRUE),Mon[mon,3],10); colnames(Yi_test)=2006:2015
Yim_test=apply(Yi_test,2,sum) 


# Y=read.csv("obs_avg.csv",header = F); Y=t(Y); colnames(Y)=1979:2003; rownames(Y)=1:365
# Yi=Y[MON[names(MON)=="JAN"],] # 31 days * 25 years
# Yim=apply(Yi,2,sum) # 1 month * 25 years 

################## 3-0 Quantile Mapping Directly ##################

# B1=B2=array(0,c(dim(Pr1)[1:2],10))
# 
# for(i in 1:dim(Pr1)[1]){
#   
#   for(j in 1:dim(Pr1)[2]){
#     print(paste("lon: ",i,"lat: ",j))
#     ST=matrix(Pr1[i,j,],Mon[mon,3],25); 
#     colnames(ST)=1979:2003
#     ST=apply(ST,2,sum)
#     
#     if(all(is.nan(ST))){
#       B1[i,j,]=B2[i,j,]=NA
#       print("All NA here!")
#       next
#     }
#     
#     ST_test=matrix(Pr1_test[i,j,],Mon[mon,3],10); 
#     colnames(ST_test)=2006:2015
#     ST_test=apply(ST_test,2,sum)
#     
#     
#     B1[i,j,]=quantile(ST,ecdf(Xim)(Xim_test),na.rm = T)
#     B2[i,j,]=Xim_test
#     for(k in 1:length(Xim_test)){
#       B2[i,j,k]=max_cond(Xim_test[k],Xim,ST)
#     }
#     
#   }
# }

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
rmse1=rmse2=rmse3=rmse4=rmse5=rep(0,10)

for(k in 1:10){
  T1=apply(Pr1_test[,,(1:Mon[mon,3])+(k-1)*Mon[mon,3]],c(1,2),sum, na.rm = TRUE) # Test truth (not available in real case) (for a given year)
  
  M=apply(Pr1,c(1,2),sum, na.rm = TRUE)/25 # Local Historical Monthly (accumulation) mean
  D1=Ci_test[k]*M/mean(Ci) #Copula; Downscale 1979 Jan, actually apply to training data (which need to try in testing data)
  D2=Qi_test[k]*M/mean(Qi) #QM Downscale 1979 Jan, actually apply to training data (which need to try in testing data)
  D3=Xim_test[k]*M/mean(Xim) #GCM
  
  
  rmse1[k]=sqrt(mean((D1-T1)^2)) 
  rmse2[k]=sqrt(mean((D2-T1)^2)) 
  rmse3[k]=sqrt(mean((D3-T1)^2)) 
  # rmse4[k]=sqrt(mean((B1[,,k]-T1)^2,na.rm = T)) 
  # rmse5[k]=sqrt(mean((B2[,,k]-T1)^2,na.rm=T))
}

plot(2006:2015,rmse1,col="black",type="o",ylim = c(0,max(c(rmse1,rmse2,rmse3))),xlab="",ylab="RMSE",
     main=paste("Downscaling on Testing Set in",Mon[mon,1]))
lines(2006:2015,rmse2,col="blue",type="o")
lines(2006:2015,rmse3,col="green",type="o")
# lines(2006:2015,rmse4,col="green",type="o")
# lines(2006:2015,rmse5,col="yellow",type="o")
legend("topright",c("NCBCSD","BCSD","GCM"),lty=1,pch=1,col=c("black","blue","green"),cex=1)

arrows(2011,5,2012,10, col='black', length=0.1, lwd=3)
text(2011,2,label="NCBCSD",col="black",cex=0.75)
arrows(2009,80,2008,70, col='blue', length=0.1, lwd=3)
text(2009,85,label="BCSD",col="blue",cex=0.75)


mean(rmse1<rmse2)

sqrt(mean((Yim_test-Xim_test)^2))
sqrt(mean((Yim_test-Qi_test)^2))
sqrt(mean((Yim_test-Ci_test)^2))







################## 5. Make Plots ##################
# 5.1 Quantile Mapping
# plot(ecdf(Z)); plot(ecdf(Xim)); plot(ecdf(Yim))
plot(sort(Xim),rank(sort(Xim))/length(Xim),type="o",xlim=c(0,max(Yim)),xlab="Preicipitation",ylab="ECDF",
     main="Quantile Mapping (Monthly accumulation, 1979 - 2003)",col="blue")
lines(sort(Yim),rank(sort(Yim))/length(Yim),type="o",col="red")
legend("bottomright",c("GCM","Obs"),lty=1,pch=1,col=c("blue","red"))

plot(sort(Xim),rank(sort(Xim))/length(Xim),type="o",xlim=c(0,max(Yim)),xlab="Preicipitation",ylab="ECDF",
     main="Quantile Mapping",col="blue")
lines(sort(Yim),rank(sort(Yim))/length(Yim),type="o",col="red")
legend("bottomright",c("GCM","Obs"),lty=1,pch=1,col=c("blue","red"))

plot(sort(Xim),rank(sort(Xim))/length(Xim),type="o",xlim=c(0,max(Yim)),xlab="Preicipitation",ylab="ECDF",
     main="Quantile Mapping (Monthly accumulation, 1979 - 2003)",col="blue")
lines(sort(Yim),rank(sort(Yim))/length(Yim),type="o",col="red")
lines(sort(Qi),rank(sort(Qi))/length(Qi),type="o",col="yellow")
legend("bottomright",c("GCM","Obs","Corrected GCM"),lty=1,pch=1,col=c("blue","red","yellow"))

plot(rank(sort(Yim))/length(Yim),rank(sort(Qi))/length(Qi))

plot(rank(Xim),rank(Qi))
plot(rank(Xim),rank(Yim))

plot(rank(Xim)/length(Xim),rank(Qi)/length(Qi),col="blue",xlab = "GCM (U)", ylab = "corrected GCM (V)")
plot(rank(Xim)/length(Xim),rank(Yim)/length(Yim),col="blue",xlab = "GCM (U)", ylab = "Observation (V)")

plot(rank(Xim)/length(Xim),rank(Ci)/length(Ci),col="blue",xlab = "GCM (U)", ylab = "Observation (V)")
plot(rank(Xim_test)/length(Xim_test),rank(Ci_test)/length(Ci_test),col="blue",xlab = "GCM (U)", ylab = "Observation (V)")



plot(sort(Yim),sort(Qi),col="blue",xlab = "Observation", ylab = "corrected GCM",main="Upper Bound")

plot(sort(Yim),sort(Xim),col="blue",xlab = "Observation", ylab = "GCM",main="Upper Bound")


# 5.2 Bias Correction
plot(1979:2003,Xim,type="o",ylim=c(0,max(Yim)),col="blue",xlab="Years",ylab="Precipitation(mm)",main="Bias Correction")
lines(1979:2003,Yim,type="o",col="red")
lines(1979:2003,Qi,type="o",col="green")
lines(1979:2003,Ci,type="o",col="black")
legend("top",c("GCM","Obs","Quantile Mapping","Nonparametric Copula"),lty=1,pch=1,col=c("blue","red","green","black"))


# 5.2 Bias Correction
plot(2006:2015,Xim_test,type="o",ylim=c(min(Yim_test)-20,max(Xim_test)+20),col="green",xlab="",ylab="Precipitation(mm)",main="Bias Correction")
lines(2006:2015,Yim_test,type="o",col="red")
lines(2006:2015,Qi_test,type="o",col="blue")
lines(2006:2015,Ci_test,type="o",col="black")
legend("topright",c("GCM","Obs","BCSD(QM)","NCBCSD"),lty=1,pch=1,col=c("green","red","blue","black"),cex=0.8)
arrows(2011,55,2012,65, col='black', length=0.1, lwd=3)
text(2011,50,label="NCBCSD",col="black",cex=0.75)
arrows(2011,150,2012,130, col='blue', length=0.1, lwd=3)
text(2011,155,label="BCSD",col="blue",cex=0.75)

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



