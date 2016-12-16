library(maptools)
data(wrld_simpl)

#============================ Plot EOF ==============================================
V5=read.table("./results_eof/ta_month_1_svd_eigenvectors_top_5.txt",sep=",")
L=read.table("./results_eof/ta_month_1_svd_eigenvalues.txt",sep=",")
V5=read.table("./results_eof/ta_month_1_svd_threshold_1_eigenvectors_top_5.txt",sep=",")
L=read.table("./results_eof/ta_month_1_svd_threshold_1_eigenvalues.txt",sep=",")
V5=read.table("./results_eof/ta_month_1_svd_threshold_10_eigenvectors_top_5.txt",sep=",")
L=read.table("./results_eof/ta_month_1_svd_threshold_10_eigenvalues.txt",sep=",")


V5=V5[,ncol(V5):1]
L=rev(L[,1])
L.sum=sum(L[L>0])
L5=L[1:5]


eof1=matrix(V5[,1],144,90,byrow = T)
eof2=matrix(V5[,2],144,90,byrow = T)
eof3=matrix(V5[,3],144,90,byrow = T)
eof4=matrix(V5[,4],144,90,byrow = T)
eof5=matrix(V5[,5],144,90,byrow = T)


# Lon=ifelse(lon<=180,lon,lon-360)
# Lon=Lon[c(73:144,1:72)]

pdf("0.pdf")
par(mfrow=c(3,2),
    #    oma = c(1,1,0,0),
    mar = c(2,2,2,2))

image(Lon,lat,eof1[c(73:144,1:72),],main=paste0(round(L5[1]/L.sum,2)*100,"%"))
plot(wrld_simpl,add=T)
contour(Lon,lat,eof1[c(73:144,1:72),],add=T,col="purple")

image(Lon,lat,eof2[c(73:144,1:72),],main=paste0(round(L5[2]/L.sum,2)*100,"%"))
plot(wrld_simpl,add=T)
contour(Lon,lat,eof2[c(73:144,1:72),],add=T,col="purple")

image(Lon,lat,eof3[c(73:144,1:72),],main=paste0(round(L5[3]/L.sum,2)*100,"%"))
plot(wrld_simpl,add=T)
contour(Lon,lat,eof3[c(73:144,1:72),],add=T,col="purple")

image(Lon,lat,eof4[c(73:144,1:72),],main=paste0(round(L5[4]/L.sum,2)*100,"%"))
plot(wrld_simpl,add=T)
contour(Lon,lat,eof4[c(73:144,1:72),],add=T,col="purple")

image(Lon,lat,eof5[c(73:144,1:72),],main=paste0(round(L5[5]/L.sum,2)*100,"%"))
plot(wrld_simpl,add=T)
contour(Lon,lat,eof5[c(73:144,1:72),],add=T,col="purple")

plot(1:10,L[1:10]/L.sum,type="o",col="blue",main="Variance Explained")
dev.off()


# get GCM lat lon
# ncin="./gcm_train/pr_day_GFDL-CM3_historical_r1i1p1_20050101-20051231.nc"
# ncin = nc_open(ncin)
# lon = ncvar_get(ncin,"lon")
# lat = ncvar_get(ncin,"lat")


#============================ Plot map ==============================================

library(ggplot2)
library(ggmap)
lon <- c(mylon,rev(mylon))
lat <- c(mylat,mylat)
df <- as.data.frame(cbind(lon,lat))

# getting the map
mapgilbert <- get_map(location = c(lon = mean(df$lon), lat = mean(df$lat)), zoom = 7,
                      maptype = "terrain", scale = 2)

# plotting the map with some points on it
ggmap(mapgilbert) +
  #geom_point(data = df, aes(x = lon, y = lat, fill = "yellow", alpha = 1), size = 4, shape = 21) +
  guides(fill=FALSE, alpha=FALSE, size=FALSE)+
  annotate('segment', x=lon[1],xend=lon[1],y=lat[1],yend=lat[2],colour=I('red'), size = 1) +
  annotate('segment', x=lon[2],xend=lon[2],y=lat[1],yend=lat[2],colour=I('red'), size = 1) +
  annotate('segment', x=lon[1],xend=lon[2],y=lat[1],yend=lat[1],colour=I('red'), size = 1) +
  annotate('segment', x=lon[1],xend=lon[2],y=lat[2],yend=lat[2],colour=I('red'), size = 1) 


#============================ Grid Label ==============================================


n=10
m=10
X_label=matrix(1:(n*m),n,m,byrow = T)
C=matrix(0,n*m,n*m)

for(i1 in 1:n){
  for(j1 in 1:m){
    for(i2 in 1:n){
      for(j2 in 1:m){
        i=(i1-1)*m+j1
        j=(i2-1)*m+j2
        C[i,j]=exp(-0.4*base::norm(as.matrix(c(i1-i2,j1-j2)),type="2"))
      }
    }
  }
}

library(cape)

par(mfrow=c(1,2))
#image(rotate.mat(C),main="S: Covariance Matrix Structure(decay with distance)")
image(rotate.mat(C),main="S: Covariance Matrix Structure")
image(rotate.mat(X_label),main="X: 2D Spatial Stations Labels", col  = terrain.colors(1))
for(i in 1:n){
  for(j in 1:m){
    count=(i-1)*m+j
    text(j/(m-1)-0.1,1-i/(n-1)+0.1,count)
  }
}

load("~/Documents/esdcybersees/results_eof/Scov.RData")
image(rotate.mat(Scov),main="Scov")
load("~/Documents/esdcybersees/results_eof/Scor.RData")
image(rotate.mat(Scor),main="Scor")





