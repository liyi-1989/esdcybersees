# Calculate the Schatten p norm
schatten_norm=function(B,p){
  res=rep(0,length(p))
  for(i in 1:length(p)){
    res[i]=(sum((svd(B)$d)^p[i]))^(1/p[i])
  }
  
  return(res)
}

# Compress the matrix A according to partition d, with Schatten p norm
norm_compress=function(A,d,p){
  n=length(d)
  A1=matrix(0,n,n)
  D=c(0,cumsum(d))
  for(i in 1:n){
    for(j in 1:n){
      idx=(D[i]+1):D[i+1]
      idy=(D[j]+1):D[j+1]
      A1[i,j]=schatten_norm(A[idx,idy],p)
    }
  }
  return(A1)
}

# Do simulation: for matrix size n by n, partition d
simu=function(n,d){
  A=matrix(runif(n*n),n,n)
  A=matrix(rnorm(n*n),n,n)
  A=toeplitz(n:1)
  A=toeplitz((-n):(-1))
  #A=t(A)%*%A
  p=(50:150)/50
  r0=r1=rep(0,length(p))
  for(i in 1:length(p)){
    A1=norm_compress(A,d,p[i])
    r0[i]=schatten_norm(A,p[i]) # norm of A
    r1[i]=schatten_norm(A1,p[i]) # norm of compressed A, which is A1
  }
  par(mfrow=c(1,1))
  plot(p,r0,col="red",type="o",xlab="p",ylab="Norm",main=paste0("Schatten p norm: ","Partition:",paste0(d,sep=",",collapse = "")," Size:",n))
  points(p,r1,col="blue",type="o")
  legend("topright",c("Norm of A","Norm of Compressed N(A)"),col=c("red","blue"),lty =1,pch=1)
  print(round(r0-r1,2))
  plot(p,r0-r1,xlab="p",ylab="Norm Difference",col="blue",type="o",main=paste0("Schatten p norm: ","Partition:",paste0(d,sep=",",collapse = "")," Size:",n))
  abline(0,0)
}

simu(9,c(3,3,3))
simu(9,c(2,7))
simu(20,c(5,5,5,5))
simu(100,c(20,30,50))
simu(100,c(50,50))
