require(quaternions)

checker <- function(x,y,...) UseMethod("checker")

checker.Q <- function(x,y,delta=1E-14) {
  print(x)
  if( norm(x-y) > delta) {
    print(y)
    stop("Error")
  }
    
}

checker.numeric <- function(x,y) {
  print(x)
  if( sum((x-y)^2) > 1E-14) {
    print(y)
    stop("Error")
  }
}

checker.complex <- function(x,y) {
  print(x)
  if( as.numeric(sum(Conj(x-y)*(x-y))) > 1E-14) {
    print(y)
    stop("Error")
  }
}

checker.character <- function(x,y) {
  print(x)
  if( !all(x==y) > 1E-14) {
    print(y)
    stop("Error")
  }
}


q1 <- Q(rnorm(25),rnorm(25),rnorm(25),rnorm(25))
q2 <- Q(rnorm(25),rnorm(25),rnorm(25),rnorm(25))
q3 <- Q(rnorm(25),rnorm(25),rnorm(25),rnorm(25))
q4 <- Q(rnorm(25),rnorm(25),rnorm(25),rnorm(25))

checker(q1-q2+q2,q1)
checker(q1-q2,-(q2-q1))
checker(q1-q1,0)
checker(q1+q2,q2+q1)
checker(Conj(q1)*Conj(q2),Conj(q2*q1))

A <- matrix(q1,nrow=5)
B <- matrix(q2,nrow=5)
C <- matrix(q3,nrow=5)
D <- matrix(q4,nrow=5)



checker( adj(A)%*%B , adj(adj(B)%*%A) )

checker(class(q1),"Q")
checker(class(A),"Q")

checker(complex2Q(1+2i),Q(1,2))

checker(complex2Q(rep(1+2i,10)),Q(rep(1,10),rep(2,10)))

checker(Q(q1),q1)

checker(Q(q1$r,q1$i,q1$j,q1$z),q1)

print(q1)

checker(solve(A, A %*% B), B,delta=1E-9)

dim(q1) <- c(5,5)
dimnames(q1)<- list(paste(1:5),paste(1:5))
checker(dim(q1),c(5,5))
checker(dimnames(q1)[[1]],paste(1:5))
checker(dimnames(q1)[[2]],paste(1:5))
names(q2)<-paste(1:25)
checker(names(q2),paste(1:25))

X <- A
checker(X,A)
X[1:4,]<-B[1:4,]
X[5,]<- B[5,]
checker(X,B)

checker( A * B / B, A)

checker(Re(A),A$r)
checker(Im(A),A$i*Qi+A$j*Qj+A$k*Qk)
checker(A,A$r+A$i*Qi+A$j*Qj+A$k*Qk)


checker( inv.Q(B) * B * A , A)

checker(Re(adj(q3)%*%q3),norm(q3)^2)

checker(adj(C),Conj(t(C)))


checker(length(A),25)


checker(c(q1,q2,q3)[26:50],q2)

checker(rep(q1,3)[26:50],q1)

checker(rep(q1,each=3),q1[rep(1:25,each=3)])











