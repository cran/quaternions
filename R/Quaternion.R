
Q <- function(...) UseMethod("Q")

Q.Q <- function(r,...) r

Q.default <- function(r=0,i=0,j=0,k=0,...) {
  S <- r+i+j+k
  S[1:length(S)]<-0
  structure(list(r=r+S,i=i+S,j=j+S,k=k+S),class="Q")
}

complex2Q <- function(x) {
  Q(Re(x),Im(x))
}

print.Q <- function(x,...) {
  a <- paste(x$r,"+",x$i,"i+",x$j,"j+",x$k,"k",sep="")
  dim(a) <- dim(x$r)
  dimnames(a) <- dimnames(x$r)
  print(a,...,quote=FALSE)
  invisible(x)
}

dim.Q <- function(x) dim(x$r)

"dim<-.Q" <- function(x,value) {
  dim(x$r) <-value
  dim(x$i) <-value
  dim(x$j) <-value
  dim(x$k) <-value
  x
}

dimnames.Q <- function(x) dimnames(x$r)

"dimnames<-.Q" <- function(x,value) {
  dimnames(x$r) <-value
  dimnames(x$i) <-value
  dimnames(x$j) <-value
  dimnames(x$k) <-value
  x
}

names.Q <- function(x) names(x$r)

"names<-.Q" <- function(x,value) {
  names(x$r) <- value
  names(x$i) <- value
  names(x$j) <- value
  names(x$k) <- value
  x
}


"[.Q" <- function(x,...) Q(x$r[...],x$i[...],x$j[...],x$k[...])
#"[[.Q" <- function(x,...) Q(x$r[[...]],x$i[[...]],x$j[[...]],x$k[[...]])
"[<-.Q" <- function(x,...,value) {
  value <- Q(value)
  x$r[...]<- value$r
  x$i[...]<- value$i
  x$j[...]<- value$j
  x$k[...]<- value$k
  x
}


"+.Q" <- function(x,y) {
  x<-Q(x)
  y<-Q(y)
  Q(x$r+y$r,x$i+y$i,x$j+y$j,x$k+y$k)
}

add.Q <- function(x,y) {
  x<-Q(x)
  y<-Q(y)
  Q(x$r+y$r,x$i+y$i,x$j+y$j,x$k+y$k)
}


"-.Q" <- function(x,y) {
  if( missing(y) )
    return(Q(-x$r,-x$i,-x$j,-x$k))
  x<-Q(x)
  y<-Q(y)
  Q(x$r-y$r,x$i-y$i,x$j-y$j,x$k-y$k)
}


"*.Q" <- function(x,y) {
  x<-Q(x)
  y<-Q(y)
  Q(x$r*y$r-x$i*y$i-x$j*y$j-x$k*y$k,
    x$r*y$i+x$i*y$r+x$j*y$k-x$k*y$j,
    x$r*y$j+x$j*y$r+x$k*y$i-x$i*y$k,
    x$r*y$k+x$k*y$r+x$i*y$j-x$j*y$i
    )
}

"/.Q" <- function(x,y) {
  if( is.numeric(y) )
    return( Q(x$r/y,x$i/y,x$j/y,x$k/y) )
  x <- Q(x)
  y <- inv.Q(Q(y))
  x*y
}


"Re<-" <-function(z,...,value) UseMethod("Re<-")
"Re<-.complex" <-function(z,...,value) Re(value)+Im(z)*1i


Re.Q <- function(z) Q(z)$r
"Re<-.Q"  <- function(z,...,value) {
  z$r <- value
  z
}

"Im<-" <-function(z,...,value) UseMethod("Im<-")
"Im<-.complex" <-function(z,...,value) Re(value)*1i+Re(z)

Im.Q <- function(z) {
  z<-Q(z);
  Q(r=0,z$i,z$j,z$k)
}

"Im<-.Q" <- function(z,...,value) {
  value <- Q(value)
  z <- Q(z)
  value <- value+z-z
  z$i <- value$i
  z$j <- value$j
  z$k <- value$k
}



Qr <- Q(1,0,0,0)
Qi <- Q(0,1,0,0)
Qj <- Q(0,0,1,0)
Qk <- Q(0,0,0,1)


inv.Q <- function(x,...) {
  Conj.Q(x)/(norm.Q(x,each=TRUE)^2)
}

norm <- function(x,...) UseMethod("norm")

norm.Q <- function(x,...,each=FALSE) {
  if( each )
    sqrt(x$r^2+x$i^2+x$j^2+x$k^2)
  else
    sqrt(sum(x$r^2+x$i^2+x$j^2+x$k^2))
}


Conj.Q <- function(z) {
  x<-z
  Q(if( is.complex(x$r) ) Conj(x$r)  else x$r ,
    if( is.complex(x$i) ) -Conj(x$i) else -(x$i),
    if( is.complex(x$j) ) -Conj(x$j) else -(x$j),
    if( is.complex(x$k) ) -Conj(x$k) else -(x$k))
}

t.Q <- function(x) {
  Q(t(x$r),t(x$i),t(x$j),t(x$k))
}

adj <- function(x,...) UseMethod("adj")

adj.numeric <- function(x,...) t(x)

adj.complex <- function(x,...) {
  Conj(t(x))
}
adj.Q <- function(x,...) {
  Conj(t(x))
}

"%*%" <- function(x,y) UseMethod("%*%")
"%*%.default" <-  function(x,y) base::"%*%"(x,y)


"%*%.Q" <- function(x,y) {
  Q(x$r %*% y$r -x$i %*% y$i -x$j %*% y$j -x$k %*% y$k,
    x$r %*% y$i +x$i %*% y$r +x$j %*% y$k -x$k %*% y$j,
    x$r %*% y$j +x$j %*% y$r +x$k %*% y$i -x$i %*% y$k,
    x$r %*% y$k +x$k %*% y$r +x$i %*% y$j -x$j %*% y$i
    )
}

matrix <- function(data,...) UseMethod("matrix")
matrix.default <- base::matrix
formals(matrix.default) <- c(formals(matrix.default),alist(...=))

matrix.Q <- function(data,...) {
  data <- Q(data)
  Q(matrix(data$r,...),
    matrix(data$i,...),
    matrix(data$j,...),
    matrix(data$k,...))
}

length.Q <- function(x) length(x$r)

c.Q      <- function(x,...) {
  kk <- lapply(list(x,...),Q)
  xx <- lapply(c(r="r",i="i",j="j",k="k"),
               function(w) {
                 unlist(lapply(kk,function(k) k[[w]]))
               }
               )
  
  structure(xx,class="Q")
}


rep.Q <- function(x,times,...) {
  Q(rep(x$r,times,...),
    rep(x$i,times,...),
    rep(x$j,times,...),
    rep(x$k,times,...))
}

solve.Q <- function(a,b,...) {
  a <- Q(a)
  b <- Q(b)
  A <- rbind(
             cbind(A$r,-A$i,-A$j,-A$k),
             cbind(A$i, A$r,-A$k, A$j),
             cbind(A$j, A$k, A$r,-A$i),
             cbind(A$k,-A$j, A$i, A$r))
  B <- rbind(t(t(b$r)),t(t(b$i)),t(t(b$j)),t(t(b$k)))
  n <- nrow(B)/4
  X <- solve(A,B,...)
  erg <- Q(X[1:n,],X[n+1:n,],X[2*n+1:n,],X[3*n+1:n,])
  dim(erg) <- dim(b)
  erg
}


