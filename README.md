# Description
Higher-Order Generalized Singular Value Decomposition for R

# Installation Instructions
```
library('devtools')
install_github('barkasn/hogsvdR')
```

# Example
```
N <- 3
nrow <- c(10,10,10)
ncol <- 10
s <- 1:N
D <- lapply(s, function(x) {matrix(rnorm(n=nrow[x]*ncol,mean = 0, sd =10),nrow[x],ncol)})
res <- hogsvd(D)

D.reconstruct <- lapply(1:N, function(n) { res$U[[n]] %*% diag(res$Sigma[[n]]) %*% t(res$V) })
```
