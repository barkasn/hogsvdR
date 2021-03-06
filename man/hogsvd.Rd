% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hogsvdR.R
\name{hogsvd}
\alias{hogsvd}
\title{Compute the Higher-order generalised singular value decomposition (HOGSVD)  of a list of matrices}
\usage{
hogsvd(D, method = "arma")
}
\arguments{
\item{D}{a list of matrices to compute the GSVD decomposition on}

\item{method}{specification of internal function to use to compute HOGSVD, 'arma' or 'rsimple'}
}
\value{
A list of U, Sigma and V. U and Sigma are lists of matrices
}
\description{
Compute the Higher-order generalised singular value decomposition (HOGSVD)  of a list of matrices
}
\examples{
# Generate 3 matrices to run example on
N <- 3
nrow <- c(10,10,10)
ncol <- 10
s <- 1:N
D <- lapply(s, function(x) {matrix(rnorm(n=nrow[x]*ncol,mean = 0, sd =10),nrow[x],ncol)})

# Perform HO GSVD on the example
res <- hogsvd(D)

# Inspect result
str(res)

# The first U matrix corresponding to D[[1]]
res$U[[1]]

# The first S diagonal matrix correspoinding to D[[1]]
res$S[[1]]

# The shared V matrix
res$V

# Reconstruct the original matrices
D.reconstruct <- lapply(1:N, function(n) { 
   res$U[[n]] \%*\% diag(res$Sigma[[n]]) \%*\% t(res$V) 
})

# Now repeat with the slow algorithm
res.slow <- hogsvd(D, method = 'rsimple')
D.reconstruct.slow <- lapply(1:N, function(n) { 
  res.slow$U[[n]] \%*\% diag(res.slow$Sigma[[n]]) \%*\% t(res.slow$V) 
})

\dontshow{
# Should return TRUE to confirm all arrays could be reconstructed correctly
all(unlist(lapply(1:N, function(n) { 
  max(D[[n]] - res$U[[n]] \%*\% diag(res$Sigma[[n]]) \%*\% t(res$V)) 
})) < 1.e-10)

# Check correct reconstruction with the slow method
all(unlist(lapply(1:N, function(n) { 
  max(D[[n]] - res.slow$U[[n]] \%*\% diag(res.slow$Sigma[[n]]) \%*\% t(res.slow$V)) 
})) < 1.e-10)

}

}
