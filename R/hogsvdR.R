#' @useDynLib hogsvdR
#' @importFrom Rcpp sourceCpp
NULL

## NULL

#' Compute the HO GSVD decomposition of a list of matrices
#' @param D a list of matrices to compute the GSVD decomposition on
#' @param method specification of internal function to use to compute HOGSVD, 'arma' or 'rsimple'
#' @return A list of U, Sigma and V. U and Sigma are lists of matrices
#' @examples 
#' N <- 3
#' nrow <- c(10,10,10)
#' ncol <- 10
#' s <- 1:N
#' D <- lapply(s, function(x) {matrix(rnorm(n=nrow[x]*ncol,mean = 0, sd =10),nrow[x],ncol)})
#' res <- hogsvd(D)
#' 
#' D.reconstruct <- lapply(1:N, function(n) { res$U[[n]] %*% diag(res$Sigma[[n]]) %*% t(res$V) })
#' @export hogsvd
hogsvd <- function(D, method = 'arma') {
  # Check all D are matrices
  if(!all(unlist(lapply(D, class)) == "matrix")) {
    stop('All elements of D have to be of class matrix'); 
  }
  
  # Check all D have same #cols
  ncols <- unlist(lapply(D, function(x) {dim(x)[2]}))
  if(!all(ncols == rep(ncols[1],length(ncols)))) {
    stop('D matrices do not have the same number of columns');
  } 
  
  if (method == 'rsimple') {
    res <- hogsvd.rsimple(D);
  } else if (method == 'arma') {
    res <- hogsvd.rArmadillo(D);
  } else {
      stop('Unknown method parameter');
  }
  
  list(U = res$U, Sigma = res$Sigma, V = res$V);
}


#' Compute the HO GSVD decomposition of a list of matrices
#' @param D a list of matrices to compute the GSVD decomposition on
#' @return A list of U, Sigma, V,  Lambda and S. U and Sigma are lists of matrices
#' @importFrom  MASS ginv
hogsvd.rsimple <- function(D) {
  # Generate named sequence along data
  N <- length(D)
  Nseq <- 1:N
  names(Nseq) <- Nseq

  # Calculate normalised S
  S <- calcNormS.R(D);
    
  # Eigen decomposition of S matrix
  eigen.dec <- eigen(S, symmetric = F)
  
  # The Lambda
  Lambda <- eigen.dec$values
  
  V <- eigen.dec$vectors
  Vinv <- MASS::ginv(eigen.dec$vectors)
  
  # Compute matrices B
  B <- lapply(D, function(x) {
    t( Vinv %*% t(x)  )
  })
  
  # Compute diagonal matrices Sigma
  Sigma <- lapply(Nseq, function(i) {
    apply(B[[i]],2,function(x) sqrt(sum(x^2)))
  })
  
  # Calculate U, the column normalised version of B
  U <- lapply(Nseq, function(i) {
    sweep(B[[i]],2,Sigma[[i]],FUN='/')
  })
  
  # Return
  list( U = U, Sigma = Sigma, V = V, Lambda = Lambda, S = S)
}

#' Calculate the normalised S matrix
#' @param D a list of matrices
calcNormS.R <- function(D) {
  N <- length(D)
  Ddim <- dim(D[[1]])
    
  # Calculate A matrices and their inverses
  A <- lapply(D, function(x) {t(x) %*% x})
  Ainv <- lapply(A, function(x) {MASS::ginv(x)})
  
  # Calculate the S matrix the mean of all the ((A[[i]] %*% Ainv[[j]]) + (A[[j]] %*% Ainv[[i]])) matrices
  S <- matrix(0, Ddim[2], Ddim[2])
  i <- 1;
  while( i <= N ) {
    j <- i + 1;
    while ( j <= N ) {
      S = S + ((A[[i]] %*% Ainv[[j]]) + (A[[j]] %*% Ainv[[i]]))
      j <- j + 1
    }
    i <- i + 1; 
  }
  
  # Normalise S
  S <- S / (N * (N - 1))

  S
}

#' Compute the HO GSVD decomposition of a list of matrices using c acceleration
#' @param D a list of matrices to compute the GSVD decomposition on
#' @return A list of U, Sigma, V,  Lambda and S. U and Sigma are lists of matrices
#' @importFrom  MASS ginv
hogsvd.rArmadillo <- function(D) {

  # Generate named sequence along data
  N <- length(D)
  Nseq <- 1:N
  names(Nseq) <- Nseq
  
  Ddim <- dim(D[[1]])

  # Calculate S in C++
  S <- calcNormS(D, Ddim[2]);

  # Eigen decomposition of S matrix
  eigen.dec <- eigen(S, symmetric = F)

  # The Lambda
  Lambda <- eigen.dec$values
  
  V <- eigen.dec$vectors
  Vinv <- MASS::ginv(eigen.dec$vectors)

  # Compute matrices B
  B <- lapply(D, function(x) {
    t( Vinv %*% t(x)  )
  })

  # Compute diagonal matrices Sigma
  Sigma <- lapply(Nseq, function(i) {
    apply(B[[i]],2,function(x) sqrt(sum(x^2)))
  })

  # Calculate U, the column normalised version of B
  U <- lapply(Nseq, function(i) {
    sweep(B[[i]],2,Sigma[[i]],FUN='/')
  })

  # Return
  list( U = U, Sigma = Sigma, V = V, Lambda = Lambda, S = S)
}
