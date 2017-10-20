# Script for benchmarking performance of hogsvd

library('hogsvdR')
library('ggplot2')

getHOGSVDruntime <- function (N, nrow, ncol, method='arma', parallel, nthreads = 1) {
  # Make some data
  s <- 1:N
  D <- lapply(s, function(x) {matrix(rnorm(n=nrow*ncol,mean = 0, sd =10),nrow,ncol)})
  
  # Perform HO GSVD on the example
  start <- Sys.time()
  res <- hogsvd(D, method=method, parallel = parallel, nthreads = nthreads)
  end <- Sys.time()
  
  (end - start)
}

testNinterval <- function (paramSeq, method, nrow = 100, ncol = 100, parallel, nthreads = 1) {
  
  timing.vsN <- unlist(lapply(paramSeq, function(N) {
    cat('.');
    getHOGSVDruntime(N, nrow, ncol, method=method, parallel = parallel, nthreads = nthreads)
  }))
  test.results <- data.frame(
    N = as.numeric(names(timing.vsN)),
    t = timing.vsN
  )
}

testNCOLinterval <- function (paramSeq, method, N = 3, nrow = 100, parallel, nthreads = 1) {
  timing.vsN <- unlist(lapply(paramSeq, function(ncol) {
    cat('.');
    getHOGSVDruntime(N, nrow, ncol, method=method, parallel = parallel, nthreads = nthreads)
  }))
  
  test.results <- data.frame(
    N = as.numeric(names(timing.vsN)),
    t = timing.vsN
  )
}

testNROWinterval <- function (paramSeq, method, N = 3, ncol = 100, parallel, nthreads = 1) {
  timing.vsN <- unlist(lapply(paramSeq, function(nrow) {
    cat('.');
    getHOGSVDruntime(N, nrow, ncol, method=method, parallel = parallel, nthreads = nthreads)
  }))
  
  test.results <- data.frame(
    N = as.numeric(names(timing.vsN)),
    t = timing.vsN
  )
}


# Run time vs N
nrep <- 1

Nseq <- seq(10, 100,10);
names(Nseq) <- Nseq

testN.rsimple.noparl <- lapply(nrep, function(x) {testNinterval(Nseq,'rsimple',100, 100, parallel = F)})
testN.rsimple.noparl <- do.call(rbind, testN.rsimple.noparl)
testN.rsimple.noparl$method = c('rsimple.noparl')

testN.arma.noparl <- lapply(nrep, function(x) {testNinterval(Nseq,'arma',100, 100, parallel = F)})
testN.arma.noparl <- do.call(rbind, testN.arma.noparl)
testN.arma.noparl$method = c('arma.noparl')

ParamSeq <- seq(10, 100,10);
names(ParamSeq) <- ParamSeq

testN.rsimple.parallel <- lapply(nrep, function(x) {testNinterval(ParamSeq,'rsimple',100, 100, parallel = T, nthreads = 4)})
testN.rsimple.parallel <- do.call(rbind, testN.rsimple.parallel)
testN.rsimple.parallel$method = c('rsimple.parallel')

testN.arma.parallel <- lapply(nrep, function(x) {testNinterval(ParamSeq,'arma',100, 100, parallel = T, nthreads = 4)})
testN.arma.parallel <- do.call(rbind, testN.arma.parallel)
testN.arma.parallel$method = c('arma.parallel')

ggplot(rbind(testN.rsimple.noparl, testN.arma.noparl,testN.rsimple.parallel, testN.arma.parallel), aes(x=N, y=t, color=method)) + geom_point() + theme_bw() + 
  ggtitle('HO GSVD Runtimes vs N (100x100 matrices)') + scale_x_continuous(name='Number of samples (N)') + scale_y_continuous(name='Time (s)');
#ggsave('runtimeVsN.png');


# Run time vs NCOLS
nrep <- 1

paramSeq <- seq(100,1000,100);
names(paramSeq) <- paramSeq;

testNcol.rsimple.noparl <- lapply(1:nrep, function(x) {testNCOLinterval(paramSeq = paramSeq, method = 'rsimple', N = 3 , nrow = 100, parallel = F)})
testNcol.rsimple.noparl <- do.call(rbind, testNcol.rsimple.noparl)
testNcol.rsimple.noparl$method = c('rsimple.noparl')

testNcol.arma.noparl <- lapply(1:nrep, function(x) {testNCOLinterval(paramSeq = paramSeq, method = 'arma', N = 3 , nrow = 100, parallel = F)})
testNcol.arma.noparl <- do.call(rbind, testNcol.arma.noparl)
testNcol.arma.noparl$method = c('arma.noparl')

testNcol.rsimple.parallel <- lapply(1:nrep, function(x) {testNCOLinterval(paramSeq = paramSeq, method = 'rsimple', N = 3 , nrow = 100, parallel =T, nthreads = 4 )})
testNcol.rsimple.parallel <- do.call(rbind, testNcol.rsimple.parallel)
testNcol.rsimple.parallel$method = c('rsimple.parallel')

testNcol.arma.parallel <- lapply(1:nrep, function(x) {testNCOLinterval(paramSeq = paramSeq, method = 'arma', N = 3 , nrow = 100, parallel = T, nthreads = 4 )})
testNcol.arma.parallel <- do.call(rbind, testNcol.arma.parallel)
testNcol.arma.parallel$method = c('arma.parallel')


ggplot(rbind(testNcol.rsimple.noparl, testNcol.arma.noparl, testNcol.rsimple.parallel, testNcol.arma.parallel), aes(x=N, y=t, color=method)) + geom_point() + theme_bw() +
  ggtitle('HO GSVD Runtimes vs ncol (3, 100xncol matrices)') + scale_x_continuous(name='Number of columns') + scale_y_continuous(name='Time (s)');

#ggsave('runtimeVsNCOL.png');


# Run time vs NROWS
nrep <- 1

paramSeq <- seq(500,2500,500);
names(paramSeq) <- paramSeq;

testNcol.rsimple.noparl <- lapply(1:nrep, function(x) {testNROWinterval(paramSeq = paramSeq, method = 'rsimple', N = 3 , ncol = 1000, parallel = F)})
testNcol.rsimple.noparl <- do.call(rbind, testNcol.rsimple.noparl)
testNcol.rsimple.noparl$method = c('rsimple.noparl')

testNcol.arma.noparl <- lapply(1:nrep, function(x) {testNROWinterval(paramSeq = paramSeq, method = 'arma', N = 3 , ncol = 1000, parallel = F)})
testNcol.arma.noparl <- do.call(rbind, testNcol.arma.noparl)
testNcol.arma.noparl$method = c('arma.noparl')

testNcol.rsimple.parallel <- lapply(1:nrep, function(x) {testNROWinterval(paramSeq = paramSeq, method = 'rsimple', N = 3 , ncol = 1000, parallel =T, nthreads = 4 )})
testNcol.rsimple.parallel <- do.call(rbind, testNcol.rsimple.parallel)
testNcol.rsimple.parallel$method = c('rsimple.parallel')

testNcol.arma.parallel <- lapply(1:nrep, function(x) {testNROWinterval(paramSeq = paramSeq, method = 'arma', N = 3 , ncol = 1000, parallel = T, nthreads = 4 )})
testNcol.arma.parallel <- do.call(rbind, testNcol.arma.parallel)
testNcol.arma.parallel$method = c('arma.parallel')

ggplot(rbind(testNcol.rsimple.noparl, testNcol.arma.noparl, testNcol.rsimple.parallel, testNcol.arma.parallel), aes(x=N, y=t, color=method)) + geom_point() + theme_bw() +
  ggtitle('HO GSVD Runtimes vs NROW (3, NROW x 1000)') + scale_x_continuous(name='Number of rows') + scale_y_continuous(name='Time (s)');

#ggsave('runtimeVsNROW.png');
