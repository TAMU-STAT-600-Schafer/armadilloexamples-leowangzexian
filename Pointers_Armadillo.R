install.packages("RcppArmadillo")
# Add Rcpp and Rcpp Armadillo library
library(Rcpp)
library(RcppArmadillo)

library(Rcpp)
library(RcppArmadillo)
sourceCpp("ArmadilloExamples.cpp")
X = matrix(rnorm(300), 30, 10)
Y = matrix(rnorm(200), 10, 20)
prodCpp = matrix_mult(X, Y)
prodR = X%*%Y
all.equal(prodCpp, prodR)

X = matrix(rnorm(30000), 300, 100)
Y = matrix(rnorm(20000), 100, 200)
library(microbenchmark)
microbenchmark(
  matrix_mult(X, Y),
  X%*%Y )

sourceCpp("ArmadilloExamples.cpp")
set.seed(20386)
X = matrix(rnorm(100), 25, 4)
beta = rep(1, 4)
Y = X %*% beta + rnorm(25, sd = 0.5)
outC = fastLm(X, Y); names(outC)
cbind(outC$coefficients,
      solve(crossprod(X), crossprod(X, Y)))

Y = rnorm(100)
normArmaV(Y, 2)
normArmaV(Y, 1)
X = matrix(rnorm(300), 30, 10)
normArmaM(X, 2)
normArmaM(X, 1)

procrustesR <- function(X, V){ 
  svdXV <- svd(X %*% V)
  U <- tcrossprod(svdXV$u, svdXV$v) 
  return(U)
}
set.seed(308723)
X <- matrix(rnorm(110), 11, 10)
V <- matrix(rnorm(30), 10, 3)
U <- procrustesR(X, V)
U_Cpp <- procrustes(X, V)
sum(abs(U_Cpp - U))

# Source two versions of timesTwo function
sourceCpp("Test_pointers.cpp")

# Copy by value
x = c(1, 3, 5)
y = timesTwo(x)
y
x

# Copy by pointer
x = c(1, 3, 5)
y = timesTwo_pointer(x)
y
x
