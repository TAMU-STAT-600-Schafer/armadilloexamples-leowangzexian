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
