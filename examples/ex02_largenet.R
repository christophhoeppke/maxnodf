###
# In this example we examine the quality and performance of maxnodf,
# when applied to larger networks.
###

# We commence this example by loading the maxnodf library and constructing
# a small test network.
library('maxnodf')
N <- 250
net <- matrix(0.0, N, N)

# The library maxnodf computes the maximal NODF value over networks of the
# same class as net. We say that networks are in the same class if 
#   (1) The number of links in both networks are eqal
#           sum(net2) = sum(net)
#   (2) The number of rows and columns of both networs are equal
#           all(dim(web3) == dim(net))
# For NODF evaluations we require a minimum of one link in every row and column.
# To ensure this we require that all entries in the first row and first column 
# are set to 1.
net[1, ] <- 1
net[, 1] <- 1

# Now we fill the matrix with randomised values
fillrate = 0.05
vals <- matrix(runif(N*N), ncol=N) 
net[vals<fillrate] <- 1.0

# To estimate the performance of maxnodf we load the library 'microbenchmark'
library('microbenchmark')
microbenchmark(maxnodf(web=net), times=100)

# On the developer test system [INTEL i7-8350U] we measured a mean runtime of
# 15.93 milliseconds.
