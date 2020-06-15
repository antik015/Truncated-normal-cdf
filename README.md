## Truncated Gaussian CDF

This is an Rcpp function I wrote to compute the cdf of a multivariate Gaussian distribution with a given mean and a covariance matrix. 
I have used an importance sampling approach to compute the cdf where the imporatnce density is product of recursively truncated univariate 
Gaussians. To use the function the user must install the "Rcpp" library and include it in the current workspace. Then sourcing the function
will do the job. First, download the file on your computer. Example set of following R commands are given below:

> install.packages("Rcpp")

> library(Rcpp)

> sourceCpp("/path/to/the/cpp/file/")

> p = cpp_tmvnorm_prob(rep(0, 2),rep(0, Inf),rep(0, 2),diag(2),5) ## Probability on the positive quadrant of a standard bivariate Gaussian

## Function parameters
The function expects five inputs listed below.

1. l = vector of lower limits
2. u = vector of upper limits
3. m = vector of means
4. S = covariance matrix
5. N = number of importance sampling iterations to compute the probability 

Note that the inputs l and u must have the same number of elements as the vector m. If one of the coordinates is unrestricted then you can
specify the lower and upper limit vectors of that coordinate as (-Inf, Inf). The covariance matrix S must be a p by p matrix where p is the
length of m.

## Computing time
The function is very fast compared to the popular ptmvnorm() function from the R package "tmvtnorm". See the plot below for a time 
comparison.

![](https://github.com/antik015/Truncated-normal-cdf/blob/master/images/tmv_time_comparison_resized.png?raw=true)

