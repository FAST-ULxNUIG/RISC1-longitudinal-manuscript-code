Test for construct fd from scores
================

``` r
source(here::here("code", "functions", "decenter_fd_around_new_mean.R"))
source(here::here("code", "functions", "construct_fd_from_scores.R"))       
library(fda)     # CRAN v5.5.1
```

    ## Loading required package: splines

    ## Loading required package: Matrix

    ## Loading required package: fds

    ## Loading required package: rainbow

    ## Loading required package: MASS

    ## Loading required package: pcaPP

    ## Loading required package: RCurl

    ## Loading required package: deSolve

    ## 
    ## Attaching package: 'fda'

    ## The following object is masked from 'package:graphics':
    ## 
    ##     matplot

``` r
library(funData) # CRAN v1.3-8
```

    ## 
    ## Attaching package: 'funData'

    ## The following object is masked from 'package:stats':
    ## 
    ##     integrate

``` r
library(MFPCA)   # CRAN v1.3-9

set.seed(1996)
# set up a grid of values:
grid_test <- seq(0, 2 * pi, length.out = 101)
# simulate bivariate functional data:
simulated_MultiFunData <- simMultiFunData(type = "split",
                                          argvals = list(grid_test, grid_test),
                                          M = 5,
                                          eFunType = "Fourier",
                                          eValType = "exponential",
                                          N = 500)


plot(simulated_MultiFunData$simData)
```

![](test_construct_fd_from_scores_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
# convert to fda objects:
fd_dim1 <- funData2fd(simulated_MultiFunData$simData[[1]])
fd_dim2 <- funData2fd(simulated_MultiFunData$simData[[2]])

par(mfrow = c(1, 2))
plot(fd_dim1); plot(fd_dim2)
```

    ## [1] "done"

![](test_construct_fd_from_scores_files/figure-gfm/unnamed-chunk-1-2.png)<!-- -->

    ## [1] "done"

``` r
# construct multivariate fd object:
mfd_coef_array <- array(data = NA, dim = c(dim(fd_dim1$coefs), 2))
mfd_coef_array[,,1] <- fd_dim1$coefs
mfd_coef_array[,,2] <- fd_dim2$coefs
mfd_obj <- fd(mfd_coef_array, basisobj = fd_dim1$basis)

par(mfrow = c(1, 2))
plot(mfd_obj)
```

![](test_construct_fd_from_scores_files/figure-gfm/unnamed-chunk-1-3.png)<!-- -->

    ## [1] "done"

``` r
# Do pca and compare with MFPCA package: ----------------------------------
pca_fd <- pca.fd(fdobj = mfd_obj, nharm = 5)


MFPCA <- MFPCA(mFData = simulated_MultiFunData$simData, M = 5,
      uniExpansions = list(list(type = "uFPCA"),
                           list(type = "uFPCA")))

par(mfrow = c(2, 2))
plot(pca_fd$harmonics)
```

    ## [1] "done"

``` r
plot(MFPCA$functions)
```

![](test_construct_fd_from_scores_files/figure-gfm/unnamed-chunk-1-4.png)<!-- -->![](test_construct_fd_from_scores_files/figure-gfm/unnamed-chunk-1-5.png)<!-- -->

``` r
# First test: -------------------------------------------------------------
# if we use the observed scores we should get back original curves:
scores <- apply(pca_fd$scores, c(1, 2), sum)
reconstruct_in_sample <- construct_fd_from_scores(pca_fd_obj = pca_fd, scores_matrix = scores, K = 5)

par(mfrow = c(2, 2))
plot(mfd_obj[1:5])
```

    ## [1] "done"

``` r
plot(reconstruct_in_sample[1:5])
```

![](test_construct_fd_from_scores_files/figure-gfm/unnamed-chunk-1-6.png)<!-- -->

    ## [1] "done"

``` r
par(mfrow = c(1, 1))
plot(mfd_obj - reconstruct_in_sample)
```

![](test_construct_fd_from_scores_files/figure-gfm/unnamed-chunk-1-7.png)<!-- -->![](test_construct_fd_from_scores_files/figure-gfm/unnamed-chunk-1-8.png)<!-- -->

    ## [1] "done"

``` r
# Looks good.


# Second test: ------------------------------------------------------------
# lets give just a single observation with 0 on every score
# we should get back mean observation.
reconstruct_0 <- construct_fd_from_scores(pca_fd_obj = pca_fd,
                                          scores_matrix = matrix(0, nrow = 1, ncol = 5),
                                          K = 5)


par(mfrow = c(2, 2))
plot(pca_fd$meanfd, col = "blue")
```

    ## [1] "done"

``` r
plot(reconstruct_0, col = "red")
```

![](test_construct_fd_from_scores_files/figure-gfm/unnamed-chunk-1-9.png)<!-- -->

    ## [1] "done"

``` r
# Third Test: -------------------------------------------------------------
# Same as above but we don't add back mean -- should get zero functions
reconstruct_0_centered <- construct_fd_from_scores(
  pca_fd_obj = pca_fd,
  scores_matrix = matrix(0, nrow = 1, ncol = 5),
  K = 5,
  add_back_mean = FALSE)


par(mfrow = c(1, 2))
plot(reconstruct_0_centered, col = "red")
```

![](test_construct_fd_from_scores_files/figure-gfm/unnamed-chunk-1-10.png)<!-- -->

    ## [1] "done"

``` r
# again, good.



# Fourth test: -------------------------------------------------------------
# Expect an erro if the dimensions of PCA and scores are don't match
testthat::expect_error(
  reconstruct_0_centered <- construct_fd_from_scores(
    pca_fd_obj = pca_fd,
    scores_matrix = matrix(rnorm(3), nrow = 1, ncol = 3),
    K = 5)
)
# and now if K agrees with scores rather than pca_fd
testthat::expect_error(
  reconstruct_0_centered <- construct_fd_from_scores(
    pca_fd_obj = pca_fd,
    scores_matrix = matrix(rnorm(3), nrow = 1, ncol = 3),
    K = 3)
)
```
