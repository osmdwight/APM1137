FA2_APM1137_GONZALES
================
Dwight Gonzales
2025-03-03

``` r
poly.neville <- function(x, y, x0) {
  
  n <- length(x)
  q <- matrix(data = 0, n, n)
  q[,1] <- y  
  
  for (i in 2:n) {
    for (j in i:n) {
      q[j,i] <- ((x0 - x[j-i+1]) * q[j,i-1] - (x0 - x[j]) * q[j-1,i-1]) / (x[j] - x[j-i+1])
    }
  }
  
  
  approx_degree_1 <- q[2,2]  
  approx_degree_2 <- q[3,3]  
  approx_degree_3 <- if (n >= 4) q[4,4] else NA 
  
  res <- list(
    'Degree 1 Approximation' = approx_degree_1,
    'Degree 2 Approximation' = approx_degree_2,
    'Degree 3 Approximation' = approx_degree_3,
    'Final Approximated Value' = q[n,n],
    'Neville Iterations Table' = q
  )
  return(res)
}

x <- c(0, 0.25, 0.5, 0.75)
y <- c(1, 1.64872, 2.71828, 4.48169)

poly.neville(x, y, 0.43)
```

    ## $`Degree 1 Approximation`
    ## [1] 2.115798
    ## 
    ## $`Degree 2 Approximation`
    ## [1] 2.376383
    ## 
    ## $`Degree 3 Approximation`
    ## [1] 2.360605
    ## 
    ## $`Final Approximated Value`
    ## [1] 2.360605
    ## 
    ## $`Neville Iterations Table`
    ##         [,1]     [,2]     [,3]     [,4]
    ## [1,] 1.00000 0.000000 0.000000 0.000000
    ## [2,] 1.64872 2.115798 0.000000 0.000000
    ## [3,] 2.71828 2.418803 2.376383 0.000000
    ## [4,] 4.48169 2.224525 2.348863 2.360605

``` r
x<- c(8.1,8.3,8.6,8.7)
y<- c(16.94410,17.56492,18.50515,18.82091)

new_interp <- function(no) {
  n <- length(x)
  mat <- cbind(x, y)
  
  A <- matrix(c(rep(0, ((n) ^ 2 - n))),
              nrow = n,
              ncol = n - 1,
              byrow = T)
  
  for (j in 1:(n - 1))
    for (i in 1:(n - j)) {
      if (j == 1)
        A[i, j] = (y[i + 1] - y[i]) / (x[i + 1] - x[i])
      else
        A[i, j] = (A[(i + 1), (j - 1)] - A[(i), (j - 1)]) / (x[i + j] - x[i])
    }
  
  nam <- as.list(seq(1, n - 1, by = 1))
  colnames(A) <- nam

  A1 <- cbind(mat, A)
  Adf <- as.data.frame(A1)
  
  cat(" Divided Difference Table \n ________________________ \n")
  print(Adf)
  
  sum = y[1]
  for (i in 1:(n - 1)) {
    prod = 1
    for (j in 1:i)
      prod = prod * (no - x[j])
    sum = sum + prod * A[1, i]
  }
  cat("\n f", n - 1, "(", no, ")= ", sum, "\n")
}

new_interp(8.4)
```

    ##  Divided Difference Table 
    ##  ________________________ 
    ##     x        y      1       2            3
    ## 1 8.1 16.94410 3.1041 0.06000 -0.002083333
    ## 2 8.3 17.56492 3.1341 0.05875  0.000000000
    ## 3 8.6 18.50515 3.1576 0.00000  0.000000000
    ## 4 8.7 18.82091 0.0000 0.00000  0.000000000
    ## 
    ##  f 3 ( 8.4 )=  17.87714

``` r
library(pracma)  
library(knitr)   

x <- c(0.1, 0.2, 0.3, 0.4)
f_x <- c(-0.62049958, -0.28398668, 0.00660995, 0.24842440)
f_prime_x <- c(3.58502082, 3.14033271, 2.66668043, 2.16529366)

target_x <- 0.25

# Hermite Interpolation
hermite_interp <- function(x, f, f_prime, x_target) {
  z <- rep(x, each = 2)  
  Q <- matrix(0, nrow = 2 * length(x), ncol = 2 * length(x))
  
  for (i in 1:length(x)) {
    Q[2*i-1,1] <- f[i]
    Q[2*i,1] <- f[i]
    Q[2*i,2] <- f_prime[i]
    if (i > 1) {
      Q[2*i-1,2] <- (Q[2*i-1,1] - Q[2*i-2,1]) / (z[2*i-1] - z[2*i-2])
    }
  }
  
  for (j in 3:(2*length(x))) {
    for (i in j:(2*length(x))) {
      Q[i,j] <- (Q[i,j-1] - Q[i-1,j-1]) / (z[i] - z[i-j+1])
    }
  }
  
  hermite_value <- Q[1,1]
  product_term <- 1
  for (j in 2:(2*length(x))) {
    product_term <- product_term * (x_target - z[j-1])
    hermite_value <- hermite_value + Q[j,j] * product_term
  }
  return(hermite_value)
}

hermite_approx <- hermite_interp(x, f_x, f_prime_x, target_x)
hermite_error <- abs(hermite_approx - (target_x * cos(target_x) - 2 * target_x^2 + 3 * target_x - 1))

# Natural Cubic Spline
natural_spline <- splinefun(x, f_x, method = "natural")
natural_approx <- natural_spline(target_x)
natural_derivative <- natural_spline(target_x, deriv = 1)
natural_error <- abs(natural_approx - (target_x * cos(target_x) - 2 * target_x^2 + 3 * target_x - 1))

clamped_spline <- splinefun(x, f_x, method = "fmm")
clamped_approx <- clamped_spline(target_x)
clamped_derivative <- clamped_spline(target_x, deriv = 1)
clamped_error <- abs(clamped_approx - (target_x * cos(target_x) - 2 * target_x^2 + 3 * target_x - 1))

results_table <- data.frame(
  Method = c("Hermite Interpolation", "Natural Cubic Spline", "Clamped Cubic Spline"),
  Approximation_f_0.25 = c(hermite_approx, natural_approx, clamped_approx),
  Approximate_f_prime_0.25 = c(NA, natural_derivative, clamped_derivative),
  Absolute_Error = c(hermite_error, natural_error, clamped_error)
)

kable(results_table)
```

| Method | Approximation_f_0.25 | Approximate_f_prime_0.25 | Absolute_Error |
|:---|---:|---:|---:|
| Hermite Interpolation | -0.1327676 | NA | 0.0000043 |
| Natural Cubic Spline | -0.1315860 | 2.908355 | 0.0011859 |
| Clamped Cubic Spline | -0.1327697 | 2.907160 | 0.0000022 |
