# Bayesian NMF Gibbs sampler
# Based on Schmidt et al. 2009, "Bayesian non-negative matrix factorization"
#
# Note:
#   Sample_truncated_gaussian is directly adapted from their code. What the
#   authors refer to as rectified gaussian is actually a truncated gaussian.
#   For sanity check we can also sample using function rtruncnorm from the
#   truncnorm package. The results should be the same.
#


if (!require("truncnorm")) install.packages("truncnorm")

initialize_matrices_bnmf <- function(X, rank, params) {
    # Randomly initializes factor matrices W and H. Prepares parameter
    # matrices alpha and beta.
    #
    # Arguments:
    #   X
    #     is our original matrix with data.
    #   rank
    #     is an integer, the factorization rank to achieve.
    #   params
    #     can be NULL or a list that can contain the following values:
    #       - init_bounds: a numeric vector with two values specifying bounds
    #                      for random sampling of the matrix elements. 
    #                      init_bounds[2] must be larger than init_bounds[1]
    #       - lambda_w: a scalar for initializing matrix alpha
    #       - lambda_h: a scalar for initializing matrix beta
    # Output:
    #   A list containing initialized matrices W, H, alpha and beta.
    #
    
    
    if(is.null(params)) params <- list()
    if(is.null(params$init_bounds)) params$init_bounds <- c(0, 0.1)
    if(is.null(params$lambda_w)) params$lambda_w <- 0.1
    if(is.null(params$lambda_h)) params$lambda_h <- 0.1
    
    
    W <- matrix(runif(nrow(X) * rank, params$init_bounds[1],
                      params$init_bounds[2]), nrow = nrow(X))
    H <- matrix(runif(rank * ncol(X), params$init_bounds[1],
                      params$init_bounds[2]), ncol = ncol(X))

    alpha <- matrix(params$lambda_w, nrow = nrow(X), ncol = rank)
    beta <- matrix(params$lambda_h, nrow = rank, ncol = ncol(X))
    
    return(list(W = W,
                H = H,
                alpha = alpha,
                beta = beta))
}

update_bnmf <- function(X, W, H, rank, alpha, beta, sigma, k, theta, x_ss) {
    # Gibbs sampler iteration as proposed by Schmidt et al.
    #
    # Arguments:
    #   X
    #     is our original matrix with data.
    #   W
    #     is the first factor matrix (basis).
    #   H
    #     is the second factor matrix (mixture).
    #   rank
    #     is an integer, the factorization rank to achieve.
    #   alpha
    #     is a matrix defining prior over W.
    #   beta
    #     is a matrix defining prior over H.
    #   sigma
    #     is a scalar, variance of the normal likelihood.
    #   k
    #     is a scalar, prior parameter for sigma.
    #   theta
    #     is a scalar, prior parameter for sigma.
    #   x_ss
    #     is a scalar, sum(X^2) / 2.
    #
    # Output:
    #   A list containing an updated version of W, H and sigma.
    #
    
    
    C <- H %*% t(H)
    D <- X %*% t(H)
    
    for(n in 1:rank) {
        
        nn <- 1:rank
        nn <- nn[-n]
        
        denom <- C[n, n] + .Machine$double.eps
        a_n <- c(D[, n] - W[, nn] %*% C[nn, n] - alpha[, n] * sigma) / denom
        
        s <- sigma / (denom)
        
        #W[, n] <- sample_truncated_gaussian(a_n, s, alpha[, n])
        W[, n] <- rtruncnorm(mean = a_n, sd = s, a = 0, n = nrow(W))
    }
    
    # Update sigma
    scale <-  (theta + x_ss + sum(W * (W %*% C - 2 * D)) / 2)
    shape <- ((nrow(X) * ncol(X)) / 2) + k + 1
    sigma <- 1 / rgamma(1, shape = shape, scale = scale)
    
    E <- t(W) %*% W
    F_ <- t(W) %*% X
    
    for(n in 1:rank) {
        
        nn <- 1:rank
        nn <- nn[-n]
        
        denom <- E[n, n] + .Machine$double.eps
        b_n <- c(F_[n, ] - E[n, nn] %*% H[nn, ] - beta[n, ] * sigma) / denom
        
        s <- sigma / (denom)
        
        #H[n, ] <- sample_truncated_gaussian(b_n, s, beta[n, ])
        H[n, ] <- rtruncnorm(mean = b_n, sd = s, a = 0, n = ncol(H))
    }
    
    return(list(W = W,
                H = H,
                sigma = sigma))
    
}

gibbs_nmf <- function(X, rank, params = NULL,
                      samples = 1000, burnin = 100, thin = 1) {
    
    # NMF Gibbs sampler as proposed by Schmidt et al. Every iteration
    # we also measure time and calculate MSE and Frobenius norm between the
    # original matrix and the matrix product of current factor matrices.
    #
    # Arguments:
    #   X
    #     is our original matrix with data.
    #   rank
    #     is an integer, the factorization rank to achieve.
    #   params
    #     can be NULL or a list that can contain the following values:
    #       - init_bounds: a numeric vector with two values specifying bounds
    #                      for random sampling of the matrix elements. 
    #                      init_bounds[1] must not be larger than init_bounds[2]
    #       - lambda_w: a scalar for initializing matrix alpha
    #       - lambda_h: a scalar for initializing matrix beta
    #       - theta: a scalar, prior parameter for sigma.
    #       - k: a scalar, prior parameter for sigma
    #       - sigma: a scalar, inital value for sigma
    #   samples
    #     is an integer, number of Gibbs samples we want to have in the end.
    #   burnin
    #     is an integer, number of burnin samples (they get discarded and are
    #     not returned).
    #   thin
    #     is an integer used for thining, we keep every thin-th sample.
    #
    # Output:
    #   A list with the results: Gibbs samples for W, H and sigma, together
    #   with performance and time for each sample.
    #   
    
    
    n <- nrow(X)
    p <- ncol(X)
    
    init <- initialize_matrices_bnmf(X, rank, params)
    
    W <- init$W
    H <- init$H
    alpha <- init$alpha
    beta <- init$beta
    
    # if not specified, initialize other parameters 
    if(is.null(params)) params <- list()
    if(is.null(params$theta)) params$theta <- 0
    if(is.null(params$k)) params$k <- 0
    
    if(is.null(params$sigma)) sigma <- 0.1
    else sigma <- params$sigma
    
    
    x_sum <- sum(X^2) / 2
    
    Wm <- list()
    Hm <- list()
    sigmam <- list()
    
    perf <- list()
    times <- c()
    
    start <- Sys.time()
    
    # Gibbs
    for(m in 1:samples) {
        
        for(i in 1:(burnin * (m == 1) + thin * (m > 1))) {
            
            tmp <- update_bnmf(X, W, H, rank, alpha, beta, sigma,
                               params$k, params$theta, x_sum)
            W <- tmp$W
            H <- tmp$H
            sigma <- tmp$sigma
            
        }
        
        time <- as.numeric(Sys.time() - start, units = "secs")
        times <- c(times, time)
        
        Wm[[m]] <- W
        Hm[[m]] <- H
        sigmam[[m]] <- sigma
        
        # Calculate current performance
        X_pred <- W %*% H
        diff <- X - X_pred
        
        perf[[m]] <- list()
        perf[[m]]$MSE <- sum(diff^2) / (nrow(diff) * ncol(diff))
        perf[[m]]$frob <- norm(diff, "F")
            
    }
    
    return(list("W" = Wm,
                "H" = Hm,
                "sigma" = sigmam,
                "perf" = perf,
                "times" = times))
}

sample_truncated_gaussian <- function(m, s, l) {
    # RANDR Random numbers from 
    #   p(x)=K*exp(-(x-m)^2/s-l'x), x>=0 
    #
    # Usage
    #   x = randr(m, s, l)
    #    
    # Copyright 2007 Mikkel N. Schmidt, ms@it.dk, www.mikkelschmidt.dk
    
    A <- (l * s - m) / sqrt(2 * s)
    a <- A > 26
    x <- rep(0, length(m))
    
    y <- runif(length(m))
    
    x[a] <- - log(y[a]) / ((l[a] * s - m[a]) / s)
    a <- !a
    
    R <- erfc(abs(A[a]))
    x[a] <- erfcinv(y[a] * R - (A[a] < 0) * (2 * y[a] + R - 2)) *sqrt(2 * s) +
        m[a] - l[a] * s
    x[!is.finite(x)] <- 0
    x[x < 0] <- 0
    
    return(x)
}


