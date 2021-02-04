# Bayesian NMF - ICM
# Based on Schmidt et al. 2009, "Bayesian non-negative matrix factorization"
#
# Note:
#   Brouwer argues in his PhD thesis that ICM has the tendency to set most
#   columns to zero, resulting in poor performance. He proposes resetting these
#   values to a small positive value like 0.1.
#


initialize_matrices_icm <- function(X, rank, params) {
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

update_icm <- function(X, W, H, rank, alpha, beta, sigma, k, theta, x_ss) {
    # ICM iteration as proposed by Schmidt et al.
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
    
    
    MINIMUM_TN <- 0.1 # ICM has the tendency to set most columns to 0's;
                      # we reset them to this value.
    
    C <- H %*% t(H)
    D <- X %*% t(H)
    
    for(n in 1:rank) {
        
        nn <- 1:rank
        nn <- nn[-n]
        
        denom <- C[n, n] + .Machine$double.eps
        a_n <- c(D[, n] - W[, nn] %*% C[nn, n] - alpha[, n] * sigma) / denom
        a_n[a_n < 0] <- MINIMUM_TN
        
        W[, n] <- a_n
    }
    
    # Update sigma
    num <-  theta + x_ss + sum(W * (W %*% C - 2 * D)) / 2
    den <- ((nrow(X) * ncol(X)) / 2) + k + 1
    sigma <- num / den
    
    E <- t(W) %*% W
    F_ <- t(W) %*% X
    
    for(n in 1:rank) {
        
        nn <- 1:rank
        nn <- nn[-n]
        
        denom <- E[n, n] + .Machine$double.eps
        b_n <- c(F_[n, ] - E[n, nn] %*% H[nn, ] - beta[n, ] * sigma) / denom
        b_n[b_n < 0] <- 0
        
        H[n, ] <- b_n
    }
    
    return(list(W = W,
                H = H,
                sigma = sigma))
    
}

icm_nmf <- function(X, rank, params = NULL, stop_conds = NULL) {
    # ICM algorithm as proposed by Schmidt et al. Every iteration
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
    #   stop_conds
    #     can be NULL or a list that contains values "type" and "cond" in one
    #     of the following ways:
    #       - type == "iter": cond is an integer specifying the maximum number
    #                         of iterations
    #       - type == "time": cond is a scalar specifying the maximum allowed
    #                         running time of the algorithm (in seconds)
    #       - type == "frob": cond is a scalar, Frobenius threshold. The
    #                         algorithm terminates when the Frobenius norm
    #                         between the original and the calculated matrix is
    #                         less than cond.
    #
    # Output:
    #   A list with the results: matrices W and H, sigma, together
    #   with performance and time for each iteration.
    #  
    
    
    n <- nrow(X)
    p <- ncol(X)
    
    init <- initialize_matrices_icm(X, rank, params)
    
    if(is.null(stop_conds)) {
        stop_conds <- list("type" = "iter", "cond" = 1000)
    }
    
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
    
    
    perf <- list()
    times <- c()
    
    start <- Sys.time()
    
    i <- 1
    
    while(TRUE) {
        
        tmp <- update_icm(X, W, H, rank, alpha, beta, sigma,
                          params$k, params$theta, x_sum)
        W <- tmp$W
        H <- tmp$H
        sigma <- tmp$sigma
        
        time <- as.numeric(Sys.time() - start, units = "secs")
        times <- c(times, time)
        
        # Calculate current performance
        X_pred <- W %*% H
        diff <- X - X_pred
        
        perf[[i]] <- list()
        perf[[i]]$MSE <- sum(diff^2) / (nrow(diff) * ncol(diff))
        perf[[i]]$frob <- norm(diff, "F")
        
        
        if(!stop_check(stop_conds, i, time, perf[[i]]$frob)) {
            break
        }
        
        i <- i + 1
    }
    
    return(list("W" = W,
                "H" = H,
                "sigma" = sigma,
                "perf" = perf,
                "times" = times))
}

stop_check <- function(stop_conds, i, t, perf) {
    # Check if we should continue with NMF iterations.
    #
    # Arguments:
    #   stop_conds
    #     is a list containing values "type" and "cond" in one of the
    #     following ways:
    #       - type == "iter": cond is an integer specifying the maximum number
    #                         of iterations
    #       - type == "time": cond is a scalar specifying the maximum allowed
    #                         running time of the algorithm (in seconds)
    #       - type == "frob": cond is a scalar, Frobenius threshold. The
    #                         algorithm terminates when the Frobenius norm
    #                         between the original and the calculated matrix is
    #                         less than cond.
    #   i
    #     is an integer, current iteration index
    #   perf
    #     is a scalar, the Frobenius norm between the original and the
    #     calculated matrix.
    # Output:
    #   TRUE if we should continue with iterations and the stopping criteria
    #   is not reached, otherwise FALSE
    #
    
    
    if(stop_conds$type == "frob" && perf < stop_conds$cond) return(FALSE)
    if(stop_conds$type == "time" && t > stop_conds$cond) return(FALSE)
    if(stop_conds$type == "iter" && i > stop_conds$cond) return(FALSE)
    
    return(TRUE)
    
}