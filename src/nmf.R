# NMF
# Based on Lee and Seung 1999, "Learning the parts of objects by non-negative
# matrix factorization"
#


initialize_matrices_nmf <- function(X, rank, params) {
    # Randomly initializes factor matrices W and H.
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
    # Output:
    #   A list containing initialized matrices W and H.
    #
    
    
    if(is.null(params)) params <- list()
    if(is.null(params$init_bounds)) params$init_bounds <- c(0, 0.1)
    
    W <- matrix(runif(nrow(X) * rank, params$init_bounds[1],
                      params$init_bounds[2]), nrow = nrow(X))
    H <- matrix(runif(rank * ncol(X), params$init_bounds[1],
                      params$init_bounds[2]), ncol = ncol(X))
    
    return(list(W = W,
                H = H))
}

nmf <- function(X, rank, params = NULL, stop_conds = NULL) {
    # Classic NMF method using multiplicative updates proposed by Lee and Seung.
    # Every iteration we also measure time and calculate MSE and Frobenius norm
    # between the original matrix and the matrix product of current factor
    # matrices. Stopping condition can be specified by the user, the default
    # value is 1000 iterations.
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
    # Output:
    #   A list with the results: matrices W and H, together with performance
    #   and time for each sample.
    #
    
    
    n <- nrow(X)
    p <- ncol(X)
    
    init <- initialize_matrices_nmf(X, rank, params)
    
    if(is.null(stop_conds)) {
        stop_conds <- list("type" = "iter", "cond" = 1000)
    }
    
    W <- init$W
    H <- init$H
    
    perf <- list()
    times <- c()
    
    start <- Sys.time()
    
    i <- 1
    
    while(TRUE) {
        
        # Update W
        XHT <- X %*% t(H)
        WHHT <- W %*% H %*% t(H)
        
        W <- W * XHT / WHHT
        
        # Update H
        WTX <- t(W) %*% X
        WTWH <- t(W) %*% W %*% H
        
        H <- H * WTX / WTWH
        
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