working_dir <- "insert/path/to/dir/with/nmf/methods"
setwd(working_dir)


source("bayesian_nmf.R")
source("nmf.R")
source("icm_nmf.R")

if (!require("ggplot2")) install.packages("ggplot2")


get_mse <- function(perf, times, method) {
    mse <- NULL
    
    for(per in perf) {
        mse <- c(mse, per[["MSE"]])
    }
    
    df <- data.frame("mse" = mse,
                     "iter" = 1:length(mse),
                     "time" = times,
                     "method" = method)
    
    df
}


set.seed(0)

n <- 50
p <- 60
r <- 3


W_true <- matrix(rexp(n * r), nrow = n)
H_true <- matrix(rexp(r * p), nrow = r)

X_true <- W_true %*% H_true + rnorm(n * p, mean=0, sd = 0.01)


res_bnmf <- gibbs_nmf(X_true, r, samples = 10000, burnin = 1)
res_nmf <- nmf(X_true, r)
res_icm <- icm_nmf(X_true, r)

df <- NULL
df <- rbind(df, get_mse(res_bnmf[["perf"]], res_bnmf$times, "gibbs"))
df <- rbind(df, get_mse(res_nmf[["perf"]], res_nmf$times, "nmf"))
df <- rbind(df, get_mse(res_icm[["perf"]], res_icm$times, "icm"))

df$method <- as.factor(df$method)

ggplot(df, aes(x = iter, y = mse, color = method)) + geom_line() + xlim(0, 100)

