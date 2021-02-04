working_dir <- "insert/path/to/dir"
saving_path <- "../fig/"

setwd(working_dir)


library(mvtnorm)
library(ggplot2)


set.seed(0)

# distr. parameters
m <- c(0, 0)
sigma_mtx <- matrix(c(1, 0.5, 0.5, 1), byrow=T, nrow=2)


# plot parameters
min_limit <- -3.5
max_limit <- 3.5
data.grid <- expand.grid(x = seq(min_limit, max_limit, length.out = 200),
                         y = seq(min_limit, max_limit, length.out = 200))

q <- cbind(data.grid,
           prob = mvtnorm::dmvnorm(data.grid, mean = m, sigma = sigma_mtx))


p <- ggplot() +
    geom_contour(data=q, aes(x=x, y=y, z=prob), alpha=0.5) +
    xlim(min_limit, max_limit) + ylim(min_limit, max_limit) +
    theme_minimal() +
    xlab("x1") + ylab("x2") +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())

tmp <- p + ggtitle("p(x1, x2)")
ggsave(plot=tmp, filename=paste0(saving_path, "bivariate.png"),
       width=5, height=4, units="in")

s <- 100 # number of gibbs samples
middle_plot_max_iter <- 3

# starting point
mu_A <- rep(-2.5, s)
mu_B <- rep(-2.5, s)

p <- p + geom_point(aes(x=mu_A[1], y = mu_B[1]))
ggsave(plot=p, filename=paste0(saving_path, "gibbs0.png"),
       width=5, height=4, units="in")

for(i in 2:s) {
    # x1 | x2
    mu_A[i] <- rnorm(1, 0.5 * mu_B[i-1], 0.75)
    
    horizontal_density <- data.frame(y1 = dnorm(q$x, 0.5 * mu_B[i-1], 0.75) + mu_B[i-1],
                                     x1 = q$x)
    
    # distr & new point
    if(i <= middle_plot_max_iter) {
        tmp <- p + geom_line(data=horizontal_density, aes(x=x1, y=y1),
                             alpha=0.5) +
            ggtitle(paste0("i = ", i - 1, ", x1 | x2"))
        ggsave(plot=tmp, filename=paste0(saving_path, "gibbs", i, "_1_1.png"),
               width=5, height=4, units="in")
        
        tmp <- tmp + geom_point(aes(x=mu_A[i], y = mu_B[i-1])) + # new point
            geom_segment(aes(x=mu_A[i-1], y = mu_B[i-1], xend=mu_A[i],
                             yend = mu_B[i-1]), alpha = 0.2) +
            ggtitle(paste0("i = ", i - 1, ", x1 | x2"))
        ggsave(plot=tmp, filename=paste0(saving_path, "gibbs", i, "_1_2.png"),
               width=5, height=4, units="in")
    }
    
    
    
    # x2 | x1
    mu_B[i] <- rnorm(1, 0.5 * mu_A[i], 0.75)
    
    vertical_density <- data.frame(y1 = q$y,
                                   x1 = dnorm(q$y, 0.5 * mu_A[i], 0.75) + mu_A[i])
    
    # distr & new point
    if(i <= middle_plot_max_iter) {
        tmp <- tmp + geom_path(data=vertical_density, aes(x=x1, y=y1),
                               alpha = 0.5) +
            ggtitle(paste0("i = ", i - 1, ", x2 | x1"))
        ggsave(plot=tmp, filename=paste0(saving_path, "gibbs", i, "_2_1.png"),
               width=5, height=4, units="in")
        
        tmp <- tmp + geom_point(aes(x=mu_A[i], y = mu_B[i])) +
            geom_segment(aes(x=mu_A[i], y = mu_B[i-1],
                             xend=mu_A[i], yend = mu_B[i]), alpha = 0.2) +
            ggtitle(paste0("i = ", i - 1, ", x2 | x1"))
        ggsave(plot=tmp, filename=paste0(saving_path, "gibbs", i, "_2_2.png"),
               width=5, height=4, units="in")
    }
    
        
    
    # final plot
    p <- p + geom_point(data=data.frame(x=mu_A[i], y = mu_B[i-1]),
                        aes(x=x, y=y)) + # new point
        geom_segment(data=data.frame(x=mu_A[i-1], y=mu_B[i-1],
                                     xend=mu_A[i], yend=mu_B[i-1]),
                     aes(x=x, y=y, xend=xend, yend=yend), alpha = 0.2) + #connect
        geom_point(data=data.frame(x=mu_A[i], y = mu_B[i]), aes(x=x, y=y)) +
        geom_segment(data=data.frame(x=mu_A[i], y=mu_B[i-1],
                                     xend=mu_A[i], yend=mu_B[i]),
                     aes(x=x, y=y, xend=xend, yend=yend), alpha = 0.2)
}

p <- p + ggtitle("i = 100")

ggsave(plot=p, filename=paste0(saving_path, "gibbs_final.png"),
       width=5, height=4, units="in")

