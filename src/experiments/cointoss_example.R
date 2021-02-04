library(ggplot2)

working_dir <- "insert/path/to/dir"
saving_path <- "./handout/fig/"

setwd(working_dir)

get_posterior <- function(n, k, a, b) {
    thetas <- seq(0, 1, by = 0.001)
    likelihood <- dbeta(thetas, 1 + k, 1 + n - k)
    prior <- dbeta(thetas, a, b)
    posterior <- dbeta(thetas, k + a, n - k + b)
    
    df <- data.frame("theta" = rep(thetas, 3),
                     "values" = c(prior, posterior, likelihood),
                     "type" = c(rep("prior", length(prior)),
                                rep("posterior", length(posterior)),
                                rep("(scaled)\n likelihood",
                                    length(likelihood))),
                     "prior" = c(rep(paste0("prior: beta(", a, ", ", b, ")"),
                                     length(thetas))))
    return(df)
}


set.seed(0)

df <- NULL

n <- 50 # number of tosses

events <- c(0, 1) # heads or tails
tosses <- sample(events, n, replace = T)

k <- sum(tosses)


# Prior 1: beta(1, 1)
a <- 1
b <- 1

df <- rbind(df, get_posterior(n, k, a, b))


# Prior 2: beta(100, 100)
a2 <- 100
b2 <- 100

df <- rbind(df, get_posterior(n, k, a2, b2))


# Prior 3: beta(1, 10)
a3 <- 1
b3 <- 10

df <- rbind(df, get_posterior(n, k, a3, b3))


df$prior <- factor(df$prior, levels = c("prior: beta(1, 1)",
                                        "prior: beta(100, 100)",
                                        "prior: beta(1, 10)"))


# Plot

p <- ggplot(df, aes(x = theta, y = values, color = type)) +
    geom_line(size=0.8) +
    facet_wrap(~prior, ncol = 1) +
    scale_color_manual(values=c("#B30000", "#2F1736", "#3374C0")) +
    theme_minimal() + theme(legend.position="bottom", legend.box="horizontal") +
    guides(color = guide_legend(title.position="top", title.hjust = 0.5))

ggsave(paste0(saving_path, "posterior_example.png"), p,
       height = 4.5, width = 3, units = "in")

