library(tidyr)
library(ggplot2)
library(ggpubr)

ALPHAS <- seq(0,1, by = 0.1)


weights_rf <- function(pi, alpha){
  if (pi == 1)
    return(sqrt(1-(alpha-1)^2)+1)
  if(pi==0)
    return(1-alpha)
  if(pi == 0.5)
    return(-sqrt(1-(alpha-1)^2)+1)
}

weights_sophie <- function(pi, alpha){
  if (pi == 1)
    return(1)
  if(pi==0.5)
    return(1-alpha)
  if(pi == 0)
    return(sqrt(1-(alpha)^2))
}


weights_lasso <- function(pi, alpha){
  return(1 - pi * alpha)
}


links <- data.frame("Pi" = rep(c(1, 0.5, 0), each = 11),
                    "alpha" = rep(ALPHAS, 3),
                    "bRF (proposition oceane)" = mapply(weights_rf, rep(c(1, 0.5, 0), each = 11), rep(ALPHAS, 3)),
                    "bRF (proposition sophie)" = mapply(weights_sophie, rep(c(1, 0.5, 0), each = 11), rep(ALPHAS, 3)),
                    "LASSO-D3S" = mapply(weights_lasso, rep(c(1, 0.5, 0), each = 11), rep(ALPHAS, 3)), check.names = F)

links %>%
  gather(key = "method", value = "weight", `bRF (proposition oceane)`, `bRF (proposition sophie)`, `LASSO-D3S`) %>%
  ggplot(aes(x=alpha, y=weight, color = factor(Pi))) + 
  geom_point() +
  facet_wrap(~method) + 
  geom_smooth() +
  theme_pubr()

