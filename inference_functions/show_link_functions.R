# this script shows how to visualize the link between TFBM priors and the weights 
# used for model estimation in weightedLASSO and weightedRF
library(tidyr)
library(ggplot2)
library(ggpubr)


ALPHAS <- seq(0,1, by = 0.001)

# the sampling weight to get picked at a decision node
weights_rf <- function(pi, alpha){
  if (pi == 1)
    return(sqrt(1-(alpha-1)^2)+1)
  if(pi==0.5)
    return(1-alpha)
  if(pi == 0)
    return(-sqrt(1-(alpha-1)^2)+1)
}

# the weighte modulating penalty strength in the LASSO
weights_lasso <- function(pi, alpha){
  return(1 - pi * alpha)
}


links <- data.frame("Pi" = rep(c(1, 0.5, 0), each = length(ALPHAS)),
                    "alpha" = rep(ALPHAS, 3),
                    "b. weightedRF" = mapply(weights_rf, rep(c(1, 0.5, 0), each = length(ALPHAS)), rep(ALPHAS, 3))/2,
                    "a. weightedLASSO" = mapply(weights_lasso, rep(c(1, 0.5, 0), each = length(ALPHAS)), rep(ALPHAS, 3)), check.names = F)

link_plot <- links %>%
  gather(key = "method", value = "weight", `b. weightedRF`, `a. weightedLASSO`) %>%
  ggplot(aes(x=alpha, y=weight, color = factor(Pi))) + 
  geom_line(size = 1.1) +
  facet_wrap(~method) + 
  theme_bw()+xlab(expression(alpha))+
  theme(strip.background = element_blank(), 
        strip.text = element_text(size = 15))+
  scale_color_brewer(palette = "Set1", name = "Pior TFBM\nscore")

ggexport(link_plot, filename = "results/link_functions.pdf", width = 7, height = 5)
