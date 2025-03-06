
library(patchwork)

alpha <- c(2,5,5,0.50)
beta <- c(5,5,2,0.50)

beta.stats <- function(alpha, beta){
  
  bmean <- alpha/(alpha+beta)
  bvariance <- (alpha*beta)/((alpha+beta)^2*(alpha+beta+1))
  bskewness <- (2*(beta-alpha)*(sqrt(alpha+beta+1)))/((alpha+beta+2)*sqrt(alpha*beta))
  bexcess.kurtosis <- (6*((alpha+beta)^2*(alpha+beta+1)-alpha*beta*(alpha+beta+2)))/
    (alpha*beta*(alpha+beta+2)*(alpha+beta+3))
  df <- data.frame(alpha, beta, bmean, bvariance, bskewness, bexcess.kurtosis)
  
  q1.fig.dat <- tibble(x = seq(-0.25, 1.25, length.out=1000))|>   # generate a grid of points
    mutate(beta.pdf = dbeta(x, alpha, beta),                      # compute the beta PDF
           norm.pdf = dnorm(x, mean = bmean, sd = sqrt(bvariance)))                         # Gaussian distribution
  
  cplot <- ggplot(data= q1.fig.dat)+                                              # specify data
    geom_line(aes(x=x, y=beta.pdf, color="Beta")) +                 # plot beta dist
    geom_line(aes(x=x, y=norm.pdf, color="Gaussian(0.2857, 0.0255)")) +  # plot gaussian dist
    geom_hline(yintercept=0)+ # plot x axis
    theme_bw()+                                                          # change theme
    xlab("x")+                                                           # label x axis
    ylab("Density")+                                                     # label y axis
    scale_color_manual("", values = c("black", "grey"))+                 # change colors
    theme(legend.position = "bottom")+ # move legend to bottom
    ggtitle(paste0("Beta(", alpha, ", ", beta, ")"))

  return(cplot)
}

plots <- list()
for(i in 1:length(alpha)){
  plots[[i]] <- beta.stats(alpha[i], beta[i])
}

full.plot <- wrap_plots(plots, ncol = 2)  
print(full.plot)
