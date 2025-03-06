
################################################################
##### Task 1 ###################################################
################################################################

library(patchwork)

beta.stats <- function(alpha, beta){
  
  mean <- alpha/(alpha+beta)
  variance <- (alpha*beta)/((alpha+beta)^2*(alpha+beta+1))
  skewness <- (2*(beta-alpha)*(sqrt(alpha+beta+1)))/((alpha+beta+2)*sqrt(alpha*beta))
  excess.kurtosis <- ((6*((alpha-beta)^2*(alpha+beta+1)-alpha*beta*(alpha+beta+2)))/
    (alpha*beta*(alpha+beta+2)*(alpha+beta+3)))
  return(tibble(alpha, beta, mean, variance, skewness, excess.kurtosis))
}

beta.plots <- function(alpha, beta){
  tib <- beta.stats(alpha, beta)
  q1.fig.dat <- tibble(x = seq(-0.25, 1.25, length.out=1000))|>   # generate a grid of points
    mutate(beta.pdf = dbeta(x, alpha, beta),                      # compute the beta PDF
           norm.pdf = dnorm(x, mean = tib$mean, sd = sqrt(tib$variance)))                         # Gaussian distribution
  
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
alpha <- c(2,5,5,0.50)
beta <- c(5,5,2,0.50)
beta.table <- tibble()

for(i in 1:length(alpha)){
  plots[[i]] <- beta.plots(alpha[i], beta[i]) + xlim(0,1) + ylim(0,3)
  btable <- beta.stats(alpha[i],beta[i])
  beta.table <- bind_rows(beta.table, btable)
}
view(beta.table)
full.plot <- wrap_plots(plots, ncol = 2)  
print(full.plot)

################################################################
##### Task 2 ###################################################
################################################################

beta.moment <- function(alpha, beta, k, centered){
  
  xmean <- alpha/(alpha+beta)
  fX.uncentered <- function(x) dbeta(x,alpha,beta)*x^k
  fX.centered   <- function(x) dbeta(x,alpha,beta)*(x-xmean)^k
  
  if(centered==F){
    kth.moment <- integrate(fX.uncentered,0,1)
  }
  else{
    kth.moment <- integrate(fX.centered,0,1)
  }
  return(kth.moment$value)
}

moment.stats <- function(alpha,beta){
  
  mean <- beta.moment(alpha, beta, 1, F)
  variance <- beta.moment(alpha, beta, 2, T)
  skewness <- beta.moment(alpha, beta, 3, T)/(variance^(3/2))
  excess.kurtosis <- (beta.moment(alpha, beta, 4, T)/(variance^2))-3
  return(tibble(alpha, beta, mean, variance, skewness, excess.kurtosis))
  
}

moment.table <- tibble()
for(i in 1:length(alpha)){
  mtable <- moment.stats(alpha[i],beta[i])
  moment.table <- bind_rows(moment.table, mtable)
}
view(moment.table)

################################################################
##### Task 3 ###################################################
################################################################

set.seed(7272)
sample.size <- 500
library(e1071)

bsample <- function(alpha, beta){
  beta.sample <- rbeta(n = sample.size,  # sample size
                       shape1 = alpha,   # alpha parameter
                       shape2 = beta)    # beta parameter
  return(beta.sample)
}

hist.plots <- function(alpha, beta){
  beta.sample <- bsample(alpha,beta)
  beta.df <- data.frame(x = beta.sample)
  hist <- ggplot(beta.df, aes(x)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30) +
    geom_density(color = "red") +
    stat_function(fun = dbeta, 
                  args = list(shape1 = alpha, shape2 = beta), 
                  color = "blue") +
    theme_bw() +
    xlab("x") +
    ylab("Density") +
    geom_hline(yintercept = 0)+
    ggtitle(paste0("Beta(", alpha, ", ", beta, ")"))

  summary <- beta.df %>%
    summarize(
      mean = mean(beta.sample),
      variance = var(beta.sample),
      skewness = skewness(beta.sample),
      excess.kurtosis = kurtosis(beta.sample))
  
  return(hist)
}

histograms <- list()
for(i in 1:length(alpha)){
  histograms[[i]] <- hist.plots(alpha[i], beta[i]) + xlim(0,1) + ylim(0,3)
}
full.hist <- wrap_plots(histograms, ncol = 2)  
print(full.hist)

################################################################
##### Task 4 ###################################################
################################################################

library(cumstats)
