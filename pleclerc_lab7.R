
################################################################
##### Task 1 ###################################################
################################################################

library(patchwork)
library(tidyverse)

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
  q1.fig.dat <- tibble(x = seq(-0.25, 1.25, length.out=1000))|>  
    mutate(beta.pdf = dbeta(x, alpha, beta))
  
  cplot <- ggplot(data= q1.fig.dat)+                
    geom_line(aes(x=x, y=beta.pdf)) +  
    geom_hline(yintercept=0) + 
    theme_bw()+  
    xlab("x")+   
    ylab("Density")+  
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
#view(beta.table)
full.plot <- wrap_plots(plots, ncol = 2)  
ggsave("task1.png", width = 10, height = 8)
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
#view(moment.table)

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
    geom_histogram(aes(y = after_stat(density)), bins = 20) +
    geom_density(color = "red") +
    stat_function(fun = dbeta, 
                  args = list(shape1 = alpha, shape2 = beta), 
                  color = "blue") +
    theme_bw() +
    xlab("x") +
    ylab("Density") +
    geom_hline(yintercept = 0)+
    ggtitle(paste0("Beta(", alpha, ", ", beta, ")"))

  summary <- beta.df |>
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
ggsave("task3.png", width = 10, height = 8)
print(full.hist)

################################################################
##### Task 4 ###################################################
################################################################

library(cumstats)

calc.summary <- function(alpha, beta){

beta.sample <- bsample(alpha,beta)
beta.df <- data.frame(x = beta.sample)

csum <- beta.df |>
  mutate(
    index = row_number(),
    cmean = cummean(x),
    cvariance = cumvar(x),
    cskewness = cumskew(x),
    ckurtosis = cumkurt(beta.sample))

return(csum)

}

csummary <- calc.summary(2,5)

sp1 <- ggplot(csummary, aes(x = index, y = cmean)) +
    geom_line(color = "blue") +
    geom_hline(yintercept = beta.stats(2,5)$mean, color="red", linetype="dashed") +
    theme_bw() +
    xlab("Sample Size") +
    ylab("Cumulative Mean")
  
sp2 <- ggplot(csummary, aes(x=index, y=cvariance)) +
    geom_line(color = "blue") +
    geom_hline(yintercept = beta.stats(2,5)$variance, color="red", linetype="dashed") +
    theme_bw() +
    xlab("Sample Size") +
    ylab("Cumulative Variance")
  
sp3 <- ggplot(csummary, aes(x = index, y = cskewness)) +
    geom_line(color = "blue") +
    geom_hline(yintercept = beta.stats(2,5)$skewness, color = "red", linetype="dashed") +
    theme_bw() +
    xlab("Sample Size") +
    ggtitle("Cumulative Skewness")
  
sp4 <- ggplot(csummary, aes(x = index, y = ckurtosis)) +
    geom_line(color = "blue") +
    geom_hline(yintercept = beta.stats(2,5)$excess.kurtosis+3, color = "red", linetype = "dashed") +
    theme_bw() +
    xlab("Sample Size") +
    ggtitle("Cumulative Kurtosis")

# I don't think I did this correctly.
for(i in 2:50){
  set.seed(7272+i)
  bstats <- calc.summary(2,5)
  sp1 <- sp1 + 
    geom_line(data = bstats, aes(x = index, y = cmean), color=i)
  sp2 <- sp2 + 
    geom_line(data = bstats, aes(x = index, y = cvariance), color=i)
  sp3 <- sp3 + 
    geom_line(data = bstats, aes(x = index, y = cskewness), color=i)
  sp4 <- sp4 + 
    geom_line(data = bstats, aes(x = index, y = ckurtosis), color=i)
}

sumplots <- list(sp1, sp2, sp3, sp4)
full.plot <- wrap_plots(sumplots, ncol = 2)
ggsave("task4.png", width = 10, height = 8)
print(full.plot) 

################################################################
##### Task Five ################################################
################################################################

library(e1071)
alpha <- 2
beta <- 5
n <- 500
stats <- tibble(mean = numeric(), variance = numeric(),
                skewness = numeric(), kurtosis = numeric())

for(i in 1:1000){
  set.seed(7272+i)
  sample <- rbeta(n, alpha, beta)
  mean.val <- mean(sample)
  var.val <- var(sample)
  skew.val <- skewness(sample)
  kurt.val <- kurtosis(sample)-3
  
  stats <- bind_rows(stats, tibble(mean = mean.val, 
                                   variance = var.val, 
                                   skewness = skew.val, 
                                   kurtosis = kurt.val))
}

p1 <- ggplot(stats, aes(x = mean)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "blue", alpha = 0.5) +
  geom_density(color = "black") + ggtitle("Sampling Distribution of Mean") +
  geom_hline(yintercept = 0)

p2 <- ggplot(stats, aes(x = variance)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "blue", alpha = 0.5) +
  geom_density(color = "black") + ggtitle("Sampling Distribution of Variance") +
  geom_hline(yintercept = 0)

p3 <- ggplot(stats, aes(x = skewness)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "blue", alpha = 0.5) +
  geom_density(color = "black") + ggtitle("Sampling Distribution of Skewness") +
  geom_hline(yintercept = 0)

p4 <- ggplot(stats, aes(x = kurtosis)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "blue", alpha = 0.5) +
  geom_density(color = "black") + ggtitle("Sampling Distribution of Kurtosis") +
  geom_hline(yintercept = 0)

library(patchwork)
(p1 | p2) / (p3 | p4)
ggsave("task5.png", width = 10, height = 8)

################################################################
##### Task Six #################################################
################################################################

death.rates <- read_csv("worldbank.csv")
head(death.rates)

death.rates <- death.rates |>
  mutate(across(-c(1:4), ~./1000))
head(death.rates)

################################################################
##### Task Seven ###############################################
################################################################

library(nleqslv)

MOM.exp <- function(data, par){
  alpha <- par[1]
  beta <- par[2]
  
  m1 <- mean(data, na.rm=T)
  EX1 <- alpha/(alpha+beta)
  
  m2 <- mean(data^2, na.rm=T)
  EX2 <- (alpha+1)*alpha/((alpha+beta+1)*(alpha+beta))
  return(c(EX1-m1, EX2-m2))
}

moms <- nleqslv(x = c(1,1),
        fn = MOM.exp,
        data=death.rates$`2022`)
(alpha.hat.mom <- moms$x[1])
(beta.hat.mom <- moms$x[2])

s = seq(0,0.025,length.out=1000)
mom.dist <- tibble(x=s) |>
  mutate(pdf = dbeta(x=s,alpha.hat.mom,beta.hat.mom))

llbeta <- function(data, par, neg=F){
  alpha <- par[1]
  beta <- par[2]
  loglik <- sum(dbeta(x=data, alpha, beta, log=T), na.rm=T)
  return(ifelse(neg, -loglik, loglik))
}

mles <- optim(par = c(1,1), 
      fn = llbeta,
      data=death.rates$`2022`,
      neg=T)
(alpha.hat.mle <- mles$par[1])
(beta.hat.mle <- mles$par[2])

mle.dist <- tibble(x=s) |>
  mutate(pdf = dbeta(x=s,alpha.hat.mle,beta.hat.mle))

ggplot() +
  geom_histogram(data=death.rates,
                 aes(x=`2022`,
                     y=after_stat(density)),
                 breaks=seq(0,0.024,0.002)) +
  ggtitle("Histogram of Death Rates") +
  geom_hline(yintercept=0) +
  theme_bw() +
  xlab("Death Rates") +
  ylab("Density") +
  geom_line(data = mom.dist, aes(x=x,y=pdf), color="red",size=1) +
  geom_line(data = mle.dist, aes(x=x,y=pdf), color="blue",size=1)
ggsave("task7.png", width = 10, height = 8)

################################################################
##### Task Eight ###############################################
################################################################

alpha <- 8
beta <- 950
n <- 266
mom.est <- tibble(alpha = numeric(), beta = numeric())
mle.est <- tibble(alpha = numeric(), beta = numeric())

for(i in 1:1000){
  set.seed(7272+i)
  sample <- rbeta(n, alpha, beta)
  
  moms <- nleqslv(x = c(1,1),
                  fn = MOM.exp,
                  data=sample)
  mom.est <- bind_rows(mom.est, 
                       tibble(alpha = moms$x[1], 
                              beta = moms$x[2]))
  
  mles <- optim(par = c(1,1), 
                fn = llbeta,
                data=sample,
                neg=T)
  mle.est <- bind_rows(mle.est, 
                       tibble(alpha = mles$par[1], 
                              beta = mles$par[2]))
}

p1 <- ggplot(mom.est, aes(x = alpha)) +
  geom_density(fill = "blue", alpha = 0.5) + 
  ggtitle("MOM Alpha Density") +
  geom_hline(yintercept = 0)
p2 <- ggplot(mom.est, aes(x = beta)) +
  geom_density(fill = "blue", alpha = 0.5) + 
  ggtitle("MOM Beta Density") +
  geom_hline(yintercept = 0)
p3 <- ggplot(mle.est, aes(x = alpha)) +
  geom_density(fill = "red", alpha = 0.5) + 
  ggtitle("MLE Alpha Density") +
  geom_hline(yintercept = 0)
p4 <- ggplot(mle.est, aes(x = beta)) +
  geom_density(fill = "red", alpha = 0.5) + 
  ggtitle("MLE Beta Density") +
  geom_hline(yintercept = 0)

library(patchwork)
(p1 | p2) / (p3 | p4)
ggsave("task8.png", width = 10, height = 8)

metrics <- function(estimate, value) {
  bias <- mean(estimate) - value
  precision <- 1 / var(estimate)
  mse <- var(estimate) + bias^2
  return(tibble(Bias = bias, 
                Precision = precision, 
                MSE = mse))
}

mom_metrics_alpha <- metrics(mom.est$alpha, alpha)
mle_metrics_alpha <- metrics(mle.est$alpha, alpha)
mom_metrics_beta <- metrics(mom.est$beta, beta)
mle_metrics_beta <- metrics(mle.est$beta, beta)

metrics.est <- tibble(
  Parameter = rep(c("Alpha", "Beta"), each = 2),
  Method = rep(c("MOM", "MLE"), 2),
  Bias = c(mom_metrics_alpha$Bias, 
           mle_metrics_alpha$Bias, 
           mom_metrics_beta$Bias, 
           mle_metrics_beta$Bias),
  Precision = c(mom_metrics_alpha$Precision, 
                mle_metrics_alpha$Precision, 
                mom_metrics_beta$Precision, 
                mle_metrics_beta$Precision),
  MSE = c(mom_metrics_alpha$MSE, 
          mle_metrics_alpha$MSE, 
          mom_metrics_beta$MSE, 
          mle_metrics_beta$MSE))

print(metrics.est)
