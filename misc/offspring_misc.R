rm(list = ls())
#######################################
###### error-in-variables models ######
#######################################
# https://www.r-bloggers.com/errors-in-variables-models-in-stan/
setwd("Littorina_offspring/misc/")
n.reps <- 3
n.repeated <- 10
n <- 50

# true covariate values
x <- runif(n, -3, 3)
y <- x + rnorm(n)  # alpha=0, beta=1, sdy=1
plot(x,y)

# random subset to perform repeat covariate measurements
which.repeated <- sample(n, n.repeated)
xsd <- 1  # measurement error
xerr <- rnorm(n + (n.repeated * (n.reps - 1)), 0, xsd)

# indx assigns measurements to sample units
indx <- c(1:n, rep(which.repeated, each = n.reps - 1))
indx <- sort(indx)
nobs <- length(indx)
xobs <- x[indx] + xerr
plot(x[indx], xobs,
     xlab = "True covariate value",
     ylab = "Observed covariate value")
abline(0, 1, lty = 2)
segments(x0 = x[indx], x1 = x[indx],
         y0 = x[indx], y1 = xobs, col = "red")
abline(v = x[which.repeated], col = "green", lty = 3)

# write the .stan file
cat("
data{
  int n;
  int nobs;
  real xobs[nobs];
  real y[n];
  int indx[nobs];
}

parameters {
  real alpha;
  real beta;
  real sigmay;
  real sigmax;
  real x[n];
}

model {
  // priors
  alpha ~ normal(0,100);
  beta ~ normal(0,100);
  sigmay ~ uniform(0,1000);
  sigmax ~ uniform(0,1000);
  
  // model structure  
  for (i in 1:nobs){
    xobs[i] ~ normal(x[indx[i]], sigmax);
  }
  for (i in 1:n){
    y[i] ~ normal(alpha + beta*x[i], sigmay);
  }
}
  ",
file = "latent_x.stan")

library(rstan)
library(modeest)
stan_d <- c("y", "xobs", "nobs", "n", "indx")
chains <- 3
iter <- 1000
thin <- 1
mod1 <- stan(file = "latent_x.stan", data = stan_d,
             chains = chains, iter = iter,
             thin = thin)
########################
###### cline plot ######
########################
library(pracma)
x <- seq(-6, 6, length.out = 101)
y2 <- abs(sigmoid(x, a = 6) - 1)

svg(filename = "figures/wild_phen_cline.svg", width = 8, height = 5)
plot(x, y2, type = "l", col = "#f1a340", lwd = 4,
     xlab = "Distance", ylab = "Phenotype")
grid()
text(x = -4, y = 0.2, labels = "CRAB", cex = 2)
text(x = 4, y = 0.2, labels = "WAVE", cex = 2)
dev.off()

# svg(filename = "figures/lab_phen_cline_1.svg", width = 8, height = 5)
# getwd()
# setwd()
pdf(file = "../../PhD/meetings/Littorina/Tja_Jan_2020/figures/lab_phen_cline_1.pdf", width = 8, height = 5)
plot(x, y2, type = "l", col = "#f1a340", lwd = 4,
     xlab = "Distance", ylab = "Phenotype")
# y1 <- sigmoid(x, a = -1)
grid()
# abline(h = 0.4, lwd=4, col = "#998ec3")
lines(x, rep(x = 0.4, 101), col = "#998ec3", lwd = 4)
text(x = -4, y = 0.2, labels = "CRAB", cex = 2)
text(x = 4, y = 0.2, labels = "WAVE", cex = 2)
legend(4, 1, legend=c("wild", "lab"),
       col=c("#f1a340", "#998ec3"), lwd = c(4,4), cex=1.25)
dev.off()

# svg(filename = "figures/lab_phen_cline_2.svg", width = 8, height = 5)
pdf(file = "../../PhD/meetings/Littorina/Tja_Jan_2020/figures/lab_phen_cline_2.pdf", width = 8, height = 5)
# y1 <- sigmoid(x, a = -5, b = 2)
# mysig = function(x) {
#   1 / (1 + exp(-x))
# }
# plot(seq(-5, 5, 0.01), mysig(seq(-5, 5, 0.01)+2), col='blue')

# y2[y2>0.91] = y2[y2>0.91] + 0.3
# y2[y2<0.19] = y2[y2<0.19] + 0.3
# plot(x, y1, type = "l", col = "#f1a340", lwd = 4,
#      xlab = "Distance", ylab = "Phenotype")
y2 <- abs(sigmoid(x, a = 6, b = 0.5) - 1)
y3 = y2+0.2
y2 = y3-0.2
length(y2)
y3[55:length(y3)] = y2[55:length(y2)]
plot(x, y3, type = "l", col = "#998ec3", lwd = 4,
     xlab = "Distance", ylab = "Phenotype", ylim = c(0,1.2))
y1 <- abs(sigmoid(x, a = 5) - 1)
var(y2)
var(y3)
var(y1)
# var(y1)
# y3 <- sigmoid(x, a = -5, b = -2)
grid()
# lines(x, y1, col = "#f1a340", lwd = 4)
lines(x, y1, col = "#f1a340", lwd = 4)
text(x = -4, y = 0.2, labels = "CRAB", cex = 2)
text(x = 4, y = 0.2, labels = "WAVE", cex = 2)
legend(4, 1.2, legend=c("wild", "lab"),
       col=c("#f1a340", "#998ec3"), lwd = c(4,4), cex=1.25)
dev.off()

pdf(file = "../../PhD/meetings/Littorina/Tja_Jan_2020/figures/lab_phen_cline_3.pdf", width = 8, height = 5)
x <- seq(-6, 6, length.out = 101)
cline_pl <- function(position,centre,w,left,right,sl,sc,sr){
  d <- position-centre
  p_x <- 1/(1+exp(0-4*(d)/w)) 
  # p_x is frequency cline as in HZAR NB p=0 for left, 1 for right (corresponding to scaling phenotypes with crab=0, wave=1)
  z_x <- left + (right-left)*p_x  
  # z_x is expected phenotype, wave-crab always positive, given scaling
  s_x <- sqrt(sl^2 + 4*p_x*(1-p_x)*sc^2 + (p_x^2)*(sr^2-sl^2))
  # unimodal variance model as in HZAR, sx is SD for wave/hybrid/crab (assumes variances are additive, unlike Cfit)
  # minusll <- -sum(dnorm(phen,z_x,s_x,log=T))
  # return(minusll)
  phen_cline = data.frame(phen_cline = z_x, sd_cline = s_x, position = position)
  return(phen_cline)
}
y2 = cline_pl(position = x, centre = 0, w = 4, left = 1.2, right = 0, sl = 0.1, sc = 0.1, sr = 0.1)[, "phen_cline"]
plot(x, y2, type = "l", col = "#998ec3", lwd = 4,
     xlab = "Distance", ylab = "Phenotype", ylim = c(0,1.2))

y1 = sigmoid(x, a = -7)
var(y1)
var(y2)
# var(y1)
# y3 <- sigmoid(x, a = -5, b = -2)
grid()
# lines(x, y1, col = "#f1a340", lwd = 4)
lines(x, y1, col = "#f1a340", lwd = 4)
text(x = -4, y = 0.2, labels = "CRAB", cex = 2)
text(x = 4, y = 0.2, labels = "WAVE", cex = 2)
legend(4, 1.2, legend=c("wild", "lab"),
       col=c("#f1a340", "#998ec3"), lwd = c(4,4), cex=1.25)
dev.off()

##########################
###### scatter plot ######
##########################
x = seq(0, to = 1, length.out = 101) + rnorm(n = 101, sd = 0.1)
var(x)
y = seq(0, to = 1, length.out = 101) + rnorm(n = 101, sd = 0.1)
var(y)
zx = seq(0, to = 1, length.out = 101)
zy = rnorm(n = 101, mean =  0.4, sd = 0.1)

pdf(file = "../../PhD/meetings/Littorina/Tja_Jan_2020/figures/lab_vs_wild_1.pdf", width = 5, height = 5)
plot(1, type="n", xlab="wild sample", ylab="lab-reared sample", xlim=c(0, 1), ylim=c(0, 1))
grid()
abline(a = 0, b = 1, lwd=4, lty='dashed', col='grey70')
abline(h = 0.4, lwd=4)
# ggplot() +
#   geom_abline(slope = 1, col="#f1a340", size = 3) +
#   geom_hline(yintercept = 0.4, col="#998ec3", size = 3) +
#   xlim(c(0,1))+
#   ylim(c(0,1)) +
#   labs(x = "wild sample", y = "lab-reared sample") +
#   theme(axis.title = element_text(size = 18),
#         axis.text = element_text(size=11),
#         axis.ticks = element_line(size = 0.7),
#         panel.background = element_blank(),
#         panel.border = element_rect(colour = "black", fill=NA, size=0.5),
#         axis.line = element_line(size = 0.2, linetype = "solid",
#                                  colour = "black"),
#         panel.grid = element_line(colour = "gray70", size = 0.2))
dev.off()

xx = seq(-1, to = 1, length.out = 101) + rnorm(n = 101, sd = 0.1)
var(x)
var(xx)
yy = seq(-1, to = 1, length.out = 101) + rnorm(n = 101, sd = 0.1)
var(y)
var(yy)
plot(xx, yy)
points(x, y, col='red')

pdf(file = "../../PhD/meetings/Littorina/Tja_Jan_2020/figures/lab_vs_wild_2.pdf", width = 5, height = 5)
plot(1, type="n", xlab="wild sample", ylab="lab-reared sample", xlim=c(0, 1), ylim=c(0, 1.2))
grid()
abline(a = 0, b = 1, lwd=4, lty='dashed', col='grey70')
abline(a = 0, b = 2.5, lwd=4)
dev.off()

pdf(file = "../../PhD/meetings/Littorina/Tja_Jan_2020/figures/lab_vs_wild_3.pdf", width = 5, height = 5)
plot(1, type="n", xlab="wild sample", ylab="lab-reared sample", xlim=c(0, 1), ylim=c(0, 1.2))
grid()
abline(a = 0, b = 1, lwd=4, lty='dashed', col='grey70')
segments(x0 = 0, y0 = 0, x1 = 0, y1 = 0.4, lwd=4)
segments(x0 = 0, y0 = 0.4, x1 = 0.6, y1 = 0.6, lwd=4)
segments(x0 = 0.6, y0 = 0.6, x1 = 1, y1 = 0.8, lwd=4)
# segments(x0 = 0.8, y0 = 0.8, x1 = 1, y1 = 0.9, lwd=4)
segments(x0 = 1, y0 = 0.8, x1 = 1, y1 = 1.2, lwd=4)
dev.off()

######################
###### stanfile ######
######################
# transformed parameters {
#   for (i in 1:N) {
#     x_meas[i] = 1 / (sigma_x[i] * sqrt(2 * pi)) * exp(- (x - mu)^2 / (2 * sigma_x[i]^2))
#   }
# }
