rm(list = ls())

.packages = c("optparse", "tidyr", "knitr", "kableExtra", "ggrepel", "bbmle", "Rmisc",
              "flextable", "officer", "ggcorrplot", "data.table", "car", "MASS", "ggplot2", "boot", "dplyr")
# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
# Load packages into session
lapply(.packages, require, character.only=TRUE)

# define island/experiment
island = "CZA"

# read data
dat_dir = paste0(island, "_off_SW/", island, "_off_final_data/")
dat_off = read.csv(file = paste0(dat_dir, island, "_off_all_phenos_main_20200406.csv"))
# remove weight typo in CZA
if (island=="CZA") {
  dat_off = dat_off[-as.integer(rownames(dat_off[dat_off$snail_ID=='C_60', ])), ]
}
# prepare data for downstream analysis
dat_off = separate(data = dat_off, col = "snail_ID", into = c("pop", "ID"), sep = "_")
dat_off[, "generation"] = 1
dat_off[which(nchar(as.character(dat_off$ID)) == 2), "generation"] = 0
dat_off$weight_cuberoot = dat_off$weight_g^(1/3)
off_phen = c("bold_score", "mean_thickness", "weight_cuberoot", "PC2")
tbl1 = dat_off[, c(off_phen, "size_mm", "generation", "pop", "ID", "sex")]
rm(dat_off)

table(tbl1$sex)
table(tbl1[tbl1$generation==0, "sex"])

# boxplot raw data
ggplot(data = tbl1, aes(x = pop, y = size_mm, fill = sex)) +
  geom_boxplot(position = position_dodge(0.9), fatten = NULL) +
  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y.., x = pop),
               width = 0.75, size = 1, position = position_dodge(0.9)) +
  labs(fill = "") +
  theme(legend.position = 'top')
tbl1[tbl1$pop=='A' & tbl1$sex=='juvenile', ]

ggplot(data = tbl1, aes(x = sex, y = size_mm, fill = sex)) +
  facet_wrap(~pop) +
  geom_boxplot(position = position_dodge(0.9), fatten = NULL) +
  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y.., x = sex),
               width = 0.75, size = 1, position = position_dodge(0.9)) +
  theme(axis.text.x = element_text(angle = 315, hjust=0))

(p3 = ggplot(data = tbl1, aes(x = pop, y = size_mm, fill = pop)) +
  facet_wrap(~sex) +
  geom_boxplot(position = position_dodge(0.9), notch = FALSE) +
  # stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y.., x = pop),
  #              width = 0.75, size = 1, position = position_dodge(0.9)) +
  labs(x = "", fill = "") +
  theme(legend.position = "null", legend.text = element_text(size = 11),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 11),
        panel.background = element_blank(),
        strip.background = element_rect(fill = "#91bfdb"),
        strip.text = element_text(size = 11),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(size = 0.2, linetype = "solid",
                                 colour = "black"),
        panel.grid = element_line(colour = "gray70", size = 0.2)))
ggsave(file=paste0(dirname(dat_dir), "/", island, "_off_final_figures/", island, "_off_boxplot_size_maturity.svg"),
       plot=p3, width=10, height=8)
ggsave(file=paste0(dirname(dat_dir), "/", island, "_off_final_figures/", island, "_off_boxplot_size_maturity.pdf"),
       plot=p3, width=10, height=8)

# prepare data for model fitting
table(tbl1$sex)
table(tbl1[tbl1$generation==0, "sex"])
tbl1 = tbl1[tbl1$sex!="missing", ]

# model 1: maturity is 0,1 for juvenile,adult
tbl1$maturity = 1
tbl1[tbl1$sex!="female" & tbl1$sex!="male", "maturity"] = 0
table(tbl1$maturity)
# remove NAs
tbl1 = tbl1[!is.na(tbl1$size_mm), ]

## Model 1a
# function for the likelihood of the data given that the logit(probability of maturity) is a linear function of size
# ‘mat’ is 0,1 for juvenile,adult
# ‘mean’ and ‘slope’ are the regression parameters
s_mat <- function(size, mat, mean, slope) {
  # multuply or divide by the slope
  logit_p <- (size - mean) / slope
  p <- exp(logit_p) / (exp(logit_p)+1)
  minusll <- -sum(dbinom(mat, 1, p, log = TRUE))
  return(minusll)
}

# can fit this to the data without clinal variation – so, overall size vs maturity relationship – here using BBMLE
tbl1$size_log = log(tbl1$size_mm)
mean(tbl1[, "size_log"])
theta.init <- list(mean = 2, slope = 2)

mle.s_mat <- mle2(s_mat, start = theta.init, data = list(size = tbl1$size_log,
                                                         mat = tbl1$maturity))
summary(mle.s_mat)
# as.vector(coef(mle.s_mat)["mean"])
AIC(mle.s_mat)
confint(mle.s_mat)

s_mat_fit <- function(x, mod) {
  # multiply or divide by the slope
  logit_p <- (x - as.vector(coef(mod)["mean"])) / as.vector(coef(mod)["slope"])
  p <- exp(logit_p) / (exp(logit_p)+1)
  # y <- rbinom(n = times, size = 1, prob = p)
  # y_rep <- dbinom(y, 1, p, log = TRUE)
  return(p)
  
}
# s_mat_fit(x = c(2.7, 3, 3), mod = mle.s_mat)

#######################
# TODO: SE fitted curve
# C <- c(1, 3)
# std.er <- sqrt(t(C) %*% vcov(mle.s_mat) %*% C)
# s_mat_fit(x = 1, mod = mle.s_mat, times = 1000)
# mean(replicate(n = 1000, rbinom(n = 1, size = 1, prob = 0.3)))
# mean(replicate(n = 10, rbinom(n = 1, size = 1, prob = 0.3)))
# boot_out <- boot(data = tbl1, statistic = s_mat_fit, R = 2000, mod = mle.s_mat)
# boot.ci(boot_out)
# boot.ci(boot_out, type = "perc", index = 10)
# boot.ci(boot_out, type = "perc", index = 50)
# dim(boot_out$t)
# boot_out$t[1:5, 1:10]
# plot(boot_out, index = 200)
# boot.ci(boot_out, type = "perc", index = 200)
#######################

f_x = seq(from = min(tbl1$size_log), to = max(tbl1$size_log), length.out = 200)
f_dat = data.frame(size_log = f_x,
                   p_fit = s_mat_fit(x = f_x, mod = mle.s_mat))
# plot(x = f_dat$size_log, y = f_dat$p_fit)

# range(tbl1$size_mm)
range(tbl1$size_log)
# breaks = c(seq(from = min(tbl1$size_mm), to = (max(tbl1$size_mm) - 0.7), by = 0.7), (max(tbl1$size_mm) + 0.7))
breaks = c(seq(from = min(tbl1$size_log), to = (max(tbl1$size_log) - 0.08), by = 0.08), (max(tbl1$size_log) + 0.08))
# tbl1$bin = cut(tbl1$size_mm, breaks, include.lowest = TRUE)
tbl1$bin = cut(tbl1$size_log, breaks, include.lowest = TRUE)
s_dat = group.CI(x = size_log ~ bin, data = tbl1)
p_dat = cbind(group.CI(x = maturity ~ bin, data = tbl1),
              size_log.mean = s_dat$size_log.mean)
p_dat[p_dat$maturity.upper > 1, "maturity.upper"] = 1
p_dat[p_dat$maturity.lower < 0, "maturity.lower"] = 0

(p_fit1 = ggplot(data = p_dat) +
  geom_errorbar(aes(x = size_log.mean, ymin = maturity.lower, ymax = maturity.upper, col = "observations"),
                width = 0.05) +
  geom_line(data = f_dat, aes(size_log, p_fit, col = "predictions"), size = 1.5) +
  geom_point(aes(x = size_log.mean, y = maturity.mean, col="observations"), size = 2.5) +
  # geom_point(data = CZ_data_bin, aes(x = mean_ratio, y = mount), col='blue', size=3.5) +
  scale_colour_manual(values=c("blue", "orange")) +
  labs(x = "ln size",
       y = "probability of maturity", col = "") +
  # scale_x_continuous(breaks = seq(-1.5,1.5,0.5)) +
  theme(legend.text = element_text(size = 18), legend.position = 'top',
        axis.title = element_text(size = 18),
        # axis.ticks = element_line(size = 2),
        axis.text = element_text(size = 13),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(size = 0.2, linetype = "solid",
                                 colour = "black"),
        panel.grid = element_line(colour = "gray70", size = 0.2)))
ggsave(file=paste0(dirname(dat_dir), "/", island, "_off_final_figures/", island, "_off_size_maturity_fit.svg"),
       plot=p_fit1, width=10, height=8)
ggsave(file=paste0(dirname(dat_dir), "/", island, "_off_final_figures/", island, "_off_size_maturity_fit.pdf"),
       plot=p_fit1, width=10, height=8)

## Model 1b
rm(list = setdiff(ls(), c("tbl1", "dat_dir", "island", "s_mat", "s_mat_fit", "mle.s_mat")))
sample_n(tbl = tbl1, size = 10)
# remove outlying samples in each sex class
# it is done across populations and not for each population separately
table(tbl1$sex)
mod1b_dt = rbindlist(lapply(names(table(tbl1$sex)), FUN = function(x) {
  m1 <- mean(tbl1[tbl1$sex==x, "size_mm"])
  var_name = tbl1[tbl1$sex==x, "size_mm"]
  outlier <- boxplot.stats(var_name)$out
  mo <- mean(outlier)
  var_name <- ifelse(var_name %in% outlier, NA, var_name)
  m2 <- mean(var_name, na.rm = T)
  cat("Mean", x, "without removing", x, "outliers:", round(m1, 2), "mm\n")
  cat("Mean", x, "if we remove", sum(is.na(var_name)), x, "outliers:", round(m2, 2), "mm\n")
  
  tbl1b = tbl1
  tbl1b[tbl1b$sex==x, "size_mm"] = var_name
  mean(tbl1[tbl1$sex==x, "size_mm"])
  mean(tbl1b[tbl1b$sex==x, "size_mm"], na.rm = TRUE)
  tbl1b = tbl1b[!is.na(tbl1b$size_mm), ]
  cat("Total mean without removing", x, "outliers:", round(mean(tbl1$size_mm), 2), "mm\n")
  cat("Total mean if we remove", sum(is.na(var_name)), x, "outliers:", round(mean(tbl1b$size_mm), 2), "mm\n")
  list_dt_nout = tbl1b[tbl1b$sex==x, ]
  return(list_dt_nout)
}))
round(mean(mod1b_dt$size_mm), 2)
table(mod1b_dt$sex)
table(mod1b_dt$maturity)

# boxplot without outliers
(p3 = ggplot(data = as.data.frame(mod1b_dt), aes(x = pop, y = size_mm, fill = pop)) +
    facet_wrap(~sex) +
    geom_boxplot(position = position_dodge(0.9), notch = FALSE) +
    labs(x = "", fill = "") +
    theme(legend.position = "null", legend.text = element_text(size = 11),
          axis.title.y = element_text(size = 14),
          axis.text = element_text(size = 11),
          panel.background = element_blank(),
          strip.background = element_rect(fill = "#91bfdb"),
          strip.text = element_text(size = 11),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.line = element_line(size = 0.2, linetype = "solid",
                                   colour = "black"),
          panel.grid = element_line(colour = "gray70", size = 0.2)))
ggsave(file=paste0(dirname(dat_dir), "/", island, "_off_final_figures/", island, "_off_boxplot_size_maturity_no_outliers.svg"),
       plot=p3, width=10, height=8)
ggsave(file=paste0(dirname(dat_dir), "/", island, "_off_final_figures/", island, "_off_boxplot_size_maturity_no_outliers.pdf"),
       plot=p3, width=10, height=8)

# fitting without outliers
rm(list = setdiff(ls(), c("mod1b_dt", "dat_dir", "island", "s_mat", "s_mat_fit", "mle.s_mat")))
tbl1 = as.data.frame(mod1b_dt)
rm(mod1b_dt)
mean(tbl1[, "size_mm"])
mean(tbl1[, "size_log"])
theta.init <- list(mean = 2, slope = 2)
mle.s_mat_nout <- mle2(s_mat, start = theta.init, data = list(size = tbl1$size_log,
                                                              mat = tbl1$maturity))
summary(mle.s_mat_nout)
# as.vector(coef(mle.s_mat)["mean"])
AIC(mle.s_mat)
AIC(mle.s_mat_nout)
confint(mle.s_mat_nout)

f_x = seq(from = min(tbl1$size_log), to = max(tbl1$size_log), length.out = 200)
f_dat = data.frame(size_log = f_x,
                   p_fit = s_mat_fit(x = f_x, mod = mle.s_mat_nout))

breaks = c(seq(from = min(tbl1$size_log), to = (max(tbl1$size_log) - 0.08), by = 0.08), (max(tbl1$size_log) + 0.08))
tbl1$bin = cut(tbl1$size_log, breaks, include.lowest = TRUE)
s_dat = group.CI(x = size_log ~ bin, data = tbl1)
p_dat = cbind(group.CI(x = maturity ~ bin, data = tbl1),
              size_log.mean = s_dat$size_log.mean)
p_dat[p_dat$maturity.upper > 1, "maturity.upper"] = 1
p_dat[p_dat$maturity.lower < 0, "maturity.lower"] = 0

(p_fit1 = ggplot(data = p_dat) +
    geom_errorbar(aes(x = size_log.mean, ymin = maturity.lower, ymax = maturity.upper, col = "observations"),
                  width = 0.05) +
    geom_line(data = f_dat, aes(size_log, p_fit, col = "predictions"), size = 1.5) +
    geom_point(aes(x = size_log.mean, y = maturity.mean, col="observations"), size = 2.5) +
    # geom_point(data = CZ_data_bin, aes(x = mean_ratio, y = mount), col='blue', size=3.5) +
    scale_colour_manual(values=c("blue", "orange")) +
    labs(x = "ln size",
         y = "probability of maturity", col = "") +
    # scale_x_continuous(breaks = seq(-1.5,1.5,0.5)) +
    theme(legend.text = element_text(size = 18), legend.position = 'top',
          axis.title = element_text(size = 18),
          # axis.ticks = element_line(size = 2),
          axis.text = element_text(size = 13),
          panel.background = element_rect(fill = "white"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.line = element_line(size = 0.2, linetype = "solid",
                                   colour = "black"),
          panel.grid = element_line(colour = "gray70", size = 0.2)))
ggsave(file=paste0(dirname(dat_dir), "/", island, "_off_final_figures/", island, "_off_size_maturity_fit_no_outliers.svg"),
       plot=p_fit1, width=10, height=8)
ggsave(file=paste0(dirname(dat_dir), "/", island, "_off_final_figures/", island, "_off_size_maturity_fit_no_outliers.pdf"),
       plot=p_fit1, width=10, height=8)

####################
# TODO: cline fit
# build that into a cline fit (cline for mean size at maturity first, constant slope of the maturity – size relationship across the cline)
# ‘c_mean’ and ‘w_mean’ are the size at maturity in crab and wave
# ‘c’ and ‘w’ are the centre and width of the cline in size at maturity (this is an old script – better to make log(w) the fitted parameter)

s_mat_cline <- function(position,size,mat,c_mean,w_mean,slope, c, w){
  d <- position-c
  
  p_x <- 1/(1+exp(0-4*(d)/w))
  
  mean <- c_mean+(w_mean-c_mean)*p_x
  logit_p <- (size-mean)*slope
  p <- exp(logit_p)/(exp(logit_p)+1)
  minusll <- -sum(dbinom(mat,1,p,log=T))
  
  return(minusll)
  
}

theta.init <- list(c_mean=1,w_mean=1,slope=2,c=80,w=20)

mle.s_mat_cline <- mle2(s_mat_cline,theta.init,data=list(position=snail$line_dist,size=log(snail$length_mm),mat=snail$mature))

summary(mle.s_mat_cline)
AIC(mle.s_mat_cline)
####################

rm(list = setdiff(ls(), c("tbl1", "mle.s_mat", "mle.s_mat_nout", "dat_dir", "island", "s_mat", "s_mat_fit")))

# model 2: maturity is 0,1,2 for juvenile,female,male
# two means two slopes
tbl1 = tbl1[!is.na(tbl1$size_mm), ]
tbl1$size_log = log(tbl1$size_mm)
sample_n(tbl = tbl1, size = 10)
tbl1$maturity = 0
tbl1[tbl1$sex=="female" , "maturity"] = 1
tbl1[tbl1$sex=="male", "maturity"] = 2
table(tbl1$maturity)
names(table(tbl1$maturity))
head(tbl1)
# write.csv(x = tbl1, file = "/Users/samuelperini/Desktop/test_mod2_size_maturity.csv", quote = FALSE, row.names = FALSE)

k_classes = names(table(tbl1$maturity))[-1]
k_pars = c("mean", "slope")
(nm_pars = as.vector(outer(k_pars, k_classes, paste, sep=".")))

parnames(twoM_twoS) = nm_pars[order(nm_pars)]

theta.init = setNames(c(2,2,1,1),  nm_pars[order(nm_pars)])
mle.s_mat_2M2S <- mle2(minuslogl = twoM_twoS, start = theta.init,
                       data = list(data = tbl1, y_col = 10, x_col = 11))
summary(mle.s_mat_2M2S)
AIC(mle.s_mat_2M2S)
round(coef(mle.s_mat_2M2S), 2)
mle.s_mat_fm_2M2S = mle.s_mat_2M2S
rm(mle.s_mat_2M2S)

# model 2: two means and one slope
k_classes = names(table(tbl1$maturity))[-1]
k_pars = "mean"
(nm_pars = c("slope", as.vector(outer(k_pars, k_classes, paste, sep="."))))

parnames(twoM_oneS) = nm_pars

theta.init = setNames(c(1,2,2),  nm_pars)
mle.s_mat_2M1S <- mle2(minuslogl = twoM_oneS, start = theta.init,
                       data = list(data = tbl1, y_col = 10, x_col = 11))
summary(mle.s_mat_2M1S)
round(coef(mle.s_mat_2M1S), 2)
mle.s_mat_fm_2M1S = mle.s_mat_2M1S
rm(mle.s_mat_2M1S)





# function for the likelihood of the data given that the logit(probability of maturity) is a linear function of size
# ‘mat’ is 0,1,2 for juvenile,female,male
# ‘mean.*’ and ‘slope.*’ are the regression parameters
s_mat_sex <- function(data, mean.f, slope.f, mean.m, slope.m) {
  
  logit_p_f <- (data[, "size_log"] - mean.f)  / slope.f
  logit_p_m <- (data[, "size_log"] - mean.m)  / slope.m
  
  p_f <- 0.5 * exp(logit_p_f) / (exp(logit_p_f)+1)
  p_m <- 0.5 * exp(logit_p_m) / (exp(logit_p_m)+1)
  p_j <- 1 - p_f - p_m
  # p_j <- (1 - p_f) * sr + (1-p_m) * (1 - sr)
  
  ll <- log(p_j)
  ll[data[, "maturity"] == 1] <- log(p_f[data[, "maturity"] == 1])
  ll[data[, "maturity"] == 2] <- log(p_m[data[, "maturity"] == 2])
  
  minusll <- -sum(ll)
  return(minusll)
}
mean(tbl1[tbl1$maturity==0, "size_log"])
mean(tbl1[tbl1$maturity==1, "size_log"])
mean(tbl1[tbl1$maturity==2, "size_log"])
theta.init <- list(mean.f = 2, slope.f = 1, mean.m = 2, slope.m = 1)
mle.s_mat_fm <- mle2(minuslogl = s_mat_sex, start = theta.init,
                     data = list(data = tbl1[, c("size_log", "maturity")]))
summary(mle.s_mat_fm)
AIC(mle.s_mat_fm)

# compare AIC between mod1a, mod1b, mod2
mle_obj = ls()[grepl(pattern = "mle.", x = ls())]
which.min(sapply(mle_obj, function(x) AIC(get(x))))

# plot fitted curves
round(coef(mle.s_mat_fm), 2)
s_mat_sex_fit <- function(var_x, mod) {
  coef_mod <- coef(mod)
  
  logit_p_f <- (var_x - coef_mod[1])  / coef_mod[2]
  logit_p_m <- (var_x - coef_mod[3])  / coef_mod[4]
  
  p_f <- data.frame(p_fit = exp(logit_p_f) / (exp(logit_p_f)+1),
                    maturity = "female")
  p_m <- data.frame(p_fit = exp(logit_p_m) / (exp(logit_p_m)+1),
                    maturity = "male")
  
  out_dt = rbind(p_f, p_m)
  return(out_dt)
}
f_x = seq(from = min(tbl1$size_log), to = max(tbl1$size_log), length.out = 200)
f_dat = cbind(size_log = rep(f_x, 2), s_mat_sex_fit(var_x = f_x, mod = mle.s_mat_fm))
# plot(x = f_dat$size_log, y = f_dat$p_fit)
sample_n(tbl = f_dat, size = 10)
with(plot(size_log, p_fit, col = maturity), data = f_dat)

breaks = c(seq(from = min(tbl1$size_log), to = (max(tbl1$size_log) - 0.08), by = 0.08), (max(tbl1$size_log) + 0.08))
tbl1$bin = cut(tbl1$size_log, breaks, include.lowest = TRUE)
s_dat = rbind(cbind(group.CI(x = size_log ~ bin, data = tbl1[tbl1$maturity!=2, ]), maturity = "female"),
              cbind(group.CI(x = size_log ~ bin, data = tbl1[tbl1$maturity!=1, ]), maturity = "male"))
p_dat = rbind(cbind(group.CI(x = maturity ~ bin, data = tbl1[tbl1$maturity!=2, ]),
                    size_log.mean = s_dat[s_dat$maturity=="female", "size_log.mean"], maturity = "female"),
              cbind(group.CI(x = maturity ~ bin, data = tbl1[tbl1$maturity!=1, ]),
                    size_log.mean = s_dat[s_dat$maturity=="male", "size_log.mean"], maturity = "male"))
p_dat[p_dat$maturity=="male", 2:4] = p_dat[p_dat$maturity=="male", 2:4] / 2
p_dat[p_dat$maturity.upper > 1, "maturity.upper"] = 1
p_dat[p_dat$maturity.lower < 0, "maturity.lower"] = 0
sample_n(tbl = p_dat, size = 10)
sample_n(tbl = f_dat, size = 10)

(p_fit2 = ggplot(data = p_dat, aes(col = maturity)) +
    # geom_errorbar(aes(x = size_log.mean, ymin = maturity.lower, ymax = maturity.upper), width = 0.05) +
    geom_line(data = f_dat, aes(size_log, p_fit), size = 1.5) +
    # geom_point(aes(x = size_log.mean, y = maturity.mean), size = 2.5) +
    scale_colour_manual(values=c("black", "red")) +
    labs(x = "ln size",
         y = "probability of maturity", col = "") +
    # scale_x_continuous(breaks = seq(-1.5,1.5,0.5)) +
    theme(legend.text = element_text(size = 18), legend.position = 'top',
          axis.title = element_text(size = 18),
          # axis.ticks = element_line(size = 2),
          axis.text = element_text(size = 13),
          panel.background = element_rect(fill = "white"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.line = element_line(size = 0.2, linetype = "solid",
                                   colour = "black"),
          panel.grid = element_line(colour = "gray70", size = 0.2)))
ggsave(file=paste0(dirname(dat_dir), "/", island, "_off_final_figures/", island, "_off_size_maturity_fm_fit.svg"),
       plot=p_fit2, width=10, height=8)
ggsave(file=paste0(dirname(dat_dir), "/", island, "_off_final_figures/", island, "_off_size_maturity_fm_fit.pdf"),
       plot=p_fit2, width=10, height=8)


# Model 3: four parameters and two generations without outliers
# The same logic as model 2, three maturity classes (juvenile, parents and offspring) and
# four parameters
rm(list = setdiff(ls(), c("tbl1", "mle.s_mat", "mle.s_mat_nout", "mle.s_mat_fm", "dat_dir",
                          "island", "s_mat_sex", "s_mat_sex_fit")))

# maturity is 0,1,2 for juvenile, adult parents and adult offspring
tbl1 = tbl1[!is.na(tbl1$size_mm), ]
tbl1$size_log = log(tbl1$size_mm)
sample_n(tbl = tbl1, size = 10)
tbl1$maturity = 1
tbl1[tbl1$sex!="female" & tbl1$sex!="male", "maturity"] = 0
tbl1$ja_maturity = tbl1$maturity

table(tbl1$ja_maturity)
tbl1$maturity = 0
tbl1[tbl1$generation==0 & tbl1$ja_maturity==1, "maturity"] = 1
tbl1[tbl1$generation==1 & tbl1$ja_maturity==1, "maturity"] = 2
table(tbl1$maturity)
names(table(tbl1$maturity))
head(tbl1)

# s_mat_sex

mean(tbl1[tbl1$maturity==0, "size_log"])
mean(tbl1[tbl1$maturity==1, "size_log"])
mean(tbl1[tbl1$maturity==2, "size_log"])
theta.init <- list(mean.f = 2, slope.f = 1, mean.m = 2, slope.m = 1)
mle.s_mat_po <- mle2(minuslogl = s_mat_sex, start = theta.init,
                     data = list(data = tbl1[, c("size_log", "maturity")]))
summary(mle.s_mat_po)
AIC(mle.s_mat_po)
round(coef(mle.s_mat_po), 2)

# compare AIC
mle_obj = ls()[grepl(pattern = "mle.", x = ls())]
which.min(sapply(mle_obj, function(x) AIC(get(x))))

# plot fitted curves
s_mat_sex_fit <- function(var_x, mod, classes) {
  coef_mod <- coef(mod)
  
  logit_p_1 <- (var_x - coef_mod[1])  / coef_mod[2]
  logit_p_2 <- (var_x - coef_mod[3])  / coef_mod[4]
  
  p_1 <- data.frame(p_fit = exp(logit_p_1) / (exp(logit_p_1)+1),
                    maturity = classes[1])
  p_2 <- data.frame(p_fit = exp(logit_p_2) / (exp(logit_p_2)+1),
                    maturity = classes[2])
  
  out_dt = rbind(p_1, p_2)
  return(out_dt)
}
f_x = seq(from = min(tbl1$size_log), to = max(tbl1$size_log), length.out = 200)
f_dat = cbind(size_log = rep(f_x, 2), s_mat_sex_fit(var_x = f_x, mod = mle.s_mat_po, classes = c("gen0", "gen1")))
sample_n(tbl = f_dat, size = 10)
with(plot(size_log, p_fit, col = maturity), data = f_dat)

m_class = c("gen0", "gen1")
breaks = c(seq(from = min(tbl1$size_log), to = (max(tbl1$size_log) - 0.08), by = 0.08), (max(tbl1$size_log) + 0.08))
tbl1$bin = cut(tbl1$size_log, breaks, include.lowest = TRUE)
s_dat = rbind(cbind(group.CI(x = size_log ~ bin, data = tbl1[tbl1$maturity!=2, ]), maturity = m_class[1]),
              cbind(group.CI(x = size_log ~ bin, data = tbl1[tbl1$maturity!=1, ]), maturity = m_class[2]))
p_dat = rbind(cbind(group.CI(x = maturity ~ bin, data = tbl1[tbl1$maturity!=2, ]),
                    size_log.mean = s_dat[s_dat$maturity==m_class[1], "size_log.mean"], maturity = m_class[1]),
              cbind(group.CI(x = maturity ~ bin, data = tbl1[tbl1$maturity!=1, ]),
                    size_log.mean = s_dat[s_dat$maturity==m_class[2], "size_log.mean"], maturity = m_class[2]))
if ("male" %in% m_class) {
  p_dat[p_dat$maturity=="male", 2:4] = p_dat[p_dat$maturity=="male", 2:4] / 2
}

p_dat[p_dat$maturity.upper > 1, "maturity.upper"] = 1
p_dat[p_dat$maturity.lower < 0, "maturity.lower"] = 0
sample_n(tbl = p_dat, size = 10)
sample_n(tbl = f_dat, size = 10)

(p_fit2 = ggplot(data = p_dat, aes(col = maturity)) +
    # geom_errorbar(aes(x = size_log.mean, ymin = maturity.lower, ymax = maturity.upper), width = 0.05) +
    geom_line(data = f_dat, aes(size_log, p_fit), size = 1.5) +
    # geom_point(aes(x = size_log.mean, y = maturity.mean), size = 2.5) +
    scale_colour_manual(values=c("black", "red")) +
    labs(x = "ln size",
         y = "probability of maturity", col = "") +
    # scale_x_continuous(breaks = seq(-1.5,1.5,0.5)) +
    theme(legend.text = element_text(size = 18), legend.position = 'top',
          axis.title = element_text(size = 18),
          # axis.ticks = element_line(size = 2),
          axis.text = element_text(size = 13),
          panel.background = element_rect(fill = "white"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.line = element_line(size = 0.2, linetype = "solid",
                                   colour = "black"),
          panel.grid = element_line(colour = "gray70", size = 0.2)))
ggsave(file=paste0(dirname(dat_dir), "/", island, "_off_final_figures/", island, "_off_size_maturity_",
                   paste(m_class, collapse = "_"), "_fit.svg"),
       plot=p_fit2, width=10, height=8)
ggsave(file=paste0(dirname(dat_dir), "/", island, "_off_final_figures/", island, "_off_size_maturity_",
                   paste(m_class, collapse = "_"), "_fit.pdf"),
       plot=p_fit2, width=10, height=8)


# colnames(tbl1)[10]
# colnames(tbl1)[11]
# y_col = 10
# x_col = 11
# data = tbl1
# mean = c(2,2)
# slope = c(1,1)
twoM_twoS <- function(data, y_col, x_col, beta) {
  k_count <- table(data[, y_col])[-1]
  logit_p <- matrix(nrow = nrow(data), ncol = length(k_count))
  p <- matrix(nrow = nrow(data), ncol = length(k_count))
  
  for (k in 1:length(k_count)) {
    # variable x should be already transformed
    logit_p[, k] <- (data[, x_col] - beta[k]) / beta[k+2]
    p[, k] <- 0.5 * exp(logit_p[, k]) / (exp(logit_p[, k])+1)
  }
  p_0 <- 1 - rowSums(p)
  # logit_p_f <- (data[, "size_log"] - mean.f)  / slope.f
  # logit_p_m <- (data[, "size_log"] - mean.m)  / slope.m
  # 
  # p_f <- 0.5 * exp(logit_p_f) / (exp(logit_p_f)+1)
  # p_m <- 0.5 * exp(logit_p_m) / (exp(logit_p_m)+1)
  # p_j <- 1 - p_f - p_m
  # p_j <- (1 - p_f) * sr + (1-p_m) * (1 - sr)
  
  ll <- log(p_0)
  for (k in 1:length(k_count)) {
    ll[data[, y_col] == k] <- log(p[data[, y_col] == k, k])
  }
  
  minusll <- -sum(ll)
  return(minusll)
}
k_classes = names(table(tbl1$maturity))[-1]
k_pars = c("mean", "slope")
(nm_pars = as.vector(outer(k_pars, k_classes, paste, sep=".")))

parnames(twoM_twoS) = nm_pars[order(nm_pars)]

theta.init = setNames(c(2,2,1,1),  nm_pars[order(nm_pars)])
mle.s_mat_2M2S <- mle2(minuslogl = twoM_twoS, start = theta.init,
                       data = list(data = tbl1, y_col = 10, x_col = 11))
summary(mle.s_mat_2M2S)
AIC(mle.s_mat_2M2S)
AIC(mle.s_mat_po)
round(coef(mle.s_mat_2M2S), 2)
mle.s_mat_po_2M2S = mle.s_mat_2M2S
rm(mle.s_mat_2M2S)


twoM_oneS <- function(data, y_col, x_col, sr, beta) {
  k_count <- table(data[, y_col])[-1]
  logit_p <- matrix(nrow = nrow(data), ncol = length(k_count))
  p <- matrix(nrow = nrow(data), ncol = length(k_count))
  
  for (k in 1:length(k_count)) {
    # variable x should be already transformed
    logit_p[, k] <- (data[, x_col] - beta[k+1]) / beta[1]
    p[, k] <- sr * exp(logit_p[, k]) / (exp(logit_p[, k])+1)
  }
  p_0 <- 1 - rowSums(p)
  # logit_p_f <- (data[, "size_log"] - mean.f)  / slope.f
  # logit_p_m <- (data[, "size_log"] - mean.m)  / slope.m
  # 
  # p_f <- 0.5 * exp(logit_p_f) / (exp(logit_p_f)+1)
  # p_m <- 0.5 * exp(logit_p_m) / (exp(logit_p_m)+1)
  # p_j <- 1 - p_f - p_m
  # p_j <- (1 - p_f) * sr + (1-p_m) * (1 - sr)
  
  ll <- log(p_0)
  for (k in 1:length(k_count)) {
    ll[data[, y_col] == k] <- log(p[data[, y_col] == k, k])
  }
  
  minusll <- -sum(ll)
  return(minusll)
}
k_classes = names(table(tbl1$maturity))[-1]
k_pars = "mean"
(nm_pars = c("slope", as.vector(outer(k_pars, k_classes, paste, sep="."))))

parnames(twoM_oneS) = nm_pars

theta.init = setNames(c(1,2,2),  nm_pars)
mle.s_mat_2M1S <- mle2(minuslogl = twoM_oneS, start = theta.init,
                       data = list(data = tbl1, y_col = 10, x_col = 11))
summary(mle.s_mat_2M1S)
AIC(mle.s_mat_2M2S)
AIC(mle.s_mat_2M1S)
round(coef(mle.s_mat_2M1S), 2)
mle.s_mat_po_2M1S = mle.s_mat_2M1S
rm(mle.s_mat_2M1S)
tbl1$po_maturity = tbl1$maturity
table(tbl1$po_maturity)

(mle_obj = ls()[grepl(pattern = "mle.s", x = ls())])
aic_mat = sapply(mle_obj, function(x) AIC(get(x)))
data.frame(AIC = aic_mat[order(aic_mat)])
# which.min(aic_mat)
(fun_obj = ls()[grepl(pattern = "twoM", x = ls())])
rm(list = setdiff(ls(), c("tbl1", mle_obj, "dat_dir", "island", fun_obj)))


## MODEL 4: 18 parameters, mean for each population
sample_n(tbl = tbl1, size = 10)
table(tbl1$ja_maturity)
table(tbl1$generation)

# for each population, define juvenile as 0 and adult as 1 in population A,
# 2 in population B and so on
table(tbl1$sex)
mat_by_pop <- function(data, pop, y_col, sex_col) {
  col_pop <- grepl(pattern = "pop", x = colnames(data))
  mat_class <- which(LETTERS == pop)
  one_pop <- data[data[, col_pop] == pop, ]
  one_pop[, y_col] = mat_class
  one_pop[one_pop[, sex_col] != "female" & one_pop[, sex_col] != "male", y_col] = 0
  return(one_pop)
}
colnames(tbl1)[10]
colnames(tbl1)[9]
tbl1 = as.data.frame(rbindlist(lapply(LETTERS[1:18], FUN = function(x) {
  mat_by_pop(data = tbl1, pop = x, y_col = 10, sex_col = 9)
})))
table(tbl1$maturity)
head(tbl1)

mean(tbl1[tbl1$maturity==0, "size_log"])
mean(tbl1[tbl1$maturity==1, "size_log"])
mean(tbl1[tbl1$maturity==2, "size_log"])
mean(tbl1[tbl1$maturity==3, "size_log"])

k_classes = names(table(tbl1$maturity))[-1]
k_pars = "mean"
(nm_pars = c("slope", as.vector(outer(k_pars, k_classes, paste, sep="."))))

parnames(twoM_oneS) = nm_pars

theta.init = setNames(c(1, rep(2, (length(nm_pars)-1))),  nm_pars)
colnames(tbl1)[10]
colnames(tbl1)[11]
mle.s_mat_2M1S <- mle2(minuslogl = twoM_oneS, start = theta.init,
                       data = list(data = tbl1, y_col = 10, x_col = 11, sr = 0.05))
summary(mle.s_mat_2M1S)
AIC(mle.s_mat_2M1S)
round(coef(mle.s_mat_2M1S), 2)

table(tbl1[tbl1$pop=="K", "maturity"])
mle2(minuslogl = twoM_oneS, start = theta.init,
     data = list(data = tbl1, y_col = 10, x_col = 11, sr = 0.05))

mle.s_mat_po_2M1S = mle.s_mat_2M1S
rm(mle.s_mat_2M1S)
tbl1$po_maturity = tbl1$maturity
table(tbl1$po_maturity)

(mle_obj = ls()[grepl(pattern = "mle.s", x = ls())])
aic_mat = sapply(mle_obj, function(x) AIC(get(x)))
data.frame(AIC = aic_mat[order(aic_mat)])
# which.min(aic_mat)
(fun_obj = ls()[grepl(pattern = "twoM", x = ls())])
rm(list = setdiff(ls(), c("tbl1", mle_obj, "dat_dir", "island", fun_obj)))






