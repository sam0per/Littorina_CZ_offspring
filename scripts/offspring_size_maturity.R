rm(list = ls())

.packages = c("optparse", "dplyr", "tidyr", "knitr", "kableExtra", "ggrepel", "bbmle", "Rmisc",
              "flextable", "officer", "ggcorrplot", "data.table", "car", "MASS", "ggplot2", "boot")
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
# remove weight type in CZA
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

tbl1 = tbl1[!is.na(tbl1$size_mm), ]

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
log(mean(tbl1[, "size_mm"]))
theta.init <- list(mean = 2, slope = 2)

mle.s_mat <- mle2(s_mat, start = theta.init, data = list(size = log(tbl1$size_mm),
                                                         mat = tbl1$maturity))
summary(mle.s_mat)
as.vector(coef(mle.s_mat)["mean"])
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
s_mat_fit(x = c(2.7, 3, 3), mod = mle.s_mat)
C <- c(1, 3)
std.er <- sqrt(t(C) %*% vcov(mle.s_mat) %*% C)
# s_mat_fit(x = 1, mod = mle.s_mat, times = 1000)

mean(replicate(n = 1000, rbinom(n = 1, size = 1, prob = 0.3)))
mean(replicate(n = 10, rbinom(n = 1, size = 1, prob = 0.3)))
# boot_out <- boot(data = tbl1, statistic = s_mat_fit, R = 2000, mod = mle.s_mat)
# boot.ci(boot_out)
# boot.ci(boot_out, type = "perc", index = 10)
# boot.ci(boot_out, type = "perc", index = 50)
# dim(boot_out$t)
# boot_out$t[1:5, 1:10]
# plot(boot_out, index = 200)
# boot.ci(boot_out, type = "perc", index = 200)

f_x = seq(from = log(min(tbl1$size_mm)), to = log(max(tbl1$size_mm)), length.out = 200)
f_dat = data.frame(size_log = f_x,
                   p_fit = s_mat_fit(x = f_x, mod = mle.s_mat))
# plot(x = f_dat$size_log, y = f_dat$p_fit)

# range(tbl1$size_mm)
tbl1$size_log = log(tbl1$size_mm)
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


# model 2: maturity is 0,1,2 for juvenile,female,male
tbl1 = tbl1[!is.na(tbl1$size_mm), ]
tbl1$log_size = log(tbl1$size_mm)
tbl1$maturity = 0
tbl1[tbl1$sex=="female" , "maturity"] = 1
tbl1[tbl1$sex=="male", "maturity"] = 2
table(tbl1$maturity)
names(table(tbl1$maturity))

head(tbl1)
# write.csv(x = tbl1, file = "/Users/samuelperini/Desktop/test_mod2_size_maturity.csv", quote = FALSE, row.names = FALSE)

# function for the likelihood of the data given that the logit(probability of maturity) is a linear function of size
# ‘mat’ is 0,1,2 for juvenile,female,male
# ‘mean’ and ‘slope’ are the regression parameters
s_mat_sex <- function(data, y_col, x_col, mean, slope) {
  k_count <- table(data[, y_col])
  logit_p <- vector(mode = "list", length = length(k_count))
  p <- vector(mode = "list", length = length(k_count))
  minusll <- vector(mode = "list", length = length(k_count))
  k=1
  for (k in seq_along(names(k_count))) {
    # variable x should be already transformed
    x <- data[data[, y_col] == names(k_count)[k], x_col]
    y <- data[data[, y_col] == names(k_count)[k], y_col]
    logit_p[[k]] <- (x - mean[k]) / slope[k]
    p[[k]] <- exp(logit_p[[k]])/(exp(logit_p[[k]])+1)
    minusll[[k]] <- -sum(dbinom(y, 1, p[[k]], log = TRUE))
  }
  # minusll <- -sum(dmultinom(x = as.vector(k_count), prob = unlist(p), log = TRUE))
  return(minusll)
}

log(mean(tbl1[tbl1$sex=="male", "size_mm"]))
log(mean(tbl1[tbl1$sex=="female", "size_mm"]))
tbl1[tbl1$sex!="male",]
log(mean(tbl1[tbl1$sex!="male", "size_mm"]))
log(mean(tbl1[tbl1$sex!="female", "size_mm"]))
theta.init <- list(mean = c(2,2,2), slope = c(2,2,2))

colnames(tbl1)[10]
mle.s_mat_m <- mle2(s_mat_sex, start = theta.init,
                    data = list(data = tbl1, y_col = 11, x_col = 10))
summary(mle.s_mat_m)
AIC(mle.s_mat_m)

mle.s_mat_f <- mle2(s_mat, start = theta.init, data = list(size = log(tbl1[tbl1$sex!="male", "size_mm"]),
                                                           mat = tbl1[tbl1$sex!="male", "maturity"]))
summary(mle.s_mat_f)
AIC(mle.s_mat_f)

f_x = seq(from = log(min(tbl1$size_mm)), to = log(max(tbl1$size_mm)), length.out = 200)
f_dat = data.frame(size_log = c(f_x, f_x), sex = c(rep("female", length(f_x)), rep("male", length(f_x))),
                   p_fit = c(s_mat_fit(x = f_x, mod = mle.s_mat_f), s_mat_fit(x = f_x, mod = mle.s_mat_m)))
# plot(x = f_dat$size_log, y = f_dat$p_fit)
sample_n(tbl = f_dat, size = 10)
with(plot(size_log, p_fit, col = sex), data = f_dat)
# range(tbl1$size_mm)
tbl1$size_log = log(tbl1$size_mm)
range(tbl1$size_log)
# breaks = c(seq(from = min(tbl1$size_mm), to = (max(tbl1$size_mm) - 0.7), by = 0.7), (max(tbl1$size_mm) + 0.7))
breaks = c(seq(from = min(tbl1$size_log), to = (max(tbl1$size_log) - 0.08), by = 0.08), (max(tbl1$size_log) + 0.08))
# tbl1$bin = cut(tbl1$size_mm, breaks, include.lowest = TRUE)
tbl1$bin = cut(tbl1$size_log, breaks, include.lowest = TRUE)
s_dat = rbind(cbind(group.CI(x = size_log ~ bin, data = tbl1[tbl1$sex!="male", ]), sex = "female"),
              cbind(group.CI(x = size_log ~ bin, data = tbl1[tbl1$sex!="female", ]), sex = "male"))
p_dat = rbind(cbind(group.CI(x = maturity ~ bin, data = tbl1[tbl1$sex!="male", ]),
                    size_log.mean = s_dat[s_dat$sex=="female", "size_log.mean"], sex = "female"),
              cbind(group.CI(x = maturity ~ bin, data = tbl1[tbl1$sex!="female", ]),
                    size_log.mean = s_dat[s_dat$sex=="male", "size_log.mean"], sex = "male"))
p_dat[p_dat$maturity.upper > 1, "maturity.upper"] = 1
p_dat[p_dat$maturity.lower < 0, "maturity.lower"] = 0
sample_n(tbl = p_dat, size = 10)
sample_n(tbl = f_dat, size = 10)

(p_fit2 = ggplot(data = p_dat, aes(col = sex)) +
    geom_errorbar(aes(x = size_log.mean, ymin = maturity.lower, ymax = maturity.upper), width = 0.05) +
    geom_line(data = f_dat, aes(size_log, p_fit), size = 1.5) +
    geom_point(aes(x = size_log.mean, y = maturity.mean), size = 2.5) +
    # geom_point(data = CZ_data_bin, aes(x = mean_ratio, y = mount), col='blue', size=3.5) +
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