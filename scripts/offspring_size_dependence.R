rm(list = ls())
str(tbl1)
library(car)
library(MASS)
library(lme4)
library(ggplot2)
one_phen = c(off_phen[3], "size_mm")
genx = c("0", "1")

mean_x0 = apply(X = tbl1[tbl1$generation==0, one_phen], MARGIN = 2, FUN = function(x) mean(x, na.rm=TRUE))
sd_x0 = apply(X = tbl1[tbl1$generation==0, one_phen], MARGIN = 2, FUN = function(x) sd(x, na.rm=TRUE))

head(tbl1[tbl1$generation==0, one_phen])
length(complete.cases(tbl1[tbl1$generation==0, one_phen]))

scaled_phen = lapply(genx, function(x) {
  one_dt = tbl1[tbl1$generation==x, c(one_phen, 'pop', 'ID', 'sex')]
  outna_dt = one_dt[complete.cases(one_dt), ]
  outna_dt[, paste0('scaled_', one_phen[1])] = (outna_dt[,1] - mean_x0[1]) / sd_x0[1]
  outna_dt[, paste0('scaled_', one_phen[2])] = (outna_dt[,2] - mean_x0[2]) / sd_x0[2]
  return(outna_dt)
})
col_phen = c(ncol(scaled_phen[[1]]) - 1, ncol(scaled_phen[[1]]))

lapply(seq_along(genx), function(x) hist(scaled_phen[[x]][, col_phen[1]]))
lapply(seq_along(genx), function(x) hist(scaled_phen[[x]][, col_phen[2]]))
lapply(seq_along(genx), function(x) qqp(scaled_phen[[x]][, col_phen[1]], "norm"))
lapply(seq_along(genx), function(x) qqp(scaled_phen[[x]][, col_phen[2]], "norm"))
lapply(seq_along(genx), function(x) qqp(scaled_phen[[x]][, col_phen[1]], "lnorm"))
lapply(seq_along(genx), function(x) qqp(scaled_phen[[x]][, col_phen[2]], "lnorm"))

basic.lm = lapply(seq_along(genx), function(x) {
  lm(scaled_phen[[x]][, col_phen[1]] ~ scaled_phen[[x]][, col_phen[2]],
     data = scaled_phen[[x]])
})
lapply(basic.lm, summary)

(prelim_plot <- lapply(seq_along(genx), function(x) {
  ggplot(scaled_phen[[x]], aes(x = scaled_phen[[x]][, col_phen[2]], y = scaled_phen[[x]][, col_phen[1]])) +
    geom_point() +
    geom_smooth(method = "lm")
}))
plot(basic.lm[[1]], which = 1)
plot(basic.lm[[2]], which = 1)
plot(basic.lm[[1]], which = 2)
plot(basic.lm[[2]], which = 2)

lapply(seq_along(genx), function(x) {
  boxplot(scaled_phen[[x]][, col_phen[1]] ~ scaled_phen[[x]][, 'pop'], data = scaled_phen[[x]])
})
