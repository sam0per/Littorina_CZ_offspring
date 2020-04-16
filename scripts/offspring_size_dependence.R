rm(list = ls())
str(tbl1)
library(car)
library(MASS)
library(lme4)
library(ggplot2)

getwd()
dat_dir = paste0(island, "_off_SW/", island, "_off_final_data/")
dat_off = read.csv(file = paste0(dat_dir, island, "_off_all_phenos_main_20200406.csv"))

dat_off = dat_off[-as.integer(rownames(dat_off[dat_off$snail_ID=='C_60', ])), ]

one_phen = c(off_phen[3], "size_mm")
genx = c("0", "1")

mean_x0 = apply(X = tbl1[tbl1$generation==0, one_phen], MARGIN = 2, FUN = function(x) mean(x, na.rm=TRUE))
sd_x0 = apply(X = tbl1[tbl1$generation==0, one_phen], MARGIN = 2, FUN = function(x) sd(x, na.rm=TRUE))

scaled_phen = lapply(genx, function(x) {
  one_dt = tbl1[tbl1$generation==x, c(one_phen, 'pop', 'ID', 'sex')]
  outna_dt = one_dt[complete.cases(one_dt), ]
  outna_dt[, paste0('scaled_', one_phen[1])] = (outna_dt[,1] - mean_x0[1]) / sd_x0[1]
  outna_dt[, paste0('scaled_', one_phen[2])] = (outna_dt[,2] - mean_x0[2]) / sd_x0[2]
  return(outna_dt)
})
col_phen = c(ncol(scaled_phen[[1]]) - 1, ncol(scaled_phen[[1]]))
colnames(scaled_phen[[1]])[col_phen]

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
plot(basic.lm[[1]], which = 1)
plot(basic.lm[[2]], which = 1)
plot(basic.lm[[1]], which = 2)
plot(basic.lm[[2]], which = 2)
lapply(basic.lm, summary)

bypop.lm = lapply(seq_along(genx), function(x) {
  lm(scaled_phen[[x]][, col_phen[1]] ~ scaled_phen[[x]][, col_phen[2]] * pop,
     data = scaled_phen[[x]])
})
lapply(bypop.lm, summary)

(prelim_plot <- lapply(seq_along(genx), function(x) {
  ggplot(scaled_phen[[x]], aes(x = scaled_phen[[x]][, col_phen[2]], y = scaled_phen[[x]][, col_phen[1]], col=pop)) +
    # facet_wrap(~pop) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE) +
    labs(y = paste0('scaled_', one_phen[1]), x = paste0('scaled_', one_phen[2]), title = paste0('generation ', genx[x]))
}))
scaled_phen[[1]][scaled_phen[[1]]$pop=='C', col_phen[1]]
scaled_phen[[1]][scaled_phen[[1]]$pop=='C', ]


lapply(seq_along(genx), function(x) {
  boxplot(scaled_phen[[x]][, col_phen[1]] ~ scaled_phen[[x]][, 'pop'], data = scaled_phen[[x]],
          xlab = "", ylab = paste0('scaled_', one_phen[1]), main = paste0('generation ', genx[x]))
})
scaled_phen[[1]]$generation = 0
scaled_phen[[2]]$generation = 1
library(data.table)
all_genx = as.data.frame(rbindlist(scaled_phen))
boxplot(all_genx[, col_phen[1]] ~ all_genx[, 'pop'] + all_genx[, 'generation'], data = all_genx,
        xlab = "", ylab = paste0('scaled_', one_phen[1]))
e <- ggplot(data = all_genx, aes(x = pop, y = all_genx[, col_phen[1]]))
e2 <- e + geom_boxplot(aes(fill = factor(generation)),position = position_dodge(0.9)) +
  scale_fill_manual(values = c("#5e3c99", "#e66101")) +
  labs(x = "", y = paste0('scaled_', one_phen[1]), fill = "generation") +
  theme(legend.position = "top",
        # plot.title = element_text(size = 19, hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size=9),
        axis.ticks = element_line(size = 0.5),
        panel.background = element_blank(),
        strip.background=element_rect(fill="#91bfdb"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(size = 0.2, linetype = "solid",
                                 colour = "black"),
        panel.grid = element_line(colour = "gray70", size = 0.2))
e2
dir.create(paste0(dirname(dat_dir), "/", island, "_off_final_figures"))
ggsave(file=paste0(dirname(dat_dir), "/", island, "_off_final_figures/", island, "_off_boxplot_scaled_", one_phen[1], ".svg"),
       plot=e2, width=10, height=8)
ggsave(file=paste0(dirname(dat_dir), "/", island, "_off_final_figures/", island, "_off_boxplot_scaled_", one_phen[1], ".pdf"),
       plot=e2, width=10, height=8)


bypop.lm = lapply(seq_along(genx), function(x) {
  lm(scaled_phen[[x]][, col_phen[1]] ~ scaled_phen[[x]][, col_phen[2]] + pop,
     data = scaled_phen[[x]])
})
lapply(bypop.lm, summary)
str(summary(bypop.lm[[1]]))
round(summary(bypop.lm[[1]])$coefficients, 3)
names(bypop.lm[[1]]$coefficients)[2] = paste0('scaled_', one_phen[2])

