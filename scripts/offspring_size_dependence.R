rm(list = ls())

.packages = c("optparse", "dplyr", "tidyr", "knitr", "kableExtra", "ggrepel",
              "flextable", "officer", "ggcorrplot", "data.table", "car", "MASS", "ggplot2")
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

# define which phenotype to analyse together with size
one_phen = c(off_phen[3], "size_mm")
genx = c("0", "1")
# scale the data
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

# lapply(seq_along(genx), function(x) hist(scaled_phen[[x]][, col_phen[1]]))
# lapply(seq_along(genx), function(x) hist(scaled_phen[[x]][, col_phen[2]]))
# lapply(seq_along(genx), function(x) qqp(scaled_phen[[x]][, col_phen[1]], "norm"))
# lapply(seq_along(genx), function(x) qqp(scaled_phen[[x]][, col_phen[2]], "norm"))
# lapply(seq_along(genx), function(x) qqp(scaled_phen[[x]][, col_phen[1]], "lnorm"))
# lapply(seq_along(genx), function(x) qqp(scaled_phen[[x]][, col_phen[2]], "lnorm"))

## model fitting
# prepare data for maturity models
scaled_phen[[1]]$generation = 0
scaled_phen[[2]]$generation = 1
lapply(seq_along(genx), function(x) {
  table(scaled_phen[[x]]$sex)
})
mod_dat = rbindlist(lapply(seq_along(genx), function(x) {
  table(scaled_phen[[x]]$sex)
  no_miss = scaled_phen[[x]][scaled_phen[[x]]$sex!="missing", ]
  no_miss$maturity = "juvenile"
  no_miss[no_miss$sex=="female" | no_miss$sex=="male", "maturity"] = "adult"
  return(no_miss)
}))
table(mod_dat[mod_dat$generation==0, "maturity"])
table(mod_dat[mod_dat$generation==1, "maturity"])
rm(list = c("tbl1", "scaled_phen"))
mod_dat = as.data.frame(mod_dat)
# change the next variable with either 0 or 1
one_gen = 1
one_dat = mod_dat[mod_dat$generation==one_gen, ]
# reminder about which phenotypes we are analysing
colnames(mod_dat)[col_phen]

basic.lm = lm(scaled_weight_cuberoot ~ scaled_size_mm, data = one_dat)

# plot(basic.lm[[1]], which = 1)
# plot(basic.lm[[2]], which = 1)
# plot(basic.lm[[1]], which = 2)
# plot(basic.lm[[2]], which = 2)
# lapply(basic.lm, summary)

bypop.lm = lm(scaled_weight_cuberoot ~ pop, data = one_dat)

add_pop.lm = lm(scaled_weight_cuberoot ~ scaled_size_mm + pop, data = one_dat)

add_mat.lm = lm(scaled_weight_cuberoot ~ scaled_size_mm + maturity, data = one_dat)

add_mat_pop.lm = lm(scaled_weight_cuberoot ~ scaled_size_mm + maturity + pop, data = one_dat)

times_pop.lm = lm(scaled_weight_cuberoot ~ scaled_size_mm * pop, data = one_dat)

times_mat.lm = lm(scaled_weight_cuberoot ~ scaled_size_mm * maturity, data = one_dat)

times_mat_pop.lm = lm(scaled_weight_cuberoot ~ scaled_size_mm * maturity * pop, data = one_dat)

aic.lm = as.data.frame(rbind(c("weight ~ size", AIC(basic.lm)),
                             c("weight ~ population", AIC(bypop.lm)),
                             c("weight ~ size + population", AIC(add_pop.lm)),
                             c("weight ~ size + maturity", AIC(add_mat.lm)),
                             c("weight ~ size + maturity + population", AIC(add_mat_pop.lm)),
                             c("weight ~ size x population", AIC(times_pop.lm)),
                             c("weight ~ size x maturity", AIC(times_mat.lm)),
                             c("weight ~ size x maturity x population", AIC(times_mat_pop.lm))))
colnames(aic.lm) = c("formula", paste0("AIC_gen", one_gen))
aic.lm
(aic.low = as.character(aic.lm[which.min(as.numeric(as.character(aic.lm[, 2]))), "formula"]))

# plot best fit (model with the lowest AIC) onto the data
best_plot <- ggplot(one_dat, aes(x = scaled_size_mm, y = scaled_weight_cuberoot)) +
  labs(y = paste0('scaled_', one_phen[1]), x = paste0('scaled_', one_phen[2]),
       title = paste0('generation ', one_gen, " - lowest AIC: ", aic.low)) +
  theme(panel.background = element_blank(),
        strip.background=element_rect(fill="#91bfdb"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(size = 0.2, linetype = "solid",
                                 colour = "black"),
        panel.grid = element_line(colour = "gray70", size = 0.2),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 11))

(fit_plot = best_plot +
    geom_smooth(method = "lm", se = TRUE, aes(col = pop, fill = pop), alpha = 0.2) +
    geom_point(aes(col = pop)))

# save figures
dir.create(paste0(dirname(dat_dir), "/", island, "_off_final_figures"))
ggsave(file=paste0(dirname(dat_dir), "/", island, "_off_final_figures/", island, "_off_best_lm_",
                   paste(one_phen, collapse = "_"), "_gen", one_gen, ".svg"),
       plot=fit_plot, width=10, height=8)
ggsave(file=paste0(dirname(dat_dir), "/", island, "_off_final_figures/", island, "_off_best_lm_",
                   paste(one_phen, collapse = "_"), "_gen", one_gen, ".pdf"),
       plot=fit_plot, width=10, height=8)

# predict weight at the overall mean size
aic.low
summary(times_pop.lm)
one_pred = data.frame(predict(times_pop.lm, data.frame(scaled_size_mm = mean(mod_dat$scaled_size_mm),
                                                       pop = unique(as.character(one_dat$pop))), interval = "confidence"),
                      pop = unique(one_dat$pop),
                      generation = one_gen)
# test for the function predict
coef(times_pop.lm)
coef(times_pop.lm)[1] + coef(times_pop.lm)[2] * mean(mod_dat$scaled_size_mm)
coef(times_pop.lm)[1] + coef(times_pop.lm)[2] * mean(mod_dat$scaled_size_mm) + coef(times_pop.lm)[3] + coef(times_pop.lm)[19] * mean(mod_dat$scaled_size_mm)
coef(times_pop.lm)[1] + coef(times_pop.lm)[2] * mean(mod_dat$scaled_size_mm) + coef(times_pop.lm)[4] + coef(times_pop.lm)[20] * mean(mod_dat$scaled_size_mm)

dir.create(paste0(dirname(dat_dir), "/", island, "_off_results"))
res_dir = paste0(dirname(dat_dir), "/", island, "_off_results")
dir.create(paste0(res_dir, "/tables"))
if (file.exists(paste0(res_dir, "/tables/", island, "_predict_mean_", one_phen[1], ".csv"))) {
  write.table(x = one_pred, file = paste0(res_dir, "/tables/", island, "_predict_mean_", one_phen[1], ".csv"),
              quote = FALSE, row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
} else {
  file.create(paste0(res_dir, "/tables/", island, "_predict_mean_", one_phen[1], ".csv"))
  write.table(x = one_pred, file = paste0(res_dir, "/tables/", island, "_predict_mean_", one_phen[1], ".csv"),
              quote = FALSE, row.names = FALSE, append = TRUE, col.names = TRUE, sep = ",")
}

# scatterplot lab vs wild predicted weight
pred_dat = read.csv(file = paste0(res_dir, "/tables/", island, "_predict_mean_", one_phen[1], ".csv"))
col_genx = c("#5e3c99", "#e66101")
ggplot(data = pred_dat[pred_dat$generation==0 & pred_dat$pop!="A", ]) +
  geom_abline(slope = 1, linetype="dashed") +
  geom_point(aes(x = fit, y = pred_dat[pred_dat$generation==1, "fit"]),
             size=3) +
  geom_errorbar(aes(x = fit,
                    ymin = pred_dat[pred_dat$generation==1, "lwr"], ymax = pred_dat[pred_dat$generation==1, "upr"]),
                size = 0.3) + 
  geom_errorbarh(aes(y = pred_dat[pred_dat$generation==1, "fit"],
                     xmin = lwr, xmax = upr),
                 size = 0.3) +
  geom_label_repel(aes(x = fit, y = pred_dat[pred_dat$generation==1, "fit"],
                       label = LETTERS[2:18]),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  labs(title = paste0("predicted ", one_phen[1]), x = "wild sample", y = "lab-reared sample") +
  theme(legend.position = "none",
        plot.title = element_text(size = 19, hjust = 0.5),
        axis.title = element_text(size = 18),
        axis.text = element_text(size=11),
        axis.ticks = element_line(size = 0.7),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(size = 0.2, linetype = "solid",
                                 colour = "black"),
        panel.grid = element_line(colour = "gray70", size = 0.2))

# boxplot wild and sample predicted weight
(pred_plot = ggplot(data = pred_dat, aes(x = pop, y = fit, col = factor(generation))) +
  geom_point(position = position_dodge(0.9), size = 3) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), size = 0.5, position = position_dodge(0.9)) +
  scale_color_manual(values = col_genx) +
  labs(x = "", y = paste0('predicted scaled ', one_phen[1]), col = "generation") +
  theme(legend.position = "top", legend.text = element_text(size = 11), legend.title = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 11),
        # axis.ticks = element_line(size = 0.5),
        panel.background = element_blank(),
        strip.background=element_rect(fill="#91bfdb"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(size = 0.2, linetype = "solid",
                                 colour = "black"),
        panel.grid = element_line(colour = "gray70", size = 0.2)))
ggsave(file=paste0(dirname(dat_dir), "/", island, "_off_final_figures/", island, "_off_boxplot_scaled_", one_phen[1], "_pred.svg"),
       plot=pred_plot, width=10, height=8)
ggsave(file=paste0(dirname(dat_dir), "/", island, "_off_final_figures/", island, "_off_boxplot_scaled_", one_phen[1], "_pred.pdf"),
       plot=pred_plot, width=10, height=8)


# summary best fit
round(summary(bypop.lm[[1]])$coefficients, 3)
names(bypop.lm[[1]]$coefficients)[2] = paste0('scaled_', one_phen[2])

## test for sampling time effect
# prepare data
mod_dat[mod_dat$generation==1, "generation"] = substr(mod_dat[mod_dat$generation==1, "ID"], start = 1, stop = 1)
table(mod_dat$generation)

## plot data
# 1 = target phenotype; 2 = size
phen_num = 2
colnames(mod_dat)[col_phen[phen_num]]
col_genx = c("#5e3c99", "#e66101")
e <- ggplot(data = mod_dat, aes(x = pop, y = mod_dat[, col_phen[phen_num]], fill = factor(generation)))
e2 <- e + geom_boxplot(position = position_dodge(0.9), fatten = NULL) +
  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y.., x = pop),
               width = 0.75, size = 1, position = position_dodge(0.9)) +
  # scale_fill_manual(values = col_genx) +
  scale_fill_manual(values = c("#5e3c99", "#e66101", "#1b9e77", "#e7298a")) +
  labs(x = "", y = paste0('scaled ', one_phen[phen_num]), fill = "generation") +
  theme(legend.position = "top", legend.text = element_text(size = 11), legend.title = element_text(size = 14),
        # plot.title = element_text(size = 19, hjust = 0.5),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 11),
        # axis.ticks = element_line(size = 0.5),
        panel.background = element_blank(),
        strip.background=element_rect(fill="#91bfdb"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(size = 0.2, linetype = "solid",
                                 colour = "black"),
        panel.grid = element_line(colour = "gray70", size = 0.2))
e2
# save plot
ggsave(file=paste0(dirname(dat_dir), "/", island, "_off_final_figures/", island, "_off_boxplot_scaled_", one_phen[1], ".svg"),
       plot=e2, width=10, height=8)
ggsave(file=paste0(dirname(dat_dir), "/", island, "_off_final_figures/", island, "_off_boxplot_scaled_", one_phen[1], ".pdf"),
       plot=e2, width=10, height=8)

