rm(list = ls())

.packages = c("optparse", "tidyr", "knitr", "kableExtra", "ggrepel", "RColorBrewer",
              "flextable", "officer", "ggcorrplot", "data.table", "car", "MASS", "ggplot2",
              "rstan", "shinystan", "dplyr")
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

# define which phenotype to analyse either together with size (thickness and weight) or
# independetly of size (boldness and shape)
off_phen
# one_phen = "bold_score"
one_phen = c(off_phen[3], "size_mm")
genx = c("0", "1")
# scale the data
if (length(one_phen) < 2) {
  mean_x0 = mean(tbl1[tbl1$generation==0, one_phen], na.rm=TRUE)
  sd_x0 = sd(tbl1[tbl1$generation==0, one_phen], na.rm=TRUE)
} else {
  mean_x0 = apply(X = tbl1[tbl1$generation==0, one_phen], MARGIN = 2, FUN = function(x) mean(x, na.rm=TRUE))
  sd_x0 = apply(X = tbl1[tbl1$generation==0, one_phen], MARGIN = 2, FUN = function(x) sd(x, na.rm=TRUE))
}

scaling <-  function(grp, trt) {
  out <- tryCatch(
    {
      one_dt <- tbl1[tbl1$generation==grp, c(trt, 'pop', 'ID', 'sex')]
      outna_dt = one_dt[complete.cases(one_dt), ]
      outna_dt[, paste0('scaled_', trt[1])] = (outna_dt[,1] - mean_x0[1]) / sd_x0[1]
      outna_dt[, paste0('scaled_', trt[2])] = (outna_dt[,2] - mean_x0[2]) / sd_x0[2]
      return(outna_dt)
    },
    # Handler when a warning occurs:
    warning = function(cond) {
      message("Here's the original warning message:")
      message(cond)
      
      # Choose a return value when such a type of condition occurs
      return(NULL)
    },
    # Handler when an error occurs:
    error = function(cond) {
      message("Here's the original error message:")
      message(cond)
      
      # Choose a return value when such a type of condition occurs
      outna_dt[, paste0('scaled_', trt[1])] = (outna_dt[,1] - mean_x0[1]) / sd_x0[1]
      return(outna_dt)
    },
    finally = {
      message(paste("Scaled phenotype(s) in group", grp, ":", paste(trt, collapse = " ")))
    }
  )
  return(out)
}
scaled_phen <- lapply(genx, FUN = function(x) scaling(grp = x, trt = one_phen))
lapply(scaled_phen, head)
head(tbl1)
# col_phen = c(ncol(scaled_phen[[1]]) - 1, ncol(scaled_phen[[1]]))

# lapply(seq_along(genx), function(x) hist(scaled_phen[[x]][, col_phen[1]]))
# lapply(seq_along(genx), function(x) hist(scaled_phen[[x]][, col_phen[2]]))
# lapply(seq_along(genx), function(x) qqp(scaled_phen[[x]][, col_phen[1]], "norm"))
# lapply(seq_along(genx), function(x) qqp(scaled_phen[[x]][, col_phen[2]], "norm"))
# lapply(seq_along(genx), function(x) qqp(scaled_phen[[x]][, col_phen[1]], "lnorm"))
# lapply(seq_along(genx), function(x) qqp(scaled_phen[[x]][, col_phen[2]], "lnorm"))

## model fitting
# prepare data
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
head(mod_dat)

# T_dat = mod_dat
# rm(list=setdiff(ls(), "T_dat"))
# W_dat = mod_dat
# head(T_dat)
# head(W_dat)
# table(T_dat$pop)
# table(W_dat$pop)
# rm(mod_dat)

# boxplot to check for outliers
head(mod_dat)
# out_phen <- one_phen
out_phen <- one_phen[2]
(p3 = ggplot(data = mod_dat, aes_string(x = "pop", y = out_phen, fill = "pop")) +
    facet_wrap(~maturity) +
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
# remove outlying samples in each maturity class
# it is done across populations and not for each population separately
table(mod_dat$maturity)
no_out <- rbindlist(lapply(names(table(mod_dat$maturity)), FUN = function(x) {
  m1 <- mean(mod_dat[mod_dat$maturity==x, out_phen])
  var_name <- mod_dat[mod_dat$maturity==x, out_phen]
  outlier <- boxplot.stats(var_name)$out
  mo <- mean(outlier)
  var_name <- ifelse(var_name %in% outlier, NA, var_name)
  m2 <- mean(var_name, na.rm = T)
  cat("Mean", x, out_phen, "without removing", x, "outliers:", round(m1, 2), "\n")
  cat("Mean", x, out_phen, "if we remove", sum(is.na(var_name)), x, "outliers:", round(m2, 2), "\n")
  
  tbl1b = mod_dat
  tbl1b[tbl1b$maturity==x, out_phen] = var_name
  mean(mod_dat[mod_dat$maturity==x, out_phen])
  mean(tbl1b[tbl1b$maturity==x, out_phen], na.rm = TRUE)
  tbl1b = tbl1b[!is.na(tbl1b[, out_phen]), ]
  cat("Total mean without removing", x, out_phen, "outliers:", round(mean(mod_dat[, out_phen]), 2), "\n")
  cat("Total mean if we remove", sum(is.na(var_name)), x, out_phen, "outliers:", round(mean(tbl1b[, out_phen]), 2), "\n")
  list_dt_nout = tbl1b[tbl1b$maturity==x, ]
  return(list_dt_nout)
}))
(p3_nout <- ggplot(data = no_out, aes_string(x = "pop", y = out_phen, fill = "pop")) +
    facet_wrap(~maturity) +
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

mod_dat <- no_out
## model all samples combined
col_genx = c("#5e3c99", "#e66101")
frmls = c(paste0("scaled_", one_phen[1], " ~ scaled_", one_phen[2], " * pop * generation"),
          paste0("scaled_", one_phen[1], " ~ scaled_", one_phen[2], " * pop + generation"),
          paste0("scaled_", one_phen[1], " ~ scaled_", one_phen[2], " * pop"),
          paste0("scaled_", one_phen[1], " ~ scaled_", one_phen[2], " * generation"),
          paste0("scaled_", one_phen[1], " ~ scaled_", one_phen[2], " * generation + pop"))
# times_pop_gen.lm = lm(scaled_weight_cuberoot ~ scaled_size_mm * pop * generation, data = mod_dat)
# times_pop_all.lm = lm(scaled_weight_cuberoot ~ scaled_size_mm * pop, data = mod_dat)
# times_gen_all.lm = lm(scaled_weight_cuberoot ~ scaled_size_mm * generation, data = mod_dat)
times_pop_gen.lm = lm(formula = frmls[1], data = mod_dat)
times_pop.gen.lm = lm(formula = frmls[2], data = mod_dat)
times_pop_all.lm = lm(formula = frmls[3], data = mod_dat)
times_gen_all.lm = lm(formula = frmls[4], data = mod_dat)
times_gen.pop.lm = lm(formula = frmls[5], data = mod_dat)
mod.lm = ls()[grepl(pattern = ".lm$", x = ls())]
aic.lm = sapply(mod.lm, function(x) AIC(get(x)))
aic.lm[order(aic.lm)]
(aic.low = frmls[which.min(c(AIC(times_pop_gen.lm),
                             AIC(times_pop.gen.lm),
                             AIC(times_pop_all.lm),
                             AIC(times_gen_all.lm),
                             AIC(times_gen.pop.lm)))])

# plot best fit (model with the lowest AIC) onto the data
best_plot <- ggplot(mod_dat, aes_string(x = paste0('scaled_', one_phen[2]), y = paste0('scaled_', one_phen[1]))) +
  labs(y = paste0('scaled_', one_phen[1]), x = paste0('scaled_', one_phen[2]),
       title = paste0("All samples combined - lowest AIC: ", aic.low)) +
  theme(panel.background = element_blank(),
        strip.background=element_rect(fill="#91bfdb"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(size = 0.2, linetype = "solid",
                                 colour = "black"),
        panel.grid = element_line(colour = "gray70", size = 0.2),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 13),
        legend.text = element_text(size = 14), legend.title = element_text(size = 14),
        title = element_text(size = 12))
(fit_plot = best_plot +
    geom_smooth(method = "lm", se = TRUE, aes(col = factor(generation), fill = factor(generation)), alpha = 0.4) +
    geom_point(aes(col = factor(generation)), alpha = 0.6) +
    scale_fill_manual(values = col_genx) +
    scale_color_manual(values = col_genx) +
    labs(col = "generation", fill = "generation"))
# save figures
dir.create(paste0(dirname(dat_dir), "/", island, "_off_final_figures"))
ggsave(file=paste0(dirname(dat_dir), "/", island, "_off_final_figures/", island, "_off_best_lm_",
                   paste(one_phen, collapse = "_"), "_all_samples.svg"),
       plot=fit_plot, width=10, height=8)
ggsave(file=paste0(dirname(dat_dir), "/", island, "_off_final_figures/", island, "_off_best_lm_",
                   paste(one_phen, collapse = "_"), "_all_samples.pdf"),
       plot=fit_plot, width=10, height=8)

## model with separate generations
rm(list = setdiff(ls(), c("mod_dat", "dat_dir", "genx", "island", "one_phen", "res_dir", "col_genx", "col_phen")))
# change the next variable with either 0 or 1
one_gen = 0
one_dat = mod_dat[mod_dat$generation==one_gen, ]
head(one_dat)
# reminder about which phenotypes we are analysing
colnames(mod_dat)[col_phen]

frmls = c(paste0("scaled_", one_phen[1], " ~ scaled_", one_phen[2]),
          paste0("scaled_", one_phen[1], " ~ pop"),
          paste0("scaled_", one_phen[1], " ~ scaled_", one_phen[2], " + pop"),
          paste0("scaled_", one_phen[1], " ~ scaled_", one_phen[2], " + maturity"),
          paste0("scaled_", one_phen[1], " ~ scaled_", one_phen[2], " + maturity + pop"),
          paste0("scaled_", one_phen[1], " ~ scaled_", one_phen[2], " * pop"),
          paste0("scaled_", one_phen[1], " ~ scaled_", one_phen[2], " * maturity"),
          paste0("scaled_", one_phen[1], " ~ scaled_", one_phen[2], " * maturity * pop"))

basic.lm = lm(formula = frmls[1], data = one_dat)

# plot(basic.lm[[1]], which = 1)
# plot(basic.lm[[2]], which = 1)
# plot(basic.lm[[1]], which = 2)
# plot(basic.lm[[2]], which = 2)
# lapply(basic.lm, summary)

bypop.lm = lm(formula = frmls[2], data = one_dat)

add_pop.lm = lm(formula = frmls[3], data = one_dat)

add_mat.lm = lm(formula = frmls[4], data = one_dat)

add_mat_pop.lm = lm(formula = frmls[5], data = one_dat)

times_pop.lm = lm(formula = frmls[6], data = one_dat)

times_mat.lm = lm(formula = frmls[7], data = one_dat)

times_mat_pop.lm = lm(formula = frmls[8], data = one_dat)

aic.lm = as.data.frame(rbind(c(frmls[1], AIC(basic.lm)),
                             c(frmls[2], AIC(bypop.lm)),
                             c(frmls[3], AIC(add_pop.lm)),
                             c(frmls[4], AIC(add_mat.lm)),
                             c(frmls[5], AIC(add_mat_pop.lm)),
                             c(frmls[6], AIC(times_pop.lm)),
                             c(frmls[7], AIC(times_mat.lm)),
                             c(frmls[8], AIC(times_mat_pop.lm))))
colnames(aic.lm) = c("formula", paste0("AIC_gen", one_gen))
aic.lm
(aic.low = as.character(aic.lm[which.min(as.numeric(as.character(aic.lm[, 2]))), "formula"]))
aic.lm[, paste0("AIC_gen", one_gen)] = round(as.numeric(as.character(aic.lm[, paste0("AIC_gen", one_gen)])), 2)
aic.lm[order(aic.lm[, paste0("AIC_gen", one_gen)]),]

col_mat = c("#1b9e77", "#e7298a")
# define which variable to use for plot colours (e.g., "pop" or "maturity")
col_mod = "maturity"
# plot best fit (model with the lowest AIC) onto the data
best_plot <- ggplot(one_dat, aes_string(x = paste0('scaled_', one_phen[2]), y = paste0('scaled_', one_phen[1]))) +
  labs(y = paste0('scaled_', one_phen[1]), x = paste0('scaled_', one_phen[2]),
       title = paste0('generation ', one_gen, " - lowest AIC: ", aic.low)) +
  theme(panel.background = element_blank(),
        strip.background=element_rect(fill="#91bfdb"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(size = 0.2, linetype = "solid",
                                 colour = "black"),
        panel.grid = element_line(colour = "gray70", size = 0.2),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 9),
        legend.text = element_text(size = 14), legend.title = element_text(size = 14),
        title = element_text(size = 13))

(fit_plot = if (col_mod=="maturity") {
  best_plot +
    # only if col_mod is "maturity"
    facet_wrap(~pop) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
    geom_point(aes_string(col = col_mod), alpha = 0.5) +
    # only if col_mod is "maturity"
    scale_color_manual(values = col_mat) + scale_fill_manual(values = col_mat)
} else {
  best_plot +
    facet_wrap(~pop) +
    geom_smooth(method = "lm", se = TRUE, aes_string(col = col_mod, fill = col_mod), alpha = 0.2) +
    geom_point(aes_string(col = col_mod), alpha = 0.5)
})
# summary(add_mat_pop.lm)
# summary(times_mat_pop.lm)
# save figures
dir.create(paste0(dirname(dat_dir), "/", island, "_off_final_figures"))
ggsave(file=paste0(dirname(dat_dir), "/", island, "_off_final_figures/", island, "_off_best_lm_",
                   paste(one_phen, collapse = "_"), "_gen", one_gen, ".svg"),
       plot=fit_plot, width=10, height=8)
ggsave(file=paste0(dirname(dat_dir), "/", island, "_off_final_figures/", island, "_off_best_lm_",
                   paste(one_phen, collapse = "_"), "_gen", one_gen, ".pdf"),
       plot=fit_plot, width=10, height=8)

# predict phenotype at the overall mean size
aic.low
head(mod_dat)
ggplot(data = mod_dat) +
  facet_wrap(~pop) +
  geom_point(aes_string(x = one_phen[2], y = one_phen[1], col=factor(mod_dat$generation)),
             size=1) +
  scale_colour_manual(values = c("#5e3c99", "#e66101")) +
  labs(col = "generation") +
  theme(legend.position = "top", legend.text = element_text(size = 12), legend.title = element_text(size = 14),
        # plot.title = element_text(size = 19, hjust = 0.5),
        axis.title = element_text(size = 14),
        axis.text = element_text(size=10),
        axis.ticks = element_line(size = 0.5),
        panel.background = element_blank(),
        strip.background=element_rect(fill="#91bfdb"), strip.text = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(size = 0.2, linetype = "solid",
                                 colour = "black"),
        panel.grid = element_line(colour = "gray70", size = 0.2))
pred_dat <- data.frame(predict(times_pop_gen.lm, data.frame(scaled_size_mm = mean(mod_dat$scaled_size_mm),
                                                            pop = mod_dat$pop,
                                                            generation = mod_dat$generation), se.fit = TRUE)[c("fit", "se.fit")],
                       pop = mod_dat$pop,
                       generation = mod_dat$generation)
head(pred_dat)

ad_dt <- mod_dat[mod_dat$maturity=="adult", ]
pred_dat <- data.frame(predict(times_pop_gen.lm, data.frame(scaled_size_mm = mean(ad_dt$scaled_size_mm),
                                                            pop = rep(unique(as.character(ad_dt$pop)), 2),
                                                            generation = c(rep(0,18), rep(1,18))), se.fit = TRUE)[c("fit", "se.fit")],
                       pop = rep(unique(as.character(ad_dt$pop)), 2),
                       generation = c(rep(0,18), rep(1,18)))
# summary(times_pop.lm)
if (grepl(pattern = "maturity", x = aic.low)) {
  one_pred = data.frame(predict(add_mat_pop.lm, data.frame(scaled_size_mm = mean(mod_dat$scaled_size_mm),
                                                           maturity = c(rep("juvenile", length(unique(as.character(one_dat$pop)))),
                                                                        rep("adult", length(unique(as.character(one_dat$pop))))),
                                                           pop = rep(unique(as.character(one_dat$pop)), 2)), interval = "confidence"),
                        maturity = c(rep("juvenile", length(unique(one_dat$pop))), rep("adult", length(unique(as.character(one_dat$pop))))),
                        pop = rep(unique(one_dat$pop), 2),
                        generation = one_gen)
} else {
  one_pred = data.frame(predict(times_pop_gen.lm, data.frame(scaled_size_mm = mean(mod_dat$scaled_size_mm),
                                                             pop = unique(as.character(one_dat$pop))), se.fit = TRUE)[c("fit", "se.fit")],
                        N = as.vector(table(one_dat$pop)),
                        pop = unique(one_dat$pop),
                        generation = one_gen)
}
one_pred$sd.fit = one_pred$se.fit * sqrt(one_pred$N)
one_coef = data.frame(parameter = names(coef(times_pop.lm)),
                      summary(times_pop.lm)$coefficients[, 1:2],
                      generation = one_gen)
# test for the function predict
coef(times_pop.lm)
coef(times_pop.lm)[1] + coef(times_pop.lm)[2] * mean(mod_dat$scaled_size_mm)
coef(times_pop.lm)[1] + coef(times_pop.lm)[2] * mean(mod_dat$scaled_size_mm) + coef(times_pop.lm)[3] + coef(times_pop.lm)[19] * mean(mod_dat$scaled_size_mm)
coef(times_pop.lm)[1] + coef(times_pop.lm)[2] * mean(mod_dat$scaled_size_mm) + coef(times_pop.lm)[4] + coef(times_pop.lm)[20] * mean(mod_dat$scaled_size_mm)
###############################

dir.create(paste0(dirname(dat_dir), "/", island, "_off_results"))
res_dir = paste0(dirname(dat_dir), "/", island, "_off_results")
dir.create(paste0(res_dir, "/tables"))
if (file.exists(paste0(res_dir, "/tables/", island, "_off_size_adjusted_", one_phen[1], ".csv"))) {
  infl = read.csv(file = paste0(res_dir, "/tables/", island, "_off_size_adjusted_", one_phen[1], ".csv"))
  if (length(table(infl$generation)) < 2) {
    write.table(x = one_pred, file = paste0(res_dir, "/tables/", island, "_off_size_adjusted_", one_phen[1], ".csv"),
                quote = FALSE, row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
    write.table(x = one_coef,
                file = paste0(res_dir, "/tables/", island, "_off_lm_coef_", one_phen[1], ".csv"),
                quote = FALSE, row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
} else {
  file.create(paste0(res_dir, "/tables/", island, "_off_size_adjusted_", one_phen[1], ".csv"))
  write.table(x = one_pred, file = paste0(res_dir, "/tables/", island, "_off_size_adjusted_", one_phen[1], ".csv"),
              quote = FALSE, row.names = FALSE, append = TRUE, col.names = TRUE, sep = ",")
  file.create(paste0(res_dir, "/tables/", island, "_off_lm_coef_", one_phen[1], ".csv"))
  write.table(x = one_coef,
              file = paste0(res_dir, "/tables/", island, "_off_lm_coef_", one_phen[1], ".csv"),
              quote = FALSE, row.names = FALSE, append = TRUE, col.names = TRUE, sep = ",")
}

# scatterplot lab vs wild size-adjusted phenotype
aic.low
rm(list=setdiff(ls(), c("times_pop.lm", "mod_dat", "aic.low", "col_genx", "dat_dir",
                        "genx", "island", "one_phen", "res_dir")))

pred_dat = read.csv(file = paste0(res_dir, "/tables/", island, "_off_size_adjusted_", one_phen[1], ".csv"))
pred_dat <- pred_dat[pred_dat$pop!="A", ]
ggplot(data = pred_dat[pred_dat$generation==0, ]) +
  geom_abline(slope = 1, linetype = "dashed") +
  geom_point(aes(x = fit, y = pred_dat[pred_dat$generation==1, "fit"]),
             size=3) +
  geom_errorbar(aes(x = fit,
                    ymin = (pred_dat[pred_dat$generation==1, "fit"] - pred_dat[pred_dat$generation==1, "se.fit"]),
                    ymax = (pred_dat[pred_dat$generation==1, "fit"] + pred_dat[pred_dat$generation==1, "se.fit"])),
                size = 0.3, width = 0.01) + 
  geom_errorbarh(aes(y = pred_dat[pred_dat$generation==1, "fit"],
                     xmin = (fit - se.fit),
                     xmax = (fit + se.fit)),
                 size = 0.3, height = 0.005) +
  geom_label_repel(aes(x = fit, y = pred_dat[pred_dat$generation==1, "fit"],
                       label = unique(pred_dat$pop)),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  labs(title = paste0("Size-adjusted ", one_phen[1]), x = "parents", y = "lab-reared adults") +
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

# boxplot wild and sample size_adjusted phenotype
head(pred_dat)
(pred_plot = ggplot(data = pred_dat, aes(x = pop, y = fit, col = factor(generation))) +
  geom_point(position = position_dodge(0.9), size = 3) +
  geom_errorbar(aes(ymin = (fit - se.fit), ymax = (fit + se.fit)), size = 0.5, position = position_dodge(0.9)) +
  scale_color_manual(values = col_genx) +
  labs(x = "", y = paste0('size-adjusted scaled ', one_phen[1]), col = "generation") +
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
ggsave(file=paste0(dirname(dat_dir), "/", island, "_off_final_figures/", island, "_off_boxplot_adjusted_", one_phen[1], ".svg"),
       plot=pred_plot, width=10, height=8)
ggsave(file=paste0(dirname(dat_dir), "/", island, "_off_final_figures/", island, "_off_boxplot_adjusted_", one_phen[1], ".pdf"),
       plot=pred_plot, width=10, height=8)
rm(pred_plot)

## boxplot of the regression coefficients of wild and lab samples
coef_dat = read.csv(file = paste0(res_dir, "/tables/", island, "_off_lm_coef_", one_phen[1], ".csv"))
head(coef_dat)
coef_dat$variance <- coef_dat$Std..Error^2

# calculate regression coefficients for each population (i.e., reference slope + population slope)
calc_regr <- function(data, par_col, grp, beta_ref) {
  grp_name <- deparse(substitute(grp))
  data[, par_col] <- as.character(data[, par_col])
  par1 <- data[data[, par_col] == beta_ref, ]
  par_betas <- data[grepl(pattern = paste0(beta_ref, ":"), x = data[, par_col]), ]
  out_regr <- matrix(nrow = nrow(par1) + nrow(par_betas), ncol = ncol(par_betas))
  colnames(out_regr) <- colnames(par_betas)
  out_regr <- as.data.frame(out_regr)
  out_regr[, grp_name] <- c(par1[, grp_name], par_betas[, grp_name])
  out_regr[1:nrow(par1), ] <- par1
  for (i in 1:nrow(par1)) {
    a_grp <- par1[i, grp_name]
    par_grp <- par_betas[par_betas[, grp_name] == a_grp, ]
    out_regr[out_regr[, grp_name] == a_grp, ][-1, par_col] <- par_grp[, par_col]
    out_regr[out_regr[, grp_name] == a_grp, ][-1, "Estimate"] <- par1[i, "Estimate"] + par_grp[, "Estimate"]
    out_regr[out_regr[, grp_name] == a_grp, ][-1, "variance"] <- par1[i, "variance"] + par_grp[, "variance"]
  }
  return(out_regr)
}
regr_dat <- calc_regr(data = coef_dat, par_col = 1, grp = generation, beta_ref = paste0("scaled_", one_phen[2]))
regr_dat$sd <- sqrt(regr_dat$variance)
regr_dat[regr_dat$generation==0, "pop"] <- unique(mod_dat[mod_dat$generation==0, "pop"])
regr_dat[regr_dat$generation==1, "pop"] <- unique(mod_dat[mod_dat$generation==1, "pop"])
regr_dat$parameter <- NULL
write.csv(x = regr_dat, file = paste0(res_dir, "/tables/", island, "_off_slopes_", one_phen[1], ".csv"),
          quote = FALSE, row.names = FALSE)

(regr_plot = ggplot(data = regr_dat, aes(x = pop, y = Estimate, col = factor(generation))) +
    geom_point(position = position_dodge(0.9), size = 3) +
    geom_errorbar(aes(ymin = (Estimate - sd), ymax = (Estimate + sd)), size = 0.5, position = position_dodge(0.9)) +
    scale_color_manual(values = col_genx) +
    labs(x = "", y = paste0(one_phen[1], " regression coef."), col = "generation") +
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
ggsave(file=paste0(dirname(dat_dir), "/", island, "_off_final_figures/", island, "_off_boxplot_slopes_", one_phen[1], ".svg"),
       plot=regr_plot, width=10, height=8)
ggsave(file=paste0(dirname(dat_dir), "/", island, "_off_final_figures/", island, "_off_boxplot_slopes_", one_phen[1], ".pdf"),
       plot=regr_plot, width=10, height=8)
rm(regr_plot)


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


## weight variation explained by thickness
# correlation of the adjusted means and regression coefficients
# using the model with the interactions between size and populations
# (i.e., best model for weight but not for thickness)
two_phen = c(off_phen[2], off_phen[3], "size_mm")
one_gen = 0
T_dat$snail_ID = paste(T_dat$pop, T_dat$ID, sep = "_")
W_dat$snail_ID = paste(W_dat$pop, W_dat$ID, sep = "_")
rm(mod_dat)
mod_dat = merge(W_dat, T_dat)
head(mod_dat)

one_dat = mod_dat[mod_dat$generation==one_gen, ]

frmls = c(paste0("scaled_", two_phen[1], " ~ scaled_", two_phen[3], " * pop"),
          paste0("scaled_", two_phen[2], " ~ scaled_", two_phen[3], " * pop"))

T_times_pop.lm = lm(formula = frmls[1], data = one_dat)
W_times_pop.lm = lm(formula = frmls[2], data = one_dat)

summary(T_times_pop.lm)
summary(W_times_pop.lm)

# correlation of slopes
regr_idx = grepl(pattern = paste0('scaled_', two_phen[3]), x = names(coef(T_times_pop.lm)))
(regr_coef = data.frame(thick_slopes = round(coef(T_times_pop.lm)[regr_idx], 2),
                        weight_slopes = round(coef(W_times_pop.lm)[regr_idx], 2)))
round(cor(regr_coef$thick_slopes, regr_coef$weight_slopes), 2)
# correlation of slopes + intercepts
(all_coef = data.frame(thickness = round(coef(T_times_pop.lm), 2),
                       weight = round(coef(W_times_pop.lm), 2)))
round(cor(all_coef$thickness, all_coef$weight), 2)

# correlation adjusted means
one_pred = rbind(data.frame(round(predict(T_times_pop.lm,
                                          data.frame(scaled_size_mm = mean(mod_dat$scaled_size_mm),
                                                     pop = unique(as.character(one_dat$pop))), interval = "confidence"), 2),
                            pop = unique(one_dat$pop),
                            generation = one_gen,
                            phenotype = two_phen[1]),
                 data.frame(round(predict(W_times_pop.lm,
                                          data.frame(scaled_size_mm = mean(mod_dat$scaled_size_mm),
                                                     pop = unique(as.character(one_dat$pop))), interval = "confidence"), 2),
                            pop = unique(one_dat$pop),
                            generation = one_gen,
                            phenotype = two_phen[2]))
one_pred = one_pred[order(one_pred$pop), ]
round(cor(one_pred[one_pred$phenotype==two_phen[1], "fit"],
          one_pred[one_pred$phenotype==two_phen[2], "fit"]), 2)

# display.brewer.all(colorblindFriendly = TRUE)
# display.brewer.pal(n = 8, name = 'Set2')
# brewer.pal(n = 8, name = 'Set2')
# display.brewer.pal(n = 8, name = 'Paired')
# brewer.pal(n = 8, name = 'Paired')

# boxplot two phenotypes adjusted means
(pred_plot = ggplot(data = one_pred, aes(x = pop, y = fit, col = factor(phenotype))) +
    geom_point(position = position_dodge(0.9), size = 3) +
    geom_errorbar(aes(ymin = lwr, ymax = upr), size = 0.5, position = position_dodge(0.9)) +
    scale_color_manual(values = c("#1F78B4", "#B3B3B3")) +
    labs(x = "", y = "Size-adjusted means", col = "", title = paste0("Generation ", one_gen)) +
    theme(legend.position = "top", legend.text = element_text(size = 11), legend.title = element_text(size = 14),
          plot.title = element_text(size = 14, hjust = 0.5),
          axis.title.y = element_text(size = 14),
          axis.text = element_text(size = 11),
          # axis.ticks = element_line(size = 0.5),
          panel.background = element_blank(),
          strip.background=element_rect(fill="#91bfdb"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.line = element_line(size = 0.2, linetype = "solid",
                                   colour = "black"),
          panel.grid = element_line(colour = "gray70", size = 0.2)))
ggsave(file=paste0(dirname(dat_dir), "/", island, "_off_final_figures/", island, "_off_boxplot_adjusted_", paste(two_phen[-3], collapse = "_"),
                   "_gen",one_gen, ".svg"),
       plot=pred_plot, width=10, height=8)
ggsave(file=paste0(dirname(dat_dir), "/", island, "_off_final_figures/", island, "_off_boxplot_adjusted_", paste(two_phen[-3], collapse = "_"),
                   "_gen",one_gen, ".pdf"),
       plot=pred_plot, width=10, height=8)

col_mod = "pop"
# plot best fit (model with the lowest AIC) onto the data
best_plot <- ggplot(one_dat, aes_string(x = paste0('scaled_', one_phen[2]), y = paste0('scaled_', one_phen[1]))) +
  labs(y = paste0('scaled_', one_phen[1]), x = paste0('scaled_', one_phen[2]),
       title = paste0('generation ', one_gen, " - lowest AIC: ", aic.low)) +
  theme(panel.background = element_blank(),
        strip.background=element_rect(fill="#91bfdb"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(size = 0.2, linetype = "solid",
                                 colour = "black"),
        panel.grid = element_line(colour = "gray70", size = 0.2),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 9),
        legend.text = element_text(size = 14), legend.title = element_text(size = 14),
        title = element_text(size = 13))

(fit_plot = best_plot +
    # only if col_mod is "maturity"
    # facet_wrap(col_mod) +
    geom_smooth(method = "lm", se = TRUE, aes_string(col = col_mod, fill = col_mod), alpha = 0.2) +
    geom_point(aes_string(col = col_mod), alpha = 0.5))



## Plasticity analysis
# Boldness
head(no_out)
table(no_out[no_out$generation==0, "maturity"])
x_meas <- no_out[no_out$generation==0 & no_out$maturity=="adult", ]
x_meas <- aggregate(x = x_meas$scaled_bold_score, by = list(pop = x_meas$pop),
                    FUN = function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE), N = length(x)))
# sum(x_meas$pop=="P")
x_meas$se <- x_meas$x[, "sd"] / sqrt(x_meas$x[, "N"])
x_meas$generation <- 0

table(no_out[no_out$generation==1, "maturity"])
y_meas <- no_out[no_out$generation==1 & no_out$maturity=="adult", ]
y_meas <- aggregate(x = y_meas$scaled_bold_score, by = list(pop = y_meas$pop),
                    FUN = function(y) c(mean = mean(y, na.rm = TRUE), sd = sd(y, na.rm = TRUE), N = length(y)))
# sum(y_meas$pop=="P")
y_meas$se <- y_meas$x[, "sd"] / sqrt(y_meas$x[, "N"])
y_meas$generation <- 1

pop_diff <- setdiff(x_meas$pop, y_meas$pop)
if (length(pop_diff) != 0) {
  x_meas <- x_meas[x_meas$pop!=pop_diff, ]
}

res_dir = paste0(island, "_off_SW/", island, "_off_results/")
dir.create(paste0(res_dir, "models"))
dir.create(paste0(res_dir, "tables"))

write.csv(x = rbind(x_meas, y_meas), file = paste0(res_dir, "tables/", island, "_xy_meas_", one_phen, ".csv"))

stanfile = "Littorina_offspring/scripts/offspring_err_in_var_model.stan"
writeLines(readLines(stanfile))
rstan_options(auto_write = TRUE)
# options(mc.cores = 4)
options(mc.cores = parallel::detectCores(logical = FALSE) - 2)
dat = list(N = nrow(x_meas), x = x_meas$x[, "mean"], sd_x = x_meas$se,
           y = y_meas$x[, "mean"], sd_y = y_meas$se)
err_in_var <- rstan::stan(file = stanfile,
                          data = dat, iter = 12000, warmup = 4000,
                          chains=4, refresh=12000,
                          control = list(stepsize = 0.01, adapt_delta = 0.99, max_treedepth = 15))



saveRDS(err_in_var, paste0(res_dir, "models/", island, "_err_in_var_", one_phen, ".rds"))
# err_in_var = readRDS(paste0(res_dir, "models/", island, "_err_in_var.rds"))

# launch_shinystan(err_in_var)
# pairs(err_in_var)
pars = c("alpha", "beta", "sigma")
# pars = c("alpha", "beta")
pdf(file = paste0(dirname(res_dir), "/", island, "_off_final_figures/", island, "_off_err_in_var_pairs_", one_phen, ".pdf"),
    width = 8, height = 8)
pairs(err_in_var, pars = pars, include = TRUE)
dev.off()
# print(err_in_var, pars=c("alpha", "beta", "sigma"), digits=3)
print(err_in_var, pars = pars, digits=3)
# stan_pars = read.csv(file = paste0(res_dir, "tables/", island, "_err_in_var_stanfit_", one_phen, "_pars.csv"))

# postd = extract(err_in_var)
# names(postd)
# dim(postd$alpha)
stbl = rstan::summary(err_in_var)
write.csv(x = stbl$summary, file = paste0(res_dir, "tables/", island, "_err_in_var_stanfit_", one_phen, ".csv"))
# stbl = read.csv(file = paste0(res_dir, "tables/", island, "_err_in_var_stanfit_", one_phen, ".csv"))
# head(stbl)
stan_x_lat <- stbl$summary[grepl(pattern = "x_lat", x = rownames(stbl$summary)), "mean"]
new_x_lat <- seq(from = min(stan_x_lat), to = max(stan_x_lat), length.out = 1000)
new_mu_yhat <- stbl$summary["alpha", "mean"] + stbl$summary["beta", "mean"] * new_x_lat
# new_mu_yhat <- stan_pars[1, "mean"] + stan_pars[2, "mean"] * new_x_lat
# plot(x = new_x_lat, y = new_mu_yhat)
# plasticty plot
(pp <- ggplot() +
  geom_abline(slope = 1, linetype="dashed") +
  geom_point(aes(x = x_meas$x[, "mean"], y = y_meas$x[, "mean"]),
             size=3) +
  geom_label_repel(aes(x = x_meas$x[, "mean"], y = y_meas$x[, "mean"],
                       label = x_meas$pop),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  geom_point(aes(x = stbl$summary[1:nrow(x_meas), "mean"], y = stbl$summary[(nrow(x_meas)+1):(nrow(x_meas)*2), "mean"]),
             size=3, col="red") +
  geom_line(aes(x = new_x_lat, y = new_mu_yhat),
            size=2, col="red") +
  labs(title = one_phen, x = "parents", y = "lab-reared adults") +
  theme(legend.position = "none",
        plot.title = element_text(size = 19, hjust = 0.5),
        axis.title = element_text(size = 18),
        axis.text = element_text(size=11),
        axis.ticks = element_line(size = 0.7),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(size = 0.2, linetype = "solid",
                                 colour = "black"),
        panel.grid = element_line(colour = "gray70", size = 0.2)))
ggsave(file=paste0(dirname(dat_dir), "/", island, "_off_final_figures/", island, "_off_plasticity_", one_phen, ".svg"),
       plot=pp, width=10, height=8)
ggsave(file=paste0(dirname(dat_dir), "/", island, "_off_final_figures/", island, "_off_plasticity_", one_phen, ".pdf"),
       plot=pp, width=10, height=8)

# diagnose parameter sigma
y_meas$pop[order(y_meas$se)]
x_meas$pop[order(x_meas$se)]
mean(y_meas$se)
mean(x_meas$se)

col_genx = c("#5e3c99", "#e66101")
(se_check <- ggplot(data = rbind(x_meas, y_meas)) +
  geom_point(aes(x = pop, y = se, col = factor(generation)), size = 3) +
  scale_color_manual(values = col_genx) +
  labs(x = "", col = "generation") +
  geom_text(aes(x = 2, y = 0.35, label = paste("mean se generation 0 (x):", round(mean(x_meas$se), 2))),
            vjust = "inward", hjust = "inward") +
  geom_text(aes(x = 16, y = 0.35, label = paste("mean se generation 1 (y):", round(mean(y_meas$se), 2))),
            vjust = "inward", hjust = "inward"))
ggsave(file=paste0(dirname(dat_dir), "/", island, "_off_final_figures/", island, "_off_se_diagnostic_", one_phen, ".pdf"),
       plot=se_check, width=10, height=8)
###################################
# Size-adjusted means after scaling
pred_dat
# rm(list = setdiff(ls(), c("island", "one_phen")))
res_dir = paste0(island, "_off_SW/", island, "_off_results/")
# pred_dat = read.csv(file = paste0(res_dir, "tables/", island, "_off_size_adjusted_", one_phen[1], ".csv"))
# col_genx = c("#5e3c99", "#e66101")

diff_pop = names(table(pred_dat$pop))[table(pred_dat$pop)<2]
mod_dat = pred_dat[pred_dat$pop!=diff_pop, ]
sample_n(tbl = mod_dat, size = 10)

stanfile = "Littorina_offspring/scripts/offspring_err_in_var_model.stan"
writeLines(readLines(stanfile))
rstan_options(auto_write = TRUE)
# options(mc.cores = 4)
options(mc.cores = parallel::detectCores(logical = FALSE) - 2)

x_meas = mod_dat[mod_dat$generation==0, ]

y_meas = mod_dat[mod_dat$generation==1, ]
# nums <- unlist(lapply(y_meas, is.numeric))
# y_meas_ppt = round(y_meas[, nums], 2)

dat = list(N = nrow(x_meas), x = x_meas$fit, sd_x = x_meas$se.fit, y = y_meas$fit, sd_y = y_meas$se.fit)
err_in_var <- rstan::stan(file = stanfile,
                          data = dat, iter = 12000, warmup = 4000,
                          chains=4, refresh=12000,
                          control = list(stepsize = 0.01, adapt_delta = 0.99, max_treedepth = 15))

dir.create(paste0(res_dir, "models"))
dir.create(paste0(res_dir, "tables"))

saveRDS(err_in_var, paste0(res_dir, "models/", island, "_err_in_var.rds"))
# err_in_var = readRDS(paste0(res_dir, "models/", island, "_err_in_var.rds"))

# launch_shinystan(err_in_var)
# pairs(err_in_var)
pars = c("alpha", "beta", "sigma")
# pars = c("alpha", "beta")
pdf(file = paste0(dirname(res_dir), "/", island, "_off_final_figures/", island, "_off_err_in_var_pairs_adj-",one_phen[1],".pdf"),
    width = 8, height = 8)
pairs(err_in_var, pars = pars, include = TRUE)
dev.off()
# print(err_in_var, pars=c("alpha", "beta", "sigma"), digits=3)
print(err_in_var, pars = pars, digits=3)
# stan_pars = read.csv(file = paste0(res_dir, "tables/", island, "_err_in_var_stanfit_adj-", one_phen[1], "_pars.csv"))

# postd = extract(err_in_var)
# names(postd)
# dim(postd$alpha)

stbl = rstan::summary(err_in_var)
write.csv(x = stbl$summary, file = paste0(res_dir, "tables/", island, "_err_in_var_stanfit_adj-", one_phen[1], ".csv"))
# xtable::xtable(stbl$summary)
# xtable::xtable(stbl)
stbl = read.csv(file = paste0(res_dir, "tables/", island, "_err_in_var_stanfit_adj-", one_phen[1], ".csv"))
head(stbl)
stan_x_lat <- stbl[grepl(pattern = "x_lat", x = stbl$X), "mean"]
new_x_lat <- seq(from = min(stan_x_lat), to = max(stan_x_lat), length.out = 1000)
new_mu_yhat <- stan_pars[1, "mean"] + stan_pars[2, "mean"] * new_x_lat

plot(x_meas$x[, 'mean'], stbl$summary[grepl(pattern = "mu_yhat", x = row.names(stbl$summary)), 'mean'], col='red', pch=19)
points(x_meas$x[, 'mean'], y_meas$x[, 'mean'])

# plasticty plot
(pp <- ggplot() +
    geom_abline(slope = 1, linetype="dashed") +
    geom_point(aes(x = x_meas$fit, y = y_meas$fit),
               size=3) +
    geom_label_repel(aes(x = x_meas$fit, y = y_meas$fit,
                         label = x_meas$pop),
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50') +
    geom_point(aes(x = stbl$summary[1:nrow(x_meas), "mean"], y = stbl$summary[(nrow(x_meas)+1):(nrow(x_meas)*2), "mean"]),
               size=3, col="red") +
    geom_line(aes(x = new_x_lat, y = new_mu_yhat),
              size=2, col="red") +
    labs(title = paste0("size-adjusted ",one_phen[1]), x = "parents", y = "lab-reared adults") +
    theme(legend.position = "none",
          plot.title = element_text(size = 19, hjust = 0.5),
          axis.title = element_text(size = 18),
          axis.text = element_text(size=11),
          axis.ticks = element_line(size = 0.7),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.line = element_line(size = 0.2, linetype = "solid",
                                   colour = "black"),
          panel.grid = element_line(colour = "gray70", size = 0.2)))
ggsave(file=paste0(dirname(dat_dir), "/", island, "_off_final_figures/", island, "_off_plasticity_adj-", one_phen[1], ".svg"),
       plot=pp, width=10, height=8)
ggsave(file=paste0(dirname(dat_dir), "/", island, "_off_final_figures/", island, "_off_plasticity_adj-", one_phen[1], ".pdf"),
       plot=pp, width=10, height=8)
