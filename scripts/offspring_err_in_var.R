library(dplyr)
library(tidyr)
library(rstan)

island = "CZA"
# cz_gen = "gen1"
cz_phen = "mean_thickness"

dat_dir = paste0(island, "_off_SW/", island, "_off_final_data/")

dat_off = read.csv(file = paste0(dat_dir, island, "_off_all_phenos_main_20200110.csv"))
colnames(dat_off)
head(dat_off)
sample_n(dat_off, size = 10)

dat_off = separate(data = dat_off, col = "snail_ID", into = c("pop", "ID"), sep = "_")
# dat_off[c(402, 403, 404, 829, 830, 831),]

# mean(scale(dat_off[, cz_phen]), na.rm = TRUE)
# mean(dat_off[, cz_phen], na.rm = TRUE)

table(nchar(as.character(dat_off$ID)))
dat_off[, "generation"] = 1
dat_off[which(nchar(as.character(dat_off$ID)) == 2), "generation"] = 0

dat_gen0 = dat_off[dat_off$generation==0, c("pop", "ID", cz_phen)]
table(dat_gen0$pop)
mean(dat_gen0[dat_gen0$pop=="A", cz_phen])
x_meas = aggregate(x = dat_gen0[, cz_phen], by = list(pop = dat_gen0$pop),
                   FUN = function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE)))
# class(x_meas)
x_meas[, paste0("scaled_", cz_phen)] = scale(x_meas$x[,"mean"])[,1]
# mean(x_meas$x)
# sd(x_meas$x)


dat_gen1 = dat_off[dat_off$generation==1, c("pop", "ID", cz_phen)]
table(dat_gen1$pop)
mean(dat_gen1[dat_gen1$pop=="B", cz_phen], na.rm = TRUE)
y_meas = aggregate(x = dat_gen1[, cz_phen], by = list(pop = dat_gen1$pop),
                   FUN = function(y) c(mean = mean(y, na.rm = TRUE), sd = sd(y, na.rm = TRUE)))
y_meas = y_meas[complete.cases(y_meas), ]
class(x_meas)
y_meas[, paste0("scaled_", cz_phen)] = scale(y_meas$x[,"mean"])[,1]
# mean(x_meas$x)
# sd(x_meas$x)

(diff_pop = setdiff(x_meas$pop, y_meas$pop))
x_meas = x_meas[x_meas$pop!=diff_pop,]

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores(logical = FALSE) - 2)

dat = list(N = nrow(x_meas), x_meas = x_meas$scaled_mean_thickness, tau = x_meas$x[, 'sd'],
           y_mean = y_meas$scaled_mean_thickness)
err_in_var = rstan::stan(file = "Littorina_offspring/scripts/offspring_err_in_var_model.stan",
                         data = dat, iter = 8000, warmup = 2000,
                         chains=4, refresh=8000)
saveRDS(err_in_var, "models/err_in_var.rds")