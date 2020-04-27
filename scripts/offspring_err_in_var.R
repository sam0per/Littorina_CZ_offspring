rm(list = ls())
.packages = c("optparse", "dplyr", "tidyr", "rstan", "shinystan", "xtable")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])

# Load packages into session
lapply(.packages, require, character.only=TRUE)


option_list = list(
  make_option(c("-i", "--island"), type="character", default=NULL,
              help="name of the island [CZA, CZD]", metavar="character"),
  make_option(c("-s", "--stanfile"), type="character", default=NULL,
              help="model written in Stan", metavar="character"),
  make_option(c("-o", "--output"), type = "character", default = "output",
              help = "prefix for output files [default: %default]", metavar = "character"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$island) | is.null(opt$stanfile)){
  print_help(opt_parser)
  stop("At least two arguments must be supplied, island name and stan file.\n", call.=FALSE)
}

# pref_out = opt$output
island = opt$island
stanfile = opt$stanfile
# stanfile = "Littorina_offspring/scripts/offspring_err_in_var_model.stan"
# cz_gen = "gen1"
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
head(mod_dat)

table(mod_dat$pop)
table(mod_dat$ID)

mod_dat = mod_dat[!grepl(pattern = "x", x = mod_dat$ID), ]
# diff_ypop = apply(X = dat_gen1[, c("pop", "ID")], MARGIN = 2,
#                   FUN = function(x) !grepl(pattern = "x", x = x))
# diff_ypop_idx = which(apply(diff_ypop, MARGIN = 1, FUN = sum)==2)
# dat_gen1 = dat_gen1[diff_ypop_idx, ]

diff_pop = setdiff(mod_dat[mod_dat$generation==0, "pop"], mod_dat[mod_dat$generation==1, "pop"])
mod_dat = mod_dat[mod_dat$pop!=diff_pop, ]
sample_n(tbl = mod_dat, size = 10)

# dat_gen0[, paste0("scaled_", cz_phen)] = (dat_gen0$mean_thickness - mean(dat_gen0$mean_thickness, na.rm = TRUE)) / sd(dat_gen0$mean_thickness, na.rm = TRUE)
# x_meas = aggregate(x = dat_gen0[, paste0("scaled_", cz_phen)], by = list(pop = dat_gen0$pop),
#                    FUN = function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE)))

x_meas = aggregate(x = mod_dat[mod_dat$generation==0, paste0('scaled_', one_phen[1])],
                   by = list(pop = mod_dat[mod_dat$generation==0, "pop"]),
                   FUN = function(x) c(mean_x = mean(x, na.rm = TRUE), sd_x = sd(x, na.rm = TRUE)))

# dat_gen1[, paste0("scaled_", cz_phen)] = (dat_gen1$mean_thickness - mean(dat_gen0$mean_thickness, na.rm = TRUE)) / sd(dat_gen0$mean_thickness, na.rm = TRUE)
# y_meas = aggregate(x = dat_gen1[, paste0("scaled_", cz_phen)], by = list(pop = dat_gen1$pop),
#                    FUN = function(y) c(mean = mean(y, na.rm = TRUE), sd = sd(y, na.rm = TRUE)))

y_meas = aggregate(x = mod_dat[mod_dat$generation==1, paste0('scaled_', one_phen[1])],
                   by = list(pop = mod_dat[mod_dat$generation==1, "pop"]),
                   FUN = function(y) c(mean_y = mean(y, na.rm = TRUE), sd_y = sd(y, na.rm = TRUE)))
# plot(x_meas$x[, 'mean_x'], y_meas$x[, 'mean_y'])

writeLines(readLines(stanfile))
rstan_options(auto_write = TRUE)
options(mc.cores = 4)
# options(mc.cores = parallel::detectCores(logical = FALSE) - 2)

dat = list(N = nrow(x_meas$x), x = x_meas$x[, 'mean'], sd_x = x_meas$x[, "sd"])
err_in_var <- rstan::stan(file = stanfile,
                          data = dat, iter = 12000, warmup =4000,
                          chains=4, refresh=12000,
                          control = list(stepsize = 0.01, adapt_delta = 0.99, max_treedepth = 15))

dir.create(paste0(island, "_off_SW/", island, "_off_results"))
res_dir = paste0(island, "_off_SW/", island, "_off_results/")
dir.create(paste0(res_dir, "models"))
dir.create(paste0(res_dir, "tables"))

saveRDS(err_in_var, paste0(res_dir, "models/", island, "_err_in_var.rds"))
# err_in_var = readRDS(paste0(res_dir, "models/", island, "_err_in_var.rds"))

# launch_shinystan(err_in_var)
# pairs(err_in_var)
pars = c("alpha", "beta", "sigma")
pairs(err_in_var, pars = pars, include = TRUE)
print(err_in_var, pars=c("alpha", "beta", "sigma"), digits=3)

# postd = extract(err_in_var)
# names(postd)
# dim(postd$alpha)

stbl = rstan::summary(err_in_var)
write.csv(x = stbl$summary, file = paste0(res_dir, "tables/", island, "_err_in_var_stanfit.csv"))
# xtable::xtable(stbl$summary)
# xtable::xtable(stbl)

plot(x_meas$x[, 'mean'], stbl$summary[grepl(pattern = "mu_yhat", x = row.names(stbl$summary)), 'mean'], col='red', pch=19)
points(x_meas$x[, 'mean'], y_meas$x[, 'mean'])



