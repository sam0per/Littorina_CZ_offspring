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
# island = "CZA"
# cz_gen = "gen1"
cz_phen = "mean_thickness"

dat_dir = paste0(island, "_off_SW/", island, "_off_final_data/")

dat_off = read.csv(file = paste0(dat_dir, island, "_off_all_phenos_main_20200110.csv"))
# colnames(dat_off)
# head(dat_off)
# sample_n(dat_off, size = 10)

dat_off = separate(data = dat_off, col = "snail_ID", into = c("pop", "ID"), sep = "_")
# dat_off[c(402, 403, 404, 829, 830, 831),]

# mean(scale(dat_off[, cz_phen]), na.rm = TRUE)
# mean(dat_off[, cz_phen], na.rm = TRUE)

# table(nchar(as.character(dat_off$ID)))
dat_off[, "generation"] = 1
dat_off[which(nchar(as.character(dat_off$ID)) == 2), "generation"] = 0

dat_gen0 = dat_off[dat_off$generation==0, c("pop", "ID", cz_phen)]

dat_gen1 = dat_off[dat_off$generation==1, c("pop", "ID", cz_phen)]
diff_ypop = apply(X = dat_gen1[, c("pop", "ID")], MARGIN = 2,
                  FUN = function(x) !grepl(pattern = "x", x = x))
diff_ypop_idx = which(apply(diff_ypop, MARGIN = 1, FUN = sum)==2)
dat_gen1 = dat_gen1[diff_ypop_idx, ]

diff_pop = setdiff(dat_gen0$pop, dat_gen1$pop)
dat_gen0 = dat_gen0[dat_gen0$pop!=diff_pop,]

dat_gen0[, paste0("scaled_", cz_phen)] = (dat_gen0$mean_thickness - mean(dat_gen0$mean_thickness, na.rm = TRUE)) / sd(dat_gen0$mean_thickness, na.rm = TRUE)
x_meas = aggregate(x = dat_gen0[, paste0("scaled_", cz_phen)], by = list(pop = dat_gen0$pop),
                   FUN = function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE)))


dat_gen1[, paste0("scaled_", cz_phen)] = (dat_gen1$mean_thickness - mean(dat_gen0$mean_thickness, na.rm = TRUE)) / sd(dat_gen0$mean_thickness, na.rm = TRUE)
y_meas = aggregate(x = dat_gen1[, paste0("scaled_", cz_phen)], by = list(pop = dat_gen1$pop),
                   FUN = function(y) c(mean = mean(y, na.rm = TRUE), sd = sd(y, na.rm = TRUE)))

# plot(x_meas$x[, 'mean'], y_meas$x[, 'mean'])

rstan_options(auto_write = TRUE)
options(mc.cores = 4)
# options(mc.cores = parallel::detectCores(logical = FALSE) - 2)

dat = list(N = nrow(x_meas$x), x = x_meas$x[, 'mean'], sd_x = x_meas$x[,"sd"],
           y = y_meas$x[, 'mean'], sd_y = y_meas$x[,"sd"])
err_in_var = rstan::stan(file = stanfile,
                         data = dat, iter = 12000, warmup =4000,
                         chains=4, refresh=12000,
                         control = list(stepsize=0.01, adapt_delta = 0.99))
dir.create(paste0(island, "_off_SW/", island, "_off_results"))
res_dir = paste0(island, "_off_SW/", island, "_off_results/")
dir.create(paste0(res_dir, "models"))
dir.create(paste0(res_dir, "tables"))

saveRDS(err_in_var, paste0(res_dir, "models/", island, "_err_in_var.rds"))
# err_in_var = readRDS(paste0(res_dir, "models/", island, "_err_in_var.rds"))

# launch_shinystan(err_in_var)
# pairs(err_in_var)
print(err_in_var, pars=c("alpha", "beta", "sigma"), digits=3)

# postd = extract(err_in_var)
# names(postd)
# dim(postd$alpha)

stbl = rstan::summary(err_in_var)
write.csv(x = stbl$summary, file = paste0(res_dir, "tables/", island, "_err_in_var_stanfit.csv"))
# xtable::xtable(stbl$summary)
# xtable::xtable(stbl)