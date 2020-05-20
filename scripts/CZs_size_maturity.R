rm(list = ls())
island = "CZD"
if (island=="CZA") {
  dat_dir = "CZA_off_SW/CZA_off_final_data"
} else {
  dat_dir = "CZD_off_SW/CZD_off_final_data"
}
odate = c("20200406", "20200511")
tbl1 = read.csv(paste0(dat_dir, "/", island, "_all_phenos_main_", odate[2], ".csv"))
# OFF = read.csv(paste0(dat_dir, "/", island, "_off_all_phenos_main_", odate[1], ".csv"))
# colnames(CZ)
# CZ <- CZ[, order(colnames(CZ))]
# OFF <- OFF[, order(colnames(OFF))]
# colnames(OFF)
# identical(colnames(OFF), colnames(CZ))
table(tbl1$sex)
tbl1 <- tbl1[!is.na(tbl1$sex), ]
tbl1$maturity = 1
tbl1[tbl1$sex!="female" & tbl1$sex!="male", "maturity"] = 0
table(tbl1$maturity)
# table(OFF$sex)

tbl1$size_log <- log(tbl1$size_mm)
tbl1 <- tbl1[!is.na(tbl1$size_log), ]
range(tbl1$size_log)
# breaks = c(seq(from = min(tbl1$size_mm), to = (max(tbl1$size_mm) - 0.7), by = 0.7), (max(tbl1$size_mm) + 0.7))
bin_size <- 0.15
breaks = c(seq(from = min(tbl1$size_log), to = (max(tbl1$size_log) - bin_size), by = bin_size), (max(tbl1$size_log) + bin_size))
# tbl1$bin = cut(tbl1$size_mm, breaks, include.lowest = TRUE)
tbl1$bin = cut(tbl1$size_log, breaks, include.lowest = TRUE)
s_dat = group.CI(x = size_log ~ bin, data = tbl1)
s_dat$size_log.upper = ifelse(test = is.na(s_dat$size_log.upper), yes = s_dat$size_log.mean, no = s_dat$size_log.upper)
s_dat$size_log.lower = ifelse(test = is.na(s_dat$size_log.lower), yes = s_dat$size_log.mean, no = s_dat$size_log.lower)

p_dat = cbind(group.CI(x = maturity ~ bin, data = tbl1),
              size_log.mean = s_dat$size_log.mean)
p_dat$maturity.upper = ifelse(test = is.na(p_dat$maturity.upper), yes = p_dat$maturity.mean, no = p_dat$maturity.upper)
p_dat$maturity.lower = ifelse(test = is.na(p_dat$maturity.lower), yes = p_dat$maturity.mean, no = p_dat$maturity.lower)
p_dat[p_dat$maturity.upper > 1, "maturity.upper"] = 1
p_dat[p_dat$maturity.lower < 0, "maturity.lower"] = 0

(p_fit1 = ggplot(data = p_dat) +
    geom_errorbar(aes(x = size_log.mean, ymin = maturity.lower, ymax = maturity.upper, col = paste0(island, " observations")),
                  width = 0.05) +
    # geom_line(data = f_dat, aes(size_log, p_fit, col = "predictions"), size = 1.5) +
    geom_point(aes(x = size_log.mean, y = maturity.mean, col = paste0(island, " observations")), size = 2.5) +
    # geom_point(data = CZ_data_bin, aes(x = mean_ratio, y = mount), col='blue', size=3.5) +
    # scale_colour_manual(values=c("blue", "orange")) +
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
dir.create(path = paste0(dirname(dat_dir), "/", island, "_off_final_figures"))
fig_dir <- paste0(dirname(dat_dir), "/", island, "_off_final_figures")
ggsave(filename = paste0(fig_dir, "/", island, "_size_maturity.pdf"), plot = p_fit1, width = 8, height = 7)

# diagnose CZA/D observations of size and maturity
(mean_size_log <- aggregate(tbl1$size_log, by = list(sex = tbl1$maturity), mean))
unique(tbl1[tbl1$maturity == 1, "sex"])
ggplot(data = tbl1, aes(fill = sex)) +
  geom_histogram(aes(x = size_log), bins = 40)
CZ_spa <- read.csv(list.files(path = paste0(dat_dir, "/Littorina_", island), pattern = "spatial", full.names = TRUE))
head(CZ_spa)
head(tbl1)
tbl1 <- merge(tbl1, CZ_spa, by = "snail_ID")
(diagn1 <- ggplot(data = tbl1, aes(col = factor(maturity))) +
  geom_point(aes(x = LCmeanDist, y = size_log), size = 2) +
  labs(col = 'maturity'))
nu_CZD <- tbl1[-which.max(tbl1$size_log), ]
(diagn1 <- ggplot(data = nu_CZD, aes(col = factor(maturity))) +
    geom_point(aes(x = LCmeanDist, y = size_log), size = 2) +
    labs(col = 'maturity'))
tbl1 <- nu_CZD
ggsave(filename = paste0(fig_dir, "/", island, "_diagn_size_maturity.pdf"), plot = diagn1, width = 8, height = 7)
