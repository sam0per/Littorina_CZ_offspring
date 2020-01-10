rm(list = ls())

.packagesdev = "thomasp85/patchwork"
.packages = c("ggplot2", "dplyr", "reshape2", "parallel", "optparse", "tidyr", "splitstackshape", "data.table", "gdata")
# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
.instdev <- basename(.packagesdev) %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
if(length(.packagesdev[!.instdev]) > 0) devtools::install_github(.packagesdev[!.instdev])
# Load packages into session
lapply(.packages, require, character.only=TRUE)
lapply(basename(.packagesdev), require, character.only=TRUE)

option_list = list(
  make_option(c("-I", "--island"), type="character", default=NULL,
              help="Island where the parents were collected: CZA or CZD", metavar="character"),
  make_option(c("-G", "--gen"), type="character", default=NULL,
              help="either gen0 for parents or gen1 for offspring", metavar="character"),
  make_option(c("-P", "--phenotype"), type="character", default=NULL,
              help="which phenotype to analyse: boldness, dissections, weight, thickness", metavar="character"))

opt_parser = OptionParser(option_list=option_list,
                          description = "",
                          epilogue = "Example: Rscript")
opt = parse_args(opt_parser)

if (is.null(opt$vcf) | is.null(opt$directory) | is.null(opt$ecotype) | is.null(opt$island)){
  print_help(opt_parser)
  stop("All arguments must be supplied.\n", call.=FALSE)
}

island = "CZA"
# cz_gen = "gen1"
cz_phen = "boldness"
if (island=="CZA") {
  dat_dir = "CZA_off_SW/CZA_off_raw_data_renamed"
} else {
  dat_dir = list.dirs(path = "CZD_off_SW/CZD_off_SW_Oct2016_rawdata")[grepl(pattern = cz_gen,
                                                                            x = list.dirs(path = "CZD_off_SW/CZD_off_SW_Oct2016_rawdata"))]
}
(phen_fl = list.files(path = dat_dir, full.names = TRUE)[grepl(pattern = cz_phen, x = list.files(path = dat_dir))])

#########################
######## boldness #######
#########################
phen_dt = lapply(1:length(phen_fl), function(x) {
  onedt = read.xls(xls = phen_fl[x], sheet = 1)
  coldt = colnames(onedt)[grepl(pattern = "test", x = colnames(onedt))]
  # return(coldt)
  onedtcol = onedt[, which(colnames(onedt) %in% coldt)]
  cbind(snail_ID = onedt[, 1], onedtcol[,1:12])
})
lapply(phen_dt, head)
lapply(phen_dt, ncol)
lapply(phen_dt, colnames)
lapply(phen_dt, str)


# function for log likelihoods of individual observations, given mean and sd, assuming normal distribution (inpu is log(obs))
LLfunction <- function(obs,x_bar,x_sd){
  LL<-dnorm(obs, mean=x_bar, sd=x_sd,log=TRUE)
  LL[obs>log(15)]<-pnorm(log(15), mean=x_bar, sd=x_sd, lower.tail=FALSE, log.p=TRUE)
  
  return(LL)
}
# this loop is not needed for CZA because it has a single measure of activity
phen_fl
dtnum = 3
bold = phen_dt[[dtnum]]
head(bold)
str(bold)
bold[,-1] = apply(bold[,-1], MARGIN = 2, FUN = function(x) as.integer(x))

for (rep in 1:3){
  bold[,paste("crawl",rep,"_dec",sep="")] <- bold[,paste("test",rep,"_min",sep="")]+bold[,paste("test",rep,"_sec",sep="")]/60
  bold[is.na(bold[,paste("crawl",rep,"_dec",sep="")]),paste("crawl",rep,"_dec",sep="")] <- 16
  bold[(bold[,paste("crawl",rep,"_dec",sep="")]),paste("crawl",rep,"_dec",sep="")==0] <- 0.1
  bold[,paste("fullyout",rep,"_dec",sep="")] <- bold[,paste("test",rep,"_min.1",sep="")]+bold[,paste("test",rep,"_sec.1",sep="")]/60
  bold[is.na(bold[,paste("fullyout",rep,"_dec",sep="")]),paste("fullyout",rep,"_dec",sep="")] <- 16
  bold[(bold[,paste("fullyout",rep,"_dec",sep="")]),paste("fullyout",rep,"_dec",sep="")==0] <- 0.1
  #bold[bold[,paste("fullyout",rep,"_dec",sep="")]<bold[,paste("bold",rep,"_dec",sep="")],paste("bold",rep,"_dec",sep="")] <- bold[bold[,paste("fullyout",rep,"_dec",sep="")]<bold[,paste("bold",rep,"_dec",sep="")],paste("fullyout",rep,"_dec",sep="")]
}
head(bold)
# function to estimate mean for boldness score, allowing for those that did not get fully out or crawl 
replace <- function(lb1,lb2,lb3){
  lmb <- (lb1+lb2+lb3)/3  # initial mean log(boldness) per snail
  # estimate of sd among trials, averaged over snails
  bold_sd <- sqrt(mean(((lb1-lmb)^2+(lb2-lmb)^2+(lb3-lmb)^2)/3, na.rm = T)) 
  
  bold_max <- rep(c(0),length(lb1))   # set initial estimate to zero for all snails
  LL_max <-  LLfunction(lb1,3,bold_sd)+LLfunction(lb2,0,bold_sd)+LLfunction(lb3,0,bold_sd)
  
  # try a sequence of values (log scale), check likelihood and save estimate if likelihood has improved
  for (bold_mean in seq(3,0,-0.05)) {
    LL_tot <- LLfunction(lb1,bold_mean,bold_sd)+LLfunction(lb2,bold_mean,bold_sd)+LLfunction(lb3,bold_mean,bold_sd)
    
    bold_max[LL_tot>LL_max] <- bold_mean
    LL_max[LL_tot>LL_max] <- LL_tot[LL_tot>LL_max]
  }
  
  # replace 16s with draws from a distribution with mean and sd 
  prob <- bold_max==3
  rep16 <- rep(0,length(bold_max[prob]))
  c <- 1
  while (c <= length(rep16)){
    draw<-rnorm(1,mean=mean(bold_max[bold_max<3]),sd=sqrt(var(bold_max[bold_max<3])))
    if (draw > log(15)){
      rep16[c]<-draw
      c<-c+1
    }
  }
  
  bold_max[prob] <- rep16
  
  # and repeat with revised mean and sd
  c <- 1
  while (c <= length(rep16)){
    draw<-rnorm(1,mean=mean(bold_max),sd=sqrt(var(bold_max)))
    if (draw > log(15)){
      rep16[c]<-draw
      c<-c+1
    }
  }
  
  bold_max[prob] <- rep16
  
  return(bold_max)
  
} # end of 'replace' function


bold$crawl_mean <- replace(log(bold$crawl1_dec),log(bold$crawl2_dec),lb3 <- log(bold$crawl3_dec)) 

bold$fullyout_mean <- replace(log(bold$fullyout1_dec),log(bold$fullyout2_dec),lb3 <- log(bold$fullyout3_dec)) 

bold$bold_score <- (bold$crawl_mean+bold$fullyout_mean)/2

head(bold)
# check distribution (if needed)
hist(bold$crawl_mean, nclass = 21)
hist(bold$fullyout_mean, nclass = 21) 
hist(bold$bold_score,nclass=21)

# write result to file
(nu_phen_fl = tools::file_path_sans_ext(phen_fl[dtnum]))
if (island=="CZD") {
  (nu_phen_dir = paste0(dirname(dirname(dirname(nu_phen_fl))), "/", island, "_off_SW_Oct2016_finaldata"))
} else {
  (nu_phen_dir = paste0(dirname(dirname(nu_phen_fl)), "/", island, "_off_final_data"))
}
dir.create(nu_phen_dir)
(cdate = as.character(as.Date(Sys.time(), "%Y-%m-%d"), "%Y%m%d"))
write.csv(bold, paste0(nu_phen_dir, "/", basename(nu_phen_fl), "_", cdate, ".csv"), row.names = FALSE)

# pc <- princomp(bold[,20:25])
# pc1 <- pc$scores[,1]
# plot(bold$score,pc1)

# rm(list = ls())
##########################
######## thickness #######
##########################
island = "CZD"
cz_gen = "gen1"
cz_phen = "thickness"
if (island=="CZA") {
  dat_dir = "CZA_off_SW/CZA_off_raw_data_renamed"
} else {
  dat_dir = list.dirs(path = "CZD_off_SW/CZD_off_SW_Oct2016_rawdata")[grepl(pattern = cz_gen,
                                                                            x = list.dirs(path = "CZD_off_SW/CZD_off_SW_Oct2016_rawdata"))]
}
(phen_fl = list.files(path = dat_dir, full.names = TRUE)[grepl(pattern = cz_phen, x = list.files(path = dat_dir))])
phen_dt = lapply(1:length(phen_fl), function(x) {
  tmp = read.xls(xls = phen_fl[x], sheet = 1)
  tmp[,2:4] = apply(tmp[, 2:4], MARGIN = 2, FUN = function(x) as.numeric(x))
  if (mean(tmp[, "thickness.1"], na.rm = TRUE) > 1) {
    tmp[,2:4] = apply(tmp[, 2:4], MARGIN = 2, FUN = function(x) x/1000)
  }
  return(tmp)
})
lapply(phen_dt, head)
lapply(phen_dt, ncol)
lapply(phen_dt, colnames)
lapply(phen_dt, summary)
lapply(phen_dt, str)

phen_dt_thick = lapply(1:length(phen_fl), function(x) {
  mutate(phen_dt[[x]], mean_thickness = rowMeans(phen_dt[[x]][, 2:4], na.rm = TRUE))
})
# write result to file
(cdate = as.character(as.Date(Sys.time(), "%Y-%m-%d"), "%Y%m%d"))
lapply(1:length(phen_fl), function(x) {
  nu_phen_fl = tools::file_path_sans_ext(phen_fl[x])
  if (island=="CZD") {
    (nu_phen_dir = paste0(dirname(dirname(dirname(nu_phen_fl))), "/", island, "_off_SW_Oct2016_finaldata"))
  } else {
    (nu_phen_dir = paste0(dirname(dirname(nu_phen_fl)), "/", island, "_off_final_data"))
  }
  dir.create(nu_phen_dir)
  write.csv(phen_dt_thick[[x]], paste0(nu_phen_dir, "/", basename(nu_phen_fl), "_", cdate, ".csv"), row.names = FALSE)
})

# rm(list = ls())
#######################
######## weight #######
# TODO:
# type in weight for CZD off generation 1 300
#######################
island = "CZD"
cz_gen = "gen1"
cz_phen = "weight"
if (island=="CZA") {
  dat_dir = "CZA_off_SW/CZA_off_raw_data_renamed"
} else {
  dat_dir = list.dirs(path = "CZD_off_SW/CZD_off_SW_Oct2016_rawdata")[grepl(pattern = cz_gen,
                                                                            x = list.dirs(path = "CZD_off_SW/CZD_off_SW_Oct2016_rawdata"))]
}
(phen_fl = list.files(path = dat_dir, full.names = TRUE)[grepl(pattern = cz_phen, x = list.files(path = dat_dir))])
phen_dt = lapply(1:length(phen_fl), function(x) {
  if (island=="CZA") {
    read.xls(xls = phen_fl[x], sheet = 1)[, 1:5]
  } else {
    read.xls(xls = phen_fl[x], sheet = 1)
  }
})
lapply(phen_dt, head)
lapply(phen_dt, ncol)
lapply(phen_dt, colnames)
lapply(phen_dt, summary)
lapply(phen_dt, str)

# write result to file
(cdate = as.character(as.Date(Sys.time(), "%Y-%m-%d"), "%Y%m%d"))
lapply(1:length(phen_fl), function(x) {
  nu_phen_fl = tools::file_path_sans_ext(phen_fl[x])
  if (island=="CZD") {
    (nu_phen_dir = paste0(dirname(dirname(dirname(nu_phen_fl))), "/", island, "_off_SW_Oct2016_finaldata"))
  } else {
    (nu_phen_dir = paste0(dirname(dirname(nu_phen_fl)), "/", island, "_off_final_data"))
  }
  dir.create(nu_phen_dir)
  write.csv(phen_dt[[x]], paste0(nu_phen_dir, "/", basename(nu_phen_fl), "_", cdate, ".csv"), row.names = FALSE)
})

# rm(list = ls())
############################
######## dissections #######
############################
island = "CZD"
cz_gen = "gen1"
cz_phen = "dissections"
if (island=="CZA") {
  dat_dir = "CZA_off_SW/CZA_off_raw_data_renamed"
} else {
  dat_dir = list.dirs(path = "CZD_off_SW/CZD_off_SW_Oct2016_rawdata")[grepl(pattern = cz_gen,
                                                                            x = list.dirs(path = "CZD_off_SW/CZD_off_SW_Oct2016_rawdata"))]
}
(phen_fl = list.files(path = dat_dir, full.names = TRUE)[grepl(pattern = cz_phen, x = list.files(path = dat_dir))])
phen_dt = lapply(1:length(phen_fl), function(x) {
  orig = read.xls(xls = phen_fl[x], sheet = 1)
  if (island=="CZD" & cz_gen=="gen0") {
    colnames(orig) = c("snail_ID", "snail_IDD", "broodpouch", "parasites", "notes")
  } else if (island=="CZD" & cz_gen=="gen1") {
    orig = orig[, 1:11]
    colnames(orig) = c("snail_ID", "broodpouch", "penis", "lengthglandrow", "lengthtip", "nglands", "nglandrows",
                       "magnif", "ciliates", "parasites", "notes")
  } else {
    if (x==1 | x==2) {
      colnames(orig) = c("snail_ID", "sexontube", "broodpouch", "penis", "parasites", "ciliates", "experimenter", "notes")
    } else {
      colnames(orig) = c("snail_ID", "broodpouch", "penis", "lengthglandrow", "lengthtip", "nglands", "experimenter", "notes")
    }
  }
  return(orig)
})
lapply(phen_dt, head)
lapply(phen_dt, ncol)
lapply(phen_dt, colnames)
lapply(phen_dt, summary)
lapply(phen_dt, str)
#### identify sex of each snail, using brood pouch and penis data ####
sex <- function(b, p, isl) {
  if (isl=="CZD" & cz_gen=="gen0") {
    if(b=="Y" & (is.na(b)==F)) y <- "female"
    if(b=="N" & (is.na(b)==F)) y <- "NA"
    if((b %in% c("Y", "N"))==F) y <- "NA"
  } else {
    if(b=="Y" & p=="N" & (is.na(b)==F)) y <- "female"
    if(b=="N" & p=="Y" & (is.na(b)==F)) y <- "male"
    if((b %in% c("Y", "N"))==F | (p %in% c("Y", "N"))==F | b==p) y<-"NA"
  }
  return(y)
}
# sex(b = "Y", p = "Y", isl = "CZD")
# island="CZD"
phen_sex = lapply(1:length(phen_fl), function(i) {
  if (island=="CZD" & cz_gen=="gen0") {
    phen_dt[[i]][, "penis"] = 0
  }
  dt_na = mutate(phen_dt[[i]], sex=apply(phen_dt[[i]][, c("broodpouch", "penis")],
                                         MARGIN = 1, FUN = function(x) sex(b = x[1], p = x[2], isl = island)))
  dt_nona = dt_na[which(dt_na[, "sex"] != "NA"), ]
  if (island=="CZD" & cz_gen=="gen0") {
    dt_nona = dt_nona[dt_nona[, "broodpouch"]=="Y", ]
    dt_nona = dt_nona[dt_nona[, "parasites"]!="Y", ]
  } else if (island=="CZD" & cz_gen=="gen1") {
    dt_nona = dt_nona[dt_nona[, "parasites"]!="Y", ]
  } else {
    if (i==1) {
      dt_nona = dt_nona[dt_nona[, "broodpouch"]=="Y", ]
      dt_nona = dt_nona[dt_nona[, "parasites"]!="Y", ]
    } else if (i==2) {
      dt_nona = dt_nona[dt_nona[, "sexontube"]=="M", ]
      dt_nona = dt_nona[dt_nona[, "parasites"]!="Y", ]
    }
  }
  return(dt_nona)
})
lapply(phen_sex, head)
lapply(phen_sex, nrow)
lapply(phen_sex, summary)

# write result to file
(cdate = as.character(as.Date(Sys.time(), "%Y-%m-%d"), "%Y%m%d"))
lapply(1:length(phen_fl), function(x) {
  nu_phen_fl = tools::file_path_sans_ext(phen_fl[x])
  if (island=="CZD") {
    (nu_phen_dir = paste0(dirname(dirname(dirname(nu_phen_fl))), "/", island, "_off_SW_Oct2016_finaldata"))
  } else {
    (nu_phen_dir = paste0(dirname(dirname(nu_phen_fl)), "/", island, "_off_final_data"))
  }
  dir.create(nu_phen_dir)
  write.csv(phen_sex[[x]], paste0(nu_phen_dir, "/", basename(nu_phen_fl), "_", cdate, ".csv"), row.names = FALSE)
})

# rm(list = ls())
#################################
######## merge phenotypes #######
#################################
island = "CZA"
phenos = c("boldness", "dissections", "thickness", "weight")
# cz_gen = "gen1"
if (island=="CZA") {
  dat_dir = "CZA_off_SW/CZA_off_final_data"
} else {
  dat_dir = "CZD_off_SW/CZD_off_SW_Oct2016_finaldata"
}
phen_fl = lapply(1:length(phenos), function(x) {
  if (island=="CZD") {
    list.files(path = dat_dir, pattern = paste0(cz_gen, "_\\d+_", phenos[x]), full.names = TRUE)
  } else {
    list.files(path = dat_dir, pattern = phenos[x], full.names = TRUE)
  }
})
phen_dt = lapply(1:length(phenos), function(x) {
  lapply(1:length(phen_fl[[x]]), function(y) {
    read.csv(file = phen_fl[[x]][y])
  })
})

phen_all = lapply(1:length(phenos), function(x) {
  phen_merged = Reduce(function(...) merge(..., all=TRUE), phen_dt[[x]])
  phen_merged[, "snail_ID"] = as.character(phen_merged[, "snail_ID"])
  return(phen_merged)
})
lapply(phen_all, head)
lapply(phen_all, colnames)
lapply(phen_all, nrow)
# head(phen_all[[4]])
phen_all_dt = Reduce(function(...) merge(..., all=TRUE, by="snail_ID"), phen_all)
head(phen_all_dt)
colnames(phen_all_dt)
colnames(phen_all_dt)[which(colnames(phen_all_dt)=="score")] = "bold_score"
sum(grepl(pattern = "snail_ID", x = colnames(phen_all_dt)))

(cdate = as.character(as.Date(Sys.time(), "%Y-%m-%d"), "%Y%m%d"))
if (exists("cz_gen")) {
  write.csv(phen_all_dt, file = paste0(dat_dir, "/", island, "_off_", cz_gen,"_all_phenos_", cdate, ".csv"), row.names = FALSE)
} else {
  write.csv(phen_all_dt, file = paste0(dat_dir, "/", island, "_off_all_phenos_", cdate, ".csv"), row.names = FALSE)
}

# foo = read.csv("CZA_off_SW/CZA_off_final_data/CZA_off_all_phenos_main_20191216.csv")
# colnames(foo)
# kcol = intersect(colnames(phen_all_dt), colnames(foo))
# mdt = phen_all_dt[, which(colnames(phen_all_dt) %in% kcol)]
# head(mdt)
# 
# (cdate = as.character(as.Date(Sys.time(), "%Y-%m-%d"), "%Y%m%d"))
# if (exists("cz_gen")) {
#   write.csv(mdt, file = paste0(dat_dir, "/", island, "_off_", cz_gen,"_all_phenos_main_", cdate, ".csv"), row.names = FALSE)
# } else {
#   write.csv(mdt, file = paste0(dat_dir, "/", island, "_off_all_phenos_main_", cdate, ".csv"), row.names = FALSE)
# }

# rm(list = ls())
#######################
######## clines #######
#######################
island = "CZA"
# cz_gen = "gen1"
if (island=="CZA") {
  dat_dir = "CZA_off_SW/CZA_off_final_data"
} else {
  dat_dir = "CZD_off_SW/CZD_off_SW_Oct2016_finaldata"
}
mdt = read.csv(list.files(path = dat_dir, full.names = TRUE)[grepl(pattern = "all_phenos_main", x = list.files(path = dat_dir))])
head(mdt)

sdt = read.csv(list.files(path = dat_dir, full.names = TRUE)[grepl(pattern = "spatial", x = list.files(path = dat_dir))])
head(sdt)

cdt = merge(mdt, sdt, by = "snail_ID")
head(cdt)
tail(cdt)
cdt$pop = substr(cdt$snail_ID, start = 1, stop = 1)
e_LCmeanDist = aggregate(cdt[, paste0(island, "_LCmeanDist")], by=list(cdt$pop), mean)
colnames(e_LCmeanDist) = c("pop", "e_LCmeanDist")
head(e_LCmeanDist)

head(mdt)
mdt$pop = substr(mdt$snail_ID, start = 1, stop = 1)
edt = merge(mdt, e_LCmeanDist, by = "pop")
edt[which(edt$snail_ID=="G_O7"), "snail_ID"] = "G_07"
edt$generation = 1
edt[grep("_[0-9][0-9]$", edt$snail_ID), "generation"] = 0
head(edt)
str(edt)

phenos = c("bold_score", "mean_thickness", "weight_g")
phenos_plot = lapply(seq_along(phenos), function(x) {
  ggplot(data = edt) +
    facet_wrap(~generation) +
    geom_point(aes(x = e_LCmeanDist, y = edt[,phenos[x]], col=factor(generation))) +
    labs(y = phenos[x]) +
    scale_color_manual(values = c("#f1a340", "#998ec3")) +
    theme(legend.position = "none",
          strip.text = element_text(size=13),
          strip.background = element_rect(fill="lightblue", colour="black", size=0.5),
          axis.title = element_text(size = 12),
          axis.text = element_text(size=11),
          axis.ticks = element_line(size = 0.7),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.line = element_line(size = 0.2, linetype = "solid",
                                   colour = "black"),
          panel.grid = element_line(colour = "gray70", size = 0.2))
})
dir.create(path = paste0(dirname(dat_dir), "/", island, "_off_final_figures"))
(cdate = as.character(as.Date(Sys.time(), "%Y-%m-%d"), "%Y%m%d"))
lapply(seq_along(phenos), function(x) {
  ggsave(filename = paste0(dirname(dat_dir), "/", island, "_off_final_figures/", island, "_off_",
                           phenos[x], "_scatter_", cdate, ".svg"),
         plot = phenos_plot[[x]], device = "svg", width = 8, height = 8)
})

weight_sort = edt[order(edt$weight_g, decreasing = TRUE),]
weight_sort[weight_sort$generation==1, c("snail_ID", "weight_g")][1:20, ]

# g0_edt = edt[grep("_[0-9][0-9]$", edt$snail_ID), ]
# head(g0_edt)
# plot(g0_edt$e_LCmeanDist, g0_edt$mean_thickness)
# plot(g0_edt$e_LCmeanDist, g0_edt$bold_score)
# plot(g0_edt$e_LCmeanDist, g0_edt$weight_g)

pcadt = read.csv(list.files(path = dat_dir, full.names = TRUE)[grepl(pattern = "PCA", x = list.files(path = dat_dir))],
                 sep = ";")
head(pcadt)
pcadt$pop = substr(pcadt$ID, start = 1, stop = 1)
pca_edt = merge(pcadt, e_LCmeanDist, by = "pop")
head(pca_edt)
pca_edt$generation = 0
# table(pca_edt$Sex)
pca_edt[pca_edt$Sex=="offspring", "generation"] = 1
table(factor(pca_edt$generation))

col_pc = c("PC1", "PC2")
pca_plot = lapply(seq_along(col_pc), function(x) {
  ggplot(data = pca_edt) +
    facet_wrap(~generation) +
    geom_point(aes(x = e_LCmeanDist, y = pca_edt[,col_pc[x]], col=factor(generation))) +
    labs(y = col_pc[x]) +
    scale_color_manual(values = c("#f1a340", "#998ec3")) +
    theme(legend.position = "none",
          strip.text = element_text(size=13),
          strip.background = element_rect(fill="lightblue", colour="black", size=0.5),
          axis.title = element_text(size = 12),
          axis.text = element_text(size=11),
          axis.ticks = element_line(size = 0.7),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.line = element_line(size = 0.2, linetype = "solid",
                                   colour = "black"),
          panel.grid = element_line(colour = "gray70", size = 0.2))
})
lapply(seq_along(col_pc), function(x) {
  ggsave(filename = paste0(dirname(dat_dir), "/", island, "_off_final_figures/", island, "_off_",
                           col_pc[x], "_scatter_", cdate, ".svg"),
         plot = pca_plot[[x]], device = "svg", width = 8, height = 8)
})

# there are no IDs starting from 200, only 300
pc2_sort = pca_edt[order(pca_edt$PC2, decreasing = TRUE),]
pc2_sort[pc2_sort$generation==1, c("ID", "PC2")][1:20, ]
tail(pc2_sort[pc2_sort$generation==1, c("ID", "PC2")])






