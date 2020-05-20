rm(list = ls())

.packagesdev = "thomasp85/patchwork"
.packages = c("ggplot2", "reshape2", "parallel", "optparse", "tidyr", "splitstackshape", "data.table", "gdata",
              "Rmisc", "ggrepel", "dplyr")
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

island = "CZD"
dat_dir <- paste0(island, "_off_SW/", island, "_off_raw_data/Littorina_", island)
# cz_gen = "gen1"
if (island=="CZA") {
  dat_dir = "CZA_off_SW/CZA_off_raw_data_renamed"
} else {
  dat_dir = list.dirs(path = "CZD_off_SW/CZD_off_SW_Oct2016_rawdata")[grepl(pattern = cz_gen,
                                                                            x = list.dirs(path = "CZD_off_SW/CZD_off_SW_Oct2016_rawdata"))]
}

#########################
######## boldness #######
#########################
cz_phen = "boldness"
(phen_fl = list.files(path = dat_dir, full.names = TRUE)[grepl(pattern = cz_phen, x = list.files(path = dat_dir))])
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
dtnum = 4
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
island = "CZA"
cz_gen = "gen1"
cz_phen = "thickness"
if (island=="CZA") {
  dat_dir = "CZA_off_SW/CZA_off_raw_data_renamed"
} else {
  dat_dir = list.dirs(path = "CZD_off_SW/CZD_off_SW_Oct2016_rawdata")[grepl(pattern = cz_gen,
                                                                            x = list.dirs(path = "CZD_off_SW/CZD_off_SW_Oct2016_rawdata"))]
}
(phen_fl = list.files(path = dat_dir, full.names = TRUE)[grepl(pattern = cz_phen, x = list.files(path = dat_dir))])

# for CZA/D wild samples #
phen_dt <- read.csv(file = phen_fl)
head(phen_dt)
phen_dt_thick <- mutate(phen_dt, mean_thickness = rowMeans(phen_dt[, 2:4], na.rm = TRUE))
(cdate = as.character(as.Date(Sys.time(), "%Y-%m-%d"), "%Y%m%d"))
nu_phen_fl = tools::file_path_sans_ext(phen_fl)
cat("Saving final", cz_phen, "data to directory... \n")
if (island=="CZD") {
  (nu_phen_dir = paste0(dirname(dirname(dirname(nu_phen_fl))), "/", island, "_off_final_data/Littorina_", island))
} else {
  (nu_phen_dir = paste0(dirname(dirname(nu_phen_fl)), "/", island, "_off_final_data"))
}
dir.create(nu_phen_dir)
write.csv(phen_dt_thick, paste0(nu_phen_dir, "/", basename(nu_phen_fl), "_", cdate, ".csv"), row.names = FALSE)
########################

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
rm(list = setdiff(ls(), c("dat_dir", "island", "nu_phen_dir")))
#######################
######## weight #######
# TODO:
# type in weight for CZD off generation 1 300
#######################
island = "CZA"
cz_gen = "gen1"
cz_phen = "weight"
if (island=="CZA") {
  dat_dir = "CZA_off_SW/CZA_off_raw_data_renamed"
} else {
  dat_dir = list.dirs(path = "CZD_off_SW/CZD_off_SW_Oct2016_rawdata")[grepl(pattern = cz_gen,
                                                                            x = list.dirs(path = "CZD_off_SW/CZD_off_SW_Oct2016_rawdata"))]
}
(phen_fl = list.files(path = dat_dir, full.names = TRUE)[grepl(pattern = cz_phen, x = list.files(path = dat_dir))])

# for CZA wild samples #
phen_dt <- read.csv(file = phen_fl)
head(phen_dt)
(cdate = as.character(as.Date(Sys.time(), "%Y-%m-%d"), "%Y%m%d"))
nu_phen_fl = tools::file_path_sans_ext(phen_fl)
cat("Saving final", cz_phen, "data to directory... \n")
if (island=="CZD") {
  (nu_phen_dir = paste0(dirname(dirname(dirname(nu_phen_fl))), "/", island, "_off_SW_Oct2016_finaldata"))
} else {
  (nu_phen_dir = paste0(dirname(dirname(nu_phen_fl)), "/", island, "_off_final_data"))
}
dir.create(nu_phen_dir)
write.csv(phen_dt, paste0(nu_phen_dir, "/", basename(nu_phen_fl), "_", cdate, ".csv"), row.names = FALSE)
########################

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
rm(list = setdiff(ls(), c("dat_dir", "island", "nu_phen_dir")))
############################
######## dissections #######
############################
island = "CZA"
cz_gen = "gen1"
cz_phen = "dissections"
if (island=="CZA") {
  dat_dir = "CZA_off_SW/CZA_off_raw_data_renamed"
} else {
  dat_dir = list.dirs(path = "CZD_off_SW/CZD_off_SW_Oct2016_rawdata")[grepl(pattern = cz_gen,
                                                                            x = list.dirs(path = "CZD_off_SW/CZD_off_SW_Oct2016_rawdata"))]
}
(phen_fl = list.files(path = dat_dir, full.names = TRUE)[grepl(pattern = cz_phen, x = list.files(path = dat_dir))])

# for CZA wild samples #
phen_dt <- read.csv(file = phen_fl)
head(phen_dt)
# go to "sex" function
########################

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
lapply(phen_dt, function(b) table(b[, 'broodpouch']))
lapply(phen_dt, function(b) table(b[, 'penis']))
lapply(phen_dt, function(b) table(b[, 'notes']))

#### identify sex of each snail, using brood pouch and penis data ####
# immatures are treated as juveniles
# table(phen_dt$notes)
sex <- function(b, p, isl) {
  if (isl=="CZD" & cz_gen=="gen0") {
    if(b=="Y" & (is.na(b)==F)) y <- "female"
    if(b=="N" & (is.na(b)==F)) y <- "NA"
    if((b %in% c("Y", "N"))==F) y <- "NA"
  } else {
    if(b=="Y" & p=="N" & (is.na(b)==F)) y <- "female"
    if(b=="N" & p=="Y" & (is.na(b)==F)) y <- "male"
    if(b=="N" & p=="small" & (is.na(b)==F)) y <- "juvenile"
    if(b=="N" & p=="N" & (is.na(b)==F) & (is.na(p)==F)) y <- "juvenile"
    if(b=="Y" & p=="Y" & (is.na(b)==F) & (is.na(p)==F)) y <- "hermaphrodite"
    if((b %in% c("Y", "N"))==F | (p %in% c("Y", "N", "small"))==F) y <- "NA"
  }
  return(y)
}
sex(b = "N", p = "N", isl = island)
sex(b = "Y", p = "N", isl = island)
sex(b = "N", p = "small", isl = island)
sex(b = "N", p = "Y", isl = island)
sex(b = "Y", p = "Y", isl = island)

# for CZA/D wild samples #
colnames(phen_dt)
dt_na = mutate(phen_dt, sex = apply(phen_dt[, c("brood", "penis")],
                                    MARGIN = 1, FUN = function(x) sex(b = x[1], p = x[2], isl = island)))
table(dt_na$sex)
dt_na <- dt_na[dt_na$sex!="NA", ]
table(dt_na$trematodes)
table(dt_na$notes)
dt_na[dt_na$notes==levels(dt_na$notes)[18], ]
dt_na$sex <- ifelse(test = dt_na$notes==levels(dt_na$notes)[8], yes = "juvenile", no = dt_na$sex)
table(dt_na$status_dissection)
setdiff(as.character(dt_na[dt_na$trematodes=="Y", "snail_ID"]), as.character(dt_na[dt_na$status_dissection=="parasitized", "snail_ID"]))
dt_nona = dt_na[which(dt_na[, "sex"] != "NA"), ]
table(dt_nona$trematodes)
table(dt_nona$status_dissection)
table(dt_nona$notes)
dt_nona[dt_nona$notes=="was indicated as male; sex changed later because there is a brood pouch photo",]
dt_nona$notes <- NULL
dt_nona = dt_nona[dt_nona[, "trematodes"]!="Y", ]
dt_nona = dt_nona[dt_nona[, "trematodes"]!="y", ]
dt_nona$sex_tube <- NULL
table(dt_nona$sex)
table(dt_nona$status_dissection)
setdiff(as.character(dt_nona[dt_nona$sex=="juvenile", "snail_ID"]),
        as.character(dt_nona[dt_nona$status_dissection=="juvenile", "snail_ID"]))
dt_nona$status_dissection <- NULL
dt_nona$sex_dissection <- NULL

(cdate = as.character(as.Date(Sys.time(), "%Y-%m-%d"), "%Y%m%d"))
nu_phen_fl = tools::file_path_sans_ext(phen_fl)
if (island=="CZD") {
  (nu_phen_dir = paste0(dirname(dirname(dirname(nu_phen_fl))), "/", island, "_off_SW_Oct2016_finaldata"))
} else {
  (nu_phen_dir = paste0(dirname(dirname(nu_phen_fl)), "/", island, "_off_final_data"))
}
dir.create(nu_phen_dir)
write.csv(dt_nona, paste0(nu_phen_dir, "/", basename(nu_phen_fl), "_", cdate, ".csv"), row.names = FALSE)
########################

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
      dt_nona = dt_nona[dt_nona[, "sexontube"]==FALSE, ]
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
lapply(phen_sex, colnames)
lapply(phen_sex, function(b) names(table(b[, 'notes'])))

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
# Change from juvenile to immature in excel using notes

# rm(list = ls())
rm(list = setdiff(ls(), c("dat_dir", "island", "nu_phen_dir")))
##########
## size ##
##########
# NO for CZA wild samples
one_phen = "size"
len_dt = read.xls(xls = "CZA_off_SW/CZA_off_raw_data_renamed/CZoff_length_size_sample1_201501.xlsx", sheet = 1)
size_dt = read.xls(xls = "CZA_off_SW/CZA_off_raw_data_renamed/CZoff_sample1_size.xlsx", sheet = 1)
if (identical(len_dt$Size.mm, size_dt$Size.mm)) {
  colnames(len_dt) = c("snail_ID", "size_mm", "notes", "date")
  write.csv(len_dt, file = "CZA_off_SW/CZA_off_final_data/CZA_off_offspring_100_size_mm.csv", quote = FALSE, row.names = FALSE)
}
#################################
######## merge phenotypes #######
# sex was added in excel for size
# and PCs dataset of generation 0
# parents
# juvenile was added in excel for
# size dataset of generation 1
# offspring 100
#################################
island = "CZA"
phenos = c("boldness", "dissections", "thickness", "weight", "size", "PCs")
dat_dir <- paste0(island, "_off_SW/", island, "_off_final_data/Littorina_", island)
# cz_gen = "gen1"
if (island=="CZA") {
  dat_dir = "CZA_off_SW/CZA_off_final_data"
} else {
  dat_dir = "CZD_off_SW/CZD_off_SW_Oct2016_finaldata"
}
phen_fl = lapply(1:length(phenos), function(x) {
  if (island=="off_D") {
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
lapply(phen_dt, length)
phen_all = lapply(1:length(phenos), function(x) {
  phen_merged = Reduce(function(...) merge(..., all=TRUE), phen_dt[[x]])
  phen_merged[, "snail_ID"] = as.character(phen_merged[, "snail_ID"])
  return(phen_merged)
})
# phen_all[[6]]
lapply(phen_all, head)
lapply(phen_all, colnames)
lapply(phen_all, nrow)
# head(phen_all[[4]])
phen_all_dt = Reduce(function(...) merge(..., all=TRUE, by="snail_ID"), phen_all)
head(phen_all_dt)
colnames(phen_all_dt)
colnames(phen_all_dt)[which(colnames(phen_all_dt)=="score")] = "bold_score"
sum(grepl(pattern = "snail", x = colnames(phen_all_dt)))
sum(grepl(pattern = "sex", x = colnames(phen_all_dt)))

# NO for CZA/D wild samples
sex_col = colnames(phen_all_dt)[(grepl(pattern = "sex", x = colnames(phen_all_dt)))]

sex_dt = phen_all_dt[, sex_col]
apply(X = sex_dt, MARGIN = 2, FUN = function(x) table(x))
# sex_dt$sex.x
sex_miss = apply(X = sex_dt, MARGIN = 2, FUN = function(x) {
  x = as.character(x)
  x[is.na(x)] <- "missing"
  return(x)
})
apply(X = sex_miss, MARGIN = 2, FUN = function(x) table(x))
# which(sex_miss[, 4]=="missing")
# which(sex_miss[, 3]=="missing")
# just one comparison between column 1 and column 2
comp_with = 1
sex_diff = setdiff(which(sex_miss[, 3]=="missing"), which(sex_miss[, comp_with]=="missing"))
sex_miss[sex_diff, 3] = sex_miss[sex_diff, comp_with]
phen_all_dt[which(sex_miss[, 3]=="missing"), "snail_ID"]
phen_all_dt[which(sex_miss[, 3]=="missing"), c("snail_ID", "broodpouch", "penis")]
phen_all_dt[duplicated(phen_all_dt$snail_ID), "snail_ID"] %in% phen_all_dt[which(sex_miss[, 3]=="missing"), "snail_ID"]
# c("R_321B", "R_56", "R_58") %in% phen_all_dt[which(sex_miss[, 3]=="missing"), "snail_ID"]
# all_miss = phen_all_dt[which(sex_miss[, 3]=="missing"), ]

phen_all_dt$sex = sex_miss[, 3]
phen_all_dt[phen_all_dt$sex=='missing', ]
length(unique(phen_all_dt$snail_ID))
(dup_id <- phen_all_dt[duplicated(as.character(phen_all_dt$snail_ID)), "snail_ID"])
phen_all_dt[phen_all_dt$snail_ID==dup_id[1], ]
phen_all_dt[401, ]
phen_all_dt <- phen_all_dt[-401, ]

phen_all_dt <- phen_all_dt[!is.na(phen_all_dt$snail_ID), ]

sum(grepl(pattern = "O", x = as.character(phen_all_dt$snail_ID)))

(cdate = as.character(as.Date(Sys.time(), "%Y-%m-%d"), "%Y%m%d"))
if (exists("cz_gen")) {
  write.csv(phen_all_dt, file = paste0(dat_dir, "/", island, "_off_", cz_gen,"_all_phenos_", cdate, ".csv"), row.names = FALSE)
} else {
  write.csv(phen_all_dt, file = paste0(dat_dir, "/", island, "_off_all_phenos_", cdate, ".csv"), row.names = FALSE)
}

# odate = "20200406"
# foo = read.csv(paste0(dat_dir, "/", island, "_off_all_phenos_main_", odate, ".csv"))
foo <- read.csv(file = "/Users/samuelperini/Documents/research/projects/1.cz_off/CZA_off_SW/CZA_off_final_data/CZA_all_phenos_main_20200511.csv")
# foo = read.csv("CZA_off_SW/CZA_off_final_data/CZA_off_all_phenos_main_20191216.csv")
# colnames(foo)
phen_all_dt$size_mm <- phen_all_dt$length_mm
(kcol = intersect(colnames(phen_all_dt), colnames(foo)))
mdt = phen_all_dt[, which(colnames(phen_all_dt) %in% kcol)]
# head(mdt)
# 
# (cdate = as.character(as.Date(Sys.time(), "%Y-%m-%d"), "%Y%m%d"))
if (exists("cz_gen")) {
  write.csv(mdt, file = paste0(dat_dir, "/", island, "_off_", cz_gen,"_all_phenos_main_", cdate, ".csv"), row.names = FALSE, quote = FALSE)
} else {
  write.csv(mdt, file = paste0(dat_dir, "/", island, "_off_all_phenos_main_", cdate, ".csv"), row.names = FALSE, quote = FALSE)
}
# for CZA wild samples
# change the name of the file manually

levels(mdt$size_mm)
is.numeric(mdt$size_mm)
# sum(dat_off$size_mm==" but that heavy...?", na.rm = TRUE)
# sum(dat_off$size_mm=="", na.rm = TRUE)
# which(dat_off$size_mm==" but that heavy...?")
# which(dat_off$size_mm=="")
# dat_off[1650, ]
# dat_off[1651, ]

rm(list = setdiff(ls(), c("cdate", "odate", "island", "dat_dir")))
odate = c("20200406", "20200511")
CZ = read.csv(paste0(dat_dir, "/", island, "_all_phenos_main_", odate[2], ".csv"))
OFF = read.csv(paste0(dat_dir, "/", island, "_off_all_phenos_main_", odate[1], ".csv"))
colnames(CZ)
CZ <- CZ[, order(colnames(CZ))]
OFF <- OFF[, order(colnames(OFF))]
colnames(OFF)
identical(colnames(OFF), colnames(CZ))
table(CZ$sex)
table(OFF$sex)
# OFF = separate(data = OFF, col = "snail_ID", into = c("pop", "ID"), sep = "_")
# OFF[, "generation"] = 1
# OFF[which(nchar(as.character(OFF$ID)) == 2), "generation"] = 0
# mdat <- rbind()

#######################
####### spatial #######
#######################
raw_dir <- paste0(island, "_off_SW/", island, "_off_raw_data")
spa_fl <- list.files(raw_dir, pattern = "spatial", full.names = TRUE)
CZ_spa <- read.csv(spa_fl[1])
head(CZ_spa)
head(CZ)
OFF_spa <- read.csv(spa_fl[2])
head(OFF_spa)
head(OFF)
OFF_spa$snail_ID <- paste(substr(OFF_spa$snail_ID, start = 1, stop = 1), substr(OFF_spa$snail_ID, start = 2, stop = 3),
                          sep = "_")
with(data = CZ_spa, plot(X, Y))
with(data = OFF_spa, plot(x, y))

library(RANN)
library(magrittr)
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
# closest <- nn2(data = CZ_spa[, c("X", "Y")], query = OFF_spa[, c("x", "y")], k = 1)
closest <- nn2(data = OFF_spa[, c("x", "y")], query = CZ_spa[, c("X", "Y")], k = 1)
closest <- sapply(closest, cbind) %>% as_tibble
closest$X <- OFF_spa[closest$nn.idx, "x"]
closest$Y <- OFF_spa[closest$nn.idx, "y"]
head(closest)
with(data = closest, plot(X, Y))
OFF_spa$LCmeanDist <- CZ_spa[closest$nn.idx, "LCmeanDist"]
head(OFF_spa)
OFF_spa = separate(data = OFF_spa, col = "snail_ID", into = c("pop", "ID"), sep = "_")
closest$pop <- OFF_spa$pop
ggplot(data = closest) +
  geom_point(aes(x = X, y = Y, col = factor(pop)))

aggregate(x = OFF_spa$LCmeanDist, by = list(pop = OFF_spa$pop), mean)


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

# mdt[which(mdt$snail_ID=="G_O7"), "snail_ID"] = "G_07"
# cdt[which(cdt$snail_ID=="G_O7"), "snail_ID"] = "G_07"
# merge(mdt, cdt, by = "snail_ID", all = TRUE)

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
library(ggplot2)
phenos_plot = lapply(seq_along(phenos), function(x) {
  ggplot(data = edt) +
    facet_wrap(~generation) +
    geom_point(aes(x = e_LCmeanDist, y = edt[,phenos[x]], col=factor(generation))) +
    labs(y = phenos[x]) +
    scale_color_manual(values = c("#f1a340", "#998ec3")) +
    theme(legend.position = "none",
          strip.text = element_text(size=13),
          strip.background = element_rect(fill="lightblue", colour="black", size=0.5),
          axis.title = element_text(size = 18),
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
ext = 'pdf'
lapply(seq_along(phenos), function(x) {
  ggsave(filename = paste0(dirname(dat_dir), "/", island, "_off_final_figures/", island, "_off_",
                           phenos[x], "_scatter_", cdate, ".", ext),
         plot = phenos_plot[[x]], device = ext, width = 8, height = 8)
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
          axis.title = element_text(size = 18),
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
                           col_pc[x], "_scatter_", cdate, ".", ext),
         plot = pca_plot[[x]], device = ext, width = 8, height = 8)
})

# there are no IDs starting from 200, only 300
pc2_sort = pca_edt[order(pca_edt$PC2, decreasing = TRUE),]
pc2_sort[pc2_sort$generation==1, c("ID", "PC2")][1:20, ]
tail(pc2_sort[pc2_sort$generation==1, c("ID", "PC2")])

# plot lab vs wild
table(edt$generation)
head(edt)

pheno_calc = function(datas, phen, grp) {
  nn_na = datas[!is.na(datas[, phen]), ]
  rawdt = aggregate(nn_na[, phen], by=list(nn_na[, "generation"], nn_na[, grp]), CI)[-1, ]
  # rawdt = aggregate(nn_na[, phen], by=list(nn_na[, "generation"], nn_na[, grp]), CI)
  # if (phen=="bold_score") {
  #   rawdt = rawdt[-1, ]
  # }
  return(rawdt)
}
phendt = pheno_calc(datas = edt, phen = "bold_score", grp = "e_LCmeanDist")
head(pca_edt)
lapply(seq_along(col_pc), function(x) {
  phendt = pheno_calc(datas = pca_edt, phen = col_pc[x], grp = "e_LCmeanDist")
})
ggplot() +
  geom_point(aes(x = phendt[phendt$Group.1==0, "x"][, "mean"], y = phendt[phendt$Group.1==1, "x"][, "mean"]))


# edt[!is.na(edt$bold_score),]
# aggregate(edt$bold_score, by=list(edt$e_LCmeanDist), CI)
# str(edt)
# table(edt$e_LCmeanDist)
# edt$e_LCmeanDist = round(edt$e_LCmeanDist, 2)
# edt[edt$e_LCmeanDist==247.71,]

phenos_plot = lapply(seq_along(phenos), function(x) {
  phendt = pheno_calc(datas = edt, phen = phenos[x], grp = "e_LCmeanDist")
  ggplot() +
    geom_abline(slope = 1, linetype="dashed") +
    geom_point(aes(x = phendt[phendt$Group.1==0, "x"][, "mean"], y = phendt[phendt$Group.1==1, "x"][, "mean"]),
               size=3) +
    geom_label_repel(aes(x = phendt[phendt$Group.1==0, "x"][, "mean"], y = phendt[phendt$Group.1==1, "x"][, "mean"],
                         label = LETTERS[2:18]),
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50') +
    labs(title = phenos[x], x = "wild sample", y = "lab-reared sample") +
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
})
phenos_plot
dir.create(path = paste0(dirname(dat_dir), "/", island, "_off_final_figures"))
(cdate = as.character(as.Date(Sys.time(), "%Y-%m-%d"), "%Y%m%d"))
ext = 'pdf'
lapply(seq_along(phenos), function(x) {
  ggsave(filename = paste0(dirname(dat_dir), "/", island, "_off_final_figures/", island, "_off_",
                           phenos[x], "_LABvsWILD_", cdate, ".", ext),
         plot = phenos_plot[[x]], device = ext, width = 8, height = 8)
})

pcs_plot = lapply(seq_along(col_pc), function(x) {
  phendt = pheno_calc(datas = pca_edt, phen = col_pc[x], grp = "e_LCmeanDist")
  ggplot() +
    geom_abline(slope = 1, linetype="dashed") +
    geom_point(aes(x = phendt[phendt$Group.1==0, "x"][, "mean"], y = phendt[phendt$Group.1==1, "x"][, "mean"]),
               size=3) +
    geom_label_repel(aes(x = phendt[phendt$Group.1==0, "x"][, "mean"], y = phendt[phendt$Group.1==1, "x"][, "mean"],
                         label = LETTERS[2:18]),
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50') +
    labs(title = col_pc[x], x = "wild sample", y = "lab-reared sample") +
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
})
pcs_plot
# dir.create(path = paste0(dirname(dat_dir), "/", island, "_off_final_figures"))
# (cdate = as.character(as.Date(Sys.time(), "%Y-%m-%d"), "%Y%m%d"))
# ext = 'pdf'
lapply(seq_along(col_pc), function(x) {
  ggsave(filename = paste0(dirname(dat_dir), "/", island, "_off_final_figures/", island, "_off_",
                           col_pc[x], "_LABvsWILD_", cdate, ".", ext),
         plot = pcs_plot[[x]], device = ext, width = 8, height = 8)
})

###########################
######### CZ data #########
###########################
island = "CZA"
dir.create(paste0(island, "_wild_SW"))
dir.create(paste0(island, "_wild_SW/", island, "_wild_raw_data"))
phenos = c("boldness", "shape", "spatial", "thickness", "weight")
# cz_gen = "gen1"
if (island=="CZA") {
  dat_dir = "CZA_off_SW/CZA_wild_raw_data"
} else {
  dat_dir = "CZD_off_SW/CZD_wild_SW_Oct2016_rawdata"
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
colnames(phen_all_dt)[which(colnames(phen_all_dt)=="log_mean_bold")] = "bold_score"
sum(grepl(pattern = "snail_ID", x = colnames(phen_all_dt)))
phen_all_dt = mutate(phen_all_dt,
                     mean_thickness = rowMeans(phen_all_dt[, grepl(pattern = "thickness", x = colnames(phen_all_dt))],
                                               na.rm = TRUE))

(cdate = as.character(as.Date(Sys.time(), "%Y-%m-%d"), "%Y%m%d"))
# getwd()
dir.create(paste0(island, "_wild_SW/", island, "_wild_final_data"))
findat_dir = paste0(island, "_wild_SW/", island, "_wild_final_data")
if (exists("cz_gen")) {
  write.csv(phen_all_dt, file = paste0(dat_dir, "/", island, "_off_", cz_gen,"_all_phenos_", cdate, ".csv"), row.names = FALSE)
} else {
  write.csv(phen_all_dt, file = paste0(findat_dir, "/", island, "_wild_all_phenos_", cdate, ".csv"), row.names = FALSE)
}

phen_all_dt = read.csv("CZA_wild_SW/CZA_wild_final_data/CZA_wild_all_phenos_20200113.csv")
summary(phen_all_dt)
phen_all_dt$LCmeanDist = round(phen_all_dt$LCmeanDist, 1)

wild_phenos = colnames(phen_all_dt)[grepl(pattern = "score|PC|mean_thick|weight", x = colnames(phen_all_dt))]
# colnames(edt)[grepl(pattern = "score|PC|mean_thick|weight", x = colnames(edt))]


edt$e_LCmeanDist = round(edt$e_LCmeanDist)
table(edt$e_LCmeanDist)
table(phen_all_dt$LCmeanDist)
intersect(phen_all_dt$LCmeanDist, edt$e_LCmeanDist)

range(edt$e_LCmeanDist)
unique(edt$pop)
unique(edt$e_LCmeanDist)
length(unique(edt$e_LCmeanDist))
dist_brk = seq(min(edt$e_LCmeanDist)-1, to = max(edt$e_LCmeanDist)+1, length.out = length(unique(edt$pop)))
length(dist_brk)
colnames(phen_all_dt)
wild_edt = phen_all_dt[, c(wild_phenos, "LCmeanDist")]

# wild_edt$dist_bin = cut(wild_edt$LCmeanDist, dist_brk)

wild_edt$dist_bin = cut(wild_edt$LCmeanDist, unique(edt$e_LCmeanDist), include.lowest = TRUE)

wild_edt_bin = wild_edt[!is.na(wild_edt$dist_bin), ]
aggregate(x = wild_edt_bin, by = list(wild_edt_bin$dist_bin), FUN = mean)

pheno_calc = function(datas, phen, grp) {
  nn_na = datas[!is.na(datas[, phen]), ]
  # rawdt = aggregate(nn_na[, phen], by=list(nn_na[, "generation"], nn_na[, grp]), CI)[-1, ]
  rawdt = aggregate(nn_na[, phen], by=list(nn_na[, grp]), CI)
  # if (phen=="bold_score") {
  #   rawdt = rawdt[-1, ]
  # }
  return(rawdt)
}
table(wild_edt_bin$dist_bin)
head(edt)
unique(edt[order(edt$e_LCmeanDist), "e_LCmeanDist"])
edt = edt[order(edt$e_LCmeanDist), ]

pheno_calc(datas = wild_edt_bin, phen = wild_phenos[5], grp = "dist_bin")
genx = 0
one_gen_edt = edt[edt$generation==genx, ]
one_gen_edt$dist_bin = cut(one_gen_edt$e_LCmeanDist, unique(edt$e_LCmeanDist), include.lowest = TRUE)
one_gen_phenos = colnames(one_gen_edt)[grepl(pattern = "score|PC|mean_thick|weight", x = colnames(one_gen_edt))]
wild_phenos = intersect(one_gen_phenos, wild_phenos)
# table(one_gen_edt$dist_bin)
# unique(one_gen_edt[, c("pop", "dist_bin")])
unique(wild_edt_bin$dist_bin)


lapply(seq_along(wild_phenos), function(x) {
  one_gen_bin = pheno_calc(datas = one_gen_edt, phen = wild_phenos[x], grp = "dist_bin")
  one_phen_bin = pheno_calc(datas = wild_edt_bin, phen = wild_phenos[x], grp = "dist_bin")
  # one_phen_bin$e_LCmeanDist = sort(unique(edt$e_LCmeanDist))
  ggplot() +
    geom_abline(slope = 1, linetype = 'dashed') +
    geom_point(aes(x = one_gen_bin[, 'x'][,'mean'], y = one_phen_bin[, 'x'][,'mean'])) +
    labs(x = "lab-reared generation 0", y = paste0("wild ", island, " sample"), title = wild_phenos[x]) +
    theme(plot.title = element_text(hjust = 0.5))
})

one_gen_pl = ggplot(data = one_gen_edt, aes(x = e_LCmeanDist, y = weight_g)) +
  geom_point()
wild_edt_pl = ggplot(data = wild_edt_bin, aes(x = LCmeanDist, y = weight_g)) +
  geom_point()
one_gen_pl + wild_edt_pl



