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
              help="which phenotype to analyse: boldness, dissections, weight, thickness", metavar="character"),
  make_option(c("-E", "--ecotype"), type="character", default=NULL,
              help="select either crab or wave", metavar="character"),
  make_option(c("-I", "--island"), type="character", default=NULL,
              help="select one or more islands: CZA, CZB and/or CZD", metavar="character"),
  make_option(c("-L", "--linkgrp"), type="integer", default=NULL,
              help="choose a linkage group for testing", metavar="integer"))

opt_parser = OptionParser(option_list=option_list,
                          description = "",
                          epilogue = "Example: Rscript")
opt = parse_args(opt_parser)

if (is.null(opt$vcf) | is.null(opt$directory) | is.null(opt$ecotype) | is.null(opt$island)){
  print_help(opt_parser)
  stop("All arguments must be supplied.\n", call.=FALSE)
}

island = "CZA"
cz_gen = "gen1"
cz_phen = "boldness"
if (island=="CZA") {
  dat_dir = "CZA_off_SW/CZA_off_raw_data_renamed"
} else {
  dat_dir = list.dirs(path = "CZD_off_SW/CZD_off_SW_Oct2016_rawdata")[grepl(pattern = cz_gen,
                                                                            x = list.dirs(path = "CZD_off_SW/CZD_off_SW_Oct2016_rawdata"))]
}
phen_fl = list.files(path = dat_dir, full.names = TRUE)[grepl(pattern = cz_phen, x = list.files(path = dat_dir))]

#########################
######## boldness #######
#########################
phen_dt = lapply(1:length(phen_fl), function(x) {
  onedt = read.xls(xls = phen_fl[x], sheet = 1)
  coldt = colnames(onedt)[grepl(pattern = "test", x = colnames(onedt))]
  onedtcol = onedt[, which(colnames(onedt) %in% coldt)]
  cbind(snail_ID = onedt[, 1], onedtcol[,1:12])
})
lapply(phen_dt, head)
lapply(phen_dt, ncol)
lapply(phen_dt, colnames)


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

bold$score <- (bold$crawl_mean+bold$fullyout_mean)/2

head(bold)
# check distribution (if needed)
hist(bold$crawl_mean, nclass = 21)
hist(bold$fullyout_mean, nclass = 21) 
hist(bold$score,nclass=21)

# write result to file
(nu_phen_fl = tools::file_path_sans_ext(phen_fl[dtnum]))
(nu_phen_dir = paste0(dirname(dirname(nu_phen_fl)), "/", island, "_off_final_data"))
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
# cz_gen = "gen1"
cz_phen = "thickness"
if (island=="CZA") {
  dat_dir = "CZA_off_SW/CZA_off_raw_data_renamed"
} else {
  dat_dir = list.dirs(path = "CZD_off_SW/CZD_off_SW_Oct2016_rawdata")[grepl(pattern = cz_gen,
                                                                            x = list.dirs(path = "CZD_off_SW/CZD_off_SW_Oct2016_rawdata"))]
}
phen_fl = list.files(path = dat_dir, full.names = TRUE)[grepl(pattern = cz_phen, x = list.files(path = dat_dir))]
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
  nu_phen_dir = paste0(dirname(dirname(nu_phen_fl)), "/", island, "_off_final_data")
  dir.create(nu_phen_dir)
  write.csv(phen_dt_thick[[x]], paste0(nu_phen_dir, "/", basename(nu_phen_fl), "_", cdate, ".csv"), row.names = FALSE)
})

# rm(list = ls())
#######################
######## weight #######
#######################
island = "CZA"
# cz_gen = "gen1"
cz_phen = "weight"
if (island=="CZA") {
  dat_dir = "CZA_off_SW/CZA_off_raw_data_renamed"
} else {
  dat_dir = list.dirs(path = "CZD_off_SW/CZD_off_SW_Oct2016_rawdata")[grepl(pattern = cz_gen,
                                                                            x = list.dirs(path = "CZD_off_SW/CZD_off_SW_Oct2016_rawdata"))]
}
phen_fl = list.files(path = dat_dir, full.names = TRUE)[grepl(pattern = cz_phen, x = list.files(path = dat_dir))]
phen_dt = lapply(1:length(phen_fl), function(x) {
  read.xls(xls = phen_fl[x], sheet = 1)[, 1:5]
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
  nu_phen_dir = paste0(dirname(dirname(nu_phen_fl)), "/", island, "_off_final_data")
  dir.create(nu_phen_dir)
  write.csv(phen_dt[[x]], paste0(nu_phen_dir, "/", basename(nu_phen_fl), "_", cdate, ".csv"), row.names = FALSE)
})

# rm(list = ls())
############################
######## dissections #######
############################
island = "CZA"
# cz_gen = "gen1"
cz_phen = "dissections"
if (island=="CZA") {
  dat_dir = "CZA_off_SW/CZA_off_raw_data_renamed"
} else {
  dat_dir = list.dirs(path = "CZD_off_SW/CZD_off_SW_Oct2016_rawdata")[grepl(pattern = cz_gen,
                                                                            x = list.dirs(path = "CZD_off_SW/CZD_off_SW_Oct2016_rawdata"))]
}
phen_fl = list.files(path = dat_dir, full.names = TRUE)[grepl(pattern = cz_phen, x = list.files(path = dat_dir))]
phen_dt = lapply(1:length(phen_fl), function(x) {
  orig = read.xls(xls = phen_fl[x], sheet = 1)
  if (x==1 | x==2) {
    colnames(orig) = c("snail_ID", "sexontube", "broodpouch", "penis", "parasites", "ciliates", "experimenter", "notes")
  } else {
    colnames(orig) = c("snail_ID", "broodpouch", "penis", "lengthglandrow", "lengthtip", "nglands", "experimenter", "notes")
  }
  return(orig)
})
lapply(phen_dt, head)
lapply(phen_dt, ncol)
lapply(phen_dt, colnames)
lapply(phen_dt, summary)
lapply(phen_dt, str)
#### identify sex of each snail, using brood pouch and penis data ####
sex <- function(b, p) {
  if(b=="Y" & p=="N" & (is.na(b)==F)) y <- "female"
  if(b=="N" & p=="Y" & (is.na(b)==F)) y <- "male"
  if((b %in% c("Y", "N"))==F | (p %in% c("Y", "N"))==F | b==p) y<-"NA"
  return(y)
}
phen_sex = lapply(1:length(phen_fl), function(i) {
  dt_na = mutate(phen_dt[[i]], sex=apply(phen_dt[[i]][, c("broodpouch", "penis")],
                                         MARGIN = 1, FUN = function(x) sex(x[1], x[2])))
  dt_nona = dt_na[which(dt_na[, "sex"] != "NA"), ]
  if (i==1) {
    dt_nona = dt_nona[dt_nona[, "broodpouch"]=="Y", ]
    dt_nona = dt_nona[dt_nona[, "parasites"]!="Y", ]
  } else if (i==2) {
    dt_nona = dt_nona[dt_nona[, "sexontube"]=="M", ]
    dt_nona = dt_nona[dt_nona[, "parasites"]!="Y", ]
  }
  return(dt_nona)
})
lapply(phen_sex, head)
lapply(phen_sex, summary)

# write result to file
(cdate = as.character(as.Date(Sys.time(), "%Y-%m-%d"), "%Y%m%d"))
lapply(1:length(phen_fl), function(x) {
  nu_phen_fl = tools::file_path_sans_ext(phen_fl[x])
  nu_phen_dir = paste0(dirname(dirname(nu_phen_fl)), "/", island, "_off_final_data")
  dir.create(nu_phen_dir)
  write.csv(phen_sex[[x]], paste0(nu_phen_dir, "/", basename(nu_phen_fl), "_", cdate, ".csv"), row.names = FALSE)
})

# rm(list = ls())
#################################
######## merge phenotypes #######
#################################
island = "CZA"
phenos = c("boldness", "dissections", "thickness", "wet_weight")
# cz_gen = "gen1"
if (island=="CZA") {
  dat_dir = "CZA_off_SW/CZA_off_final_data"
} else {
  dat_dir = list.dirs(path = "CZD_off_SW/CZD_off_SW_Oct2016_finaldata")[grepl(pattern = cz_gen,
                                                                              x = list.dirs(path = "CZD_off_SW/CZD_off_SW_Oct2016_finaldata"))]
}
phen_fl = lapply(1:length(phenos), function(x) {
  list.files(path = dat_dir, pattern = phenos[x], full.names = TRUE)
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
head(phen_all[[4]])
phen_all_dt = Reduce(function(...) merge(..., all=TRUE, by="snail_ID"), phen_all)
colnames(phen_all_dt)
sum(grepl(pattern = "snail_ID", x = colnames(phen_all_dt)))

(cdate = as.character(as.Date(Sys.time(), "%Y-%m-%d"), "%Y%m%d"))
write.csv(phen_all_dt, file = paste0(dat_dir, "/CZA_off_all_phenos_", cdate, ".csv"), row.names = FALSE)


