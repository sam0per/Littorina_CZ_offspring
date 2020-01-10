rm (list=ls())

CZoff_P_space = read.csv("CZA_off_SW/CZA_off_raw_data_renamed/CZoff_spatial_2014.csv")
head(CZoff_P_space)
tail(CZoff_P_space)
colnames(CZoff_P_space) = c("snail_ID", "X", "Y", "Z", "Env")
# CZoff_P_shape = read.csv("results/CZ-off-shape-sex.csv")
# head(CZoff_P_shape)
# str(CZoff_P_shape)

island = "CZA"
if (island=="CZA") {
  dat_dir = "CZA_off_SW/CZA_off_final_data"
} else {
  dat_dir = "CZD_off_SW/CZD_off_SW_Oct2016_finaldata"
}
mat_dir = "/Users/samuelperini/Documents/research/projects/2.mating/data/"
CZA_space = read.csv(paste0(mat_dir, island, "_spatial_LCP_20190121.csv"))
head(CZA_space)
# CZA_bez = read.table("data/Littorina_CZ_final_data/CZA_bezierPath_20150617.csv", sep=",", stringsAsFactors=T, header=T)
# head(CZA_bez)
# length(unique(CZA_bez$pathX))
# CZE = read.csv("data/Littorina_CZoff_raw_data_renamed/CZE_Snail_habitat_SN_2m.csv")

library(ggplot2)
library(tidyverse)
ggplot(CZA_space, aes(X,Y)) +
  geom_point() +
  geom_point(data = CZoff_P_space, aes(X,Y,col='off'))

range(CZoff_P_space$X)
CZA_space_red = CZA_space[CZA_space$x>=min(CZoff_P_space$x) & CZA_space$x<=max(CZoff_P_space$x) &
                            CZA_space$y>=min(CZoff_P_space$y) & CZA_space$y<=max(CZoff_P_space$y) &
                            !is.na(CZA_space$DistAlongBezPath),]
CZA_bez_red = CZA_bez[CZA_bez$pathX>=min(CZoff_P_space$x)-5 & CZA_bez$pathX<=max(CZoff_P_space$x)+3 &
                        CZA_bez$pathY>=min(CZoff_P_space$y)-5 & CZA_bez$pathY<=max(CZoff_P_space$y),]

CZoff_P_space = CZoff_P_space %>%
  separate(snail_ID,into = c("pop","ID"),sep = 1)

ggplot(CZA_space_red,aes(x,y,fill='CZA')) +
  geom_point() +
  geom_point(data = CZoff_P_space,aes(x,y,col=pop)) +
  #scale_colour_manual(values = c('black','red')) +
  labs(col='',fill='')

ggplot(CZA_space_red,aes(BezX,BezY,fill='CZA')) +
  geom_point() +
  geom_point(data = CZoff_P_space,aes(x,y,col=pop)) +
  #scale_colour_manual(values = c('black','red')) +
  labs(col='',fill='')

ggplot(CZA_bez_red,aes(pathX,pathY,fill='CZA')) +
  geom_point() +
  geom_point(data = CZoff_P_space,aes(x,y,col=pop)) +
  #scale_colour_manual(values = c('black','red')) +
  labs(col='',fill='')

# CZoff_P_space$YX = CZoff_P_space$y - CZoff_P_space$x
# CZA_bez_red$YX = CZA_bez_red$pathY - CZA_bez_red$pathX

head(CZoff_P_space)
head(CZA_space)
# head(CZA_bez_red)


library(RANN)
# closest <- nn2(CZA_bez_red[,1:2],CZoff_P_space[, 2:3], 1)
closest <- nn2(CZA_space[,3:5], CZoff_P_space[, 2:4], 1)
sapply(closest, head)
# CZA_space[338,]
# CZoff_P_space[338,]
CZoff_P_space[, paste0(island, "_LCmeanDist")] = CZA_space[closest$nn.idx, "LCmeanDist"]
head(CZoff_P_space)
CZoff_P_space$snail_ID = paste0(substr(x = CZoff_P_space$snail_ID, start = 1, stop = 1), "_",
                                substr(x = CZoff_P_space$snail_ID, start = 2, stop = 3))
(cdate = as.character(as.Date(Sys.time(), "%Y-%m-%d"), "%Y%m%d"))
if (exists("cz_gen")) {
  write.csv(CZoff_P_space, file = paste0(dat_dir, "/", island, "_off_", cz_gen,"_all_phenos_", cdate, ".csv"), row.names = FALSE)
} else {
  write.csv(CZoff_P_space, file = paste0(dat_dir, "/", island, "_off_spatial_", cdate, ".csv"), row.names = FALSE)
}


# CZoff_P_space$closest.id <- closest$nn.idx
# CZoff_bez = CZA_bez_red[CZoff_P_space$closest.id,]
# rownames(CZoff_bez) = seq(1:nrow(CZoff_bez))
# head(CZoff_bez)
# 
# CZoff_P_space = merge(CZoff_P_space,CZoff_bez,by = "row.names",all = TRUE)
# head(CZoff_P_space)
CZoff_P_space %>%
  separate(snail_ID,into = c("pop","ID"),sep = 1) %>% select(-one_of("Row.names")) %>%
  ggplot() +
  geom_point(data = CZA_bez_red,aes(pathX,pathY,fill='CZA'),alpha=0.1) +
  geom_point(aes(pathX,pathY,col=pop),size=2.5) +
  #geom_point(aes(x,y,col=pop)) +
  labs(fill='')


########################################################################
########################################################################
########################################################################

BezX = lapply(CZoff_P_space$x,function(x) min(abs(x - CZA_bez_red$pathX)))
CZoff_P_space$pathX = CZoff_P_space$x - unlist(BezX)
CZoff_P_space$pathX2 = CZoff_P_space$x + unlist(BezX)

BezX = lapply(CZoff_P_space$x,function(x) min(abs(x - CZA_space_red$x)))
CZoff_P_space$x = CZoff_P_space$x - unlist(BezX)
head(CZoff_P_space)
head(CZA_space_red)
CZoff_P_space$X2 = CZoff_P_space$x + unlist(BezX)

CZoff_P_space %>% left_join(CZA_space_red, by = "x") %>% mutate(x = ifelse(is.na(y.y),X2,x)) %>%
  select(-(4:29)) %>%
  inner_join(CZA_space_red, by = "x") %>%
  ggplot() +
  geom_point(aes(BezX,BezY,col=pop),size=2.5) +
  geom_point(data = CZoff_P_space,aes(x,y,col=pop)) +
  geom_point(data = CZA_bez_red,aes(pathX,pathY,fill='CZA'),alpha=0.1) +
  labs(fill='')


lapply(CZoff_P_space,function(x) print(x[[1]]))

BezY = lapply(CZoff_P_space$y,function(y) min(abs(y - CZA_bez_red$pathY)))
CZoff_P_space$BezY = unlist(BezY) - CZoff_P_space$y

ggplot(CZoff_P_space,aes(BezX,BezY,col='bez')) +
  geom_point() +
  geom_point(aes(x,y,col='coords'))

range(CZA_space_red$x)
range(CZoff_P_space$x)
x.breaks = c(seq(min(as.integer(CZoff_P_space$x)),max(as.integer(CZoff_P_space$x)),1),2050)
CZA_space_red$x.bin = cut(CZA_space_red$x,x.breaks)

range(CZA_space_red$y)
range(CZoff_P_space$y)
y.breaks = c(seq(min(as.integer(CZoff_P_space$y)),max(as.integer(CZoff_P_space$y)),4),2049)
CZA_space_red$y.bin = cut(CZA_space_red$y,y.breaks)


CZA_space_bin = CZA_space_red %>%
  group_by(x.bin,y.bin) %>%
  summarise(x.mean=mean(x),y.mean=mean(y),count=n(),dist.mean=mean(DistAlongBezPath))


#coors$snail_ID=as.character(coors$snail_ID)
CZoff_P_space = CZoff_P_space %>% inner_join(CZoff_P_shape)

ggplot(CZA_space_red,aes(x,y,col='CZA')) +
  geom_point() +
  geom_point(data = CZoff_P_space,aes(x,y,col='OFF')) +
  scale_colour_manual(values = c('black','red')) +
  labs(col='')

CZoff_P_space$x.bin = cut(CZoff_P_space$x,x.breaks)
CZoff_P_space$y.bin = cut(CZoff_P_space$y,y.breaks)
CZoff_P_bin = CZoff_P_space %>%
  separate(snail_ID,into = c("pop","ID"),sep = 1) %>%
  group_by(pop,x.bin,y.bin,sex) %>%
  summarise(x.mean=mean(x),y.mean=mean(y),count=n(),PC1.mean=mean(PC1),PC2.mean=mean(PC2),CS.mean=mean(CS))

ggplot(CZA_space_bin,aes(x.mean,y.mean,col='CZA')) +
  geom_point() +
  geom_point(data = CZoff_P_bin,aes(x.mean,y.mean,col='OFF')) +
  scale_colour_manual(values = c('black','red')) +
  labs(col='')


CZoff_P_dist = CZA_space_bin %>% inner_join(CZoff_P_bin,by = c("x.bin","y.bin"))
head(CZoff_P_dist)

ggplot(CZoff_P_dist,aes(dist.mean,CS.mean,col=pop,shape=sex))+
  geom_point(size=2)

ggplot(CZoff_P_dist,aes(x.mean.y,y.mean.y,col=pop,shape=sex))+
  geom_point(size=2)



########################################################################
########################################################################
########################################################################




# based on CZ_spatial.R

#setwd("D:/Dropbox/Littorina_CZ/Littorina_CZ_data/Littorina_CZ_scripts")
#zone = "CZE"
zone = "CZoff"
getwd()
########################
#### GET SNAIL DATA ####

# get file names
apertures_file = list.files("data/Littorina_CZoff_raw_data_renamed", pattern = "apertures")
boldness_file = list.files("data/Littorina_CZoff_raw_data_renamed", pattern = "boldness")
dissections_file = list.files("data/Littorina_CZoff_raw_data_renamed", pattern = "dissections")
egg_counts_file = list.files("data/Littorina_CZoff_raw_data_renamed", pattern = "egg_counts")
egg_measurements_file = list.files("data/Littorina_CZoff_raw_data_renamed", pattern = "measurements")
foot_area_measurement_file = list.files("data/Littorina_CZoff_raw_data_renamed", pattern = "foot_area_measurement")
length_file = list.files("data/Littorina_CZoff_raw_data_renamed", pattern = "length")
mating_trials_file = list.files("data/Littorina_CZoff_raw_data_renamed", pattern = "mating_trials")
SnailDataLitPath_file = list.files("data/Littorina_CZoff_raw_data_renamed", pattern = "SnailDataLitPath")
surroundings_file = list.files(paste("../Littorina_CZ_", zone, "_final_data", sep=""), pattern = "surroundings")
thickness_colour_scars_file = list.files(paste("../Littorina_CZ_", zone, "_final_data", sep=""), pattern = "thickness_colour_scars")
wet_weight_file = list.files(paste("../Littorina_CZ_", zone, "_final_data", sep=""), pattern = "wet_weight")
shape_file = list.files(paste("../Littorina_CZ_", zone, "_final_data", sep=""), pattern = "shape")

# get phenotypic data
apertures = read.table(paste("../Littorina_CZ_", zone, "_final_data/", apertures_file, sep=""), header=T, sep=",", stringsAsFactors=F)
boldness = read.table(paste("../Littorina_CZ_", zone, "_final_data/", boldness_file, sep=""), header=T, sep=",", stringsAsFactors=F)
dissections = read.table(paste("../Littorina_CZ_", zone, "_final_data/", dissections_file, sep=""), header=T, sep=",", stringsAsFactors=F)
egg_counts = read.table(paste("../Littorina_CZ_", zone, "_final_data/", egg_counts_file, sep=""), header=T, sep=",", stringsAsFactors=F)
egg_measurements = read.table(paste("../Littorina_CZ_", zone, "_final_data/", egg_measurements_file, sep=""), header=T, sep=",", stringsAsFactors=F)
len = read.table(paste("../Littorina_CZ_", zone, "_final_data/", length_file, sep=""), header=T, sep=",", stringsAsFactors=F)
mating_trials = read.table(paste("../Littorina_CZ_", zone, "_final_data/", mating_trials_file, sep=""), header=T, sep=",", stringsAsFactors=F)
SnailDataLitPath = read.table(paste("../Littorina_CZ_", zone, "_final_data/", SnailDataLitPath_file, sep=""), header=T, sep=",", stringsAsFactors=F)
surroundings = read.table(paste("../Littorina_CZ_", zone, "_final_data/", surroundings_file, sep=""), header=T, sep=",", stringsAsFactors=F)
thickness_colour_scars = read.table(paste("../Littorina_CZ_", zone, "_final_data/", thickness_colour_scars_file, sep=""), header=T, sep=",", stringsAsFactors=F)
wet_weight = read.table(paste("../Littorina_CZ_", zone, "_final_data/", wet_weight_file, sep=""), header=T, sep=",", stringsAsFactors=F)
shape = read.table(paste("../Littorina_CZ_", zone, "_final_data/", shape_file, sep=""), header=T, sep=",", stringsAsFactors=F)


# make list of all available data frames
list1 = list()
x = 1
for (item in c("apertures", "boldness", "dissections", "egg_counts", "egg_measurements",
  "len", "SnailDataLitPath", "surroundings", "thickness_colour_scars", "wet_weight", "shape")){
  if (exists(item)==T){    
       list1[[x]] = get(item)
       x = x+1
  }
}

# merge all available data frames
merged = Reduce(function(...) merge(..., by="snail_ID", all=T), list1)

# reduce data frame by keeping only variables shared among all sites
keep = names(merged)[names(merged) %in% c("snail_ID",
                                          "aperture",
                                          "thickness1_mm", 
                                          "thickness2_mm",
                                          "thickness3_mm",
                                          "scars",
                                          "ridges",
                                          "colour_category",
                                          "pattern_category",
                                          "x",
                                          "y",
                                          "z",
                                          "DistAlongPath",
                                          "ShortestDistToBezierPath",
                                          "BezX",
                                          "BezY",
                                          "BezZ",
                                          "ResidX",
                                          "ResidY",
                                          "ResidZ",
                                          "length_mm",
                                          "BP_rock_1m",
                                          "BP_barnacle_1m",
                                          "BP_fucus_1m",
                                          "BP_rock_2m",
                                          "BP_barnacle_2m",
                                          "BP_fucus_2m",
                                          "BP_rock_4m",
                                          "BP_barnacle_4m",
                                          "BP_fucus_4m",
                                          "SN_rock_1m",
                                          "SN_barnacle_1m",
                                          "SN_fucus_1m",
                                          "SN_rock_2m",
                                          "SN_barnacle_2m",
                                          "SN_fucus_2m",
                                          "SN_rock_4m",
                                          "SN_barnacle_4m",
                                          "SN_fucus_4m",
                                          "foot_cm2_1",
                                          "foot_cm2_2",
                                          "foot_cm2_3",
                                          "X1_egg_normal",
                                          "X2_early_veliger_normal",
                                          "X3_veliger_normal",
                                          "X4_hatched_normal",
                                          "X5_egg_misdev",
                                          "X6_early_veliger_misdev",
                                          "X7_veliger_misdev",
                                          "X8_hatched_misdev",
                                          "X9_egg_twins",
                                          "X10_early_veliger_twins",
                                          "X11_veliger_twins",
                                          "ciliates",
                                          "brood" ,
                                          "penis",
                                          "trematodes",
                                          "ciliates_diss",
                                          "MLE_bold_mean_ln",
                                          "weight_g",
                                          "egg_width_av",
                                          "early_veliger_width_av",
                                          "veliger_width_av",
                                          "hatched_width_av",
                                          "RW1",
                                          "RW2",
                                          "PC1",
                                          "PC2",
                                          "CS")]
                                          
merged = merged[, keep]
merged$thickness_av = (merged$thickness1_mm + merged$thickness2_mm + merged$thickness3_mm ) / 3
merged$X12_total_offspring = apply(merged[, c("X1_egg_normal", "X2_early_veliger_normal", "X3_veliger_normal", "X4_hatched_normal", "X5_egg_misdev", "X6_early_veliger_misdev", "X7_veliger_misdev", "X8_hatched_misdev", "X9_egg_twins", "X10_early_veliger_twins", "X11_veliger_twins")],
                             1, function (x) sum(x, na.rm=T))
mature = merged[(merged$brood=="Y" | merged$penis=="Y") & is.na(merged$brood)==F & is.na(merged$penis)==F, ]
mature = mature[(mature$snail_ID %in% c("CZA486", "CZA509"))==F, ] # exclude snails with strange aperture/length ratio


##########################
#### GET SPATIAL DATA ####

spatial_file = list.files(paste("../Littorina_CZ_", zone, "_final_data", sep=""), pattern = "spatial")
spatial = read.table(paste("../Littorina_CZ_", zone, "_final_data/", spatial_file, sep=""), header=T, sep=",", stringsAsFactors = F)
mature = merge(mature, spatial, all=F, by="snail_ID")

# make habitat df: R/S: rock/stone; B/N: barnacles/none; F/N: Fucus/none
habitat_categories = c("RBF", "RBN", "RNF", "RNN", "SBF", "SBN", "SNF", "SNN")
habitat = spatial[spatial$type %in% habitat_categories, ]
habitat_sep = data.frame(do.call('rbind', strsplit(as.character(habitat$type), "", fixed=T)))
habitat$rock = habitat_sep$X1
habitat$barnacle = habitat_sep$X2
habitat$fucus = habitat_sep$X3

# plot habitat
habitat$col = NA
habitat$col[habitat$rock=="S"] = rgb(211, 242, 198, maxColorValue=255)
habitat$col[habitat$rock=="R"] = rgb(22, 162, 135, maxColorValue=255)
plot(habitat$x, habitat$y, col=habitat$col, pch=16, xlab="x coordinate (m)", ylab="y coordinate (m)")


##########################
#### COLOUR BY TRAIT #####
rbPal = colorRampPalette(c('blue', 'yellow'))
mature$col_aperture = rbPal(10)[as.numeric(cut(mature$aperture, breaks=10))]
mature$col_aperture_corr = rbPal(10)[as.numeric(cut(mature$aperture/(mature$length_mm^2), breaks=10))]
mature$col_len = rbPal(10)[as.numeric(cut(mature$length_mm, breaks=10))]
mature$col_weight = rbPal(10)[as.numeric(cut(mature$weight_g, breaks=10))]
mature$col_thickness = rbPal(10)[as.numeric(cut(mature$thickness_av, breaks=10))]
mature$col_boldness = rbPal(10)[as.numeric(cut(mature$MLE_bold_mean_ln, breaks=10))]
mature$col_RW1 = rbPal(10)[as.numeric(cut(mature$PC1, breaks=10))]


################################
### PLOT FOR ALL MATURE INDS ###

points(mature$x, mature$y, col=mature$col_aperture, pch=16)
points(mature$x, mature$y, col=mature$col_aperture_corr, pch=16)
points(mature$x, mature$y, col=mature$col_len, pch=16)
points(mature$x, mature$y, col=mature$col_weight, pch=16)
points(mature$x, mature$y, col=mature$col_thickness, pch=16)
points(mature$x, mature$y, col=mature$col_boldness, pch=16)
points(mature$x, mature$y, col=mature$col_veliger_width_av, pch=16)
points(mature$x, mature$y, col=mature$col_ciliates, pch=16)
points(mature$x, mature$y, col=mature$col_total_offspring, pch=16)


## ~~~ TEMP CZD
use = mature[mature$y>2005 & is.na(mature$thickness1_mm)==F & is.na(mature$weight_g)==F &
                          is.na(mature$PC1)==F & is.na(mature$colour_category)==F, ]

# remove some females from W to get a total of 120 snails
rem = c("CZE466", "CZE481", "CZE425", "CZE495", "CZE487",
        "CZE480", "CZE438", "CZE492", "CZE462", "CZE411", "CZE444")

use = use[(use$snail_ID%in%rem)==F, ]
dim(use)
dim(use[use$brood=="Y",])
dim(use[use$penis=="Y",])
points(use$x[use$brood=="Y"], use$y[use$brood=="Y"], col="blue", pch=16, cex=0.6)
points(use$x[use$penis=="Y"], use$y[use$penis=="Y"], col="red", pch=16, cex=0.6)

write.table(use, "CZE_extract.txt", append=F, quote=F, row.names=F, col.names=T)
### ~ temp end
         

#############################
### PLOT FOR FEMALES ONLY ###

females = mature[mature$brood=="Y", ]
females$ciliates[females$ciliates=="no"] = 0
females$ciliates[females$ciliates=="yes"] = 1
females$ciliates = as.numeric(females$ciliates)
females$col_ciliates = rbPal(2)[as.numeric(cut(females$ciliates, breaks=2))]
females$col_total_offspring = rbPal(10)[as.numeric(cut(females$X12_total_offspring, breaks=20))]

plot(habitat$x, habitat$y, col="white", pch=16, xlab="x coordinate (m)", ylab="y coordinate (m)")
points(females$x, females$y, col=females$col_aperture, pch=16)
points(females$x, females$y, col=females$col_aperture_corr, pch=16)
points(females$x, females$y, col=females$col_len, pch=16)
points(females$x, females$y, col=females$col_weight, pch=16)
points(females$x, females$y, col=females$col_thickness, pch=16)
points(females$x, females$y, col=females$col_boldness, pch=16)
points(females$x, females$y, col=females$col_veliger_width_av, pch=16)
points(females$x, females$y, col=females$col_ciliates, pch=16)
points(females$x, females$y, col=females$col_total_offspring, pch=16)
points(females$x, females$y, col=females$col_RW1, pch=16)


# get Bezier path & points
bezier_file = list.files(paste("../Littorina_CZ_", zone, "_final_data", sep=""), pattern = "bezierPath")
bezier = read.table(paste("../Littorina_CZ_", zone, "_final_data/", bezier_file, sep=""), header=T, sep=",")
ecopoints_file = list.files(paste("../Littorina_CZ_", zone, "_final_data", sep=""), pattern = "ECOdataLitPath")
ecopoints = read.table(paste("../Littorina_CZ_", zone, "_final_data/", ecopoints_file, sep=""), header=T, sep=",")
points(bezier$pathX, bezier$pathY, col="green")

### TEMP
a = females[females$DistAlongPath<20 & (is.na(females$DistAlongPath)==F), ]
a = females[females$DistAlongPath>160 & (is.na(females$DistAlongPath)==F), ]

points(mature$x, mature$y, col="yellow", pch=16)
points(a$x, a$y, col="red", pch=16)
centre = females[females$DistAlongPath>70 & females$DistAlongPath < 95,]
points(centre$BezX, centre$BezY, pch=16, col="blue")




#################
#### 3D PLOT ####

# plot
library("rgl", lib.loc="C:/Users/Anja/Documents/R/win-library/3.0")
minHeight = min(ecopoints$z) - 7
maxHeight = max(ecopoints$z) + 7
plot3d(ecopoints$x, ecopoints$y, ecopoints$z, zlim=c(minHeight, maxHeight)) # ecopoints
plot3d(mature$x, mature$y, mature$z, add=T, col=mature$col_rock) # snails
plot3d(bezier$pathX, bezier$pathY, bezier$pathZ, add=T, col="blue") # path
