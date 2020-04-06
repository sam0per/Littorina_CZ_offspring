rm(list = ls())

.packagesdev = "thomasp85/patchwork"
.packages = c("ggplot2", "dplyr", "reshape2", "geomorph", "tidyr", "data.table", "ggrepel", "abind", "RANN")
# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
.instdev <- basename(.packagesdev) %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
if(length(.packagesdev[!.instdev]) > 0) devtools::install_github(.packagesdev[!.instdev])
# Load packages into session
lapply(.packages, require, character.only=TRUE)
lapply(basename(.packagesdev), require, character.only=TRUE)

# install.packages('ape')
# install.packages('rgl')
# library(ape)
# library(rgl)
# install.packages('geomorph')
# library(geomorph)
# library(dplyr)
# library(abind)

island = "CZA"
fin_dir = paste0(island, "_off_SW/", island, "_off_final_data/")
morpho_dir = paste0(fin_dir, "morphometrics")
cohort = "male.TPS"
######################################################
## Match female IDs before and after the experiment ##
# Run commented commands for PCA of a single cohort ##
######################################################
F_before = readland.tps(file = list.files(path = morpho_dir, pattern = cohort, full.names = TRUE), specID = "imageID")
F_after = readland.tps(file = list.files(path = morpho_dir, pattern = "after", full.names = TRUE), specID = "imageID")
# F_after = readland.tps(file = list.files(path = morpho_dir, pattern = cohort, full.names = TRUE), specID = "imageID")

# strsplit(dimnames(F_before)[[3]][1], split = "_")[[1]][3]
F_before_df = rbindlist(lapply(1:dim(F_before)[3], function(x) {
  snID = strsplit(dimnames(F_before)[[3]][x], split = "_")[[1]][3]
  if (grepl(pattern = "offspring", x = cohort)) {
    data.frame(snail_ID = paste0(substr(snID, start = 1, stop = 1), "_", substr(snID, start = 2, stop = 4)),
               size_mm = F_before[1,2,x] - F_before[9,2,x], row.names = NULL)
  } else {
    data.frame(snail_ID = paste0(substr(snID, start = 1, stop = 1), "_", substr(snID, start = 2, stop = 3)),
               size_mm = F_before[1,2,x] - F_before[9,2,x], row.names = NULL)
  }
}))
F_after_df = rbindlist(lapply(1:dim(F_after)[3], function(x) {
  snID = strsplit(dimnames(F_after)[[3]][x], split = "_")[[1]][3]
  data.frame(snail_ID = paste0(substr(snID, start = 1, stop = 1), "_", substr(snID, start = 2, stop = 3)),
             size_mm = F_after[1,2,x] - F_after[9,2,x], row.names = NULL)
}))
cohort = tools::file_path_sans_ext(cohort)
write.csv(F_before_df, file = paste0(fin_dir, island, "_off_", cohort, "_size_mm.csv"), row.names = FALSE)
write.csv(F_after_df, file = paste0(fin_dir, island, "_off_female_after_size_mm.csv"), row.names = FALSE)
sort(x = F_before_df$size_mm)
sort(x = F_after_df$size_mm)

FF_lm = abind(F_before, F_after)
# FF_lm = F_after
FF_pop = sapply(1:dim(FF_lm)[3], function(x) {
  substr(strsplit(dimnames(FF_lm)[[3]][x], split = "_")[[1]][3], start = 1, stop = 1)
})
pop = "A"
FF_idx = which(FF_pop==pop)
dim(FF_lm)
plotAllSpecimens(FF_lm[,,FF_idx])
plotOutliers(FF_lm[,,FF_idx])
# plotAllSpecimens(FF_lm)
# plotOutliers(FF_lm)

lm.gpa <- gpagen(FF_lm[,,FF_idx]) # GPA-alignment with outliers
# M_CZoff_G52_1.png looks fine
# M_CZoff_P58_1.png looks small crab-like
# outliers = c("M_CZoff_P58_1.png")
# outliers = c("M_CZoff_I18_1.PNG", "M_CZoff_O13_1.PNG", "M_CZoff_G25_1.PNG", "M_CZoff_F18_1.PNG")
# outliers = c("M_CZoff_L303_1.PNG", "M_CZoff_NXX2_1.PNG", "M_CZoff_L301_1.PNG", "M_CZoff_K301_1.PNG",
#              "M_CZoff_K320_1.PNG", "M_CZoff_J320_1.PNG", "M_CZoff_NXX3_1.PNG", "M_CZoff_N324_1.PNG",
#              "M_CZoff_O314_1.PNG", "M_CZoff_I324_1.PNG", "M_CZoff_K307_1.PNG")
# outliers = c("M_CZoff_K202_1.png", "M_CZoff_L202_1.png", "M_CZoff_L204_1.png", "M_CZoff_K206_1.png",
#              "M_CZoff_E206_1.png")
# plotOutliers(FF_lm[,,!(dimnames(FF_lm)[[3]] %in% outliers)])

# lm.gpa <- gpagen(FF_lm[,,!(dimnames(FF_lm)[[3]] %in% outliers)]) # GPA-alignment without outliers
plotAllSpecimens(lm.gpa$coords)
plotOutliers(lm.gpa$coords) # identify outliers by shape
# outliers.sh = c("M_CZoff_C56_1.png","M_CZoff_R54_1.png", "M_CZoff_Q67_1.png", "M_CZoff_D58_1.png",
#                 "M_CZoff_O63_1.png", "M_CZoff_Q66_1.png")
# outliers.sh = c("M_CZoff_G23_1.PNG")
# outliers.sh = c("M_CZoff_K308_1.PNG", "M_CZoff_K308_1.PNG", "M_CZoff_E321_1.PNG", "M_CZoff_I329_1.PNG",
#                 "M_CZoff_O303_1.PNG", "M_CZoff_N305_1.PNG")
# outliers.sh = c("M_CZoff_M201_1.png", "M_CZoff_L203_1.png", "M_CZoff_M209_1.png", "M_CZoff_O202_1.png")
# plotOutliers(lm.gpa$coords[,,!(dimnames(lm.gpa$coords)[[3]] %in% outliers.sh)])

# check after excluding using index given by plotOutliers
# plotTangentSpace(lm.gpa$coords[,,!(dimnames(lm.gpa$coords)[[3]] %in% outliers.sh)], warpgrids = T)

# get scores
# scores <- as.data.frame(plotTangentSpace(lm.gpa$coords[,,!(dimnames(lm.gpa$coords)[[3]] %in% outliers.sh)])[[2]])
# plotTangentSpace(lm.gpa$coords, warpgrids = TRUE)

# get scores
scores <- as.data.frame(plotTangentSpace(lm.gpa$coords)[[2]])

out <- cbind(scores[1:2],lm.gpa$Csize)
# out <- cbind(scores[1:2],lm.gpa$Csize[(!names(lm.gpa$Csize) %in% outliers.sh)])
head(out)
colnames(out)[3] <- "CS"
cat = strsplit(rownames(scores,1,3), "_")
class <- matrix(unlist(cat), ncol=4, byrow=T)
out$pop = substr(class[,3], start = 1, stop = 1)
out$snail_ID = as.integer(substr(class[,3], start = 2, stop = 3))
# if (grepl(pattern = "offspring", x = cohort)) {
#   out$snail_ID = paste0(out$pop, "_", substr(class[,3], start = 2, stop = 4))
# } else {
#   out$snail_ID = paste0(out$pop, "_", substr(class[,3], start = 2, stop = 3))
# }
head(out)
str(out)
write.csv(out, file = paste0(fin_dir, island, "_off_female_after_PCs.csv"), row.names = FALSE, quote = FALSE) # output first two axes
write.csv(out, file = paste0(fin_dir, island, "_off_", cohort, "_PCs.csv"), row.names = FALSE, quote = FALSE)

out_before = out[out$snail_ID < 49, ]
out_after = out[out$snail_ID > 49, ]
closest <- nn2(out_before[, 1:2], out_after[, 1:2], k = 3, searchtype = "standard")
# closest <- sapply(closest, cbind) %>% as_tibble
afterID = 1
beforeID = 1
eq_idx = which(closest$nn.idx[,1]==afterID)
dist_dt = data.frame(snail_ID_after = out_after[eq_idx, "snail_ID"],
                     snail_ID_bef = out_before[beforeID, "snail_ID"],
                     dist = closest$nn.dists[eq_idx, 1])
after_nm = row.names(out_after[which(out_after$snail_ID == dist_dt[which.min(dist_dt$dist), 1]), ])
before_nm = row.names(out_before[beforeID, ])

plot(out$CS,out$PC1)
plot(out$CS,out$PC2)
head(out)
######################################################
######################################################




Male <- readland.tps("CZ-offspring/data/Littorina_CZoff_raw_data_renamed/CZoff_photos_shell_male_201407.txt", specID = "imageID")
head(Male[,,2])
dim(Male)
dimnames(Male)
male.cat = strsplit(dimnames(Male)[[3]], "_")
male.ID <- matrix(unlist(male.cat), ncol=4, byrow=T)[,3]
#class <- cbind(rownames(scores,1,3), class)
data.frame(snail_ID=male.ID,sex="male")


Female <- readland.tps("CZ-offspring/data/Littorina_CZoff_raw_data_renamed/CZoff_photos_shell_females_201407.txt", specID = "imageID")
Female[,,1]
dim(Female)

female.cat = strsplit(dimnames(Female)[[3]], "_")
female.ID <- matrix(unlist(female.cat), ncol=4, byrow=T)[,3]
#class <- cbind(rownames(scores,1,3), class)
femmale = data.frame(snail_ID=c(female.ID,male.ID),sex=c(rep("female",269),rep("male",172)))
head(femmale)


#Off <- readland.tps("Littorina_CZ_CZoff_pictures shell offspring3_august2015.tps", specID = "imageID")
#dim(Off)

#Tot_land <- abind(Female, Male, Off)
FM_lm = abind(Female,Male)
library('abind')
dim(FM_lm)
FM_lm[,,1]
strsplit(dimnames(FM_lm)[[3]], "_")

plotAllSpecimens(FM_lm)

plotOutliers(FM_lm) # identify outliers by size
# M_CZoff_I18_1.PNG is the aperture picture
# M_CZoff_G25_1.PNG looks normal with a light bump along the aperture
# M_CZoff_F18_1.PNG is the aperture picture
# M_CZoff_C02_1.PNG looks normal
# M_CZoff_C08_1.PNG looks normal
# M_CZoff_O13_1.PNG is the aperture picture
# CZoff_M_F04_1.JPG bad positioned
# CZoff_M_G18_1.JPG slightly bad positioned
outliers=c("M_CZoff_I18_1.PNG","M_CZoff_G25_1.PNG","M_CZoff_F18_1.PNG","M_CZoff_C02_1.PNG","M_CZoff_C08_1.PNG",
           "M_CZoff_O13_1.PNG","CZoff_M_F04_1.JPG","CZoff_M_G18_1.JPG")



lm.gpa<-gpagen(FM_lm[,,!(dimnames(FM_lm)[[3]] %in% outliers)]) #GPA-alignment without outliers

plotAllSpecimens(lm.gpa$coords)

plotOutliers(lm.gpa$coords) # identify outliers by shape
outliers.sh = c("CZoff_M_M19_1.PNG","CZoff_M_C17_1.PNG","CZoff_M_E12_1.PNG","CZoff_M_D11_1.PNG","CZoff_M_G08_1.PNG",
                "CZoff_M_N14_1.PNG","CZoff_M_O11_1.PNG","CZoff_M_R02_1.PNG","CZoff_M_H13_1.PNG")

# check after excluding using index given by plotOutliers
plotTangentSpace(lm.gpa$coords[,,!(dimnames(lm.gpa$coords)[[3]] %in% outliers.sh)],warpgrids = T)


# get scores
scores <- as.data.frame(plotTangentSpace(lm.gpa$coords[,,!(dimnames(lm.gpa$coords)[[3]] %in% outliers.sh)])[[2]])

out <- cbind(scores[1:2],lm.gpa$Csize[(!names(lm.gpa$Csize) %in% outliers.sh)])
head(out)
colnames(out)[3] <- "CS"
cat = strsplit(rownames(scores,1,3), "_")
class <- matrix(unlist(cat), ncol=4, byrow=T)
#class <- cbind(rownames(scores,1,3), class)
out$snail_ID = class[,3]
head(out)
str(out)
out = femmale %>% inner_join(out)

#out$snail_ID <- paste(substring(rownames(scores,1,3)),substring(rownames(scores),5,nchar(rownames(scores))-6),sep="")

#plots and head to check
plot(out$CS,out$PC1)
plot(out$CS,out$PC2)
head(out)

write.csv(out,"results/CZ-off-shape-sex.csv",row.names = F,quote = F) # output first two axes


# plotting warp grids

findMeanSpec(lm.gpa$coords)


PCA <- plotTangentSpace(lm.gpa$coords)
PC <- PCA$pc.scores[,1]
preds <- shape.predictor(lm.gpa$coords, x=PC, Intercept = FALSE,
                         method = "PLS",
                         pred1 = -0.1, pred2 = 0.1) # using PLS plot as a guide
M <- mshape(lm.gpa$coords)

plotRefToTarget(M, preds$pred1) # shape at left on PC1
plotRefToTarget(M, preds$pred2) # shape at right on PC1


###############################################################
###############################################################
###############################################################



f <- read.csv("L_females.csv", sep = ";", header = F)
female_sample <- grep(pattern = "IMAGE", x = f$V1, value = T)
female_sample

o <- read.csv("L_offspring.csv", sep = ";", header = F)
off_sample <- grep(pattern = "IMAGE", x = o$V1, value = T)
off_sample

m <- read.csv("L_male.csv", sep = ";", header = F)
male_sample <- grep(pattern = "IMAGE", x = m$V1, value = T)
male_sample

sample_list <- c(female_sample, male_sample, off_sample)
write.table(sample_list, "Snail_ID.txt", sep = " ")

Y.gpa <- gpagen(Tot_land)                              # GPA-alignment
ID <- read.csv("Snail_ID.csv", sep = ";")
PCA <- plotTangentSpace(Y.gpa$coords)
summary(PCA)
xlab <- paste("Principal Component 1 ", "(", round(PCA$pc.summary$importance[2,1]*100, 1), "%)", sep="")
xlab
ylab <- paste("Principal Component 2 ", "(", round(PCA$pc.summary$importance[2,2]*100, 1), "%)", sep="")
ylab
PCA_scores <- as_data_frame(PCA$pc.scores) 
PCA_scores <- subset(PCA_scores, select =  c("PC1", "PC2")) 
PCA_scores$Sex <- ID$Sex
PCA_scores$ID <- ID$ID
write.table(PCA_scores, "PCA_scores.txt", sep = "\t")
ggplot(PCA_scores, aes(x=PCA_scores$PC1, y=PCA_scores$PC2, colour=PCA_scores$Sex)) +
  geom_point()
plot(PCA_scores$PC1, PCA_scores$PC2, pch=21, cex=0.8, bg=PCA_scores$Sex, xlab=xlab, ylab=ylab)
legend("topright", legend= unique(PCA_scores$Sex), pch=19,  col=unique(PCA_scores$Sex), cex = 0.8)
