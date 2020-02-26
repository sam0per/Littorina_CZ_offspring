rm(list = ls())

install.packages('ape')
install.packages('rgl')
library(ape)
library(rgl)
install.packages('geomorph')
library(geomorph)
library(dplyr)
library(abind)

Male <- readland.tps("data/Littorina_CZoff_raw_data_renamed/CZoff_photos_shell_male_201407.txt", specID = "imageID")
head(Male[,,2])
dim(Male)
dimnames(Male)
male.cat = strsplit(dimnames(Male)[[3]], "_")
male.ID <- matrix(unlist(male.cat), ncol=4, byrow=T)[,3]
#class <- cbind(rownames(scores,1,3), class)
data.frame(snail_ID=male.ID,sex="male")


Female <- readland.tps("data/Littorina_CZoff_raw_data_renamed/CZoff_photos_shell_females_201407.txt", specID = "imageID")
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
