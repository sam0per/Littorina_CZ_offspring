rm (list=ls())

#install.packages("bbmle")
library(bbmle)

#start at line 95 unless it is necessary to modify data or line (but need ANGsize, ANGshape, ANGxyzp)

#get data
setwd("C:/Users/Roger/Dropbox/Littorina_CZ/Littorina_CZ_data/Littorina_CZ_ANG_final_data")
ANGxyzp <- read.csv('ANG_spatial.csv', header = T, col.names=c("snail_ID","X","Y","Z","Env"))
ANGshape <- read.csv('Littorina_CZ_ANG_morphometrics/ANG_RWscores_20150517.CSV')
ANGsize <- read.csv('ANG_length_20150606.csv')

snail <- ANGxyzp[substr(ANGxyzp$snail_ID,1,3)=="ANG",]
snail <- merge(snail, ANGshape, by="snail_ID", all=T)
snail <- merge(snail, ANGsize, by="snail_ID", all=T)


#remove snails with missing values and re-centre
snail <- na.omit(snail)
snail$X <- snail$X-2000
snail$Y <- snail$Y-2000
snail$Z <- snail$Z-200


#######################
#apply centre line
#######################

#plot X and Y
plot(snail$X,snail$Y)

#method for making a set of straight-line segments to approximate centre line

#set break points and view (adjust to look good!)
x_lim <- c(-31,5,11,13,16,17,21,30,39,49,63,75)
y_lim <- c(-62,-20,-5,3,3,5,2.5,-1,-1,15,10,10)

plot(snail$X,snail$Y)
lines(x_lim,y_lim,col="blue")

#find snail displacement from this line (distance) and nearest point on line (line_x,line_y)
distance <- rep(0,length(snail$X))
seg <- rep(0,length(snail$X))
snail$line_x <- rep(0,length(snail$X))
snail$line_y <- rep(0,length(snail$X))
m <- (y_lim[-1]-y_lim[-length(y_lim)])/(x_lim[-1]-x_lim[-length(x_lim)])

for (n in 1:length(snail$X)) {
  dist <- abs(snail$Y[n]-(m*snail$X[n])-y_lim[-length(y_lim)]+(m*x_lim[-length(x_lim)]))/sqrt(1+(m*m))
  dist_l <- sqrt((snail$X[n]-x_lim[-length(x_lim)])^2+(snail$Y[n]-y_lim[-length(y_lim)])^2)
  dist_r <- sqrt((snail$X[n]-x_lim[-1])^2+(snail$Y[n]-y_lim[-1])^2) 
  #find closest point on each line
  k <- ((y_lim[-1]-y_lim[-length(y_lim)]) * (snail$X[n]-x_lim[-length(x_lim)]) - (x_lim[-1]-x_lim[-length(x_lim)]) * (snail$Y[n]-y_lim[-length(y_lim)])) / ((y_lim[-1]-y_lim[-length(y_lim)])^2 + (x_lim[-1]-x_lim[-length(x_lim)])^2)
  line_x <- (snail$X[n] - k * (y_lim[-1]-y_lim[-length(y_lim)]))
  line_y <- (snail$Y[n] - k * (x_lim[-1]-x_lim[-length(x_lim)]))
  #convert to distance within segment
  dist[line_x<x_lim[-length(x_lim)]] <- dist_l[line_x<x_lim[-length(x_lim)]]
  dist[line_x>x_lim[-1]] <- dist_r[line_x>x_lim[-1]]
  #find distance and segment
  distance[n] <- min(dist)
  seg[n] <- which(dist==min(dist))
  snail$line_x[n] <- (snail$X[n] - k * (y_lim[-1]-y_lim[-length(y_lim)]))[dist==min(dist)]
  snail$line_y[n] <- (snail$Y[n] - k * (x_lim[-1]-x_lim[-length(x_lim)]))[dist==min(dist)]
}
#output mean square displacement (can fiddle with x-lim and y_lim to try to improve)
sum(distance*distance)/length(snail$X)

snail$dist_from_line <- distance

#plot to check
plot(snail$line_x, snail$X)
plot(snail$line_y, snail$Y)

#find distance along line
seg_len <- sqrt((x_lim[-1]-x_lim[-length(x_lim)])^2+(y_lim[-1]-y_lim[-length(y_lim)])^2)
seg_dist <- sqrt((snail$line_x-x_lim[seg])^2+(snail$line_y-y_lim[seg])^2)
for (n in 1:length(snail$X)) {
  snail$line_dist[n] <- sum(seg_len[1:seg[n]])+seg_dist[n]-seg_len[seg[n]]
}

#plots to check
plot(snail$X,snail$line_dist)
plot(snail$X,snail$line_x)
plot(snail$Y,snail$line_y)

plot(snail$X,snail$Y,col=seg)
lines(x_lim,y_lim,col="blue")

#output
setwd("~/R projects/ANG")
write.csv(snail,"ANG_spatial_with_line.csv")

################
# end of centre line routine
################

#read input with line distance
setwd("~/R projects/ANG")
snail <- read.csv("ANG_spatial_with_line.csv")

#add dissection data
setwd("~/Roger/R projects/ANG")
ANGdiss <- read.csv('ANG_dissections_20150617.csv')
snail <- merge(snail, ANGdiss, by="snail_ID", all=T)

#add colour data
ANGcolour <- read.csv("ANG_colour.csv")
snail <- merge(snail, ANGcolour, by="snail_ID", all=T)

##################
# size at maturity analyses (NB before standardisation of size)
##################

plot(snail$line_dist,snail$length_mm)

snail$mature <- rep(1,length(snail$snail_ID))
snail$mature[snail$sex_dissection=='juvenile'] <- 0
snail$mature[is.na(snail$sex_dissection)] <- NA

plot(snail$length_mm,snail$mature)

size_mat <- glm(snail$mature ~ log(snail$length_mm), family=binomial(logit), na.action=na.exclude)
summary(size_mat)
plot_data <- na.omit(cbind(log(snail$length_mm),fitted(size_mat),snail$mature))
points(plot_data[,c(1)], plot_data[,c(2)],col="red")

size_mat_crab <- glm(snail$mature[snail$line_dist<70] ~ log(snail$length_mm[snail$line_dist<70]), family=binomial(logit), na.action=na.exclude)
summary(size_mat_crab)
plot_data <- na.omit(cbind(log(snail$length_mm[snail$line_dist<70]),fitted(size_mat_crab),snail$mature[snail$line_dist<70]))
plot(plot_data[,c(1)], plot_data[,c(3)],col="black")
points(plot_data[,c(1)], plot_data[,c(2)],col="green")

size_mat_wave <- glm(snail$mature[snail$line_dist>90] ~ log(snail$length_mm[snail$line_dist>90]), family=binomial(logit), na.action=na.exclude)
summary(size_mat_wave)
plot_data <- na.omit(cbind(log(snail$length_mm[snail$line_dist>90]),fitted(size_mat_wave),snail$mature[snail$line_dist>90]))
plot(plot_data[,c(1)], plot_data[,c(3)],col="black")
points(plot_data[,c(1)], plot_data[,c(2)],col="blue")

snail$ecotype <- "crab"
snail$ecotype[snail$line_dist>90] <- "wave"
snail$ecotype[snail$line_dist<90 & snail$line_dist>70] <- "hybrid"

lsize_mat_eco <- glm(snail$mature ~ log(snail$length_mm)*snail$ecotype, family=binomial(logit), na.action=na.exclude)
summary(lsize_mat_eco)

# linear function that is more easily incorporated in a cline fit

#need to get rid of NAs first

snail <- snail[(is.na(snail$length_mm) == F & is.na(snail$mature) == F), ]

s_mat <- function(size,mat,mean,slope){
  logit_p <- (size-mean)*slope
  p <- exp(logit_p)/(exp(logit_p)+1)
  minusll <- -sum(dbinom(mat,1,p,log=T))

  return(minusll)

}

theta.init <- list(mean=5,slope=5)

mle.s_mat <- mle2(s_mat,theta.init,data=list(size=log(snail$length_mm),mat=snail$mature))

summary(mle.s_mat)
AIC(mle.s_mat)

# build that into a cline fit (cline for mean size at maturity first, constant slope)

s_mat_cline <- function(position,size,mat,c_mean,w_mean,slope, c, w){
  d <- position-c
  
  p_x <- 1-exp(d/w)/(1+exp(d/w)) # p_x is frequency cline as in Cfit NB p=1 for left, 0 for right
  
  mean <- w_mean+(c_mean-w_mean)*p_x
  logit_p <- (size-mean)*slope
  p <- exp(logit_p)/(exp(logit_p)+1)
  minusll <- -sum(dbinom(mat,1,p,log=T))
  
  return(minusll)
  
}

theta.init <- list(c_mean=1,w_mean=1,slope=2,c=80,w=20)

mle.s_mat_cline <- mle2(s_mat_cline,theta.init,data=list(position=snail$line_dist,size=log(snail$length_mm),mat=snail$mature))

summary(mle.s_mat_cline)
AIC(mle.s_mat_cline)

# try adding variation in the slope, Cfit style 

s_mat_cline_slope <- function(position,size,mat,c_mean,w_mean,c_slope,h_slope,w_slope, c, w){
  d <- position-c
  
  p_x <- 1-exp(d/w)/(1+exp(d/w)) # p_x is frequency cline as in Cfit NB p=1 for left, 0 for right
  
  mean <- w_mean+(c_mean-w_mean)*p_x
  slope <- (p_x^2)*c_slope + 2*p_x*(1-p_x)*h_slope + ((1-p_x)^2)*w_slope
    
  logit_p <- (size-mean)*slope
  p <- exp(logit_p)/(exp(logit_p)+1)
  minusll <- -sum(dbinom(mat,1,p,log=T))
  
  return(minusll)
  
}

theta.init <- list(c_mean=2,w_mean=1,c_slope=5,h_slope=5,w_slope=5,c=90,w=10)

mle.s_mat_cline_slope <- mle2(s_mat_cline_slope,theta.init,data=list(position=snail$line_dist,size=log(snail$length_mm),mat=snail$mature))

summary(mle.s_mat_cline_slope)
AIC(mle.s_mat_cline_slope)

# 2 slope version

s_mat_cline_2slope <- function(position,size,mat,c_mean,w_mean,c_slope,w_slope, c, w){
  d <- position-c
  
  p_x <- 1-exp(d/w)/(1+exp(d/w)) # p_x is frequency cline as in Cfit NB p=1 for left, 0 for right
  
  mean <- w_mean+(c_mean-w_mean)*p_x
  slope <- p_x*c_slope + (1-p_x)*w_slope
  
  logit_p <- (size-mean)*slope
  p <- exp(logit_p)/(exp(logit_p)+1)
  minusll <- -sum(dbinom(mat,1,p,log=T))
  
  return(minusll)
  
}

theta.init <- list(c_mean=2,w_mean=1,c_slope=5,w_slope=5,c=90,w=10)

mle.s_mat_cline_2slope <- mle2(s_mat_cline_2slope,theta.init,data=list(position=snail$line_dist,size=log(snail$length_mm),mat=snail$mature))

summary(mle.s_mat_cline_2slope)
AIC(mle.s_mat_cline_2slope)

#plot result
dist_plot <- 1:155
plot(snail$line_dist[snail$mature==1],log(snail$length_mm[snail$mature==1]),ylim=c(0.8,3))
points(snail$line_dist[snail$mature==0],log(snail$length_mm[snail$mature==0]),pch=16)
lines(dist_plot,0.85+(1.82-0.85)*(1-exp((dist_plot-91.8)/3.56)/(1+exp((dist_plot-91.8)/3.56))))


######################
# sort environmental variables
######################


#finding 'rockiness' for each snail position
#get topo points,recentre and remove topo points with missing values and re-centre


#ANGtopo <- ANGxyzp[substr(ANGxyzp$snail_ID,1,1) %in% 1:9 | ANGxyzp$X > 2045, ] # no topo at wave end so use positions of snails and assume all rock
ANGtopo <- ANGxyzp[substr(ANGxyzp$snail_ID,1,1) %in% 1:9 , ] # use only topo points, not snail points (because they have no env info for ANG)

ANGtopo$X <- ANGtopo$X-2000
ANGtopo$Y <- ANGtopo$Y-2000
ANGtopo$Z <- ANGtopo$Z-200
ANGtopo <- na.omit(ANGtopo)

plot(ANGtopo$X,ANGtopo$Y)

#calculate distances
snailxyz <- as.matrix(snail[,c("X","Y","Z")])
distarray <- apply(snailxyz,1,function(x) sqrt((x[1]-ANGtopo$X)^2+(x[2]-ANGtopo$Y)^2+(x[3]-ANGtopo$Z)))

#set maximum distance over which weight is calculated and find weights
dmax <- 2
wt <- (1-(1/dmax)*distarray)
wt[wt<0] <- NA

wtdf<-data.frame(wt)
rm(distarray,wt)

#use weights to assign rockiness
wtdfR <- wtdf[substr(ANGtopo$Env,1,1)=="R" ,]
wtdfN <- wtdf[substr(ANGtopo$Env,1,1)=="S" ,]

sumwtN <-colSums(as.matrix(wtdfN),na.rm=TRUE)
sumwtR <-colSums(as.matrix(wtdfR),na.rm=TRUE)

snail$Rockiness <- (sumwtR-sumwtN)/(sumwtR+sumwtN)

plot(snail$line_dist,snail$Rockiness)

rbPal = colorRampPalette(c('yellow', 'blue'))
snail$col_rock = rbPal(10)[as.numeric(cut(snail$Rockiness, breaks=10))]
plot(snail$X, snail$Y, xaxt='n', yaxt='n', xlab="", ylab="", col=snail$col_rock, pch=16, cex=0.9)


#use weights to assign Barnacleness
wtdfB <- wtdf[substr(ANGtopo$Env,2,2)=="B",]
wtdfN <- wtdf[substr(ANGtopo$Env,2,2)=="N",]

sumwtN <-colSums(as.matrix(wtdfN),na.rm=TRUE)
sumwtB <-colSums(as.matrix(wtdfB),na.rm=TRUE)

snail$Barnacleness <- (sumwtB-sumwtN)/(sumwtB+sumwtN)

#use weights to assign Fucusness
wtdfF <- wtdf[substr(ANGtopo$Env,3,3)=="F",]
wtdfN <- wtdf[substr(ANGtopo$Env,3,3)=="N",]

sumwtN <-colSums(as.matrix(wtdfN),na.rm=TRUE)
sumwtF <-colSums(as.matrix(wtdfF),na.rm=TRUE)

snail$Fucusness <- (sumwtF-sumwtN)/(sumwtF+sumwtN)

plot(snail$line_dist,snail$Barnacleness)
plot(snail$line_dist,snail$Fucusness)

#recover space, if necessary
rm(wtdf,wtdfN,wtdfR,wtdfF,wtdfB)

# combine into a PC
pca <- princomp(~snail[,22]+snail[,24]+snail[,25],na.action=na.exclude)
plot(snail$line_dist,pca$scores[,1])
plot(snail$line_dist,pca$scores[,2])
plot(snail$line_dist,pca$scores[,3])
plot(snail$Z,pca$scores[,1])
plot(snail$Z,pca$scores[,2])
plot(snail$Z,pca$scores[,3])

snail$PC1 <- pca$scores[,1]


###################
# clines for CS and RW1 and habitat
###################

# get dissection information (if not in already)
setwd("C:/Users/bo1rkb/Dropbox/Littorina_CZ/Littorina_CZ_data/Littorina_CZ_ANG_final_data")
ANGdiss <- read.csv('ANG_dissections_20150723.csv')
snail <- merge(snail, ANGdiss, by="snail_ID", all=T)

#######################
# calculate centroid size (if needed)
#######################

setwd("~/R projects/ANG")
ANGtps <- read.csv('ANGtps.csv')
ANGtps <- merge(ANGtps, ANGsize, by="snail_ID", all=T)

scale <- sqrt((ANGtps$X1-ANGtps$X9)^2+(ANGtps$Y1-ANGtps$Y9)^2)/ANGtps$length_mm

centroid_x <- rowSums(ANGtps[ ,c(seq(2,30, by=2))])/15
centroid_y <- rowSums(ANGtps[ ,c(seq(3,31, by=2))])/15

# centroid size is sqrt of sum of squared distances of landmarks from centroid
ANGtps$CS <- sqrt(rowSums(((ANGtps[ ,c(seq(2,30, by=2))]-centroid_x)/scale)^2+((ANGtps[ ,c(seq(3,31, by=2))]-centroid_y)/scale)^2))

snail <- merge(snail, ANGtps[ ,c(1,33)], by="snail_ID", all=T)

#snail <- na.omit(snail)


male <- snail[snail$sex_dissection=="male" & is.na(snail$CS)==F, ]
male <- male[substr(row.names(male),1,2)!="NA",]
female <- snail[snail$sex_dissection=="female" & is.na(snail$CS)==F, ]
female <- female[substr(row.names(female),1,2)!="NA",]

adult <- snail[(snail$sex_dissection=="male" | snail$sex_dissection=="female") & is.na(snail$CS)==F, ]
adult <- adult[substr(row.names(adult),1,2)!="NA",]

# start by putting on standard scale 0,1 with Wave>crab, sexes separate

male$size <- 1-(log(male$CS)-min(log(male$CS),na.rm=T))/(max(log(male$CS),na.rm=T)-min(log(male$CS),na.rm=T)) # '1-' to reverse direction
plot(male$line_dist,male$size)
female$size <- 1-(log(female$CS)-min(log(female$CS),na.rm=T))/(max(log(female$CS),na.rm=T)-min(log(female$CS),na.rm=T))
points(female$line_dist,female$size,pch=16)

# can scale all adults together (preferred)
adult$size <- 1-(log(adult$CS)-min(log(adult$CS),na.rm=T))/(max(log(adult$CS),na.rm=T)-min(log(adult$CS),na.rm=T))
# or males and females separately
adult$size_sep[adult$sex_dissection=="male"] <- 1-(log(adult$CS[adult$sex_dissection=="male"])-min(log(adult$CS[adult$sex_dissection=="male"]),na.rm=T))/(max(log(adult$CS[adult$sex_dissection=="male"]),na.rm=T)-min(log(adult$CS[adult$sex_dissection=="male"]),na.rm=T))
adult$size_sep[adult$sex_dissection=="female"] <- 1-(log(adult$CS[adult$sex_dissection=="female"])-min(log(adult$CS[adult$sex_dissection=="female"]),na.rm=T))/(max(log(adult$CS[adult$sex_dissection=="female"]),na.rm=T)-min(log(adult$CS[adult$sex_dissection=="female"]),na.rm=T))

# fit a simple cline for each sex (now using HZAR formulation)

cline <- function(phen,position,centre,w,crab,wave,sc,sh,sw){
  
  d <- position-centre
  
  p_x <- 1/(1+exp(0-4*(d)/w)) 
  
  # p_x is frequency cline as in HZAR NB p=0 for left, 1 for right (corresponding to scaling phenotypes with crab=0, wave=1)
  
  z_x <- crab + (wave-crab)*p_x  
  # z_x is expected phenotype, wave-crab always positive, given scaling
  
  s_x <- sqrt(sc^2 + 4*p_x*(1-p_x)*sh^2 + (p_x^2)*(sw^2-sc^2))
  # unimodal variance model as in HZAR, sx is SD for wave/hybrid/crab (assumes variances are additive, unlike Cfit)
  
  minusll <- -sum(dnorm(phen,z_x,s_x,log=T))
  
  return(minusll)
  
}

theta.init <- list(centre=90,w=10,crab=0.8,wave=0.2,sc=0.2,sh=0.2,sw=0.2)

mle.cline.m <- mle2(cline, theta.init, data=list(position=male$line_dist,phen=male$size))
summary(mle.cline.m)
AIC(mle.cline.m)

mle.cline.f <- mle2(cline, theta.init, data=list(position=female$line_dist,phen=female$size))
summary(mle.cline.f)
AIC(mle.cline.f)

mle.cline.a <- mle2(cline, theta.init, data=list(position=adult$line_dist,phen=adult$size))
summary(mle.cline.a)
AIC(mle.cline.a)

# add centre lines to plot

pars.m <- coef(mle.cline.m) # get coefficients for fitted model
pars.f <- coef(mle.cline.f)
pos <- c(0:160) # ordered distance vector

dist.m <- pos-pars.m[1]
dist.f <- pos-pars.f[1]

freq.m <- 1/(1+exp(0-4*(dist.m)/pars.m[2]))
freq.f <- 1/(1+exp(0-4*(dist.f)/pars.f[2]))
mean.m <- pars.m[4]+(pars.m[3]-pars.m[4])*freq.m  # fitted phenotype mean - male
mean.f <- pars.f[4]+(pars.f[3]-pars.f[4])*freq.f  # fitted phenotype mean - female

lines(pos, mean.m, col="blue")
lines(pos, mean.f, col="red")

# fits allowing for male-female mean difference at wave end

cline.wmf <- function(phen,position,sex,centre,w,crab,wave.m,wave.f,sc,sh,sw){
  
  d <- position-centre
  
  p_x <- 1/(1+exp(0-4*(d)/w)) 
  
  # p_x is frequency cline as in HZAR NB p=0 for left, 1 for right (corresponding to scaling phenotypes with crab=0, wave=1)
  
  z_x <- crab + (wave.m-crab)*p_x  
  # z_x is expected phenotype, wave-crab always positive, given scaling 
  z_x[sex=="female"] <- crab + (wave.f-crab)*p_x[sex=="female"]
  
  
  s_x <- sqrt(sc^2 + 4*p_x*(1-p_x)*sh^2 + (p_x^2)*(sw^2-sc^2))
  # unimodal variance model as in HZAR, sx is SD for wave/hybrid/crab (assumes variances are additive, unlike Cfit)
  
  minusll <- -sum(dnorm(phen,z_x,s_x,log=T))
  
  return(minusll)
  
}

theta.init <- list(centre=90,w=10,crab=0.8,wave.m=0.2,wave.f=0.2,sc=0.2,sh=0.2,sw=0.2)

mle.cline.wmf <- mle2(cline.wmf, theta.init, data=list(position=adult$line_dist,phen=adult$size, sex=adult$sex_dissection))
summary(mle.cline.wmf)
AIC(mle.cline.wmf)


# adding sex difference in variance at wave end
cline.wmfv <- function(phen,position,sex,centre,w,crab,wave.m,wave.f,sc,sh,sw.m,sw.f){
  
  d <- position-centre
  
  p_x <- 1/(1+exp(0-4*(d)/w)) 
  
  # p_x is frequency cline as in HZAR NB p=0 for left, 1 for right (corresponding to scaling phenotypes with crab=0, wave=1)
  
  z_x <- crab + (wave.m-crab)*p_x  
  # z_x is expected phenotype, wave-crab always positive, given scaling  
  z_x[sex=="female"] <- crab + (wave.f-crab)*p_x[sex=="female"]
  
  
  s_x <- sqrt(sc^2 + 4*p_x*(1-p_x)*sh^2 + (p_x^2)*(sw.m^2-sc^2))
  s_x[sex=="female"] <- sqrt(sc^2 + 4*p_x[sex=="female"]*(1-p_x[sex=="female"])*sh^2 + (p_x[sex=="female"]^2)*(sw.f^2-sc^2))
  
  
  minusll <- -sum(dnorm(phen,z_x,s_x,log=T))
  
  return(minusll)
  
}

theta.init <- list(centre=90,w=10,crab=0.8,wave.m=0.2,wave.f=0.2,sc=0.2,sh=0.2,sw.m=0.2,sw.f=0.2)

mle.cline.wmfv <- mle2(cline.wmfv, theta.init, data=list(position=adult$line_dist,phen=adult$size, sex=adult$sex_dissection))
summary(mle.cline.wmfv)
AIC(mle.cline.wmfv)

# adding crab end difference
cline.wcmf <- function(phen,position,sex,centre,w,crab.m,crab.f,wave.m,wave.f,sc,sh,sw){
  
  d <- position-centre
  
  p_x <- 1/(1+exp(0-4*(d)/w)) 
  
  # p_x is frequency cline as in HZAR NB p=0 for left, 1 for right (corresponding to scaling phenotypes with crab=0, wave=1)
  
  z_x <- crab.m + (wave.m-crab.m)*p_x  
  # z_x is expected phenotype, wave-crab always positive, given scaling  
  z_x[sex=="female"] <- crab.f + (wave.f-crab.f)*p_x[sex=="female"]
  
  
  s_x <- sqrt(sc^2 + 4*p_x*(1-p_x)*sh^2 + (p_x^2)*(sw^2-sc^2))
  
  minusll <- -sum(dnorm(phen,z_x,s_x,log=T))
  
  return(minusll)
  
}

theta.init <- list(centre=90,w=10,crab.m=0.8,crab.f=0.8,wave.m=0.2,wave.f=0.2,sc=0.2,sh=0.2,sw=0.2)

mle.cline.wcmf <- mle2(cline.wcmf, theta.init, data=list(position=adult$line_dist,phen=adult$size, sex=adult$sex_dissection))
summary(mle.cline.wcmf)
AIC(mle.cline.wcmf)

# replace with a single offset
cline.mf <- function(phen,position,sex,centre,w,crab,wave,z_s,sc,sh,sw){
  
  d <- position-centre
  
  p_x <- 1/(1+exp(0-4*(d)/w)) 
  
  # p_x is frequency cline as in HZAR NB p=0 for left, 1 for right (corresponding to scaling phenotypes with crab=0, wave=1)
  
  z_x <- carb + (wave-crab)*p_x  
  z_x[sex=="female"] <- z_x[sex=="female"]+z_s
  # z_x is expected phenotype, crab-wave always positive, given scaling
  
  s_x <- sqrt(sc^2 + 4*p_x*(1-p_x)*sh^2 + (p_x^2)*(sw^2-sc^2))
  
  minusll <- -sum(dnorm(phen,z_x,s_x,log=T))
  
  return(minusll)
  
}

theta.init <- list(centre=90,w=10,crab=0.8,wave=0.2,z_s=0.1,sc=0.2,sh=0.2,sw=0.2)

mle.cline.mf <- mle2(cline.mf, theta.init, data=list(position=adult$line_dist,phen=adult$size, sex=adult$sex_dissection))
summary(mle.cline.mf)
AIC(mle.cline.mf)

# replace with a clinal offset (preferred model)
cline.sex <- function(phen,position,sex,centre,w,crab,wave,zs_c,zs_w,sc,sh,sw){
  
  d <- position-centre
  
  p_x <- 1/(1+exp(0-4*(d)/w)) 
  
  # p_x is frequency cline as in HZAR NB p=0 for left, 1 for right (corresponding to scaling phenotypes with crab=0, wave=1)
  
  z_x <- crab + (wave-crab)*p_x  
  z_x[sex=="female"] <- z_x[sex=="female"]+zs_c + (zs_w-zs_c)*p_x[sex=="female"]
  # female offset depends on position in cline
  
  s_x <- sqrt(sc^2 + 4*p_x*(1-p_x)*sh^2 + (p_x^2)*(sw^2-sc^2))
  
  minusll <- -sum(dnorm(phen,z_x,s_x,log=T))
  
  return(minusll)
  
}

theta.init <- list(centre=90,w=10,crab=0.8,wave=0.2,zs_c=0.1,zs_w=0.1,sc=0.2,sh=0.2,sw=0.2)

mle.cline.sex <- mle2(cline.sex, theta.init, data=list(position=adult$line_dist,phen=adult$size, sex=adult$sex_dissection))
summary(mle.cline.sex)
AIC(mle.cline.sex)

#plot this model

pars <- coef(mle.cline.sex) # get coefficients for fitted model
pos <- c(0:160) # ordered distance vector

dist <- pos-pars[1]

freq <- 1/(1+exp(0-4*(dist)/pars[2]))
mean.m <- pars[3]+(pars[4]-pars[3])*freq  # fitted phenotype mean - male
mean.f <- pars[3]+(pars[4]-pars[3])*freq+pars[5]+(pars[6]-pars[5])*freq  # fitted phenotype mean - female

plot(adult$line_dist,adult$size,col=adult$sex_dissection,xlab="Distance(m)",ylab="1-(scaled centroid size)")
lines(pos, mean.m, col="green")
lines(pos, mean.f, col="black")

# and check residuals
freq_ld <- 1/(1+exp(0-4*(adult$line_dist-pars[1])/pars[2]))
adult$fit <- pars[3]+(pars[4]-pars[3])*freq_ld  # fitted phenotype mean - male by line distance
adult$fit[adult$sex=="female"] <- pars[3]+(pars[4]-pars[3])*freq_ld[adult$sex=="female"]+pars[5]+(pars[6]-pars[5])*freq_ld[adult$sex=="female"]  # fitted phenotype mean - female

adult$sd <- sqrt(pars[7]+4*freq_ld*(1-freq_ld)*pars[8]+freq_ld^2*(pars[9]-pars[7]))

adult$st_res <- (adult$size-adult$fit)/adult$sd

plot(adult$line_dist,adult$st_res)
points(adult$line_dist, fitted(loess(adult$st_res~adult$line_dist)),col="blue")
hist(adult$st_res)


#########################
# testing for effects of environmental variables

adult$hos <- adult$Z-mean(adult$Z) #to express height on shore relative to mean hos

plot(adult$hos,adult$st_res)
plot(adult$Rockiness,adult$st_res)
plot(adult$Barnacleness,adult$st_res)
plot(adult$Fucusness,adult$st_res)

# be is the regression for the env variable, differing between ends (be_c and be_w) and changing clinally

cline.sex.env <- function(phen,position,sex,env,centre,w,crab,wave,zs_c,zs_w,be_c,be_w,sc,sh,sw){
  
  d <- position-centre
  
  p_x <- 1/(1+exp(0-4*(d)/w)) 
  
  # p_x is frequency cline as in HZAR NB p=0 for left, 1 for right (corresponding to scaling phenotypes with crab=0, wave=1)
  
  z_x <- crab + (wave-crab)*p_x  
  z_x[sex=="female"] <- z_x[sex=="female"]+zs_c + (zs_w-zs_c)*p_x[sex=="female"]
  z_x <- z_x+env*(be_c+(be_w-be_c)*p_x)
  # z_x is expected phenotype, crab-wave always positive, given scaling
  # female offset depends on position in cline
  # env regression depends on position in cline
  
  s_x <- sqrt(sc^2 + 4*p_x*(1-p_x)*sh^2 + (p_x^2)*(sw^2-sc^2))
  
  minusll <- -sum(dnorm(phen,z_x,s_x,log=T))
  
  return(minusll)
  
}

theta.init <- list(centre=90,w=10,crab=0.8,wave=0.2,zs_c=0.1,zs_w=0.1,be_c=0,be_w=0,sc=0.2,sh=0.2,sw=0.2)

# for height on shore
mle.cline.sex.env <- mle2(cline.sex.env, theta.init, data=list(position=adult$line_dist,phen=adult$size, sex=adult$sex_dissection, env=adult$hos))
summary(mle.cline.sex.env)
AIC(mle.cline.sex.env)

#plot this model (effectively corrected for hos)

pars <- coef(mle.cline.sex.env) # get coefficients for fitted model
pos <- c(0:160) # ordered distance vector

dist <- pos-pars[1]

freq <- 1/(1+exp(0-4*(dist)/pars[2]))
mean.m <- pars[3]+(pars[4]-pars[3])*freq  # fitted phenotype mean - male
mean.f <- pars[3]+(pars[4]-pars[3])*freq+pars[5]+(pars[6]-pars[5])*freq  # fitted phenotype mean - female

plot(adult$line_dist,adult$size,col=adult$sex_dissection,xlab="Distance(m)",ylab="1-(scaled centroid size)",pch=19)
lines(pos, mean.m, col="green", lwd=2)
lines(pos, mean.f, col="black", lwd=2)

# and get residuals for association analysis (not standardised in this case)
freq_ld <- 1/(1+exp(0-4*(adult$line_dist-pars[1])/pars[2]))
adult$fit <- pars[3]+(pars[4]-pars[3])*freq_ld  # fitted phenotype mean - male by line distance
adult$fit[adult$sex=="female"] <- pars[3]+(pars[4]-pars[3])*freq_ld[adult$sex=="female"]+pars[5]+(pars[6]-pars[5])*freq_ld[adult$sex=="female"]  # fitted phenotype mean - female

#adult$sd <- sqrt(pars[7]+4*freq_ld*(1-freq_ld)*pars[8]+freq_ld^2*(pars[9]-pars[7]))

adult$CS_res <- (adult$size-adult$fit)

# function for cline with sex, height and one other variable
cline.sex.envh <- function(phen,position,sex,height,env,centre,w,crab,wave,zs_c,zs_w,bh_c,bh_w,be_c,be_w,sc,sh,sw){
  
  d <- position-centre
  
  p_x <- 1/(1+exp(0-4*(d)/w)) 
  
  # p_x is frequency cline as in HZAR NB p=0 for left, 1 for right (corresponding to scaling phenotypes with crab=0, wave=1)
  
  z_x <- crab + (wave-crab)*p_x  
  z_x[sex=="female"] <- z_x[sex=="female"]+zs_c + (zs_w-zs_c)*p_x
  z_x <- z_x+height*(bh_c+(bh_w-bh_c)*p_x)
  z_x <- z_x+env*(be_c+(be_w-be_c)*p_x)
  # z_x is expected phenotype, crab-wave always positive, given scaling
  # female offset depends on position in cline
  # height and env regressions depend on position in cline
  
  s_x <- sqrt(sc^2 + 4*p_x*(1-p_x)*sh^2 + (p_x^2)*(sw^2-sc^2))
  
  minusll <- -sum(dnorm(phen,z_x,s_x,log=T))
  
  return(minusll)
  
}

#########################
#### remaining CS functions NOT yet updated 2_2_17
#########################

# for rockiness
mle.cline.sex.env <- mle2(cline.sex.env, theta.init, data=list(position=adult$line_dist,phen=adult$size, sex=adult$sex_dissection, env=adult$Rockiness))
summary(mle.cline.sex.env)
AIC(mle.cline.sex.env)

# for barnacleness
mle.cline.sex.env <- mle2(cline.sex.env, theta.init, data=list(position=adult$line_dist,phen=adult$size, sex=adult$sex_dissection, env=adult$Barnacleness))
summary(mle.cline.sex.env)
AIC(mle.cline.sex.env)

# for fucusness
mle.cline.sex.env <- mle2(cline.sex.env, theta.init, data=list(position=adult$line_dist,phen=adult$size, sex=adult$sex_dissection, env=adult$Fucusness))
summary(mle.cline.sex.env)
AIC(mle.cline.sex.env)

plot(adult$Z,adult$Rockiness)

# now with height and one other env variable



theta.init <- list(centre=90,slope=-0.02,crab=0.8,wave=0.2,zs_c=0.1,zs_w=0.1,bh_c=0,bh_w=0,be_c=0,be_w=0,sc=0.2,sh=0.2,sw=0.2)

# for rockiness
mle.cline.sex.envh <- mle2(cline.sex.envh, theta.init, data=list(position=adult$line_dist,phen=adult$size, sex=adult$sex_dissection, height=adult$Z, env=adult$Rockiness))
summary(mle.cline.sex.envh)
AIC(mle.cline.sex.envh)

# for barnacleness
mle.cline.sex.envh <- mle2(cline.sex.envh, theta.init, data=list(position=adult$line_dist,phen=adult$size, sex=adult$sex_dissection, height=adult$Z, env=adult$Barnacleness))
summary(mle.cline.sex.envh)
AIC(mle.cline.sex.envh)

# for fucusness
mle.cline.sex.envh <- mle2(cline.sex.envh, theta.init, data=list(position=adult$line_dist,phen=adult$size, sex=adult$sex_dissection, height=adult$Z, env=adult$Fucusness))
summary(mle.cline.sex.envh)
AIC(mle.cline.sex.envh)

# now with height plus rockiness plus barnacleness

cline.sex.env3 <- function(phen,position,sex,height,rness,bness,centre,slope,crab,wave,zs_c,zs_w,bh_c,bh_w,br_c,br_w,bb_c,bb_w,sc,sh,sw){
  
  d <- position-centre
  
  p_x <- exp(slope*d)/(1+exp(slope*d)) 
  
  # p_x is frequency cline as in Cfit NB p=0 for right, 1 for left (corresponding to scaling phenotypes with crab=1, wave=0)
  # frequency slope is then negative
  
  z_x <- wave + (crab-wave)*p_x  
  z_x[sex=="female"] <- z_x[sex=="female"]+zs_w + (zs_c-zs_w)*p_x
  z_x <- z_x+height*(bh_w+(bh_c-bh_w)*p_x)
  z_x <- z_x+rness*(br_w+(br_c-br_w)*p_x)
  z_x <- z_x+bness*(bb_w+(bb_c-bb_w)*p_x)
  # z_x is expected phenotype, crab-wave always positive, given scaling
  # female offset depends on position in cline
  # height and env regressions depend on position in cline
  
  s_x <- sqrt((p_x^2)*sc^2 + 2*p_x*(1-p_x)*sh^2 + ((1-p_x)^2)*sw^2)
  # unimodal variance model as in Cfit, sx is SD for wave/hybrid/crab (but I assume variances are additive)
  
  minusll <- -sum(dnorm(phen,z_x,s_x,log=T))
  
  return(minusll)
  
}

theta.init <- list(centre=90,slope=-0.02,crab=0.8,wave=0.2,zs_c=0.1,zs_w=0.1,bh_c=0,bh_w=0,br_c=0,br_w=0,bb_c=0,bb_w=0,sc=0.2,sh=0.2,sw=0.2)

# for rockiness
mle.cline.sex.env3 <- mle2(cline.sex.env3, theta.init, data=list(position=adult$line_dist,phen=adult$size, sex=adult$sex_dissection, height=adult$Z, rness=adult$Rockiness, bness=adult$Barnacleness))
summary(mle.cline.sex.env3)
AIC(mle.cline.sex.env3)

#try confidence intervals on this fit (slow because many dimensions!)
confint(mle.cline.sex.env3)

####################################################################
#
#  RW1 - shape
#
####################################################################

#first, rescale to 0-1 with higher values for wave
adult$shape <- (adult$RW1-min(adult$RW1))/(max(adult$RW1)-min(adult$RW1))
plot(adult$line_dist,adult$shape)
plot(adult$size,adult$shape)

# cline with sex and size (using the cline.sex.env function - be is the size effect)

theta.init <- list(centre=90,w=10,crab=0.3,wave=0.8,zs_c=0.1,zs_w=0.1,be_c=0,be_w=0,sc=0.2,sh=0.2,sw=0.2)

mle.shape.sex.size <- mle2(cline.sex.env, theta.init, data=list(position=adult$line_dist,phen=adult$shape, sex=adult$sex_dissection, env=adult$size))
summary(mle.shape.sex.size)
AIC(mle.shape.sex.size)

#compare to the fit without sex just by fixing zs_c=zs_w=0
theta.init <- list(centre=90,w=10,crab=0.3,wave=0.8,zs_c=0,zs_w=0,be_c=0,be_w=0,sc=0.2,sh=0.2,sw=0.2)

mle.shape.sex.size <- mle2(cline.sex.env, theta.init, 
                           fixed=list(zs_c=0,zs_w=0),
                           data=list(position=adult$line_dist,phen=adult$shape, sex=adult$sex_dissection, env=adult$size))
summary(mle.shape.sex.size)
AIC(mle.shape.sex.size)

#consider the fit with size and hos, not sex, using the cline.sex.envh function

theta.init <- list(centre=90,w=10,crab=0.3,wave=0.8,zs_c=0,zs_w=0,bh_c=0,bh_w=0,be_c=0,be_w=0,sc=0.2,sh=0.2,sw=0.2)

mle.shape.size.hos <- mle2(cline.sex.envh, theta.init, 
                           fixed=list(zs_c=0,zs_w=0),
                           data=list(position=adult$line_dist,phen=adult$shape, sex=adult$sex_dissection, env=adult$size, height=adult$hos))
summary(mle.shape.size.hos)
AIC(mle.shape.size.hos)

#check that size effect is significant
theta.init <- list(centre=90,w=10,crab=0.3,wave=0.8,zs_c=0,zs_w=0,be_c=0,be_w=0,sc=0.2,sh=0.2,sw=0.2)

mle.shape <- mle2(cline.sex.env, theta.init, 
                           fixed=list(zs_c=0,zs_w=0,be_c=0,be_w=0),
                           data=list(position=adult$line_dist,phen=adult$shape, sex=adult$sex_dissection, env=adult$size))
summary(mle.shape)
AIC(mle.shape)

#check model without variance change in centre
cline.sex.env.2v <- function(phen,position,sex,env,centre,w,crab,wave,zs_c,zs_w,be_c,be_w,sc,sw){
  
  d <- position-centre
  
  p_x <- 1/(1+exp(0-4*(d)/w)) 
  
  # p_x is frequency cline as in HZAR NB p=0 for left, 1 for right (corresponding to scaling phenotypes with crab=0, wave=1)
  
  z_x <- crab + (wave-crab)*p_x  
  z_x[sex=="female"] <- z_x[sex=="female"]+zs_c + (zs_w-zs_c)*p_x[sex=="female"]
  z_x <- z_x+env*(be_c+(be_w-be_c)*p_x)
  # z_x is expected phenotype, crab-wave always positive, given scaling
  # female offset depends on position in cline
  # env regression depends on position in cline
  
  s_x <- sqrt(sc^2 +  (p_x)*(sw^2-sc^2))
  
  minusll <- -sum(dnorm(phen,z_x,s_x,log=T))
  
  return(minusll)
  
}

theta.init <- list(centre=90,w=10,crab=0.3,wave=0.8,zs_c=0,zs_w=0,be_c=0,be_w=0,sc=0.2,sw=0.2)

mle.shape.size.2v <- mle2(cline.sex.env.2v, theta.init, 
                           fixed=list(zs_c=0,zs_w=0),
                           data=list(position=adult$line_dist,phen=adult$shape, sex=adult$sex_dissection, env=adult$size))
summary(mle.shape.size.2v)
AIC(mle.shape.size.2v)

# or with constant variance
cline.sex.env.1v <- function(phen,position,sex,env,centre,w,crab,wave,zs_c,zs_w,be_c,be_w,s_x){
  
  d <- position-centre
  
  p_x <- 1/(1+exp(0-4*(d)/w)) 
  
  # p_x is frequency cline as in HZAR NB p=0 for left, 1 for right (corresponding to scaling phenotypes with crab=0, wave=1)
  
  z_x <- crab + (wave-crab)*p_x  
  z_x[sex=="female"] <- z_x[sex=="female"]+zs_c + (zs_w-zs_c)*p_x[sex=="female"]
  z_x <- z_x+env*(be_c+(be_w-be_c)*p_x)
  # z_x is expected phenotype, crab-wave always positive, given scaling
  # female offset depends on position in cline
  # env regression depends on position in cline
  
  
  minusll <- -sum(dnorm(phen,z_x,s_x,log=T))
  
  return(minusll)
  
}

theta.init <- list(centre=90,w=10,crab=0.3,wave=0.8,zs_c=0,zs_w=0,be_c=0,be_w=0,s_x=0.1)

mle.shape.size.1v <- mle2(cline.sex.env.1v, theta.init, 
                          fixed=list(zs_c=0,zs_w=0),
                          data=list(position=adult$line_dist,phen=adult$shape, sex=adult$sex_dissection, env=adult$size))
summary(mle.shape.size.1v)
AIC(mle.shape.size.1v)

# try to visualise this by plotting shape-size relationships either side of the cline and cline for size-corrected shape

plot(adult$size[adult$line_dist<70],adult$shape[adult$line_dist<70],xlim=c(0,1),ylim=c(0,1))
points(adult$size[adult$line_dist>75],adult$shape[adult$line_dist>75],col="blue")
points(adult$size[adult$line_dist>100],adult$shape[adult$line_dist>100],col="blue",pch=19)

#plot the cline model

pars <- coef(mle.shape.size.1v) # get coefficients for fitted model
pos <- c(0:160) # ordered distance vector

dist <- pos-pars[1]

freq <- 1/(1+exp(0-4*(dist)/pars[2]))
b_pos <- pars[7]+(pars[8]-pars[7])*freq
mean <- pars[3]+(pars[4]-pars[3])*freq 
mean <- mean+b_pos*0.5 # fitted phenotype mean for size 0.5


distx <- adult$line_dist-pars[1]

freqx <- 1/(1+exp(0-4*(distx)/pars[2]))
b_x <- pars[7]+(pars[8]-pars[7])*freqx
corr_shape <- adult$shape-b_x*(adult$size-0.5)
mean_x <- pars[3]+(pars[4]-pars[3])*freqx+b_x*0.5

plot(adult$line_dist,corr_shape,col=adult$sex_dissection, pch=19, xlab="Distance(m)", ylab="Scaled shape (at scaled size = 0.5)")
lines(pos, mean, col="black",lwd=2)

# save residuals for association analysis (after size correction and cline fit)
adult$RW1_res <- corr_shape-mean_x

###############
## fit cline for environmental PC1
###############

cline.1v <- function(phen,position,centre,w,crab,wave,s_x){
  
  d <- position-centre
  
  p_x <- 1/(1+exp(0-4*(d)/w)) 
  
  # p_x is frequency cline as in HZAR NB p=0 for left, 1 for right (corresponding to scaling phenotypes with crab=0, wave=1)
  
  z_x <- crab + (wave-crab)*p_x  
  # z_x is expected phenotype, wave-crab always positive, given scaling
  
  minusll <- -sum(dnorm(phen,z_x,s_x,log=T))
  
  return(minusll)
  
}
theta.init <- list(centre=70,w=1,crab=-1.5,wave=0.5,s_x=0.1)

mle.cline.PC1 <- mle2(cline.1v, theta.init, data=list(position=adult$line_dist,phen=adult$PC1))
summary(mle.cline.PC1)
AIC(mle.cline.PC1)

plot(adult$line_dist,adult$PC1,pch=19,xlab="Distance(m)",ylab="Environment (PC1)")

###############
# no longer in use 2_2_17
#
# new function for size effect on shape (or similar)
# with separate regressions fitted only using 
# points far from the cline
#
###############



cline.cov <- function(phen,position,cov,centre,slope,crab,wave,sc,sh,sw){
  
  d <- position-centre
  
  # NB slope is negative!

  be_w <- cov(phen[(d*slope/4)<(-2)],cov[(d*slope/4)<(-2)])/var(cov[(d*slope/4)<(-2)])#regression on crab side >2 cline widths from centre

  be_c <- cov(phen[(d*slope/4)>2],cov[(d*slope/4)>2])/var(cov[(d*slope/4)>2])#regression on wave side >2 cline widths from centre
  
  p_x <- exp(slope*d)/(1+exp(slope*d)) 
  
  # p_x is frequency cline as in Cfit NB p=0 for right, 1 for left (corresponding to scaling phenotypes with crab=1, wave=0)
  # frequency slope is then negative
  
  be <- (be_w+(be_c-be_w)*p_x) # dependence of phenotype on covariate at local position
  phen_corr <- phen - (cov-mean(cov))*be # adjust phenotype to mean of covariate
  
  z_x <- wave + (crab-wave)*p_x  
  
  # z_x is expected phenotype (corrected for covariate), crab-wave always positive, given scaling
  
  s_x <- sqrt((p_x^2)*sc^2 + 2*p_x*(1-p_x)*sh^2 + ((1-p_x)^2)*sw^2)
  # unimodal variance model as in Cfit, sx is SD for wave/hybrid/crab (but I assume variances are additive)
  
  minusll <- -sum(dnorm(phen_corr,z_x,s_x,log=T))
  
  return(minusll)
  
}

theta.init <- list(centre=80,slope=-0.6,crab=0.8,wave=0.3,sc=0.2,sh=0.2,sw=0.2)

# for RW1 with cov=size
mle.RW1.size <- mle2(cline.cov, theta.init, data=list(position=adult$line_dist,phen=adult$shape, cov=adult$size))
summary(mle.RW1.size)
AIC(mle.RW1.size)
#print fitted regressions
pars <- coef(mle.RW1.size)
sc_dist <- (adult$line_dist-pars[1])*pars[2]/4
print("Regression in crab end")
b_crab <- cov(adult$shape[sc_dist>2],adult$size[sc_dist>2])/var(adult$size[sc_dist>2])
b_crab
print("Regression in wave end")
b_wave <- cov(adult$shape[sc_dist<(-2)],adult$size[sc_dist<(-2)])/var(adult$size[sc_dist<(-2)])
b_wave

#plot fit
p <- exp(pars[2]*(adult$line_dist-pars[1]))/(1+exp(pars[2]*(adult$line_dist-pars[1])))
b <- (b_wave+(b_crab-b_wave)*p) 
shape_corr <- adult$shape - (adult$size-mean(adult$size))*b 
plot(adult$line_dist,shape_corr)
fit <- pars[4]+(pars[3]-pars[4])*p
points(adult$line_dist,fit,col="blue")



#######################
#### colour clines ####
#######################

# function to fit a phenotype frequency, such as shell colour coded as present=1, absent=0
cline.f <- function(shell,position,centre,width,crab,wave){
  
  d <- position-centre
  
  p_x <- 1/(1+exp(0-4*(d/width)))  #exp(slope*d)/(1+exp(slope*d)) 
  
  z_x <- crab + (wave-crab)*p_x  
  # z_x is expected phenotype, wave-crab always positive, given scaling
  
  minusll <- -sum(dbinom(shell,1,z_x,log=T))
  
  return(minusll)
  
}



shell_col <- c("Beige","Dark_beige","Banded","Black","Other")

# make 1,0 variable for presence/absence  of common types then plot
adult$Banded <- as.numeric(grepl("Banded",adult$Colour))
adult$Black <- as.numeric(adult$Colour=="Black")
adult$Beige <- as.numeric(adult$Colour=="Beige")
adult$Dark_beige <- as.numeric(adult$Colour=="Dark_beige")
adult$Other <- as.numeric((adult$Banded+adult$Black+adult$Beige+adult$Dark_beige)==0)

for (sc in shell_col){
  plot(adult$line_dist,adult[,sc], main=sc, xlab="Distance(m)", ylab="Colour class frequency")

  # add moving average frequency
  for (j in 1:length(adult$Colour)){adult$ma[j]<-mean(adult[adult$line_dist<adult$line_dist[j]+5 & adult$line_dist>adult$line_dist[j]-5,sc],na.rm=T)}
  points(adult$line_dist,adult$ma, col="green", pch=20)
}

#fit clines 
theta.init <- list(centre=95,width=5,crab=0.8,wave=0.05)

mle.cline.f <- mle2(cline.f, theta.init, method="L-BFGS-B",
                    upper=list(centre=150,width=0.5,crab=0.999,wave=0.999),
                    lower=list(centre=1,width=100,crab=0.001,wave=0.001),
                    control=list(parscale=abs(unlist(theta.init))), 
                    data=list(position=adult$line_dist,shell=adult$Beige))
summary(mle.cline.f)
AIC(mle.cline.f)


dist <- 1:150
ci <- as.numeric(coef(mle.cline.f)[1])
wi <- as.numeric(coef(mle.cline.f)[2])
pc <- as.numeric(coef(mle.cline.f)[3])
pw <- as.numeric(coef(mle.cline.f)[4])
cfit <- pc+(pw-pc)*1/(1+exp(0-4*(dist-ci)/wi))
sx <- (4*(pw-pc)*exp(4*(dist-ci)/wi))/(wi*(1+exp(4*(dist-ci)/wi))^2) # slope at dist
lines(dist,cfit,col="green",lwd=2)

# output for association analysis

write.csv(adult[,c("snail_ID","sex_dissection","line_dist","hos","size","CS_res","shape","RW1_res","Banded","Black","Beige","Dark_beige")],
          file="ANG_assoc_phen.csv")


###################
# unmodified old stuff below here
###################

plot(ANG$length.mm,ANG$RW1,col=ANG$sex)
legend('topright', pch = 19, col = c(1:8), legend = c('female', 'juvenile', 'male'),cex=0.6)

#Fit cline to sizes
#first set up a tanh function
#pars is the parameter list - centre c, width w, parent means p1,p2
cline <- function(c,w,p1,p2,d){
  cdist <- (2*(d-c)/w)
  y <- (p1+((p2-p1)*(0.5-0.5*(tanh(cdist)))))
  return(y)
}
#starting values
c <- 0
w <- 5
p1 <- 3
p2 <- 12

test <- cline(st,ANG$dist)
plot(ANG$dist,test)


clinefit <- nls(length.mm ~ cline(c,w,p1,p2,dist),start=list(c=0,w=5,p1=3,p2=12),data=ANG)

clinefit

plot(ANG$dist,fitted(clinefit))

plot(ANG$dist,ANG$length.mm)
lines(ANG$dist,fitted(clinefit),col="red")

clineres <- residuals(clinefit)

#cline fitting with sex offset
ANG$female <- as.numeric(ANG$sex)
#female 1, juvenile 2, male 3
ANGnj <- ANG[ANG$sex!="juvenile",]
clinesex <- function(a,c,w,p1,p2,d,fem){
  cdist <- (2*(d-c)/w)
  y <- a*fem+(p1+((p2-p1)*(0.5-0.5*(tanh(cdist)))))
  return(y)
}
#starting values
a <- -5
c <- 0
w <- 5
p1 <- 3
p2 <- 12


clinefit.sex <- nls(length.mm ~ clinesex(a,c,w,p1,p2,dist,female),start=list(a=-5,c=0,w=5,p1=3,p2=12),data=ANGnj)

clinesexres <- residuals(clinefit.sex)
plot(ANGnj$dist,fitted(clinefit.sex))
ANGnj$fit <- fitted(clinefit.sex)

plot(ANGnj$dist[ANGnj$sex=="male"],ANGnj$length.mm[ANGnj$sex=="male"],col="blue",xlab="Distance",ylab="Length")
points(ANGnj$dist[ANGnj$sex=="female"],ANGnj$length.mm[ANGnj$sex=="female"],col="red")
points(ANGnj$dist[ANGnj$sex=="male"],ANGnj$fit[ANGnj$sex=="male"],col="blue",pch=16)
points(ANGnj$dist[ANGnj$sex=="female"],ANGnj$fit[ANGnj$sex=="female"],col="red",pch=16)

dev.print(pdf, file="ClineBySex.pdf", width=14, height=5)




#add environmental variable to cline fit residuals

cline.env <- lm(clineres~z+sex+Rockiness+Fucusness+Barnacleness, data=ANG)
cline.env
anova(cline.env)
drop1(cline.env)
cline.env <- lm(clineres~z+sex+Rockiness+Fucusness, data=ANG)

boxplot(clineres~sex,data=ANG)
plot(ANG$z,clineres)
plot(ANG$Rockiness,clineres)
plot(ANG$Fucusness,clineres)

clinesex.env <- lm(clinesexres~z+Rockiness+Fucusness+Barnacleness, data=ANGnj)
clinesex.env
anova(clinesex.env)

hist(clinesexres)


#stuff for using sp or spatstat
coords <- ANGxyzp[c(4,5)]
phens <- ANGxyzp[c(2,3,6)]
snail <- ANGxyzp[c(1)]

snail.sp <- SpatialPoints(coords)
plot(snail.sp, pch=1)
snail.spdf <- SpatialPointsDataFrame(snail.sp, phens)
spplot(snail.spdf, "sex")
bubble(snail.spdf, "length.mm")
spplot(snail.spdf, "z")

area <- convexhull.xy(coords)
plot(area)
area <- ripras(coords)


#############################################
#cline fitting, using mle and incorporating sex etc



ANG$female <- as.numeric(ANG$sex)
ANGnj <- ANG[ANG$sex!="juvenile",]
ANGnj <- na.omit(ANGnj)
plot(ANGnj$dist,ANGnj$RW1, col=ANGnj$col_rock, pch=16, cex=0.9,xlab="Distance (m)",ylab="RW1")

ANGnj$z <- ANGnj$z-mean(ANGnj$z) #rescale height to zero mean

#females 0, males 1
ANGnj$sexi <- ANGnj$female
ANGnj$sexi[ANGnj$sexi==1] <- 0
ANGnj$sexi[ANGnj$sexi==3] <- 1

#write out this data set
write.table(ANGnj, file="ANG_nojuv.txt", quote=F, row.names=F, append=F)


# 1 centres, 3 variances, plus z, plus sex

cline1c3s.z.sex <- function(phen,position,height,sex,c1,slope,diff,lower,sl,sh,sr,bz,bs){
  
  d <- position-c1
    
  p_x <- exp(slope*d)/(1+exp(slope*d)) # p_x is frequency cline as in Cfit NB p=0 for right, 1 for left
  
  # to ensure frequency slope is negative
    
  z_x <- bs*sex+bz*height+lower+diff*p_x  # z_x is expected phenotype
  z_x[diff<0] <- bs*sex+bz*height+lower-diff*p_x
  
  s_x <- (p_x^2)*sl + 2*p_x*(1-p_x)*sh + ((1-p_x)^2)*sr # unimodal variance model as in Cfit, sx is SD for wave/hybrid/crab
  
  minusll <- -sum(dnorm(phen,z_x,s_x,log=T))
  
  return(minusll)
  
}

theta.init <- list(c1=1,slope=-0.15,lower=1.5,diff=2.5,sl=0.5,sh=0.5,sr=0.5,bz=-0.5,bs=-0.1)

mle.cline1c3s.z.sex <- mle2(cline1c3s.z.sex,theta.init,data=list(position=ANGnj$dist,phen=log(ANGnj$length.mm),height=ANGnj$z,sex=ANGnj$sexi))

summary(mle.cline1c3s.z.sex)
AIC(mle.cline1c3s.z.sex)

pars <- coef(mle.cline1c3s.z.sex)
(pars[5]/4+pars[6]/2+pars[7]/4)^2-(pars[5]/2+pars[7]/2)^2

## fitting clines with fixed centres to check coincidence

theta.init <- list(c1=-3.6,slope=0,wave=2.49,crab=1.5,sw=0.25,sh=0.25,sc=0.25,bz=0,bs=0)

for (fix.c1 in seq(-5,5,by=0.5)) {
  
 mle.cline1c3s.z.sex <- mle2(cline1c3s.z.sex,theta.init,fixed=list(c1=fix.c1),data=list(position=ANGnj$dist,phen=log(ANGnj$length.mm),height=ANGnj$z,sex=ANGnj$sexi))
 print(fix.c1)
 print(logLik(mle.cline1c3s.z.sex))

}
############################################################################
# 1 centres, 2 variances (no central elevation), plus z, plus sex

cline1c2s.z.sex <- function(phen,position,height,sex,c1,slope,wave,crab,sw,sc,bz,bs){
  
  d <- position-c1
  
  p_x <- exp(slope*d)/(1+exp(slope*d)) # p_x is frequency cline as in Cfit NB p=0 for wave, 1 for crab
  
  z_x <- bs*sex+bz*height+crab+(wave-crab)*p_x  # z_x is expected phenotype 
  
  s_x <- (p_x)*sw + (1-p_x)*sc # unimodal variance model as in Cfit, no central elevation, sx is SD for wave/crab
  
  minusll <- -sum(dnorm(phen,z_x,s_x,log=T))
  
  return(minusll)
  
}

theta.init <- list(c1=-10,slope=0.1,wave=2.5,crab=1.5,sw=0.5,sc=0.5,bz=0,bs=0)

mle.cline1c2s.z.sex <- mle2(cline1c2s.z.sex,theta.init,data=list(position=ANGnj$dist,phen=log(ANGnj$length.mm),height=ANGnj$z,sex=ANGnj$sexi))

summary(mle.cline1c2s.z.sex)
AIC(mle.cline1c2s.z.sex)


################################################################################

# 1 centres, 3 variances, plus z, plus sex - asymmatrical

acline1c3s.z.sex <- function(phen,position,height,sex,c1,slope.left,slope.right,diff,lower,sl,sh,sr,bz,bs){
  
  d <- position-c1
  
  p_x <- exp(slope.left*d)/(1+exp(slope.left*d)) # p_x is frequency cline as in Cfit NB p=0 for right, 1 for left
  p_x[d>0] <- exp(slope.right*d[d>0])/(1+exp(slope.right*d[d>0]))
  
  # to ensure frequency slope is negative
  
  z_x <- bs*sex+bz*height+lower+diff*p_x  # z_x is expected phenotype
  z_x[diff<0] <- bs*sex+bz*height+lower-diff*p_x
  
  s_x <- (p_x^2)*sl + 2*p_x*(1-p_x)*sh + ((1-p_x)^2)*sr # unimodal variance model as in Cfit, sx is SD for wave/hybrid/crab
  
  minusll <- -sum(dnorm(phen,z_x,s_x,log=T))
  
  return(minusll)
  
}

theta.init <- list(c1=-3,slope.left=-0.15,slope.right=-0.15,lower=-0.08,diff=0.13,sl=0.05,sh=0.05,sr=0.05,bz=0.02,bs=0.02)

mle.acline1c3s.z.sex <- mle2(acline1c3s.z.sex,theta.init,data=list(position=ANGnj$dist,phen=ANGnj$RW1,height=ANGnj$z,sex=ANGnj$sexi))

summary(mle.acline1c3s.z.sex)
AIC(mle.acline1c3s.z.sex)

############################################################################
# 1 centres, symmetrical step, 3 variances, plus z, plus sex

## NB must have deccreasing frequency cline

sstep1c3s.z.sex <- function(phen,position,height,sex,c1,slope,lower,diff,sl,sh,sr,bz,bs,cd,ef){
  
  d <- position-c1
  
  p_x <- exp(slope*d)/(1+exp(slope*d)) # p_x is frequency cline in centre as in Cfit NB p=0 for wave, 1 for crab, therefore negative slope if crab on left
  p_x[d>cd] <- (exp(slope*cd)/(1+exp(slope*cd)))*exp(-((-slope*ef)*(d[d>cd]-cd))/(1+exp(slope*cd))) #right tail, cd=distance from centre to start of tail, ef=slope in tail
  p_x[d<(0-cd)] <- 1-(1-exp(-slope*cd)/(1+exp(-slope*cd)))*exp(-((slope*ef)*(d[d<(0-cd)]+cd))/(1+exp(slope*cd)))  #left tail
  
  z_x <- bs*sex+bz*height+lower+diff*p_x  # z_x is expected phenotype
  z_x[diff<0] <- bs*sex+bz*height+lower-diff*p_x # to ensure freq cline has negative slope
  
  s_x <- (p_x^2)*sl + 2*p_x*(1-p_x)*sh + ((1-p_x)^2)*sr # unimodal variance model as in Cfit, sx is SD for wave/hybrid/crab
  
  minusll <- -sum(dnorm(phen,z_x,s_x,log=T))
  
  return(minusll)
  
}

theta.init <- list(c1=1,slope=-0.15,lower=1.5,diff=1,sl=0.2,sh=0.4,sr=0.2,bz=-0.5,bs=-0.1,cd=10,ef=1)

mle.sstep1c3s.z.sex <- mle2(sstep1c3s.z.sex,theta.init,data=list(position=ANGnj$dist,phen=log(ANGnj$length.mm),height=ANGnj$z,sex=ANGnj$sexi))

summary(mle.sstep1c3s.z.sex)
AIC(mle.sstep1c3s.z.sex)


pars <- coef(mle.sstep1c3s.z.sex)
(pars[5]/4+pars[6]/2+pars[7]/4)^2-(pars[5]/2+pars[7]/2)^2


#####################

# plot fit (1 centre)

pars <- coef(mle.cline1c3s.z.sex) # get coefficients for fitted model
pos <- c(-80:80) # ordered distance vector

dist <- pos-pars[1]

freq <- exp(pars[2]*dist)/(1+exp(pars[2]*dist))
freq <- 1-freq
meanF <- pars[3]+(pars[4]-pars[3])*freq  # fitted phenotype mean - female
meanM <- pars[9]+pars[3]+(pars[4]-pars[3])*freq  # fitted phenotype mean - male

SD <- (freq^2)*pars[7] + 2*freq*(1-freq)*pars[6] + ((1-freq)^2)*pars[5] # unimodal variance model as in Cfit, sx is SD for wave/hybrid/crab
mean.max <- pmax(meanF,meanM)
mean.min <- pmin(meanF,meanM)
plusSD <- mean.max+SD
minusSD <- mean.min-SD


plot(ANGnj$dist,ANGnj$RW1,xlab="Distance (m)",ylab="Shape - RW1")
lines(pos, meanM, col="blue")
lines(pos, meanF, col="red")
lines(pos, plusSD, col="green")
lines(pos, minusSD, col="green")

############################
#plot step function

sstep.freq <- function(position,c1,slope,cd,ef){
  
  d <- position-c1
  
  p_x <- exp(slope*d)/(1+exp(slope*d)) # p_x is frequency cline in centre as in Cfit NB p=0 for wave, 1 for crab
  p_x[d>cd] <- (exp(slope*cd)/(1+exp(slope*cd)))*exp(-((-slope*ef)*(d[d>cd]-cd))/(1+exp(slope*cd))) #right tail, cd=distance from centre to start of tail, ef=slope in tail
  p_x[d<(0-cd)] <- 1-(1-exp(-slope*cd)/(1+exp(-slope*cd)))*exp(-((slope*ef)*(d[d<(0-cd)]+cd))/(1+exp(slope*cd)))  #left tail
  
  
  return(p_x)
  
}

pp<- sstep.freq(ANGnj$dist,0.6,-0.12,11.17,3.43)
plot(ANGnj$dist,pp,xlab="Distance (m)", ylab="Frequency (RW1 fit)")
plot(ANGnj$dist,log(pp/(1-pp)))



############################
#variances and covariances along the shore

ANGnj$dist_cut <- cut(ANGnj$dist,c(-100,-40,-30,-20,-10,0,5,10,20,30,40,50,80))
ANGnj$log_length <- log(ANGnj$length.mm)

RW1_var <- aggregate(ANGnj$RW1, by=list(ANGnj$dist_cut), FUN="var",na.rm=T)
dist_mean <- aggregate(ANGnj$dist, by=list(ANGnj$dist_cut), FUN="mean")

tmp <- split(ANGnj,ANGnj$dist_cut)
RW1_len_cov <- as.data.frame(lapply(tmp,function(x) cov(x[,8],x[,13],use="pairwise")))

plot(dist_mean[,2], RW1_var[,2], xlab="Shore position", ylab="RW1 variance")
plot(dist_mean[,2], -RW1_len_cov[1,], xlab="Shore position", ylab="RW1_length covariance")
lines(dist_mean[,2], -RW1_len_cov[1,],col="blue")
