########################################################################################
############################## Step 0: Load data #######################################
########################################################################################

#clear any already loaded workspace
rm(list=ls()) 

library(tidyr)
library(ggplot2)
library(plyr)
library(scales)
library(caret)


#Set the working directory
setwd("C:/Users/Vanil/Documents/R/Rscripts/Dakar") ###DO UPDATE TO YOUR OWN WORKING DIRECTORY

# Load Sample data
Sample.Data <-read.csv("SampleData.csv", header = TRUE)

# In order to spot probems, look at dataset summary:
summary(Sample.Data)


# Make sure data is integer or numeric
str(Sample.Data)
Sample.Data$ID <- as.numeric(Sample.Data$ID)
Sample.Data$NAME_4 <- as.numeric(Sample.Data$NAME_4)
Sample.Data$Contamination <- as.numeric(Sample.Data$Contamination)
Sample.Data$Repeat <- as.numeric(Sample.Data$Repeat)
Sample.Data$Date <- as.numeric(Sample.Data$Date)
Sample.Data$Time <- as.numeric(Sample.Data$Time)
#And check again:
str(Sample.Data)


########################################################################################
############################## Step 1: ESDA ############################################
########################################################################################


##### Descriptive stats #####
library(Hmisc)
#describe() provides min, max, mean and median
describe(Sample.Data$LogTTC)
describe(Sample.Data$TTC)
describe(Sample.Data$TLF)
describe(Sample.Data$FC)
sapply(Sample.Data, sd, na.rm=TRUE) #get all standard deviations


##### Density plots #####

# Density plot of TLF
ggplot(Sample.Data) + geom_density(aes(x=TLF)) + 
  labs(x="TLF real-time reading", y = "Density", title = "Density plot of TLF")


# Density plot of TTC
ggplot(Sample.Data) + geom_density(aes(x=TTC)) + 
  labs(x="Thermotolerant coliforms (in CFU/100mL)", y = "Density", title = "Density plot of TTC")

# Log10 to restore symmetry in the TTC distribution
ggplot(Sample.Data) + geom_density(aes(x=LogTTC))+ 
  labs(x="Log(TTC)", y = "Density", title = "Density plot of Log(TTC)")


##### Boxplots #####

# Box plot for y = TLF and x = TTC grouped like on a log10 scale: 0, 1-9, 10-99, 100-999, 1000+ :
# Grouping TTC values by magnitude
Sample.Data$TTCBreaks <- cut(Sample.Data$TTC, breaks = c(-Inf, 0, 1, 10, 100, 1000, Inf),right = FALSE)
Sample.Data$TTCBreaks <- revalue(Sample.Data$TTCBreaks, c("[-Inf,0)"= "0", "[0,1)" = "0-1", "[1,10)"= "1-10", "[10,100)"= "10-100","[100,1e+03)"="100-1000", "[1e+03, Inf)"= "1000+"))
#Plotting with ggplot
p <- ggplot(Sample.Data, aes(x = TTCBreaks ,  y = TLF))
p + geom_boxplot() + labs(x="Thermotolerant coliforms (in CFU/100mL)", 
                          y = "Tryptophan-like fluorescence (in QSU)",
                          title = "Distribution of TLF values by Thermotolerant coliform results") 


# Boxplot of TLF against CDOM
Sample.Data$CDOMBreaks <- cut(Sample.Data$CDOM, breaks = c(-Inf, 0, 50, 100, 150, Inf),right = FALSE)
Sample.Data$CDOMBreaks <- revalue(Sample.Data$CDOMBreaks, c("[-Inf,0)"= "0", "[0,50)" = "0-50", "[50,100)"= "50-100", "[100,150)"= "100-150","[150, Inf)"= "150+"))
p2 <- ggplot(Sample.Data, aes(x = CDOMBreaks ,  y = TLF))
p2 + geom_boxplot() + labs(x="CDOM real-time reading", 
                           y = "TLF real-time reading",
                           title = "Distribution of TLF values by CDOM results")


# Plot relationship between TLF and LogTTC (and CDOM) by source type:
ggplot(data = Sample.Data) +  
  geom_point(mapping = aes(x = TLF, y = LogTTC, colour = CDOM)) +
  facet_wrap(~ Type, nrow = 2)




##### Correlation matrix #####

#We select all relevant and non-redundant variables for our correlation matrix. Non-parametric correlation
# using Spearman rho is chosen over Pearson coefficient due to the distribution of the data.
res <- cor(Sample.Data[,7:40], method = "spearman", use = "complete.obs")
round(res,2)

# rcorr Computes a matrix of Pearson's r or Spearman's rho rank correlation coefficients for all possible 
# pairs of columns of a matrix.
library(Hmisc)
res2 <- rcorr(as.matrix(Sample.Data[,7:40]))

# Using flattenCorrMatrix we aggregate
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
flattenCorrMatrix(res2$r, res2$P)
flatflat <- flattenCorrMatrix(res2$r, res2$P)

#Export the correlation matrix as a table
write.csv(flatflat, "CorrelationMatrix.csv")

library(corrplot)
symnum(res,abbr.colnames = TRUE)
corrplot(res, type = "upper", order = "original", 
         tl.col = "black", tl.srt = 45)

# Insignificant correlations are leaved blank
corrplot(res2$r, type="upper", order="original", 
         p.mat = res2$P, sig.level = 0.01, insig = "blank", tl.col="black")


library("PerformanceAnalytics")
my_data <- Sample.Data[, 7:40]
png("CorrelationChart.png", width = 1000, height = 900)
chart.Correlation(my_data, histogram=TRUE, pch=19)
dev.off()

#Correlation
cor(Sample.Data$LogTTC, Sample.Data$TLF, method = "spearman")



##### Subset analysis #####

# 1. Pre-Rain vs Post-Rain
ggplot(data = Sample.Data) +  
  geom_point(mapping = aes(x = TLF, y = LogTTC, colour = CDOM)) +
  facet_grid(. ~ Sample.Data$Rain) + labs(x = "TLF reading", 
                                          y = "Log(TTC)",
                                          title = "Distribution of TTC values by TLF, before and after the first rain")

PreRain <-read.csv("PreRainSample.csv", header = TRUE)
#describe() provides min, max, mean and median
describe(PreRain$LogTTC)
describe(PreRain$TTC)
describe(PreRain$TLF)
describe(PreRain$FC)
sapply(PreRain, sd, na.rm=TRUE) #get all standard deviations
PostRain <-read.csv("PostRainSample.csv", header = TRUE)
describe(PostRain$LogTTC)
describe(PostRain$TTC)
describe(PostRain$TLF)
describe(PostRain$FC)
sapply(PostRain, sd, na.rm=TRUE) #get all standard deviations

# 2. DOC

#DOC / CDOM plot
ggplot(data = Sample.Data, mapping = aes(x = TLF, y = CDOM)) +  
  geom_point(mapping = aes(colour = DOC)) +
  geom_smooth()

#Plot of DOC against CDOM with TLF as colour
ggplot(data = Sample.Data, mapping = aes(x = DOC, y = CDOM)) +  
  geom_point(mapping = aes(colour = TLF)) +
  geom_smooth()


# 3. Filtered Samples

# TLF filtered vs unfiltered
ggplot(data = Sample.Data, mapping = aes(x = TLF, y = TLF_filtered)) +  
  geom_point(mapping = aes(colour = TLF)) +
  geom_smooth()

#CDOM filtered vs unfiltered
ggplot(data = Sample.Data, mapping = aes(x = CDOM, y = CDOM_filtered)) +  
geom_point(mapping = aes(colour = CDOM)) +
  geom_smooth()


# 4. Nitrates & phosphates
# Correlation matrix shows no correlation exists between nitrates or phosphates and other variables.
#plotting Nitrates vs TLF, with CDOM as the colour
ggplot(data = Sample.Data, mapping = aes(x = Nitrates, y = TLF)) +
  geom_point(mapping = aes(colour = CDOM)) +
  geom_smooth()

ggplot(data = Sample.Data, mapping = aes(x = Phosphates, y = TLF)) +
  geom_point(mapping = aes(colour = CDOM)) +
  geom_smooth()


########################################################################################
########################## Step 2: Stepwise regression #################################
########################################################################################


#Check if there are NAs
library(Amelia)
library(mlbench)
missmap(Sample.Data, col=c("blue", "red"), legend=FALSE)
#Remove columns with missing data so we can work with all observations for the regression
Sample.Data <- Sample.Data[,1:40]
Sample.Data <- Sample.Data[,-18]
missmap(Sample.Data, col=c("blue", "red"), legend=FALSE)

# Stepwise Regression using original variables and Cross Validation in backward stepwise regression. 
# Contamination is response, other variables are predictors.

#Scaling the data except the response variable and categorical

scaledData = as.data.frame(scale(Sample.Data[,c(7,10,11,12,13,14,15,16,17,18,19,20,21,22,23,39)]))
scaledData <- cbind(Sample.Data[,c(6,8,9,24,25,26,27)],scaledData) # Add column 2 back in
colnames(scaledData)[1] <- "Contamination"
colnames(scaledData)[2] <- "Type"
colnames(scaledData)[3] <- "Rain"
colnames(scaledData)[4] <- "Sanitation"
colnames(scaledData)[5] <- "SepticTank"
colnames(scaledData)[6] <- "Soakpit"
colnames(scaledData)[7] <- "Latrines"

# Perform 5 fold Cross validation
ctrl <- trainControl(method = "repeatedcv", number = 5, repeats = 5)

lmFit_Step <- train(Contamination ~ ., data = scaledData, "lmStepAIC", scope = 
                      list(lower = Contamination~1, upper = Contamination~.), direction = "backward",trControl=ctrl)

# Step:  AIC=-156.1
# .outcome ~ Type + x + Temperature + Salinity + DistanceToCemetery + 
#   DistanceToIndustry + DistanceToLandfill
# 
# Df Sum of Sq    RSS     AIC
# <none>                            16.452 -156.10
# - x                   1   0.36932 16.822 -155.95
# - DistanceToIndustry  1   0.56585 17.018 -154.82
# - Temperature         1   0.57487 17.027 -154.77
# - Salinity            1   0.63317 17.086 -154.44
# - DistanceToLandfill  1   0.91549 17.368 -152.85
# - DistanceToCemetery  1   1.16038 17.613 -151.49
# - Type                1   1.64034 18.093 -148.88

#Fitting a new model with these 7 variables
mod_Step = glm(Contamination ~ Type + x + Temperature + Salinity + DistanceToCemetery + 
                 DistanceToIndustry + DistanceToLandfill, data = scaledData)
summary(mod_Step)


# Test model's predictive power.
# tests based on https://www.datacamp.com/community/tutorials/logistic-regression-R
glm.probs <- predict(mod_Step,type = "response")
glm.probs[1:5]
# Ween does the model accurately predict contamination when it's actually contaminated?
glm.pred <- ifelse(glm.probs > 0.5, "1", "0")
attach(scaledData)
table(glm.pred,Contamination)
mean(glm.pred == Contamination)

# We don't run the usual plot summary like Residuals vs Fitted, Normal QQ etc because
# they are not relevant to LOGISTIC regression and would be really hard to interpret.
# Instead we use a binned plot.
library(arm)
old.par <- par(no.readonly = TRUE)
x <- predict(mod_Step)
y <- resid(mod_Step)

png("BinnedPlot.png", width = 762, height = 462)
binnedplot(x,y)
dev.off()
par(old.par)
mod_Step$deviance
par(mfrow=c(1,1))




########################################################################################
############################ Step 3: Interpolation #####################################
########################################################################################

## This method follows the tutorial developed by Guy Lansley and James Cheshire for the CDRC
## available at https://data.cdrc.ac.uk/tutorial/aa5491c9-cbac-4026-97c9-f9168462f4ac/81eada9f-18a9-4311-86e4-b68b9ead7be9

# load libraries
library(sp)
library(rgdal)
library(rgeos)
library(tmap)
library(spatstat)  
library(maptools)
library(raster)
library(rasterVis)
library(gstat)

# Load study area shapefile
Study.Area <- readOGR(".", "SelectedRegions")

# load sample point shapefile
Sample.Points <- readOGR(".", "SamplePoints")




##### Autocorreation tests #####


# First, explore spatial autocorrelation patterns to determine 
# which interpolation method to adopt

library(ape)

# Moran's I
sample.dists <- as.matrix(dist(cbind(Sample.Data$x, Sample.Data$y)))
sample.dists.inv <- 1/sample.dists
diag(sample.dists.inv) <- 0
sample.dists.inv[1:5, 1:5]
sample.dists.inv[is.infinite(sample.dists.inv)] <- 0
Moran.I(Sample.Data$TLF, sample.dists.inv)
Moran.I(Sample.Data$TTC, sample.dists.inv)

# Mantel test
station.dists <- dist(cbind(Sample.Data$x, Sample.Data$y))
sample.dists2 <- dist(Sample.Data$TTC)
sample.dists3 <- dist(Sample.Data$TLF)
as.matrix(station.dists)[1:5, 1:5]
as.matrix(sample.dists2)[1:5, 1:5]
mantel.rtest(station.dists, sample.dists2, nrepet = 9999)
mantel.rtest(station.dists, sample.dists3, nrepet = 9999)

# semivariogram
v <- variogram(LogTTC~1, Sample.Points)
fitv <- fit.variogram(v, vgm(model="Nug"))
plot(v, model=fitv, as.table=TRUE)



gModel <- gstat(NULL, id="G0yKrig",formula= LogTTC~1,locations=Sample.Points, model=fitv)

# Create a tessellated surface
dat.pp <- as(dirichlet(as.ppp(Sample.Points)), "SpatialPolygons")
dat.pp <- as(dat.pp,"SpatialPolygons")

# Project to UTM 28N
proj4string(dat.pp) <- CRS("+init=EPSG:32628")
proj4string(Sample.Points) <- CRS("+init=EPSG:32628")

# Assign to each polygon the data from Sample.Points 
int.Z <- over(dat.pp,Sample.Points, fn=mean) 

# Create a SpatialPolygonsDataFrame from the thiessen polygons
thiessen <- SpatialPolygonsDataFrame(dat.pp, int.Z)


# Map  thiessen polygons and Sample.Points
tm_shape(Study.Area) + tm_fill(alpha=.3, col = "grey") +
  tm_shape(thiessen) +  tm_borders(alpha=.5, col = "black") +
  tm_shape(Sample.Points) + tm_dots(col = "blue", scale = 0.5)

library(raster)
# Crop the Thiessen polygon by study area shape
thiessen.crop <-crop(thiessen, Study.Area)

# Map cropped thiessen polygons and Sample.Points
tm_shape(Study.Area) + tm_fill(alpha=.3, col = "grey") +
  tm_shape(thiessen.crop) +  tm_borders(alpha=.5, col = "black") +
  tm_shape(Sample.Points) + tm_dots(col = "blue", scale = 0.5)

# Map TTC count across thiessen polygons

png("InterpoTTC.png", width = 762, height = 462)
tm_shape(thiessen.crop) + tm_fill(col = "LogTTC", style = "quantile", palette = "Reds", title = "Thermotolerant coliform count (Log10) across the study area") + tm_borders(alpha=.3, col = "black") +
  tm_shape(Sample.Points) + tm_dots(col = "black", scale = 0.5) +
  tm_layout(legend.position = c("left", "top"),  legend.text.size = 1, legend.title.size = 2, frame = FALSE)
dev.off()

library(gstat)
library(xts)

# define sample grid based on the extent of the Sample.Points file
grid <-spsample(Sample.Points, type = 'regular', n = 10000)

# runs the idw for the Price variable of Sample.Points
idw <- idw(Sample.Points$LogTTC ~ 1, Sample.Points, newdata= grid)
idw.output = as.data.frame(idw)
names(idw.output)[1:3] <- c("long", "lat", "prediction")

# Create a spatial points data frame
spg <- idw.output
coordinates(spg) <- ~ long + lat

# Coerce to SpatialPixelsDataFrame
gridded(spg) <- TRUE
# coerce to raster
raster_idw <- raster(spg)

# Set projection to UTM 28N
projection(raster_idw) <- CRS("+init=EPSG:32628")

library(tmap)


tm_shape(raster_idw) + tm_raster("prediction", style = "quantile", n = 100, palette = "Reds",legend.show = FALSE) +
  tm_shape(Study.Area) + tm_borders(alpha=.5) + 
  tm_shape(Sample.Points) + tm_bubbles(size = "LogTTC", col = "LogTTC", palette = "Blues", style = "quantile", legend.size.show = FALSE, title.col = "LogTTC") +
  tm_layout(legend.position = c("left", "top"),  legend.text.size = 1.1, legend.title.size = 1.4, frame = FALSE, legend.bg.color = "white", legend.bg.alpha = 0.5)


# Mask the raster by the STudy area shape
masked_idw <- mask(raster_idw, Study.Area)

# Plot the masked raster
tm_shape(masked_idw) + tm_raster("prediction", style = "quantile", n = 100, legend.show = FALSE) +
  tm_shape(Sample.Points) + tm_bubbles(size = "LogTTC", col = "LogTTC", palette = "Blues", style = "quantile", legend.size.show = FALSE, title.col = "LogTTC") +
  tm_layout(legend.position = c("left", "top"),  legend.text.size = 1, legend.title.size = 1.2, frame = TRUE) +
  tm_shape(Study.Area) +tm_borders(alpha=1, col = "black", group = "NAME_4")

# Export the IDW LogTTC plot
png("INTERPOLogTTC.png", width = 762, height = 462)
tm_shape(masked_idw) + tm_raster("prediction", style = "quantile", n = 100, legend.show = FALSE) +
  tm_shape(Sample.Points) + tm_bubbles(size = "LogTTC", col = "LogTTC", palette = "Blues", style = "quantile", legend.size.show = FALSE, title.col = "LogTTC") +
  tm_layout(legend.position = c("left", "top"),  legend.text.size = 1, legend.title.size = 1.2, frame = TRUE) +
  tm_shape(Study.Area) +tm_borders(alpha=1, col = "black", group = "NAME_4")
dev.off()



##### Now do the same for TLF #####

v <- variogram(TLF~1, Sample.Points)
fitv <- fit.variogram(v, vgm(model="Nug"))
plot(v, model=fitv, as.table=TRUE)

# Export the semivariogram
png("AutocorrSemivariogramTLF.png", width = 762, height = 462)
plot(v, model=fitv, as.table=TRUE)
dev.off()

gModel <- gstat(NULL, id="G0yKrig",formula= TLF~1,locations=Sample.Points, model=fitv)

# Create a tessellated surface
dat.pp <- as(dirichlet(as.ppp(Sample.Points)), "SpatialPolygons")
dat.pp <- as(dat.pp,"SpatialPolygons")

# Project to UTM 28N
proj4string(dat.pp) <- CRS("+init=EPSG:32628")
proj4string(Sample.Points) <- CRS("+init=EPSG:32628")

# Assign to each polygon the data from Sample.Points 
int.Z <- over(dat.pp,Sample.Points, fn=mean) 

# Create a SpatialPolygonsDataFrame from the thiessen polygons
thiessen <- SpatialPolygonsDataFrame(dat.pp, int.Z)


# Map  thiessen polygons and Sample.Points
tm_shape(Study.Area) + tm_fill(alpha=.3, col = "grey") +
  tm_shape(thiessen) +  tm_borders(alpha=.5, col = "black") +
  tm_shape(Sample.Points) + tm_dots(col = "blue", scale = 0.5)

library(raster)
# Crop the Thiessen polygon by study area shape
thiessen.crop <-crop(thiessen, Study.Area)

# Map cropped thiessen polygons and Sample.Points
tm_shape(Study.Area) + tm_fill(alpha=.3, col = "grey") +
  tm_shape(thiessen.crop) +  tm_borders(alpha=.5, col = "black") +
  tm_shape(Sample.Points) + tm_dots(col = "blue", scale = 0.5)

#Map TLF across polygons
tm_shape(thiessen.crop) + tm_fill(col = "TLF", style = "quantile", palette = "Reds", title = "TLF readings acoss the study area") + tm_borders(alpha=.3, col = "black") +
  tm_shape(Sample.Points) + tm_dots(col = "black", scale = 0.5) +
  tm_layout(legend.position = c("left", "top"),  legend.text.size = 0.50, legend.title.size = 0.8, frame = FALSE)


library(gstat)
library(xts)

# define sample grid based on the extent of the Sample.Points file
grid <-spsample(Sample.Points, type = 'regular', n = 10000)

# runs the idw for the Price variable of Sample.Points
idw <- idw(Sample.Points$TLF ~ 1, Sample.Points, newdata= grid)
idw.output = as.data.frame(idw)
names(idw.output)[1:3] <- c("long", "lat", "prediction")

# Create a spatial points data frame
spg <- idw.output
coordinates(spg) <- ~ long + lat

# Coerce to SpatialPixelsDataFrame
gridded(spg) <- TRUE
# coerce to raster
raster_idw <- raster(spg)

# Set projection to UTM 28N
projection(raster_idw) <- CRS("+init=EPSG:32628")

library(tmap)


tm_shape(raster_idw) + tm_raster("prediction", style = "quantile", n = 100, palette = "Reds",legend.show = FALSE) +
  tm_shape(Study.Area) + tm_borders(alpha=.5) + 
  tm_shape(Sample.Points) + tm_bubbles(size = "TLF", col = "TLF", palette = "Blues", style = "quantile", legend.size.show = FALSE, title.col = "TLF") +
  tm_layout(legend.position = c("left", "top"),  legend.text.size = 1.1, legend.title.size = 1.4, frame = FALSE, legend.bg.color = "white", legend.bg.alpha = 0.5)


# Mask the raster by the STudy area shape
masked_idw <- mask(raster_idw, Study.Area)

# Plot the masked raster
tm_shape(masked_idw) + tm_raster("prediction", style = "quantile", n = 100, legend.show = FALSE) +
  tm_shape(Sample.Points) + tm_bubbles(size = "TLF", col = "TLF", palette = "Blues", style = "quantile", legend.size.show = FALSE, title.col = "TLF") +
  tm_layout(legend.position = c("left", "top"),  legend.text.size = 1, legend.title.size = 1.2, frame = TRUE) +
  tm_shape(Study.Area) +tm_borders(alpha=1, col = "black", group = "NAME_4")

# Export the IDW TLF plot
png("IDWTLF.png", width = 762, height = 462)
tm_shape(masked_idw) + tm_raster("prediction", style = "quantile", n = 100, legend.show = FALSE) +
  tm_shape(Sample.Points) + tm_bubbles(size = "TLF", col = "TLF", palette = "Blues", style = "quantile", legend.size.show = FALSE, title.col = "TLF") +
  tm_layout(legend.position = c("left", "top"),  legend.text.size = 1, legend.title.size = 1.2, frame = TRUE) +
  tm_shape(Study.Area) +tm_borders(alpha=1, col = "black", group = "NAME_4")
dev.off()


########################################################################################
######################## Step 4: Hierarchical CLustering ###############################
########################################################################################


# https://www.r-bloggers.com/how-to-perform-hierarchical-clustering-using-r/ 
# From https://www.r-bloggers.com/clustering-mixed-data-types-in-r/
library(RColorBrewer)
library(ade4)
library(fastcluster)
library(plyr)
library(JLutils)
library(questionr)
library(purrr)
library(dplyr) # for data cleaning
library(cluster) # for gower similarity and pam
library(Rtsne) # for t-SNE plot
library(ggplot2) # for visualization
library(dendextend)
library(rgdal)
library(factoextra)


#Define a dissimilarity matrix using Gower distance
gower_dist <- daisy(Sample.Data,
                    metric = "gower",
                    stand=TRUE)
summary(gower_dist)
gower_mat <- as.matrix(gower_dist)
Sample.Data[
  which(gower_mat == min(gower_mat[gower_mat != min(gower_mat)]),
        arr.ind = TRUE)[1, ], ]
Sample.Data[
  which(gower_mat == max(gower_mat[gower_mat != max(gower_mat)]),
        arr.ind = TRUE)[1, ], ]


##### Testing what method is best
hc2 <- agnes(Sample.Data, method = "complete")
# Agglomerative coefficient
hc2$ac
# vector of methods to compare
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")
# Determbine which function to best compute the coefficient
ac <- function(x) {
  agnes(Sample.Data, method = x)$ac
}
map_dbl(m, ac)      


# Buid dendogram
res.hc <- hclust(gower_dist, method = "ward.D2")
plot(res.hc, cex = 0.6, hang = -1, xlab= "Gower distance", ylab = "Clustering height", main = "Dendogram of All Sample points", sub = "", axes = FALSE)


#Cut the dendogram: Identify jumps in the dendogram's inertia
inertie <- sort(res.hc$height, decreasing = TRUE)
plot(inertie[1:20], type = "s", xlab = "Number of classes", ylab = "Inertia", main="All Sample Points")
points(c(3,5), inertie[c(3,5)], col = c("green3", "red3", "blue3", "yellow3"), cex = 2, 
       lwd = 3)


plot(res.hc, labels = FALSE, main = "Partition of All Sample Points in 3 or 5 classes", xlab = "", ylab = "", 
     sub = "", axes = FALSE, hang = -1)
rect.hclust(res.hc, 3, border = "green3")
rect.hclust(res.hc, 5, border = "red3")





# VALIDATION from http://www.sthda.com/english/wiki/print.php?id=241#silhouette-plot-for-hierarchical-clustering
# With 3 clusters

png("cluster3.png", width = 762, height = 462)
res.hc <- eclust(gower_mat, "hclust", k = 3,
                 method = "ward.D2", graph = FALSE) 
fviz_silhouette(res.hc, print.summary = TRUE)
dev.off()

# With 5 clusters

png("cluster5.png", width = 762, height = 462)
res.hc <- eclust(gower_mat, "hclust", k = 5,
                 method = "ward.D2", graph = FALSE) 
fviz_silhouette(res.hc, print.summary = TRUE)
dev.off()



#Tree cut accordingly and showed in colours
A2Rplot(res.hc, k = 3, boxes = FALSE, col.up = "gray50", 
        col.down = brewer.pal(3, "Dark2"), show.labels = FALSE, main = "Clustering of All Samples")


hclusters <- cutree(res.hc, k=3)
Sample.Data <- Sample.Data %>%
  mutate(cluster = hclusters)
#write.csv(conflicts, file="conflictsclustered.csv")



map<-readOGR("communedarrondissement.shp")
#ggplot() +  geom_polygon(data=map, aes(x=long, y=lat, group=group))
# map both polys and points together
#ggplot()+geom_polygon(data=map, aes(x=long, y=lat, group=group)) +  
#  geom_point(data=conflicts, aes(x=longitude, y=latitude), color="red")


map2 <- fortify(map, region = "GCNOM")

#Nicer map
Sample.Data$cluster <-  as.numeric(Sample.Data$cluster)
cols <- c("purple", "darkblue","blue","green", "yellow", "orange", "red")


Sample.Data$cluster <- as.character(Sample.Data$cluster)

#Nicer map of the clusters
ggplot() +
  geom_polygon(data=map2, aes(x=long, y=lat, group=group), fill="#bdbdbd", 
               colour="grey90", alpha=1)+
  labs(x="", y="", title="All Samples")+ #labels
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), # get rid of x ticks/text
        axis.ticks.x = element_blank(),axis.text.x = element_blank(), # get rid of y ticks/text
        plot.title = element_text(lineheight=.8, face="bold", vjust=1))+ # make title bold and add space
  geom_point(aes(x=x, y=y, color=cluster), data=Sample.Data, alpha=1, size=1.5, color="grey20")+# to get outline
  geom_point(aes(x=x, y=y, color=cluster), data=Sample.Data, alpha=1, size=1.5)+
  scale_colour_brewer(palette = "Set1")+
  coord_equal(ratio=1) # square plot to avoid the distortion



# Cluster 1 investigation (playing with the variable displayed)
cluster1 <- filter(Sample.Data, cluster == 1)
ggplot() +  
  geom_polygon(data=map2, aes(x=long, y=lat, group=group), fill="#bdbdbd", 
               colour="grey90", alpha=1)+
  labs(x="", y="", title="All Samples: cluster 1")+ #labels
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), # get rid of x ticks/text
        axis.ticks.x = element_blank(),axis.text.x = element_blank(), # get rid of y ticks/text
        plot.title = element_text(lineheight=.8, face="bold", vjust=1))+ # make title bold and add space
  geom_point(aes(x=x, y=y, color=TLF), data=cluster1, alpha=1, size=1.5, color="grey20")+# to get outline
  geom_point(aes(x=x, y=y, color=TLF), data=cluster1, alpha=1, size=1.5)+
  scale_colour_brewer(palette = "Set1")+  
  scale_colour_gradientn("Source Type", 
                         colours=c('#fff5f0','#fee0d2','#fcbba1','#fc9272','#fb6a4a','#ef3b2c','#cb181d','#a50f15','#67000d'))+ # change color scale
  coord_equal(ratio=1) # square plot to avoid the distortion

# AVec le cluster2
cluster2 <- filter(Sample.Data, cluster == 2)
ggplot() +  
  geom_polygon(data=map2, aes(x=long, y=lat, group=group), fill="#bdbdbd", 
               colour="grey90", alpha=1)+
  labs(x="", y="", title="All Samples: cluster 2")+ #labels
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), # get rid of x ticks/text
        axis.ticks.x = element_blank(),axis.text.x = element_blank(), # get rid of y ticks/text
        plot.title = element_text(lineheight=.8, face="bold", vjust=1))+ # make title bold and add space
  geom_point(aes(x=x, y=y, color=PopDensity), data=cluster2, alpha=1, size=1.5, color="grey20")+# to get outline
  geom_point(aes(x=x, y=y, color=PopDensity), data=cluster2, alpha=1, size=1.5)+
  scale_colour_brewer(palette = "Set1")+  
  scale_colour_gradientn("Source Type", 
                         colours=c('#fff5f0','#fee0d2','#fcbba1','#fc9272','#fb6a4a','#ef3b2c','#cb181d','#a50f15','#67000d'))+ # change color scale
  coord_equal(ratio=1) # square plot to avoid the distortion

# Cluster 3
cluster3 <- filter(Sample.Data, cluster == 3)
ggplot() +  
  geom_polygon(data=map2, aes(x=long, y=lat, group=group), fill="#bdbdbd", 
               colour="grey90", alpha=1)+
  labs(x="", y="", title="All Samples: cluster 3")+ #labels
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), # get rid of x ticks/text
        axis.ticks.x = element_blank(),axis.text.x = element_blank(), # get rid of y ticks/text
        plot.title = element_text(lineheight=.8, face="bold", vjust=1))+ # make title bold and add space
  geom_point(aes(x=x, y=y, color=Type), data=cluster3, alpha=1, size=1.5, color="grey20")+# to get outline
  geom_point(aes(x=x, y=y, color=Type), data=cluster3, alpha=1, size=1.5)+
  scale_colour_brewer(palette = "Set1")+  
  scale_colour_gradientn("Source Type", 
                         colours=c('#fff5f0','#fee0d2','#fcbba1','#fc9272','#fb6a4a','#ef3b2c','#cb181d','#a50f15','#67000d'))+ # change color scale
  coord_equal(ratio=1) # square plot to avoid the distortion