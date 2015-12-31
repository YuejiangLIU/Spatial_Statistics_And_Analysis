# Code by M.Schleiss and A.Berne, EPFL-LTE
  
library(gstat)
library(sp)
library(hydroGOF)

# Define the path to the working directory.
# For Windows: use double backslashes ("Z:\\...\\")
# For Linux: use simple slashes ("Z:/.../")
# Make sure your path does not contain any blank spaces.
path <- "D:\\Ma3\\Spatial_statistics_and_analysis\\Part_2\\Ex4\\"

# Graphical parameters (do not change!)
jet.colors <- c("#00007F","blue","#007FFF","cyan","#7FFF7F","yellow","#FF7F00","red","#7F0000")
col.palette <- colorRampPalette(jet.colors)
cuts  <- seq(400,2600,200)
ncuts <- length(cuts)

# Read the prediction dataset (do not change!)
file <- paste(path,"precipitation_2011_prediction.txt",sep="")
data <- read.table(file,header=TRUE)
coordX <- data[["x"]]  		# vector with x coordinates [km]
coordY <- data[["y"]]			# vector with y coordinates [km]
coordZ <- data[["z"]]			# vector with z coordinates [km]
precip <- data[["precip"]]		# vector with total annual precipitation [mm]
coordinates(data) <- ~x+y+z

# Map with total annual precipitation (do not change!)
main <- "Total annual precipitation (prediction dataset) [mm]"
xlab <- "East [km]"
ylab <- "North [km]"
file <- paste(path,"map_precip_prediction.png",sep="")
png(file,width=1288,height=800,res=200)
print(spplot(data["precip"],scales=list(draw=TRUE),xlab=xlab,ylab=ylab,main=main,key.space="right",cuts=cuts,region=TRUE,col.regions=col.palette(ncuts)))
dev.off()


# Q1: Compute the isotropic sample variogram of the precipitation values and fit a spherical model on it
vario.iso <- variogram(precip~1,data)
iso_mod <- vgm(psill=120000,model="Sph",range=80,nugget=10000)
iso_fit <- fit.variogram(object=vario.iso,model=iso_mod)

main <- "Isotropic sample variogram of the precipitation"
xlab <- "Distance lag [km]"
ylab <- expression(paste("Semivariance [",mm^2,"]"))
xlim <- c(0,150)
ylim <- c(0,max(vario.iso$gamma)*1.1)
file <- paste(path,"isotropic_sample_variogram_precip.png",sep="")
png(file,width=1200,height=1000,res=200)
par(mai=c(0.95,0.95,0.5,0.3))
plot(vario.iso,iso_fit,main=main,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,lwd=1.5,cex=0.8,cex.lab=1.25,cex.axis=1.1)
dev.off()

# Q2: Compute the variogram map (2D variogram) of the precipitation values
vario.map <- variogram(precip~1,data,cutoff=90,width=10,map=TRUE)

# Uncomment the following lines in order to plot the variogram map.
vario.map <- as.data.frame(vario.map)[,c("map.dx","map.dy","map.var1")]

names(vario.map) <- c("dx","dy","gamma")
coordinates(vario.map) <- ~dx+dy
gridded(vario.map) <- TRUE
main <- expression(paste("Semivariance [",mm^2,"]"))
xlab <- "dx [km]"
ylab <- "dy [km]"
file <- paste(path,"vario_map.png",sep="")
png(file,width=1200,height=1200,res=200)
print(spplot(vario.map["gamma"],scales=list(draw=TRUE),xlab=xlab,ylab=ylab,main=main,col.regions=col.palette(200)))
dev.off()


# Q3: Compute the variograms along the direction of min/max variability
vario.min <- variogram(precip~1,data,alpha=50,tol.hor=15)
min_mod <- vgm(psill=100000,model="Sph",range=100,nugget=10000)
min_fit <- fit.variogram(object=vario.min,model=min_mod)
main <- "50 degree sample variogram of the precipitation"
xlab <- "Distance lag [km]"
ylab <- expression(paste("Semivariance [",mm^2,"]"))
xlim <- c(0,150)
ylim <- c(0,max(vario.min$gamma)*1.1)
file <- paste(path,"50_sample_variogram_precip.png",sep="")
png(file,width=1200,height=1000,res=200)
par(mai=c(0.95,0.95,0.5,0.3))
plot(vario.min,min_fit,main=main,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,lwd=1.5,cex=0.8,cex.lab=1.25,cex.axis=1.1)
dev.off()

vario.max <- variogram(precip~1,data,alpha=140,tol.hor=15)
max_mod <- vgm(psill=160000,model="Sph",range=70,nugget=10000)
max_fit <- fit.variogram(object=vario.max,model=max_mod)
main <- "140 degree sample variogram of the precipitation"
xlab <- "Distance lag [km]"
ylab <- expression(paste("Semivariance [",mm^2,"]"))
xlim <- c(0,150)
ylim <- c(0,max(vario.max$gamma)*1.1)
file <- paste(path,"140_sample_variogram_precip.png",sep="")
png(file,width=1200,height=1000,res=200)
par(mai=c(0.95,0.95,0.5,0.3))
plot(vario.max,max_fit,main=main,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,lwd=1.5,cex=0.8,cex.lab=1.25,cex.axis=1.1)
dev.off()


# Q4: Fit a spherical variogram model on both directional variograms
min_fit
max_fit

# Q5: Define a variogram with geometric anisotropy
vgm.anis <- vgm(psill=min_fit$psill[2],model="Sph",range=min_fit$range[2],nugget=min_fit$psill[1],anis=c(50,max_fit$range[2]/min_fit$range[2]))
anis_fit <- fit.variogram(object=vario.iso,model=vgm.anis)
main <- "aniso sample variogram of the precipitation"
xlab <- "Distance lag [km]"
ylab <- expression(paste("Semivariance [",mm^2,"]"))
xlim <- c(0,150)
ylim <- c(0,max(vario.iso$gamma)*1.1)
file <- paste(path,"aniso_sample_variogram_precip.png",sep="")
png(file,width=1200,height=1000,res=200)
par(mai=c(0.95,0.95,0.5,0.3))
plot(vario.iso,anis_fit,main=main,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,lwd=1.5,cex=0.8,cex.lab=1.25,cex.axis=1.1)
dev.off()

# Read the validation dataset (do not change!)
file <- paste(path,"precipitation_2011_validation.txt",sep="")
newdata <- read.table(file,header=TRUE)
validationX <- newdata[["x"]]			# vector with x coordinates of validation points [km]
validationY <- newdata[["y"]]			# vector with y coordinates of validation points [km]
validationZ <- newdata[["z"]]			# vector with z coordinates of validation points [km]
validation.precip <- newdata[["precip"]]	# vector with total annual precipitation for validation [mm]
coordinates(newdata) <- ~x+y+z

# Map with total annual precipitation for the validation dataset (do not change!)
main <- "Total annual precipitation (validation dataset) [mm]"
xlab <- "East [km]"
ylab <- "North [km]"
file <- paste(path,"map_precip_validation.png",sep="")
png(file,width=1288,height=800,res=200)
print(spplot(newdata["precip"],scales=list(draw=TRUE),xlab=xlab,ylab=ylab,main=main,key.space="right",cuts=cuts,region=TRUE,col.regions=col.palette(ncuts)))
dev.off()


# Q6: Interpolate the precipitation values at the locations of the validation dataset
# (1) using inverse distance weighting (IDW)
# (2) using simple kriging and an isotropic variogram (SK.iso)
# (3) using ordinary kriging and an isotropic variogram (OK.iso)
# (4) using ordinary kriging and an anisotropic variogram (OK.anis)

IDW <- krige(formula=precip~1,locations=data,newdata=newdata)
SK.iso <- krige(formula=precip~1,locations=data,newdata=newdata,model=iso_mod,beta=mean(precip,na.rm=T))
OK.iso <- krige(formula=precip~1,locations=data,newdata=newdata,model=iso_mod)
OK.anis <- krige(formula=precip~1,locations=data,newdata=newdata,model=vgm.anis)

predictions.IDW <- IDW[[1]]
predictions.SK.iso  <- SK.iso[[1]]
predictions.OK.iso  <- OK.iso[[1]]
predictions.OK.anis <- OK.anis[[1]]


# Q7: Compare the predictions with the measurements. Compute the bias and the rmse.
mean(predictions.IDW-validation.precip)
rmse(predictions.IDW,validation.precip)
mean(predictions.SK.iso-validation.precip)
rmse(predictions.SK.iso,validation.precip)
mean(predictions.OK.iso-validation.precip)
rmse(predictions.OK.iso,validation.precip)
mean(predictions.OK.anis-validation.precip)
rmse(predictions.OK.anis,validation.precip)

# Read the 1-km digital elevation model for Switzerland (do not change!)
file <- paste(path,"DEM_Switzerland_1km.txt",sep="")
DEM  <- read.table(file,header=TRUE)
DEM.coordX <- DEM[["x"]]	# vector with x coordinates [km]
DEM.coordY <- DEM[["y"]]	# vector with y coordinates [km]
DEM.coordZ <- DEM[["z"]]	# vector with z coordinates [km]
coordinates(DEM) <- ~x+y+z

# Q8: Interpolate the precipitation values at the locations given by the DEM
# (1) using inverse distance weighting (IDW)
# (2) using ordinary kriging and an isotropic variogram (OK.iso)
# (3) using ordinary kriging and an anisotropic variogram (OK.anis)

DEM.IDW <- krige(formula=precip~1,locations=data,newdata=DEM)
DEM.OK.iso  <- krige(formula=precip~1,locations=data,newdata=DEM,model=iso_mod)
DEM.OK.anis <- krige(formula=precip~1,locations=data,newdata=DEM,model=vgm.anis)

DEM.predictions.IDW <- DEM.IDW[[1]]
DEM.predictions.OK.iso  <- DEM.OK.iso[[1]]
DEM.predictions.OK.anis <- DEM.OK.anis[[1]]


# Uncomment the following lines in order to plot the interpolated precipitation values
map <- data.frame(DEM.coordX,DEM.coordY,DEM.predictions.IDW,DEM.predictions.OK.iso,DEM.predictions.OK.anis)
names(map) <- c("x","y","DEM.IDW","DEM.OK.iso","DEM.OK.anis")
coordinates(map) <- ~x+y
gridded(map) <- TRUE
pts = list("sp.points",coordinates(data),pch=3,lwd=1.5,cex=0.15,col="black") 
main <- "Interpolated total annual precipitation [mm]"
xlab <- "East [km]"
ylab <- "North [km]"
file <- paste(path,"map_interpolated_precip.png",sep="")
png(file,width=2560,height=1600,res=200)
print(spplot(map[,c("DEM.OK.iso","DEM.OK.anis","DEM.IDW")],xlab=xlab,ylab=ylab,main=main,col.regions=col.palette(200),sp.layout=list(pts),scales=list(draw=TRUE)))
dev.off()