# Spatial statistics and analysis 2013, Ex5
# Code by M.Schleiss and A.Berne, EPFL-LTE

library(gstat)
library(sp)
library(hydroGOF)

# Define the path to the working directory.
# For Windows: use double backslashes ("Z:\\...\\")
# For Linux: use simple slashes ("Z:/.../")
# Make sure your path does not contain any blank spaces.
path <- "D:\\Ma3\\Spatial_statistics_and_analysis\\Part_2\\Ex5\\"

# Graphical parameters (do not change!)
jet.colors  <- c("#00007F","blue","#007FFF","cyan","#7FFF7F","yellow","#FF7F00","red","#7F0000")
col.palette <- colorRampPalette(jet.colors)
cuts <- seq(-8,16,2)
ncuts <- length(cuts)
boundaries <- c(seq(20,60,20),seq(90,270,30))
nclass <- 6

# Read the "temperature 2011" prediction dataset (do not change!)
file <- paste(path,"temperatures_2011_prediction.txt",sep="")
data <- read.table(file,header=TRUE)
coordX <- data[["x"]]		# vector with x coordinates [km]
coordY <- data[["y"]]		# vector with y coordinates [km]
coordZ <- data[["z"]]		# vector with z coordinates [km]
temp   <- data[["temp"]]	# vector with average annual temperatures for 2011 [째C]
coordinates(data) <- ~x+y+z
pts <- list("sp.points",coordinates(data),pch=3,lwd=2,col="black") 

# Plot the map with the average annual temperature for 2011  (do not change!)
main <- "Average annual temperature [째C]"
xlab <- "East [km]"
ylab <- "North [km]"
file <- paste(path,"map_temperatures.png",sep="")
png(file,width=1288,height=800,res=200)
print(spplot(data["temp"],xlab=xlab,ylab=ylab,key.space="right",cuts=cuts,region=TRUE,col.regions=col.palette(ncuts),main=main,scales=list(draw=TRUE)))
dev.off()

# Q1: Plot the empirical pdf of the average annual temperature for 2011
breaks <- seq(min(temp),max(temp),length=nclass+1)
main <- "Empirical probability density of temperature"
xlab <- "Average annual temperature [째C]"
ylab <- "Density [-]"
file <- paste(path,"pdf_temperature.png",sep="")
png(file,width=1200,height=800,res=200)
hist(temp,breaks=breaks,main=main,xlab=xlab,ylab=ylab,freq=FALSE)
dev.off()

# Q2: Plot the temperature w.r.t. the altitude
main <- "Temperature vs altitude"
xlab <- "Altitude [km]"
ylab <- "Average annual temperature [째C]"
file <- paste(path,"temperature_vs_altitude.png",sep="")
png(file,width=1200,height=800,res=200)
plot(coordZ,temp,main=main,xlab=xlab,ylab=ylab)
#dev.off()

# Q3: Fit a linear model for the temperature vs the altitude
reg <- lm(temp~coordZ)
abline(reg)
dev.off()


# Q4: Compute the isotropic sample variogram of the temperature (without trend)
# Fit a spherical model on it.
vario.temp <- variogram(temp~1,data,boundaries=boundaries,cressie=TRUE)
main <- "Isotropic sample variogram of temperature (Cressie)"
xlab <- "Distance lag [m]"
ylab <- "Semivariance [캜^2]"
xlim <- c(0,300)
ylim <- c(0,max(vario.temp$gamma)*1.1)
file <- paste(path,"Cressie_variogram_temp.png",sep="")
png(file,width=1200,height=1000,res=200)
par(mai=c(0.95,0.95,0.5,0.3))
temp_vgm <- vgm(psill=15,model="Sph",range=200,nugget=0.5)
temp_fit <- fit.variogram(object=vario.temp,model=temp_vgm)
plot(vario.temp,temp_fit,main=main,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,lwd=1.5,cex=0.8,cex.lab=1.25,cex.axis=1.1)
dev.off()



# Q5: Compute the isotropic sample variogram of the temperature with a linear trend w.r.t. the altitude
# Fit a spherical model on it.
vario.temp_alt<- variogram(temp~coordZ,data,boundaries=boundaries,cressie=TRUE)
main <- "Variogram of temperature vs altitude (Cressie)"
xlab <- "Distance lag [m]"
ylab <- "Semivariance [캜^2]"
xlim <- c(0,300)
ylim <- c(0,max(vario.temp_alt$gamma)*1.1)
file <- paste(path,"Cressie_variogram_temp_altitude.png",sep="")
png(file,width=1200,height=1000,res=200)
par(mai=c(0.95,0.95,0.5,0.3))
temp_alt_vgm <- vgm(psill=0.9,model="Sph",range=200,nugget=0.5)
temp_alt_fit <- fit.variogram(object=vario.temp_alt,model=temp_alt_vgm)
plot(vario.temp_alt,temp_alt_fit,main=main,xlab=xlab,ylab=ylab)
dev.off()

# Read the validation dataset (do not change!)
file <- paste(path,"temperatures_2011_validation.txt",sep="")
newdata <- read.table(file,header=TRUE)
validationX <- newdata[["x"]]			# vector with x coordinates of validation points [km]
validationY <- newdata[["y"]]			# vector with y coordinates of validation points [km]
validationZ <- newdata[["z"]]			# vector with z coordinates of validation points [km]
validation.temp <- newdata[["temp"]]		# vector with average annual temperature for 2011 [째C]
coordinates(newdata) <- ~x+y+z

# Q6: Interpolate the temperature values at the locations defined in the validation set
# (1) using ordinary kriging (OK)
# (2) using a linear model (LM) for the temperature vs the altitude
# (3) using universal kriging (UK) and a linear trend for the temperature vs the altitude

OK.iso <- krige(formula=temp~1,locations=data,newdata=newdata,model=temp_vgm)
UK.iso <- krige(formula=temp~z,locations=data,newdata=newdata,model=temp_alt_vgm)

predictions.OK.iso  <- OK.iso[[1]]
predictions.LM.iso <- reg$coefficients[1] + reg$coefficients[2]*newdata$z
predictions.UK.iso  <- UK.iso[[1]]

mean(predictions.OK.iso-validation.temp)
rmse(predictions.OK.iso,validation.temp)
mean(predictions.LM.iso-validation.temp)
rmse(predictions.LM.iso,validation.temp)
mean(predictions.UK.iso-validation.temp)
rmse(predictions.UK.iso,validation.temp)

# Compare the interpolated values to the validation data. 
# Compute the bias and the rmse.
# ???

# Read the 1-km DEM for Switzerland (do not change!)
file <- paste(path,"DEM_Switzerland_1km.txt",sep="")
DEM  <- read.table(file,header=TRUE)
DEM.coordX <- DEM[["x"]]	# vector with x coordinates [km]
DEM.coordY <- DEM[["y"]]	# vector with y coordinates [km]
DEM.coordZ <- DEM[["z"]]	# vector with z coordinates [km]
coordinates(DEM) <- ~x+y+z

# Q7: Interpolate the temperature values at the locations of the DEM:
# (1) using ordinary kriging (OK)
# (2) using universal kriging (UK) and a linear trend w.r.t. the altitude
DEM.OK.iso  <- krige(formula=temp~1,locations=data,newdata=DEM,model=temp_vgm)
DEM.UK.iso <- krige(formula=temp~z,locations=data,newdata=DEM,model=temp_alt_vgm)
predictions.DEM.OK.iso  <- DEM.OK.iso[[1]]
predictions.DEM.UK.iso  <- DEM.UK.iso[[1]]

# Plot the map with the predicted temperature values
# Hint: see the code for exercise series 4 
map <- data.frame(DEM.coordX,DEM.coordY,predictions.DEM.OK.iso,predictions.DEM.UK.iso)
names(map) <- c("x","y","DEM.OK.iso","DEM.UK.iso")
coordinates(map) <- ~x+y
gridded(map) <- TRUE
pts = list("sp.points",coordinates(data),pch=2,lwd=1.5,cex=0.15,col="black") 
main <- "Interpolated average annual temperature [째C]"
xlab <- "East [km]"
ylab <- "North [km]"
file <- paste(path,"map_interpolated_temp.png",sep="")
png(file,width=2560,height=1000,res=200)
print(spplot(map[,c("DEM.OK.iso","DEM.UK.iso")],xlab=xlab,ylab=ylab,main=main,col.regions=col.palette(200),sp.layout=list(pts),scales=list(draw=TRUE)))
dev.off()
