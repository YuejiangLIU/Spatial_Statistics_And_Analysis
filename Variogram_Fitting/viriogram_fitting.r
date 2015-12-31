# Spatial statistics and analysis 2013, Ex3
# Code by M.Schleiss and A.Berne, EPFL-LTE

library(stats)
library(gstat)
library(sp)

# Define the path to the working directory.
# For Windows: use double backslashes ("Z:\\...\\")
# For Linux: use simple slashes ("Z:/.../")
# Make sure your path does not contain any blank spaces.
path <- "D:\\Ma3\\Spatial_statistics_and_analysis\\Part_2\\Ex3\\"

# Graphical parameters (do not change!)
nclass <- 10
boundaries <- seq(250,1500,250)

# Read the dataset and extract the metal concentrations (do not change!)
data(meuse)
coordinates(meuse) <- ~x+y
coordX  <- coordinates(meuse)[,1]	# vector of Y coordinates (in m)
coordY  <- coordinates(meuse)[,2]	# vector of X coordinates (in m)
coordZ  <- meuse[["elev"]]		# vector of relative elevations (in m)
cadmium <- meuse[["cadmium"]]		# vector of cadmium concentrations (in ppm)
copper  <- meuse[["copper"]]		# vector of copper concentrations (in ppm)
lead    <- meuse[["lead"]]		# vector of lead concentrations (in ppm)
zinc    <- meuse[["zinc"]]		# vector of zinc concentrations (in ppm)
om      <- meuse[["om"]]		# vector of organic matter (in percentage)


# Q1: Plot the empirical probability density function of the cadmium and lead concentrations
# cadmium
breaks <- seq(min(cadmium),max(cadmium),length=nclass+1)
main <- "Empirical probability density of cadmium"
file <- paste(path,"pdf_cadmium.png",sep="")
png(file,width=1200,height=1000,res=200)
hist(cadmium,breaks=breaks,main=main,xlab="Concentration [ppm]",ylab="Density [-]",freq=FALSE)
dev.off()

# lead
breaks <- seq(min(lead),max(lead),length=nclass+1)
main <- "Empirical probability density of lead"
file <- paste(path,"pdf_lead.png",sep="")
png(file,width=1200,height=1000,res=200)
hist(lead,breaks=breaks,main=main,xlab="Concentration [ppm]",ylab="Density [-]",freq=FALSE)
dev.off()

# Q2: Apply a square root transform on the values of cadmium and lead.
# Show the pdf of the transformed concentration values.
# Do the same using a log transform.
# cadmium sqrt
cadmium_sqrt <- sqrt(cadmium)
breaks <- seq(min(cadmium_sqrt),max(cadmium_sqrt),length=nclass+1)
main <- "Empirical probability density of the square root of cadmium"
file <- paste(path,"pdf_cadmium_sqrt.png",sep="")
png(file,width=1200,height=1000,res=200)
hist(cadmium_sqrt,breaks=breaks,main=main,xlab="Concentration [ppm]",ylab="Density [-]",freq=FALSE)
dev.off()

# lead sqrt
lead_sqrt <- sqrt(lead)
breaks <- seq(min(lead_sqrt),max(lead_sqrt),length=nclass+1)
main <- "Empirical probability density of the squrare root of lead"
file <- paste(path,"pdf_lead_sqrt.png",sep="")
png(file,width=1200,height=1000,res=200)
hist(lead_sqrt,breaks=breaks,main=main,xlab="Concentration [ppm]",ylab="Density [-]",freq=FALSE)
dev.off()

# cadmium log
cadmium_log <- log(cadmium)
breaks <- seq(min(cadmium_log),max(cadmium_log),length=nclass+1)
main <- "Empirical probability density of the log of cadmium"
file <- paste(path,"pdf_cadmium_log.png",sep="")
png(file,width=1200,height=1000,res=200)
hist(cadmium_log,breaks=breaks,main=main,xlab="Concentration [ppm]",ylab="Density [-]",freq=FALSE)
dev.off()

# lead sqrt
lead_log <- log(lead)
breaks <- seq(min(lead_log),max(lead_log),length=nclass+1)
main <- "Empirical probability density of the log of lead"
file <- paste(path,"pdf_lead_log.png",sep="")
png(file,width=1200,height=1000,res=200)
hist(lead_log,breaks=breaks,main=main,xlab="Concentration [ppm]",ylab="Density [-]",freq=FALSE)
dev.off()


# Q3: Apply the Shapiro-Wilk test:
# cadmium
shapiro.test(cadmium)
shapiro.test(cadmium_sqrt)
shapiro.test(cadmium_log)
# lead
shapiro.test(lead)
shapiro.test(lead_sqrt)
shapiro.test(lead_log)


# Q4: Compute the isotropic sample variogram of log(cadmium) using the provided boundaries.
# fit a spherical variogram model on it.
boundaries <- seq(250,1500,250)
vario.cadmium_log <- variogram(cadmium_log~1,meuse,boundaries=boundaries)
# modelling and fitting
cad_mod <- vgm(psill=1.9,model="Sph",range=1100,nugget=0.5)
cad_fit <- fit.variogram(object=vario.cadmium_log,model=cad_mod)

main <- "Isotropic sample variogram of the log of cadmium"
xlab <- "Distance lag [m]"
ylab <- "Semivariance [ppm^2]"
xlim <- c(0,1500)
ylim <- c(0,max(vario.cadmium_log$gamma)+0.2)
file <- paste(path,"250m_isotropic_sample_variogram_cadmium_log.png",sep="")
png(file,width=1200,height=1000,res=200)
par(mai=c(0.95,0.95,0.5,0.3))
plot(vario.cadmium_log,cad_fit,main=main,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,lwd=1.5,cex=0.8,cex.lab=1.25,cex.axis=1.1)

dev.off()

# Q5: Compute the isotropic sample variogram of log(lead) using the provided boundaries.
# Fit both a spherical and an exponential variogram model on it.
boundaries <- seq(250,1500,250)
vario.lead_log <- variogram(lead_log~1,meuse,boundaries=boundaries)
# modelling and fitting
lead_sph <- vgm(psill=1.9,model="Sph",range=1100,nugget=0.5)
lead_fit_sph <- fit.variogram(object=vario.lead_log,model=lead_sph)
lead_exp <- vgm(psill=1.9,model="Exp",range=1100,nugget=0.5)
lead_fit_exp <- fit.variogram(object=vario.lead_log,model=lead_exp)

main <- "Isotropic sample variogram of the log of lead (sph)"
xlab <- "Distance lag [m]"
ylab <- "Semivariance [ppm^2]"
xlim <- c(0,1500)
ylim <- c(0,max(vario.lead_log$gamma)+0.1)
file <- paste(path,"250m_isotropic_sample_variogram_lead_sph.png",sep="")
png(file,width=1200,height=1000,res=200)
par(mai=c(0.95,0.95,0.5,0.3))
plot(vario.lead_log,lead_fit_sph,col=c("blue","black"),main=main,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,lwd=1.5,cex=0.8,cex.lab=1.25,cex.axis=1.1)
dev.off()

main <- "Isotropic sample variogram of the log of lead (exp)"
xlab <- "Distance lag [m]"
ylab <- "Semivariance [ppm^2]"
xlim <- c(0,1500)
ylim <- c(0,max(vario.lead_log$gamma)+0.1)
file <- paste(path,"250m_isotropic_sample_variogram_lead_exp.png",sep="")
png(file,width=1200,height=1000,res=200)
par(mai=c(0.95,0.95,0.5,0.3))
plot(vario.lead_log,lead_fit_exp,col=c("blue","black"),main=main,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,lwd=1.5,cex=0.8,cex.lab=1.25,cex.axis=1.1)
dev.off()


