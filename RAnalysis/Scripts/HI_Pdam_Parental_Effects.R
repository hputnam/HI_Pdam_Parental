#Examination of Trans-generational Acclimatization in Pocillopora damicornis exposed to Ocean Acifiication
#Data published in X
#Title: 
#Contact: Hollie Putnam hollieputnam@gmail.com
#Supported by: NSF Ocean Sciencs Postdoctoral Research Fellowship (NSF OCE PRF-1323822) and NSF EPSCOR (NSF EPS-0903833)
#last modified 20171113
#See Readme file for details on data files and metadata

rm(list=ls()) # removes all prior objects

#R Version: R version 3.3.1
#RStudio Version: 1.0.44
######Read in required libraries#####
library(car) #version 2.1-4 Date: 2016-11-30 Depends: R (>= 3.2.0) Imports:MASS, mgcv, nnet, pbkrtest (>= 0.3-2), quantreg, grDevices, utils, stats, graphics, Rcpp
library(ggplot2) #version 2.2.1 Date/Publication: 2016-12-30 Depends: R (>= 3.1) Imports: digest, grid, gtable (>= 0.1.1), MASS, plyr (>= 1.7.1),reshape2, scales (>= 0.4.1), stats, tibble, lazyeval
library(gridExtra) #version: 2.2.1 Date/Publication: 2016-02-29 Depends: R(>= 2.5.0) Imports: gtable, grid, grDevices, graphics, utils
library(lsmeans)  #version: 2.26-3 Date: 2017-05-09 Depends: estimability, methods, R (>= 3.0) Imports: graphics, stats, utils, nlme, coda (>= 0.17), multcomp, plyr,mvtnorm, xtable (>= 1.8-2)
library(multcomp) #version: 1.4-6 Date: 2016-07-14 Depends: stats, graphics, mvtnorm (>= 1.0-3), survival (>= 2.39-4), TH.data (>= 1.0-2)
library(nlme) #version: 3.1-131 Date: 2017-02-06 Depends: R (>= 3.0.2) Imports: graphics, stats, utils, lattice
library(plotrix) #version: 3.6-5 Date: 2017-05-09 Depends: NA Imports: grDevices, graphics, stats, utils
library(plyr) #version: 1.8.4 Date/Publication: 2016-06-08 Depends: R (>= 3.1.0) Imports: Rcpp (>= 0.11.0)
library(reshape) #version: 3.3.1 Date/Publication: 2016-06-24  Depends: R (>= 3.3.1)
library(seacarb) #version: 3.2 Date/Publication: 2017-06-19 Depends: R (>= 2.10), oce, gsw Imports: NA
library(grid) #version: 3.3.1 Date/Publication: 2016-06-24  Depends: R (>= 3.3.1)
library(xtable) #version 1.8-2 Date/Publication: 2016-01-08 Depends: R (>= 2.10.0)
library(lme4) #version: 1.1-13 Date/Publication: 2017-04-19 Depends: R (>= 3.0.2), Matrix (>= 1.1.1), methods, stats Imports: graphics, grid, splines, utils, parallel, MASS, lattice, nlme(>= 3.1-123), minqa (>= 1.1.15), nloptr (>= 1.0.4)
library(blmeco) #version: 1.1 Date/Publication: 2015-08-22 Depends: R (>= 3.0.0), stats, MASS Imports: MuMIn, arm, lme4
library(MuMIn) #version: 1.15.6 Date/Publication: 2016-01-07 Depends: R (>= 3.0.0) Imports: graphics, methods, Matrix, stats, stats4

#####Required Data files#####
#Light_Calibration_Data.csv
#Temperature_Calibration_Data.csv
#Acclimation_Data.csv
#Adult_Tank_Light.csv
#Adult_Tank_Temp.csv
#Adult_Tank_NBS_pH.csv
#august.larval.release.data.csv
#CRM_TA_Data.csv
#Daily_Temp_pH_Sal.csv
#Field_Temp.csv
#july.larval.release.data.csv
#june.larval.release.data.csv
#Larval_Data_M0.csv
#Larval_Data_M1.csv
#Month1_Larval_Size.csv
#Month1_Tank_Light.csv
#Month1_Tank_NBS_pH.csv
#Month1_Tank_Temp.csv
#Month6_Larval_Size.csv
#Month6_Tank_Light.csv
#Month6_Tank_NBS_pH.csv
#Month6_Tank_Temp.csv
#~/MyProjects/HI_Pdam_Parental/RAnalysis/Data/pH_Calibration_Files/
#~/MyProjects/HI_Pdam_Parental/RAnalysis/Data/TA/
#TA_mass_data.csv
#RM_Release_Data.csv
#20140701_20140731_Tides_Kaneohe.txt

#############################################################
setwd("~/MyProjects/HI_Pdam_Parental/RAnalysis/Data") #set working directory
mainDir<-'~/MyProjects/HI_Pdam_Parental/RAnalysis/' #set main directory
#############################################################

##### CALIBRATIONS #####
#Light Logger Calibrations
#Tank2 Odyssey PAR Logger Serial Number = 2484
#Tank3 Odyssey PAR Logger Serial Number = 2485
#Tank4 Odyssey PAR Logger Serial Number = 2486
#Tank5 Odyssey PAR Logger Serial Number = 2487
Cal.L.data <- read.csv("Light_Calibration_Data.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
#Convert light from counts to instantaneous PAR
Cal.L.data$Licor.quanta <- (Cal.L.data$Licor.mol.m2*10^6)/(15*60) #convert 15 min integrated data to instantaneous µmol m-2 s-1
Cal.L.data$Tank2.quanta <- Cal.L.data$Tank2/(15*60) #convert 15 min integrated data to instantaneous µmol m-2 s-1
Cal.L.data$Tank3.quanta <- Cal.L.data$Tank3/(15*60) #convert 15 min integrated data to instantaneous µmol m-2 s-1
Cal.L.data$Tank4.quanta <- Cal.L.data$Tank4/(15*60) #convert 15 min integrated data to instantaneous µmol m-2 s-1
Cal.L.data$Tank5.quanta <- Cal.L.data$Tank5/(15*60) #convert 15 min integrated data to instantaneous µmol m-2 s-1
#Run linear models of loggers against standard and extract model coeffecients to be applied to raw data
L2.lm <- coef(lm(Licor.quanta ~ Tank2.quanta, data=Cal.L.data)) #extract model coefficients
L3.lm <- coef(lm(Licor.quanta ~ Tank3.quanta, data=Cal.L.data)) #extract model coefficients
L4.lm <- coef(lm(Licor.quanta ~ Tank4.quanta, data=Cal.L.data)) #extract model coefficients
L5.lm <- coef(lm(Licor.quanta ~ Tank5.quanta, data=Cal.L.data)) #extract model coefficients

#Temperature Logger Calibrations
#Tank1 Hobo Water Temp Pro Logger Serial Number = 10487931
#Tank2 Hobo Water Temp Pro Logger Serial Number = 10487932
#Tank3 Hobo Water Temp Pro Logger Serial Number = 10487933
#Tank4 Hobo Water Temp Pro Logger Serial Number = 10487934
#Tank5 Hobo Water Temp Pro Logger Serial Number = 10487935
Cal.T.data <- read.csv("Temperature_Calibration_Data.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
#Run linear models of loggers against standard and extract model coeffecients to be applied to raw data
T2.lm <- coef(lm(Tank1 ~ Tank2, data=Cal.T.data)) #extract model coefficients
T3.lm <- coef(lm(Tank1 ~ Tank3, data=Cal.T.data)) #extract model coefficients
T4.lm <- coef(lm(Tank1 ~ Tank4, data=Cal.T.data)) #extract model coefficients
T5.lm <- coef(lm(Tank1 ~ Tank5, data=Cal.T.data)) #extract model coefficients

##### FIELD DATA #####
##Data available in Putnam et al 2016 Evolutionary Applications and plotted there. Also provided here for viewing.
##Field Temperature Data (28March14 - 01April14)
##load field collection/acclimation temp data
# Field.data <- read.csv("Field_Temp.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
# mydate2 <- strptime(Field.data$Date.Time, format="%m/%d/%y %H:%M") #convert date format to characters
# quarterhours2 <- format(as.POSIXct(mydate2) ,format = "%H:%M") #set time as every 15 min
# Field.data <- cbind(Field.data, quarterhours2) #make a dataframe out of data and new times
# Field.data #View data
# min(Field.data$Temperature) #view minimum
# max(Field.data$Temperature) #view maximum
# quarterly.temp.mean2 <- aggregate(Temperature ~ quarterhours2, data=Field.data, mean, na.rm=TRUE) #calculate mean of temperature for every 15 min interval
# quarterly.temp.se2 <- aggregate(Temperature ~ quarterhours2, data=Field.data, std.error, na.rm=TRUE)  #calculate standard error of the mean of temperature for every 15 min interval
# Field.temp.N <- sum(!is.na(Field.data$Temperature)) #Count sample size without NA values
# field.temp.means <- data.frame(quarterly.temp.mean2, quarterly.temp.se2$Temperature) #combine mean and standard error results
# colnames(field.temp.means) <- c("Time", "mean", "se")  #rename columns to describe contents
# 
# #Plot average diurnal cycle of temp data
# FigFieldTemp <- ggplot(field.temp.means) + #Plot average diurnal cycle of temperature data
#   geom_point(aes(x = Time, y = mean), colour="black") + #Plot points using time as the x axis, light as the Y axis and black dots
#   geom_errorbar(aes(x=Time, ymax=mean+se, ymin=mean-se), position=position_dodge(0.9), data=field.temp.means) + #set values for standard error bars and offset on the X axis for clarity
#   scale_x_discrete(breaks=c("0:00", "01:00", "02:00", "03:00", "04:00", "05:00", "06:00", "07:00", "08:00", "09:00", "10:00", "11:00", "12:00","13:00", "14:00", "15:00", "16:00", "17:00", "18:00","19:00", "20:00", "21:00", "22:00", "23:00")) + #set discrete breaks on the X axis
#   ggtitle("F Field Acclimation") + #Label the graph with the main title
#   ylim(23.5,29) + #Set Y axis limits
#   xlab("Time") + #Label the X Axis
#   ylab("Temperature (°C)") + #Label the Y Axis
#   theme_bw() + #Set the background color
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
#         axis.line = element_line(color = 'black'), #Set the axes color
#         panel.border = element_blank(), #Set the border
#         panel.grid.major = element_blank(), #Set the major gridlines
#         panel.grid.minor = element_blank(), #Set the minor gridlines
#         plot.background=element_blank(), #Set the plot background
#         plot.title=element_text(hjust=0)) #Justify the title to the top left
# FigFieldTemp #View figure

##### TANK ACCLIMATION DATA #####
Acclim.data <- read.csv("Acclimation_Data.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
Acclim.data[Acclim.data == 0] <- NA #replace 0 with NA
light <-(Acclim.data$Tank2.Light)/(15*60) #Assign light column in dataframe and convert to units of µmol m-2 s-1
light <-(light*L2.lm[2])+L2.lm[1] #Apply the cross calibration of the odyssey light to Licor cosine sensor standard 192SA cosine sensor
temp <-Acclim.data$Tank2.Temp #Assign temperature column in dataframe
temp <-(temp*T2.lm[2])+T2.lm[1] #Apply the cross calibration of temperature to standard logger #1
Acc <- data.frame(Acclim.data$Date.Time, light, temp) #combine light and temperature data
mydate <- strptime(Acc$Acclim.data.Date.Time, format="%m/%d/%y %H:%M") #convert date format to characters
quarterhours <- format(as.POSIXct(mydate) ,format = "%H:%M") #set time as every 15 min
Acc.data <- cbind(Acc, quarterhours) #make a dataframe out of data and new times
Acc.data #View data
range(na.omit(Acc.data$temp)) #view data range
plot(Acc.data$Acclim.data.Date.Time, Acc.data$temp) #plot temperature by time

#Tank Temperature Data for Acclimation Period (02April14 - 05May14)
quarterly.temp.mean <- aggregate(temp ~ quarterhours, data=Acc.data, mean, na.rm=TRUE) #calculate mean of temperature for every 15 min interval
quarterly.temp.se <- aggregate(temp ~ quarterhours, data=Acc.data, std.error, na.rm=TRUE) #calculate standard error of the mean of temperature for every 15 min interval
Acc.temp.N <- sum(!is.na(Acc.data$temp)) #Count sample size without NA values
Acc.temp.N #View data
temp.means <- data.frame(quarterly.temp.mean, quarterly.temp.se$temp) #combine mean and standard error results
colnames(temp.means) <- c("Time", "mean", "se")  #rename columns to describe contents
range(na.omit(temp.means$mean)) #range of average temp levels

Fig2 <- ggplot(temp.means) + #Plot average diurnal cycle of temperature data
  geom_point(aes(x = Time, y = mean), colour="black", cex=0.8) + #Plot points using time as the x axis, light as the Y axis and black dots
  geom_errorbar(aes(x=Time, ymax=mean+se, ymin=mean-se), position=position_dodge(0.9), data=temp.means, col="black", width=0) + #set values for standard error bars and offset on the X axis for clarity
  scale_x_discrete(breaks=c("01:00", "06:00", "12:00", "18:00", "23:00")) + #set discrete breaks on the X axis
  ggtitle("ACCLIMATION\n(a)") + #Label the graph with the main title
  ylim(23.5,29) + #Set Y axis limits
  xlab("Time") + #Label the X Axis
  ylab("Temperature (°C)") + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), angle = 90, vjust = 0.5, hjust=1), #set margins on labels
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        plot.title=element_text(hjust=0)) #Justify the title to the top left
Fig2 #View figure

#Tank Light Data for Acclimation Period (02April14 - 05May14)
Acc.light.N <- sum(!is.na(Acc.data$light)) #Count sample size
Acc.light.N #View data
quarterly.light.mean <- aggregate(light ~ quarterhours, data=Acc.data, mean, na.rm=TRUE) #calculate mean of light for every 15 min interval
quarterly.light.mean <- rbind(c("05:45","NA"), quarterly.light.mean) #combine data 
quarterly.light.mean[56,c(1,2)] <- c("19:30", "NA") #add row of time to complete data frame
quarterly.light.mean[57,c(1,2)] <- c("19:45", "NA") #add row of time to complete data frame
quarterly.light.mean[58,c(1,2)] <- c("20:00", "NA") #add row of time to complete data frame
quarterly.light.mean[59,c(1,2)] <- c("20:15", "NA") #add row of time to complete data frame
quarterly.light.mean$light <- as.numeric(quarterly.light.mean$light)
quarterly.light.se <- aggregate(light ~ quarterhours, data=Acc.data, std.error, na.rm=TRUE) #calculate standard error of the mean of light for every 15 min interval
quarterly.light.se <- rbind(c("05:45","NA"), quarterly.light.se) #combine data 
quarterly.light.se[56,c(1,2)] <- c("19:30", "NA") #add row of time to complete data frame
quarterly.light.se[57,c(1,2)] <- c("19:45", "NA") #add row of time to complete data frame
quarterly.light.se[58,c(1,2)] <- c("20:00", "NA") #add row of time to complete data frame
quarterly.light.se[59,c(1,2)] <- c("20:15", "NA") #add row of time to complete data frame
quarterly.light.se$light <- as.numeric(quarterly.light.se$light) #set column as numeric
light.means <- data.frame(quarterly.light.mean, quarterly.light.se$light) #combine mean and standard error results
colnames(light.means) <- c("Time", "mean", "se")  #rename columns to describe contents
range(na.omit(light.means$mean)) #range of average light levels

Fig3 <- ggplot(light.means) + #Plot average diurnal cycle of light data
  geom_point(aes(x = Time, y = mean), colour="black", cex=0.8) + #Plot points using time as the x axis, light as the Y axis and black points
  geom_errorbar(aes(x=Time, ymax=mean+se, ymin=mean-se), data=light.means, col="black", width=0) + #set values for standard error bars
  scale_x_discrete(breaks=c("0:00", "06:00", "12:00", "18:00")) + #set discrete breaks on the X axis
  ylim(0,300) + #set Y axis limits
  ggtitle("(e)") + #Label the graph with the main title
  xlab("Time") + #Label the X Axis
  ylab(bquote('Irradiance ('*mu~'mol' ~photons ~ m^-2~s^-1*')')) + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), angle = 90, vjust = 0.5, hjust=1), #set margins on labels
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        plot.title=element_text(hjust=0)) #Justify the title to the top left
Fig3 #View figure

##### TANK EXPERIMENTAL DATA #####
#Tank Temperature Data for Adult Exposure Experimental Period (06May14 - 17Aug14)
tank.data <- read.csv("Adult_Tank_Temp.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
mydate.tanks <- strptime(tank.data$Date.Time, format="%m/%d/%y %H:%M") #convert date format to characters

tank.data$Tank3.cal <-(tank.data$Tank3*T3.lm[2])+T3.lm[1] #Apply the cross calibration of temperature to standard logger #1
tank.data$Tank4.cal <-(tank.data$Tank4*T4.lm[2])+T4.lm[1] #Apply the cross calibration of temperature to standard logger #1
tank.data$Tank5.cal <-(tank.data$Tank5*T5.lm[2])+T5.lm[1] #Apply the cross calibration of temperature to standard logger #1
tank.data$Tank5.cal[is.na(tank.data$Tank5.cal)] <- tank.data$Tank3.cal[is.na(tank.data$Tank5.cal)] #merge tanks 3 and 5 into one column so the data from first 3 days when corals were in tank 3 is now showing in tank 5
tank.tempdata <-data.frame(mydate.tanks, tank.data$Tank4.cal, tank.data$Tank5.cal) #make a dataframe of temperature and time
colnames(tank.tempdata) <- c("Date.Time", "Tank4", "Tank5")
Tank4.temp.N <- sum(!is.na(tank.tempdata$Tank4)) #Count sample size
Tank5.temp.N <- sum(!is.na(tank.tempdata$Tank5)) #Count sample size
range(na.omit(tank.tempdata$Tank4)) #range of average ambient temp levels
range(na.omit(tank.tempdata$Tank5)) #range of average high temp levels

Fig4 <- ggplot(tank.tempdata, aes(Date.Time)) + #plot tank temperature data
  geom_line(aes(y = Tank4, colour="Ambient")) + #plot Temperature data as a line on the Y axis with date as the X axis 
  geom_line(aes(y = Tank5, colour="High")) + #plot Temperature data as a line on the Y axis with date as the X axis 
  scale_colour_manual("Treatment", values = c("grey","black")) + #add colors for treatments
  xlab("Date") + #Label the X Axis
  ylab("Temperature °C") + #Label the Y Axis
  ggtitle("") + #label the main title
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), angle = 90, vjust = 0.5, hjust=1), #set margins on labels
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        plot.title=element_text(hjust=0)) #Justify the title to the top left
Fig4 #View figure

tanks.T.mean <- mean(rbind(tank.tempdata[,2], tank.tempdata[,3]), na.rm=T) #Calculate the grand average of both tanks
tanks.T.mean #View data
tanks.T.se <- sd(rbind(tank.tempdata[,2], tank.tempdata[,3]), na.rm=T)/sqrt(Tank4.temp.N+Tank5.temp.N) #Calculate the overal standard error of both tanks
tanks.T.se #View data

#Plotting diurnal cycles
tank.time <- format(as.POSIXct(mydate.tanks) ,format = "%H:%M") #Format time into only hours and minutes
tank.temperatures <- cbind(tank.tempdata, tank.time) #create a dataframe 
colnames(tank.temperatures) <- c("Date", "Tank4", "Tank5", "Time") #Rename columns to describe contents
tank.temperatures #View Data

quarterly.tank.temp.mean4 <- aggregate(Tank4 ~ Time, data=tank.temperatures, mean, na.rm=TRUE) #calculate mean of temperature for every 15 min interval
quarterly.tank.temp.se4 <- aggregate(Tank4 ~ Time, data=tank.temperatures, std.error, na.rm=TRUE)  #calculate standard error of the mean of temperature for every 15 min interval
quarterly.tank.temp.mean5 <- aggregate(Tank5 ~ Time, data=tank.temperatures, mean, na.rm=TRUE) #calculate mean of temperature for every 15 min interval
quarterly.tank.temp.se5 <- aggregate(Tank5 ~ Time, data=tank.temperatures, std.error, na.rm=TRUE)  #calculate standard error of the mean of temperature for every 15 min interval
tank.temp.means <- data.frame(quarterly.tank.temp.mean4, quarterly.tank.temp.se4$Tank4, quarterly.tank.temp.mean5$Tank5, quarterly.tank.temp.se5$Tank5) #combine mean and standard error results
colnames(tank.temp.means) <- c("Time", "Tank4.mean", "Tank4.se", "Tank5.mean", "Tank5.se")  #Rename columns to describe contents
range(na.omit(tank.temp.means$Tank4.mean)) #range of average ambient temp levels
range(na.omit(tank.temp.means$Tank5.mean)) #range of average high temp levels

Fig5 <- ggplot(tank.temp.means, aes(Time)) + # plot mean temp by tank
  geom_point(aes(y = Tank4.mean, colour="Ambient"), cex=0.8) + #plot points
  geom_errorbar(aes(x=Time, ymax=Tank4.mean+Tank4.se, ymin=Tank4.mean-Tank4.se), position=position_dodge(0.9), data=tank.temp.means, col="darkgrey", width=0) + #set values for standard error bars and offset on the X axis for clarity
  geom_point(aes(y = Tank5.mean, colour="High"), cex=0.8) + #plot points
  geom_errorbar(aes(x=Time, ymax=Tank5.mean+Tank5.se, ymin=Tank5.mean-Tank5.se), position=position_dodge(0.9), data=tank.temp.means, col="black", width=0) + #set values for standard error bars and offset on the X axis for clarity
  scale_colour_manual("Treatment", values = c("grey","black")) +
  scale_x_discrete(breaks=c("01:00", "06:00", "12:00", "18:00", "23:00")) + #set discrete breaks on the X axis
  ggtitle("ADULT EXPOSURE\n(b)") + #Label graphic title
  ylim(23.5,29) + #Set Y axis limits
  xlab("Time") + #Label the X Axis
  ylab("Temperature (°C)") + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), angle = 90, vjust = 0.5, hjust=1), #set margins on labels
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        plot.title=element_text(hjust=0), #Justify the title to the top left
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        legend.position='none') #remove legend background
Fig5 #View figure

#Tank Light Data for Adult Exposure Experimental Period (06May14 - 17Aug14)
tank.light.data <- read.csv("Adult_Tank_Light.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
tank.light.data[tank.light.data == 0] <- NA
mydate.tanks <- strptime(tank.light.data$Date.Time, format="%m/%d/%y %H:%M") #convert date format to characters
#Convert into instantaneous PAR
tank.light.data$Tank3.quanta <-(tank.light.data$Tank3)/(15*60) #Assign light column in dataframe and convert to units of µmol m-2 s-1
tank.light.data$Tank4.quanta <-(tank.light.data$Tank4)/(15*60) #Assign light column in dataframe and convert to units of µmol m-2 s-1
tank.light.data$Tank5.quanta <-(tank.light.data$Tank5)/(15*60) #Assign light column in dataframe and convert to units of µmol m-2 s-1
#apply calibration models
tank.light.data$Tank3.cal <-(tank.light.data$Tank3.quanta*L3.lm[2])+L3.lm[1] #Apply the cross calibration of temperature to standard logger #1
tank.light.data$Tank4.cal <-(tank.light.data$Tank4.quanta*L4.lm[2])+L4.lm[1] #Apply the cross calibration of temperature to standard logger #1
tank.light.data$Tank5.cal <-(tank.light.data$Tank5.quanta*L5.lm[2])+L5.lm[1] #Apply the cross calibration of temperature to standard logger #1
tank.light.data$Tank5.cal[is.na(tank.light.data$Tank5.cal)] <- tank.light.data$Tank3.cal[is.na(tank.light.data$Tank5.cal)] #merge tanks 3 and 5 into one column so the data from first 3 days when corals were in tank 3 is now showing in tank 5
tank.lightdata <-data.frame(mydate.tanks, tank.light.data$Tank4.cal, tank.light.data$Tank5.cal) #make a dataframe of temperature and time
colnames(tank.lightdata) <- c("Date.Time", "Tank4", "Tank5") #rename columns
range(na.omit(tank.lightdata$Tank4)) #range ambient light levels
range(na.omit(tank.lightdata$Tank5)) #range high light levels

Tank4.light.N <- sum(!is.na(tank.lightdata$Tank4)) #Count sample size
Tank5.light.N <- sum(!is.na(tank.lightdata$Tank5)) #Count sample size

Fig6 <- ggplot(tank.lightdata, aes(Date.Time)) + #plot tank light data
  geom_line(aes(y = Tank4, colour="Ambient")) + #plot light data as a line on the Y axis with date as the X axis 
  geom_line(aes(y = Tank5, colour="High")) + #plot light data as a line on the Y axis with date as the X axis 
  scale_colour_manual("Treatment", values = c("grey","black")) + #add colors for treatments
  xlab("Date") + #Label the X Axis
  ylab(bquote('Irradiance ('*mu~'mol' ~photons ~ m^-2~s^-1*')')) + #Label the Y Axis
  ggtitle("") + #label the main title
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background =element_blank(), #Set the plot background
        legend.key = element_blank(), #Set plot legend key
        plot.title=element_text(hjust=0)) #Justify the title to the top left
Fig6 #View figure

tanks.L.mean <- mean(rbind(tank.lightdata[,2], tank.lightdata[,3]), na.rm=T) #Calculate the grand average of both tanks
tanks.L.mean #View data
tanks.L.se <- sd(rbind(tank.lightdata[,2], tank.lightdata[,3]), na.rm=T)/sqrt(Tank4.light.N+Tank5.light.N) #Calculate the overal standard error of both tanks
tanks.L.se #View data

#Plotting diurnal cycles
tank.time <- format(as.POSIXct(mydate.tanks) ,format = "%H:%M") #Format time into only hours and minutes
tank.lights <- cbind(tank.lightdata, tank.time) #create a dataframe 
colnames(tank.lights) <- c("Date", "Tank4", "Tank5", "Time") #Rename columns to describe contents
tank.lights #View Data

quarterly.tank.light.mean4 <- aggregate(Tank4 ~ Time, data=tank.lights, mean, na.rm=TRUE) #calculate mean of temperature for every 15 min interval
quarterly.tank.light.se4 <- aggregate(Tank4 ~ Time, data=tank.lights, std.error, na.rm=TRUE)  #calculate standard error of the mean of temperature for every 15 min interval
quarterly.tank.light.mean5 <- aggregate(Tank5 ~ Time, data=tank.lights, mean, na.rm=TRUE) #calculate mean of temperature for every 15 min interval
quarterly.tank.light.mean5[58,1] <- "20:00" #add empty row 
quarterly.tank.light.mean5[59,1] <- "20:15" #add empty row
quarterly.tank.light.se5 <- aggregate(Tank5 ~ Time, data=tank.lights, std.error, na.rm=TRUE)  #calculate standard error of the mean of temperature for every 15 min interval
quarterly.tank.light.se5[58,1] <- "20:00" #add empty row
quarterly.tank.light.se5[59,1] <- "20:15" #add empty row
tank.light.means <- data.frame(quarterly.tank.light.mean4, quarterly.tank.light.se4$Tank4, quarterly.tank.light.mean5$Tank5, quarterly.tank.light.se5$Tank5) #combine mean and standard error results
colnames(tank.light.means) <- c("Time", "Tank4.mean", "Tank4.se", "Tank5.mean", "Tank5.se")  #Rename columns to describe contents
range(na.omit(tank.light.means$Tank4.mean)) #range of average ambient light levels
range(na.omit(tank.light.means$Tank5.mean)) #range of average high light levels

Fig7 <- ggplot(tank.light.means, aes(Time)) + # plot mean temp by tank
  geom_point(aes(y = Tank4.mean, colour="Ambient"), cex=0.8) + #plot points
  geom_errorbar(aes(x=Time, ymax=Tank4.mean+Tank4.se, ymin=Tank4.mean-Tank4.se), position=position_dodge(0.9), data=tank.light.means, col="darkgrey", width=0) + #set values for standard error bars and offset on the X axis for clarity
  geom_point(aes(y = Tank5.mean, colour="High"), cex=0.8) + #plot points
  geom_errorbar(aes(x=Time, ymax=Tank5.mean+Tank5.se, ymin=Tank5.mean-Tank5.se), position=position_dodge(0.9), data=tank.light.means, col="black", width=0) + #set values for standard error bars and offset on the X axis for clarity
  scale_colour_manual("Treatment", values = c("grey","black")) +
  scale_x_discrete(breaks=c("0:00", "06:00", "12:00", "18:00")) + #set discrete breaks on the X axis
  ylim(0,300) + #set Y axis limits
  ggtitle("(f)") + #Label graphic title
  xlab("Time") + #Label the X Axis
  ylab(bquote('Irradiance ('*mu~'mol' ~photons ~ m^-2~s^-1*')')) + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), angle = 90, vjust = 0.5, hjust=1), #set margins on labels
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        plot.title=element_text(hjust=0), #Justify the title to the top left
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        legend.position='none') #remove legend background
Fig7 #View figure

#Tank Temperature Data for Larval Month 1 Exposure Period (12July - 25Aug14)
M1.tank.data <- read.csv("Month1_Tank_Temp.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
mydate.tanks <- strptime(M1.tank.data$Date.Time, format="%m/%d/%y %H:%M") #convert date format to characters
M1.tank.data$Tank4.cal <-(M1.tank.data$Tank4*T4.lm[2])+T4.lm[1] #Apply the cross calibration of temperature to standard logger #1
M1.tank.data$Tank5.cal <-(M1.tank.data$Tank5*T5.lm[2])+T5.lm[1] #Apply the cross calibration of temperature to standard logger #1
M1.tank.tempdata <-data.frame(mydate.tanks, M1.tank.data$Tank4.cal, M1.tank.data$Tank5.cal) #make a dataframe of temperature and time
colnames(M1.tank.tempdata) <- c("Date.Time", "Tank4", "Tank5")
Tank4.M1.temp.N <- sum(!is.na(M1.tank.tempdata$Tank4)) #Count sample size
Tank5.M1.temp.N <- sum(!is.na(M1.tank.tempdata$Tank5)) #Count sample size

Fig8 <- ggplot(M1.tank.tempdata, aes(Date.Time)) + #plot tank temperature data
  geom_line(aes(y = Tank4, colour="Ambient")) + #plot Temperature data as a line on the Y axis with date as the X axis 
  geom_line(aes(y = Tank5, colour="High")) + #plot Temperature data as a line on the Y axis with date as the X axis 
  scale_colour_manual("Treatment", values = c("grey","black")) + #add colors for treatments
  xlab("Date") + #Label the X Axis
  ylab("Temperature °C") + #Label the Y Axis
  ggtitle("") + #label the main title
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background =element_blank(), #Set the plot background
        legend.key = element_blank(), #Set plot legend key
        plot.title=element_text(hjust=0)) #Justify the title to the top left
Fig8 #View figure

M1.tanks.T.mean <- mean(rbind(M1.tank.tempdata[,2], M1.tank.tempdata[,3]), na.rm=T) #Calculate the grand average of both tanks
M1.tanks.T.mean #View data
M1.tanks.T.se <- sd(rbind(M1.tank.tempdata[,2], M1.tank.tempdata[,3]), na.rm=T)/sqrt(Tank4.M1.temp.N+Tank5.M1.temp.N) #Calculate the overal standard error of both tanks
M1.tanks.T.se #View data

#Plotting diurnal cycles
M1.tank.time <- format(as.POSIXct(mydate.tanks) ,format = "%H:%M") #Format time into only hours and minutes
M1.tank.temperatures <- cbind(M1.tank.tempdata, M1.tank.time) #create a dataframe 
colnames(M1.tank.temperatures) <- c("Date", "Tank4", "Tank5", "Time") #Rename columns to describe contents
M1.tank.temperatures #View Data
quarterly.M1.tank.temp.mean4 <- aggregate(Tank4 ~ Time, data=M1.tank.temperatures, mean, na.rm=TRUE) #calculate mean of temperature for every 15 min interval
quarterly.M1.tank.temp.se4 <- aggregate(Tank4 ~ Time, data=M1.tank.temperatures, std.error, na.rm=TRUE)  #calculate standard error of the mean of temperature for every 15 min interval
quarterly.M1.tank.temp.mean5 <- aggregate(Tank5 ~ Time, data=M1.tank.temperatures, mean, na.rm=TRUE) #calculate mean of temperature for every 15 min interval
quarterly.M1.tank.temp.se5 <- aggregate(Tank5 ~ Time, data=M1.tank.temperatures, std.error, na.rm=TRUE)  #calculate standard error of the mean of temperature for every 15 min interval
M1.tank.temp.means <- data.frame(quarterly.M1.tank.temp.mean4, quarterly.M1.tank.temp.se4$Tank4, quarterly.M1.tank.temp.mean5$Tank5, quarterly.M1.tank.temp.se5$Tank5) #combine mean and standard error results
colnames(M1.tank.temp.means) <- c("Time", "Tank4.mean", "Tank4.se", "Tank5.mean", "Tank5.se")  #Rename columns to describe contents
range(na.omit(M1.tank.temp.means$Tank4.mean)) #range of average amb temp levels
range(na.omit(M1.tank.temp.means$Tank5.mean)) #range of average high temp levels

Fig9 <- ggplot(M1.tank.temp.means, aes(Time)) + # plot mean temp by tank
  geom_point(aes(y = Tank4.mean, colour="Ambient"), cex=0.8) + #plot points
  geom_errorbar(aes(x=Time, ymax=Tank4.mean+Tank4.se, ymin=Tank4.mean-Tank4.se), position=position_dodge(0.9), data=M1.tank.temp.means, col="darkgrey", width=0) + #set values for standard error bars and offset on the X axis for clarity
  geom_point(aes(y = Tank5.mean, colour="High"), cex=0.8) + #plot points
  geom_errorbar(aes(x=Time, ymax=Tank5.mean+Tank5.se, ymin=Tank5.mean-Tank5.se), position=position_dodge(0.9), data=M1.tank.temp.means, col="black", width=0) + #set values for standard error bars and offset on the X axis for clarity
  scale_colour_manual("Treatment", values = c("grey","black")) +
  scale_x_discrete(breaks=c("01:00", "06:00", "12:00", "18:00", "23:00")) + #set discrete breaks on the X axis
  ggtitle("LARVAL MONTH 1\n(c)") + #Label graphic title
  ylim(23.5,29) + #Set Y axis limits
  xlab("Time") + #Label the X Axis
  ylab("Temperature (°C)") + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), angle = 90, vjust = 0.5, hjust=1), #set margins on labels
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        plot.title=element_text(hjust=0), #Justify the title to the top left
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        legend.position='none') #remove legend background
Fig9 #View figure

#Tank Light Data for Larval Month 1 Exposure Period (12July - 25Aug14)
M1.tank.light.data <- read.csv("Month1_Tank_Light.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
M1.tank.light.data[M1.tank.light.data == 0] <- NA
mydate.tanks <- strptime(M1.tank.light.data$Date.Time, format="%m/%d/%y %H:%M") #convert date format to characters
M1.tank.light.data$Tank4.quanta <-(M1.tank.light.data$Tank4)/(15*60) #Assign light column in dataframe and convert to units of µmol m-2 s-1
M1.tank.light.data$Tank5.quanta <-(M1.tank.light.data$Tank5)/(15*60) #Assign light column in dataframe and convert to units of µmol m-2 s-1
M1.tank.light.data$Tank4.cal <-(M1.tank.light.data$Tank4.quanta*L4.lm[2])+L4.lm[1] #Apply the cross calibration of temperature to standard logger #1
M1.tank.light.data$Tank5.cal <-(M1.tank.light.data$Tank5.quanta*L5.lm[2])+L5.lm[1] #Apply the cross calibration of temperature to standard logger #1
M1.tank.lightdata <-data.frame(mydate.tanks, M1.tank.light.data$Tank4.cal, M1.tank.light.data$Tank5.cal) #make a dataframe of temperature and time
colnames(M1.tank.lightdata) <- c("Date.Time", "Tank4", "Tank5")
Tank4.M1.light.N <- sum(!is.na(M1.tank.lightdata$Tank4)) #Count sample size
Tank5.M1.light.N <- sum(!is.na(M1.tank.lightdata$Tank5)) #Count sample size

Fig10 <- ggplot(M1.tank.lightdata, aes(Date.Time)) + #plot tank light data
  geom_line(aes(y = Tank4, colour="Ambient")) + #plot light data as a line on the Y axis with date as the X axis 
  geom_line(aes(y = Tank5, colour="High")) + #plot light data as a line on the Y axis with date as the X axis 
  scale_colour_manual("Treatment", values = c("grey","black")) + #add colors for treatments
  xlab("Date") + #Label the X Axis
  ylab(bquote('Irradiance ('*mu~'mol' ~photons ~ m^-2~s^-1*')')) + #Label the Y Axis
  ggtitle("") + #label the main title
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background =element_blank(), #Set the plot background
        legend.key = element_blank(), #Set plot legend key
        plot.title=element_text(hjust=0)) #Justify the title to the top left
Fig10 #View figure

M1.tanks.L.mean <- mean(rbind(M1.tank.lightdata[,2], M1.tank.lightdata[,3]), na.rm=T) #Calculate the grand average of both tanks
M1.tanks.L.mean #View data
M1.tanks.L.se <- sd(rbind(M1.tank.lightdata[,2], M1.tank.lightdata[,3]), na.rm=T)/sqrt(Tank4.M1.light.N+Tank5.M1.light.N) #Calculate the overal standard error of both tanks
M1.tanks.L.se #View data

#Plotting diurnal cycles
M1.tank.time <- format(as.POSIXct(mydate.tanks) ,format = "%H:%M") #Format time into only hours and minutes
M1.tank.lights <- cbind(M1.tank.lightdata, M1.tank.time) #create a dataframe 
colnames(M1.tank.lights) <- c("Date", "Tank4", "Tank5", "Time") #Rename columns to describe contents
M1.tank.lights #View Data

quarterly.M1.tank.light.mean4 <- aggregate(Tank4 ~ Time, data=M1.tank.lights, mean, na.rm=TRUE) #calculate mean for every 15 min interval
quarterly.M1.tank.light.mean4 <- rbind(c("06:00","NA"), quarterly.M1.tank.light.mean4)
quarterly.M1.tank.light.mean4 <- rbind(c("05:45","NA"), quarterly.M1.tank.light.mean4)
quarterly.M1.tank.light.mean4[56,c(1,2)] <- c("19:30", "NA") #add empty row
quarterly.M1.tank.light.mean4[57,c(1,2)] <- c("19:45", "NA") #add empty row
quarterly.M1.tank.light.mean4[58,c(1,2)] <- c("20:00", "NA") #add empty row
quarterly.M1.tank.light.mean4[59,c(1,2)] <- c("20:15", "NA") #add empty row
quarterly.M1.tank.light.mean4$Tank4 <- as.numeric(quarterly.M1.tank.light.mean4$Tank4) #set as numeric
quarterly.M1.tank.light.se4 <- aggregate(Tank4 ~ Time, data=M1.tank.lights, std.error, na.rm=TRUE)  #calculate standard error of the mean for every 15 min interval
quarterly.M1.tank.light.se4 <- rbind(c("06:00","NA"), quarterly.M1.tank.light.se4) #combine rows 
quarterly.M1.tank.light.se4 <- rbind(c("05:45","NA"), quarterly.M1.tank.light.se4) #combine rows
quarterly.M1.tank.light.se4[56,c(1,2)] <- c("19:30", "NA") #add empty row
quarterly.M1.tank.light.se4[57,c(1,2)] <- c("19:45", "NA") #add empty row
quarterly.M1.tank.light.se4[58,c(1,2)] <- c("20:00", "NA") #add empty row
quarterly.M1.tank.light.se4[59,c(1,2)] <- c("20:15", "NA") #add empty row
quarterly.M1.tank.light.se4$Tank4 <- as.numeric(quarterly.M1.tank.light.se4$Tank4) #set as numeric
quarterly.M1.tank.light.mean5 <- aggregate(Tank5 ~ Time, data=M1.tank.lights, mean, na.rm=TRUE) #calculate mean for every 15 min interval
quarterly.M1.tank.light.mean5 <- rbind(c("05:45","NA"), quarterly.M1.tank.light.mean5) #combine rows
quarterly.M1.tank.light.mean5[57,c(1,2)] <- c("19:45", "NA") #add empty row
quarterly.M1.tank.light.mean5[58,c(1,2)] <- c("20:00", "NA") #add empty row
quarterly.M1.tank.light.mean5[59,c(1,2)] <- c("20:15", "NA") #add empty row
quarterly.M1.tank.light.mean5$Tank5 <- as.numeric(quarterly.M1.tank.light.mean5$Tank5) #set as numeric
quarterly.M1.tank.light.se5 <- aggregate(Tank5 ~ Time, data=M1.tank.lights, std.error, na.rm=TRUE)  #calculate standard error of the mean  for every 15 min interval
quarterly.M1.tank.light.se5 <- rbind(c("05:45","NA"), quarterly.M1.tank.light.se5) #combine rows
quarterly.M1.tank.light.se5[57,c(1,2)] <- c("19:45", "NA") #add empty row
quarterly.M1.tank.light.se5[58,c(1,2)] <- c("20:00", "NA") #add empty row
quarterly.M1.tank.light.se5[59,c(1,2)] <- c("20:15", "NA") #add empty row
quarterly.M1.tank.light.se5$Tank5 <- as.numeric(quarterly.M1.tank.light.se5$Tank5) #set as numeric
M1.tank.light.means <-data.frame(quarterly.M1.tank.light.mean4,quarterly.M1.tank.light.se4$Tank4,quarterly.M1.tank.light.mean5$Tank5,quarterly.M1.tank.light.se5$Tank5) #make a dataframe 
colnames(M1.tank.light.means) <- c("Time", "Tank4.mean", "Tank4.se", "Tank5.mean", "Tank5.se")  #Rename columns to describe contents
range(na.omit(M1.tank.light.means$Tank4.mean)) #range of average ambient light levels
range(na.omit(M1.tank.light.means$Tank5.mean)) #minimum of average high light levels

Fig11 <- ggplot(M1.tank.light.means, aes(Time)) + # plot mean light by tank
  geom_point(aes(y = Tank4.mean, colour="Ambient"), cex=0.8) + #plot points
  geom_errorbar(aes(x=Time, ymax=Tank4.mean+Tank4.se, ymin=Tank4.mean-Tank4.se), position=position_dodge(0.9), data=M1.tank.light.means, col="darkgrey", width=0) + #set values for standard error bars and offset on the X axis for clarity
  geom_point(aes(y = Tank5.mean, colour="High"), cex=0.8) + #plot points
  geom_errorbar(aes(x=Time, ymax=Tank5.mean+Tank5.se, ymin=Tank5.mean-Tank5.se), position=position_dodge(0.9), data=M1.tank.light.means, col="black", width=0) + #set values for standard error bars and offset on the X axis for clarity
  scale_colour_manual("Treatment", values = c("grey","black")) +
  scale_x_discrete(breaks=c("0:00", "06:00", "12:00", "18:00")) + #set discrete breaks on the X axis
  ylim(0,300) + #set Y axis limits
  ggtitle("(g)") + #Label graphic title
  xlab("Time") + #Label the X Axis
  ylab(bquote('Irradiance ('*mu~'mol' ~photons ~ m^-2~s^-1*')')) + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), angle = 90, vjust = 0.5, hjust=1), #set margins on labels
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        plot.title=element_text(hjust=0), #Justify the title to the top left
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        legend.position='none') #remove legend background
Fig11 #View figure

#Tank Temperature Data for Larval Month 6 Exposure Period (12July - 28January15)
M6.tank.data <- read.csv("Month6_Tank_Temp.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
mydate.tanks <- strptime(M6.tank.data$Date.Time, format="%m/%d/%y %H:%M") #convert date format to characters
M6.tank.data$Tank4.cal <-(M6.tank.data$Tank4*T4.lm[2])+T4.lm[1] #Apply the cross calibration of temperature to standard logger #1
M6.tank.data$Tank5.cal <-(M6.tank.data$Tank5*T5.lm[2])+T5.lm[1] #Apply the cross calibration of temperature to standard logger #1
M6.tank.tempdata <-data.frame(mydate.tanks, M6.tank.data$Tank4.cal, M6.tank.data$Tank5.cal) #make a dataframe of temperature and time
colnames(M6.tank.tempdata) <- c("Date.Time", "Tank4", "Tank5")
Tank4.M6.temp.N <- sum(!is.na(M6.tank.tempdata$Tank4)) #Count sample size
Tank5.M6.temp.N <- sum(!is.na(M6.tank.tempdata$Tank5)) #Count sample size

Fig12 <- ggplot(M6.tank.tempdata, aes(Date.Time)) + #plot tank temperature data
  geom_line(aes(y = Tank4, colour="Ambient")) + #plot Temperature data as a line on the Y axis with date as the X axis 
  geom_line(aes(y = Tank5, colour="High")) + #plot Temperature data as a line on the Y axis with date as the X axis 
  scale_colour_manual("Treatment", values = c("grey","black")) + #add colors for treatments
  xlab("Date") + #Label the X Axis
  ylab("Temperature °C") + #Label the Y Axis
  ggtitle("") + #label the main title
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background =element_blank(), #Set the plot background
        legend.key = element_blank(), #Set plot legend key
        plot.title=element_text(hjust=0)) #Justify the title to the top left
Fig12 #View figure

M6.tanks.T.mean <- mean(rbind(M6.tank.tempdata[,2], M6.tank.tempdata[,3]), na.rm=T) #Calculate the grand average of both tanks
M6.tanks.T.mean #View data
M6.tanks.T.se <- sd(rbind(M6.tank.tempdata[,2], M6.tank.tempdata[,3]), na.rm=T)/sqrt(Tank4.M6.temp.N+Tank5.M6.temp.N) #Calculate the overal standard error of both tanks
M6.tanks.T.se #View data

#Plotting diurnal cycles
M6.tank.time <- format(as.POSIXct(mydate.tanks) ,format = "%H:%M") #Format time into only hours and minutes
M6.tank.temperatures <- cbind(M6.tank.tempdata, M6.tank.time) #create a dataframe 
colnames(M6.tank.temperatures) <- c("Date", "Tank4", "Tank5", "Time") #Rename columns to describe contents
M6.tank.temperatures #View Data

quarterly.M6.tank.temp.mean4 <- aggregate(Tank4 ~ Time, data=M6.tank.temperatures, mean, na.rm=TRUE) #calculate mean of temperature for every 15 min interval
quarterly.M6.tank.temp.se4 <- aggregate(Tank4 ~ Time, data=M6.tank.temperatures, std.error, na.rm=TRUE)  #calculate standard error of the mean of temperature for every 15 min interval
quarterly.M6.tank.temp.mean5 <- aggregate(Tank5 ~ Time, data=M6.tank.temperatures, mean, na.rm=TRUE) #calculate mean of temperature for every 15 min interval
quarterly.M6.tank.temp.se5 <- aggregate(Tank5 ~ Time, data=M6.tank.temperatures, std.error, na.rm=TRUE)  #calculate standard error of the mean of temperature for every 15 min interval
M6.tank.temp.means <- data.frame(quarterly.M6.tank.temp.mean4, quarterly.M6.tank.temp.se4$Tank4, quarterly.M6.tank.temp.mean5$Tank5, quarterly.M6.tank.temp.se5$Tank5) #combine mean and standard error results
colnames(M6.tank.temp.means) <- c("Time", "Tank4.mean", "Tank4.se", "Tank5.mean", "Tank5.se")  #Rename columns to describe contents
range(na.omit(M6.tank.temp.means$Tank4.mean)) #range of average amb temp levels
range(na.omit(M6.tank.temp.means$Tank5.mean)) #range of average high temp levels

Fig13 <- ggplot(M6.tank.temp.means, aes(Time)) + # plot mean temp by tank
  geom_point(aes(y = Tank4.mean, colour="Ambient"), cex=0.8) + #plot points
  geom_errorbar(aes(x=Time, ymax=Tank4.mean+Tank4.se, ymin=Tank4.mean-Tank4.se), position=position_dodge(0.9), data=M6.tank.temp.means, col="darkgrey", width=0) + #set values for standard error bars and offset on the X axis for clarity
  geom_point(aes(y = Tank5.mean, colour="High"), cex=0.8) + #plot points
  geom_errorbar(aes(x=Time, ymax=Tank5.mean+Tank5.se, ymin=Tank5.mean-Tank5.se), position=position_dodge(0.9), data=M6.tank.temp.means, col="black", width=0) + #set values for standard error bars and offset on the X axis for clarity
  scale_colour_manual("Treatment", values = c("grey","black")) +
  scale_x_discrete(breaks=c("01:00", "06:00", "12:00", "18:00", "23:00")) + #set discrete breaks on the X axis
  ggtitle("LARVAL MONTH 6\n(d)") + #Label graphic title
  ylim(23.5,29) + #Set Y axis limits
  xlab("Time") + #Label the X Axis
  ylab("Temperature (°C)") + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), angle = 90, vjust = 0.5, hjust=1), #set margins on labels
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        plot.title=element_text(hjust=0), #Justify the title to the top left
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        legend.position='none') #remove legend background
Fig13 #View figure

#Tank Light Data for Larval Month 6 Exposure Period (12July - 25Aug14)
M6.tank.light.data <- read.csv("Month6_Tank_Light.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
M6.tank.light.data[M6.tank.light.data == 0] <- NA
mydate.tanks <- strptime(M6.tank.light.data$Date.Time, format="%m/%d/%y %H:%M") #convert date format to characters
M6.tank.light.data$Tank4.quanta <-(M6.tank.light.data$Tank4)/(15*60) #Assign light column in dataframe and convert to units of µmol m-2 s-1
M6.tank.light.data$Tank5.quanta <-(M6.tank.light.data$Tank5)/(15*60) #Assign light column in dataframe and convert to units of µmol m-2 s-1
M6.tank.light.data$Tank4.cal <-(M6.tank.light.data$Tank4.quanta*L4.lm[2])+L4.lm[1] #Apply the cross calibration of temperature to standard logger #1
M6.tank.light.data$Tank5.cal <-(M6.tank.light.data$Tank5.quanta*L5.lm[2])+L5.lm[1] #Apply the cross calibration of temperature to standard logger #1
M6.tank.lightdata <-data.frame(mydate.tanks, M6.tank.light.data$Tank4.cal, M6.tank.light.data$Tank5.cal) #make a dataframe of temperature and time
colnames(M6.tank.lightdata) <- c("Date.Time", "Tank4", "Tank5")
Tank4.M6.light.N <- sum(!is.na(M6.tank.lightdata$Tank4)) #Count sample size
Tank5.M6.light.N <- sum(!is.na(M6.tank.lightdata$Tank5)) #Count sample size

Fig14 <- ggplot(M6.tank.lightdata, aes(Date.Time)) + #plot tank light data
  geom_line(aes(y = Tank4, colour="Ambient")) + #plot light data as a line on the Y axis with date as the X axis 
  geom_line(aes(y = Tank5, colour="High")) + #plot light data as a line on the Y axis with date as the X axis 
  scale_colour_manual("Treatment", values = c("grey","black")) + #add colors for treatments
  xlab("Date") + #Label the X Axis
  ylab(bquote('Irradiance ('*mu~'mol' ~photons ~ m^-2~s^-1*')')) + #Label the Y Axis
  ggtitle("") + #label the main title
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background =element_blank(), #Set the plot background
        legend.key = element_blank(), #Set plot legend key
        plot.title=element_text(hjust=0)) #Justify the title to the top left
Fig14 #View figure

M6.tanks.L.mean <- mean(rbind(M6.tank.lightdata[,2], M6.tank.lightdata[,3]), na.rm=T) #Calculate the grand average of both tanks
M6.tanks.L.mean #View data
M6.tanks.L.se <- sd(rbind(M6.tank.lightdata[,2], M6.tank.lightdata[,3]), na.rm=T)/sqrt(Tank4.M6.light.N+Tank5.M6.light.N) #Calculate the overal standard error of both tanks
M6.tanks.L.se #View data

#Plotting diurnal cycles
M6.tank.time <- format(as.POSIXct(mydate.tanks) ,format = "%H:%M") #Format time into only hours and minutes
M6.tank.lights <- cbind(M6.tank.lightdata, M6.tank.time) #create a dataframe 
colnames(M6.tank.lights) <- c("Date", "Tank4", "Tank5", "Time") #Rename columns to describe contents
M6.tank.lights #View Data

quarterly.M6.tank.light.mean4 <- aggregate(Tank4 ~ Time, data=M6.tank.lights, mean, na.rm=TRUE) #calculate mean for every 15 min interval
quarterly.M6.tank.light.mean4 <- rbind(c("06:00","NA"), quarterly.M6.tank.light.mean4) #combine rows 
quarterly.M6.tank.light.mean4 <- rbind(c("05:45","NA"), quarterly.M6.tank.light.mean4) #combine rows 
quarterly.M6.tank.light.mean4[56,c(1,2)] <- c("19:30", "NA") #add empty row
quarterly.M6.tank.light.mean4[57,c(1,2)] <- c("19:45", "NA") #add empty row
quarterly.M6.tank.light.mean4[58,c(1,2)] <- c("20:00", "NA") #add empty row
quarterly.M6.tank.light.mean4[59,c(1,2)] <- c("20:15", "NA") #add empty row
quarterly.M6.tank.light.mean4$Tank4 <- as.numeric(quarterly.M6.tank.light.mean4$Tank4) #set as numeric
quarterly.M6.tank.light.se4 <- aggregate(Tank4 ~ Time, data=M6.tank.lights, std.error, na.rm=TRUE)  #calculate standard error of the mean for every 15 min interval
quarterly.M6.tank.light.se4 <- rbind(c("06:00","NA"), quarterly.M6.tank.light.se4) #combine rows 
quarterly.M6.tank.light.se4 <- rbind(c("05:45","NA"), quarterly.M6.tank.light.se4) #combine rows 
quarterly.M6.tank.light.se4[56,c(1,2)] <- c("19:30", "NA") #add empty row
quarterly.M6.tank.light.se4[57,c(1,2)] <- c("19:45", "NA") #add empty row
quarterly.M6.tank.light.se4[58,c(1,2)] <- c("20:00", "NA") #add empty row
quarterly.M6.tank.light.se4[59,c(1,2)] <- c("20:15", "NA") #add empty row
quarterly.M6.tank.light.se4$Tank4 <- as.numeric(quarterly.M6.tank.light.se4$Tank4) #set as numeric
quarterly.M6.tank.light.mean5 <- aggregate(Tank5 ~ Time, data=M6.tank.lights, mean, na.rm=TRUE) #calculate mean  for every 15 min interval
quarterly.M6.tank.light.mean5 <- rbind(c("05:45","NA"), quarterly.M6.tank.light.mean5) #combine rows 
quarterly.M6.tank.light.mean5[57,c(1,2)] <- c("19:45", "NA") #add empty row
quarterly.M6.tank.light.mean5[58,c(1,2)] <- c("20:00", "NA") #add empty row
quarterly.M6.tank.light.mean5[59,c(1,2)] <- c("20:15", "NA") #add empty row
quarterly.M6.tank.light.mean5$Tank5 <- as.numeric(quarterly.M6.tank.light.mean5$Tank5) #set as numeric
quarterly.M6.tank.light.se5 <- aggregate(Tank5 ~ Time, data=M6.tank.lights, std.error, na.rm=TRUE)  #calculate standard error of the mean for every 15 min interval
quarterly.M6.tank.light.se5 <- rbind(c("05:45","NA"), quarterly.M6.tank.light.se5) #combine rows 
quarterly.M6.tank.light.se5[57,c(1,2)] <- c("19:45", "NA") #add empty row
quarterly.M6.tank.light.se5[58,c(1,2)] <- c("20:00", "NA") #add empty row
quarterly.M6.tank.light.se5[59,c(1,2)] <- c("20:15", "NA") #add empty row
quarterly.M6.tank.light.se5$Tank5 <- as.numeric(quarterly.M6.tank.light.se5$Tank5) #set as numeric
M6.tank.light.means <-data.frame(quarterly.M6.tank.light.mean4,quarterly.M6.tank.light.se4$Tank4,quarterly.M6.tank.light.mean5$Tank5,quarterly.M6.tank.light.se5$Tank5) #make a dataframe
colnames(M6.tank.light.means) <- c("Time", "Tank4.mean", "Tank4.se", "Tank5.mean", "Tank5.se")  #Rename columns to describe contents
range(na.omit(M6.tank.light.means$Tank4.mean)) #range of average ambient light levels
range(na.omit(M6.tank.light.means$Tank5.mean)) #minimum of average high light levels

Fig15 <- ggplot(M6.tank.light.means, aes(Time)) + # plot mean light by tank
  geom_point(aes(y = Tank4.mean, colour="Ambient"), cex=0.8) + #plot points
  geom_errorbar(aes(x=Time, ymax=Tank4.mean+Tank4.se, ymin=Tank4.mean-Tank4.se), position=position_dodge(0.9), data=M6.tank.light.means, col="darkgrey", width=0) + #set values for standard error bars and offset on the X axis for clarity
  geom_point(aes(y = Tank5.mean, colour="High"), cex=0.8) + #plot points
  geom_errorbar(aes(x=Time, ymax=Tank5.mean+Tank5.se, ymin=Tank5.mean-Tank5.se), position=position_dodge(0.9), data=M6.tank.light.means, col="black", width=0) + #set values for standard error bars and offset on the X axis for clarity
  scale_colour_manual("Treatment", values = c("grey","black")) +
  scale_x_discrete(breaks=c("0:00", "06:00", "12:00", "18:00")) + #set discrete breaks on the X axis
  ylim(0,300) + #set Y axis limits
  ggtitle("(h)") + #Label graphic title
  xlab("Time") + #Label the X Axis
  ylab(bquote('Irradiance ('*mu~'mol' ~photons ~ m^-2~s^-1*')')) + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), angle = 90, vjust = 0.5, hjust=1), #set margins on labels
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        plot.title=element_text(hjust=0), #Justify the title to the top left
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        legend.position='none') #remove legend background
Fig15 #View figure

##### DISCRETE pH CALCULATIONS #####
path <-("~/MyProjects/HI_Pdam_Parental/RAnalysis/Data/pH_Calibration_Files/")
file.names<-list.files(path = path, pattern = "csv$") #list all the file names in the folder to get only get the csv files
pH.cals <- data.frame(matrix(NA, nrow=length(file.names), ncol=3, dimnames=list(file.names,c("Date", "Intercept", "Slope")))) #generate a 3 column dataframe with specific column names

for(i in 1:length(file.names)) { # for every file in list start at the first and run this following function
  Calib.Data <-read.table(file.path(path,file.names[i]), header=TRUE, sep=",", na.string="NA", as.is=TRUE) #reads in the data files
  model <-lm(mVTris ~ TTris, data=Calib.Data) #runs a linear regression of mV as a function of temperature
  coe <- coef(model) #extracts the coeffecients
  pH.cals[i,2:3] <- coe #inserts them in the dataframe
  pH.cals[i,1] <- substr(file.names[i],1,8) #stores the file name in the Date column
}
colnames(pH.cals) <- c("Calib.Date",  "Intercept",  "Slope") #rename columns
pH.cals #view data

#constants for use in pH calculation 
R <- 8.31447215 #gas constant in J mol-1 K-1 
F <-96485.339924 #Faraday constant in coulombs mol-1

#read in probe measurements of pH, temperature, and salinity from tanks
daily <- read.csv("Daily_Temp_pH_Sal.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA

#merge with Seawater chemistry file
SW.chem <- merge(pH.cals, daily, by="Calib.Date")

mvTris <- SW.chem$Temperature*SW.chem$Slope+SW.chem$Intercept #calculate the mV of the tris standard using the temperature mv relationships in the measured standard curves 
STris<-34.5 #salinity of the Tris
phTris<- (11911.08-18.2499*STris-0.039336*STris^2)*(1/(SW.chem$Temperature+273.15))-366.27059+ 0.53993607*STris+0.00016329*STris^2+(64.52243-0.084041*STris)*log(SW.chem$Temperature+273.15)-0.11149858*(SW.chem$Temperature+273.15) #calculate the pH of the tris (Dickson A. G., Sabine C. L. and Christian J. R., SOP 6a)
SW.chem$pH.Total<-phTris+(mvTris/1000-SW.chem$pH.MV/1000)/(R*(SW.chem$Temperature+273.15)*log(10)/F) #calculate the pH on the total scale (Dickson A. G., Sabine C. L. and Christian J. R., SOP 6a)

##### DISCRETE TA CALCULATIONS #####
setwd(file.path(mainDir, 'Data'))
massfile<-"TA_mass_data.csv" # name of your file with masses
path<-"~/MyProjects/HI_Pdam_Parental/RAnalysis/Data/TA" #the location of all your titration files
Sample.Info <- read.csv("TA_mass_data.csv", header=T, sep=",", na.string="NA", as.is=T) #load data
Mass<-read.csv(massfile, header=T, sep=",", na.string="NA", as.is=T, row.names=1)  #load Sample Info Data

# Select the mV for pH=3 and pH=3.5 based on probe calibration
pH35<-mean(Sample.Info$pH35, na.rm=T) #take the average mV reading for pH 3.5 across all samples
pH3<-mean(Sample.Info$pH3, na.rm=T) #take the average mV reading for pH 3.0 across all samples

#find all the titration data files
file.names<-list.files(path=path) #list all the file names in your data and sample directory
file.names <- file.names[grep("[.]csv", file.names)] # select only get the csv files

#create an empty dataframe to put the TA values in
nrow<-length(file.names) #set number of rows to the number of samples
TA <- matrix(nrow = nrow, ncol = 3) #set the dimensions of the dataframe
rownames(TA)<-file.names #identify row names
colnames(TA)<-c('Sample.ID','Mass','TA.Mes') #identify column names

setwd(file.path(mainDir, 'Data/TA')) # set working directory to where the data are
#run a for loop to bring in the titration files one at a time and calculate TA
for(i in 1: length(file.names)) {
  Data<-read.table(file.names[i], header=F, sep=",", na.string="NA",as.is=T) #read in each data file
  Data<-Data[-1:-6,] #remove the rows with characters
  
  # everything was brought in as a character because of the second line, converts back to numeric
  Data$Temperature<-as.numeric(Data[,7]) #convert to numeric and assign temperature column
  Data$Signal<-as.numeric(Data[,3]) #convert to numeric and assign mV column
  Data$Volume<-as.numeric(Data[,1]) #convert to numeric and assign volumn of titrant column
  
  #name of the file without .csv
  name<-unlist(strsplit(file.names[i], split='.', fixed=TRUE))[1]
  
  #identifies the indices of values between pH 3 and 3.5 
  mV<-which(Data$Signal<pH3 & Data$Signal>pH35) 
  
  #density of your titrant: specific to each bottle
  d1<-100*(-0.00000331*mean(Data$Temperature[mV], na.rm=T)^2-0.0001401*mean(Data$Temperature[mV], na.rm=T)+1.02933)/1000
  d2<-100*(-0.00000350*mean(Data$Temperature[mV], na.rm=T)^2-0.0001319*mean(Data$Temperature[mV], na.rm=T)+1.02907)/1000
  d3<-100*(-0.00000379*mean(Data$Temperature[mV], na.rm=T)^2-0.00012043*mean(Data$Temperature[mV], na.rm=T)+1.0296876)/1000
  
d <- if(Mass[name,4] =="d1") {
    d1                              #if density function = d1 use d1
  } else if(Mass[name,4] =="d2") {
    d2                              #if density function = d2 use d2
  } else if(Mass[name,4] =="d3")
    d3                              #if density function = d3 use d3
  
  #concentration of your titrant: specific to each bottle
  c<-Mass[name,3]

  #Salinity of your samples: changed with every sample
  s<-Mass[name,2]
  
  #mass of sample in g: changed with every sample
  mass<-Mass[name,1]

  #Calculate TA
  #at function is based on code in saecarb package by Steeve Comeau, Heloise Lavigne and Jean-Pierre Gattuso
  TA[i,1]<-name #add sample name to data output
  TA[i,2]<-mass #add mass to data output
  TA[i,3]<-10000000*at(S=s,T=mean(Data$Temperature[mV], na.rm=T), C=c, d=d, pHTris=NULL, ETris=NULL, weight=mass, E=Data$Signal[mV], volume=Data$Volume[mV]) #add TA to data output
}

TA <- data.frame(TA) #make a dataframe from the TA results

#load CRM standard Info
setwd(file.path(mainDir, 'Data')) # set working directory to where the data are
CRMs <- read.csv("CRM_TA_Data.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
Refs <- merge(CRMs, TA, by="Sample.ID") #merge the TA calculations with the Reference metadata
Refs$TA.Mes <- as.numeric(paste(Refs$TA.Mes)) #set valuse as numeric
Refs$Per.Off <- 100*((Refs$TA.Mes-Refs$CRM.TA)/Refs$CRM.TA) #calculate the percent difference of the TA from the CRM
Refs$TA.Corr <- Refs$CRM.TA-Refs$TA.Mes #correct TA for offset from CRM
Refs <- Refs[order(Refs$Date, abs(Refs$TA.Corr) ), ] #sort by id and reverse of abs(value)
Refs <- Refs[ !duplicated(Refs$Date), ]  # take the first row within each id
Refs #view data
CRM.res <- mean(Refs$Per.Off, na.rm=T) #calculate the average % difference of TA from CRM values over the course of the experiment
CRM.res #view the resolution of TA assay according to tests againsts Dickson CRMs

##### SEAWATER CHEMISTRY ANALYSIS FOR DISCRETE MEASUREMENTS#####
#Seawater chemistry table from simultaneous TA, pH, temperature and salinity measurements
#merge calculated pH and daily measures with TA data and run seacarb
SW.chem$Sample.ID <- paste(SW.chem$Date, SW.chem$Tank, sep='_') #generate new row with concatenated sample id
SW.chem <- merge(SW.chem,TA, by="Sample.ID", all = TRUE, sort = T) #merge seawater chemistry with total alkalinity
SW.chem <- na.omit(SW.chem) #remove NA
SW.chem <- merge(SW.chem, Refs[c("Date", "TA.Corr")], by="Date", all = F, sort = F) #merge seawater chemistry with total alkalinity
SW.chem$TA.Mes <- as.numeric(paste(SW.chem$TA.Mes)) #set as numeric
SW.chem$TA.Corr <- as.numeric(paste(SW.chem$TA.Corr)) #set as numeric
SW.chem$Corrected.TA <- SW.chem$TA.Mes - SW.chem$TA.Corr #correct for offset from CRM
SW.chem <- na.omit(SW.chem) #remove NA

#Calculate CO2 parameters using seacarb
carb.ouptput <- carb(flag=8, var1=SW.chem$pH.Total, var2=SW.chem$Corrected.TA/1000000, S= SW.chem$Salinity, T=SW.chem$Temperature, P=0, Pt=0, Sit=0, pHscale="T", kf="pf", k1k2="l", ks="d") #calculate seawater chemistry parameters using seacarb
carb.ouptput$ALK <- carb.ouptput$ALK*1000000 #convert to µmol kg-1
carb.ouptput$CO2 <- carb.ouptput$CO2*1000000 #convert to µmol kg-1
carb.ouptput$HCO3 <- carb.ouptput$HCO3*1000000 #convert to µmol kg-1
carb.ouptput$CO3 <- carb.ouptput$CO3*1000000 #convert to µmol kg-1
carb.ouptput$DIC <- carb.ouptput$DIC*1000000 #convert to µmol kg-1
carb.ouptput <- carb.ouptput[,-c(1,4,5,8,10:13,19)] #subset variables of interest
carb.ouptput <- cbind(SW.chem$Date,  SW.chem$Tank,  SW.chem$Treatment, SW.chem$Period1,SW.chem$Period2, SW.chem$Period3, carb.ouptput) #combine the sample information with the seacarb output
colnames(carb.ouptput) <- c("Date",  "Tank",  "Treatment",	"Period1", "Period2", "Period3",	"Salinity",	"Temperature", "pH",	"CO2",	"pCO2","HCO3",	"CO3",	"DIC", "TA",	"Aragonite.Sat") #Rename columns to describe contents
write.table(carb.ouptput, "~/MyProjects/HI_Pdam_Parental/RAnalysis/Output/Seawater_chemistry_table_Output_All.csv", sep=",", row.names = FALSE) #save data

carb.ouptput.Acc <- subset(carb.ouptput, Period1 == "Acc") #subset data
carb.ouptput.Adult <- subset(carb.ouptput, Period1 == "Adult") #subset data
carb.ouptput.Adult <- carb.ouptput.Adult[,-c(1,5:6)] #subset data
carb.ouptput.M1 <- subset(carb.ouptput, Period2 == "Month1") #subset data
carb.ouptput.M1 <- carb.ouptput.M1[,-c(1,4,6)] #subset data
carb.ouptput.M6 <- subset(carb.ouptput, Period3 == "Month6") #subset data
carb.ouptput.M6 <- carb.ouptput.M6[,-c(1,4,5)] #subset data

carbo.melted.Adult <- melt(carb.ouptput.Adult) #reshape the dataframe to more easily summarize all output parameters
mean.carb.output.Adult <-ddply(carbo.melted.Adult, .(Treatment, variable), summarize, #For each subset of a data frame, apply function then combine results into a data frame.
                         N = length(na.omit(value)), #number of records
                         mean = (mean(value)),       #take the average of the parameters (variables) summarized by treatments
                         sem = (sd(value)/sqrt(N))) #calculate the SEM as the sd/sqrt of the count or data length
mean.carb.output.Adult # display mean and sem 
mean.carb.output.Adult <- mean.carb.output.Adult[with(mean.carb.output.Adult, order(variable)), ] #order the data by the variables
adult.carb.table <- mean.carb.output.Adult[,-c(3)] #remove column
adult.carb.table <- reshape(adult.carb.table, direction="wide", timevar="variable", idvar="Treatment") #reshape data
adult.carb.table$N <- c(mean.carb.output.Adult[1,3],mean.carb.output.Adult[2,3]) #include sample size
#create an empty dataframe 
adult.chem.table <- matrix(nrow = 2, ncol = 1) #set the dimensions of the dataframe
colnames(adult.chem.table)<-c("Treatment") #identify column names
adult.chem.table <- data.frame(adult.chem.table) #change to dataframe
adult.chem.table$Treatment <- adult.carb.table$Treatment #add treatment info
adult.chem.table$N <- adult.carb.table$N #add sample size
adult.chem.table$Temperature <- paste(round(adult.carb.table$mean.Temperature, digits=2), round(adult.carb.table$sem.Temperature, digits=2), sep=' ± ') #add combined mean and sem with ± separating them
adult.chem.table$Salinity <- paste(round(adult.carb.table$mean.Salinity, digits=1), round(adult.carb.table$sem.Salinity, digits=1), sep=' ± ') #add combined mean and sem with ± separating them
adult.chem.table$Total.Alkalinity <- paste(round(adult.carb.table$mean.TA, digits=0), round(adult.carb.table$sem.TA, digits=0), sep=' ± ') #add combined mean and sem with ± separating them
adult.chem.table$pH <- paste(round(adult.carb.table$mean.pH, digits=2), round(adult.carb.table$sem.pH, digits=2), sep=' ± ') #add combined mean and sem with ± separating them
adult.chem.table$pCO2 <- paste(round(adult.carb.table$mean.pCO2, digits=0), round(adult.carb.table$sem.pCO2, digits=0), sep=' ± ') #add combined mean and sem with ± separating them
adult.chem.table$CO2 <- paste(round(adult.carb.table$mean.CO2, digits=0), round(adult.carb.table$sem.CO2, digits=0), sep=' ± ') #add combined mean and sem with ± separating them
adult.chem.table$HCO3 <- paste(round(adult.carb.table$mean.HCO3, digits=0), round(adult.carb.table$sem.HCO3, digits=0), sep=' ± ') #add combined mean and sem with ± separating them
adult.chem.table$CO3 <- paste(round(adult.carb.table$mean.CO3, digits=0), round(adult.carb.table$sem.CO3, digits=0), sep=' ± ') #add combined mean and sem with ± separating them
adult.chem.table$DIC <- paste(round(adult.carb.table$mean.DIC, digits=0), round(adult.carb.table$sem.DIC, digits=0), sep=' ± ') #add combined mean and sem with ± separating them
adult.chem.table$Arag.Sat <- paste(round(adult.carb.table$mean.Aragonite.Sat, digits=1), round(adult.carb.table$sem.Aragonite.Sat, digits=1), sep=' ± ') #add combined mean and sem with ± separating them

carbo.melted.M1 <- melt(carb.ouptput.M1) #reshape the dataframe to more easily summarize all output parameters
mean.carb.output.M1 <-ddply(carbo.melted.M1, .(Treatment, variable), summarize, #For each subset of a data frame, apply function then combine results into a data frame.
                               N = length(na.omit(value)), #number of records
                               mean = (mean(value)),       #take the average of the parameters (variables) summarized by treatments
                               sem = (sd(value)/sqrt(N))) #calculate the SEM as the sd/sqrt of the count or data length
mean.carb.output.M1 # display mean and sem 
mean.carb.output.M1 <- mean.carb.output.M1[with(mean.carb.output.M1, order(variable)), ] #order the data by the variables
M1.carb.table <- mean.carb.output.M1[,-c(3)] #remove sample size
M1.carb.table <- reshape(M1.carb.table, direction="wide", timevar="variable", idvar="Treatment")
M1.carb.table$N <- c(mean.carb.output.M1[1,3],mean.carb.output.M1[2,3]) #include sample size
#create an empty dataframe 
M1.chem.table <- matrix(nrow = 2, ncol = 1) #set the dimensions of the dataframe
colnames(M1.chem.table)<-c("Treatment") #identify column names
M1.chem.table <- data.frame(M1.chem.table) #change to dataframe
M1.chem.table$Treatment <- M1.carb.table$Treatment #add treatment info
M1.chem.table$N <- M1.carb.table$N #add sample size
M1.chem.table$Temperature <- paste(round(M1.carb.table$mean.Temperature, digits=2), round(M1.carb.table$sem.Temperature, digits=2), sep=' ± ') #add combined mean and sem with ± separating them
M1.chem.table$Salinity <- paste(round(M1.carb.table$mean.Salinity, digits=1), round(M1.carb.table$sem.Salinity, digits=1), sep=' ± ') #add combined mean and sem with ± separating them
M1.chem.table$Total.Alkalinity <- paste(round(M1.carb.table$mean.TA, digits=0), round(M1.carb.table$sem.TA, digits=0), sep=' ± ') #add combined mean and sem with ± separating them
M1.chem.table$pH <- paste(round(M1.carb.table$mean.pH, digits=2), round(M1.carb.table$sem.pH, digits=2), sep=' ± ') #add combined mean and sem with ± separating them
M1.chem.table$pCO2 <- paste(round(M1.carb.table$mean.pCO2, digits=0), round(M1.carb.table$sem.pCO2, digits=0), sep=' ± ') #add combined mean and sem with ± separating them
M1.chem.table$CO2 <- paste(round(M1.carb.table$mean.CO2, digits=0), round(M1.carb.table$sem.CO2, digits=0), sep=' ± ') #add combined mean and sem with ± separating them
M1.chem.table$HCO3 <- paste(round(M1.carb.table$mean.HCO3, digits=0), round(M1.carb.table$sem.HCO3, digits=0), sep=' ± ') #add combined mean and sem with ± separating them
M1.chem.table$CO3 <- paste(round(M1.carb.table$mean.CO3, digits=0), round(M1.carb.table$sem.CO3, digits=0), sep=' ± ') #add combined mean and sem with ± separating them
M1.chem.table$DIC <- paste(round(M1.carb.table$mean.DIC, digits=0), round(M1.carb.table$sem.DIC, digits=0), sep=' ± ') #add combined mean and sem with ± separating them
M1.chem.table$Arag.Sat <- paste(round(M1.carb.table$mean.Aragonite.Sat, digits=1), round(M1.carb.table$sem.Aragonite.Sat, digits=1), sep=' ± ') #add combined mean and sem with ± separating them

carbo.melted.M6 <- melt(carb.ouptput.M6) #reshape the dataframe to more easily summarize all output parameters
mean.carb.output.M6 <-ddply(carbo.melted.M6, .(Treatment, variable), summarize, #For each subset of a data frame, apply function then combine results into a data frame.
                            N = length(na.omit(value)), #number of records
                            mean = (mean(value)),       #take the average of the parameters (variables) summarized by treatments
                            sem = (sd(value)/sqrt(N))) #calculate the SEM as the sd/sqrt of the count or data length
mean.carb.output.M6 # display mean and sem 
mean.carb.output.M6 <- mean.carb.output.M6[with(mean.carb.output.M6, order(variable)), ] #order the data by the variables
M6.carb.table <- mean.carb.output.M6[,-c(3)] #remove sample size
M6.carb.table <- reshape(M6.carb.table, direction="wide", timevar="variable", idvar="Treatment")
M6.carb.table$N <- c(mean.carb.output.M6[1,3],mean.carb.output.M6[2,3]) #include sample size
#create an empty dataframe 
M6.chem.table <- matrix(nrow = 2, ncol = 1) #set the dimensions of the dataframe
colnames(M6.chem.table)<-c("Treatment") #identify column names
M6.chem.table <- data.frame(M6.chem.table) #change to dataframe
M6.chem.table$Treatment <- M6.carb.table$Treatment #add treatment info
M6.chem.table$N <- M6.carb.table$N #add sample size
M6.chem.table$Temperature <- paste(round(M6.carb.table$mean.Temperature, digits=2), round(M6.carb.table$sem.Temperature, digits=2), sep=' ± ') #add combined mean and sem with ± separating them
M6.chem.table$Salinity <- paste(round(M6.carb.table$mean.Salinity, digits=1), round(M6.carb.table$sem.Salinity, digits=1), sep=' ± ') #add combined mean and sem with ± separating them
M6.chem.table$Total.Alkalinity <- paste(round(M6.carb.table$mean.TA, digits=0), round(M6.carb.table$sem.TA, digits=0), sep=' ± ') #add combined mean and sem with ± separating them
M6.chem.table$pH <- paste(round(M6.carb.table$mean.pH, digits=2), round(M6.carb.table$sem.pH, digits=2), sep=' ± ') #add combined mean and sem with ± separating them
M6.chem.table$pCO2 <- paste(round(M6.carb.table$mean.pCO2, digits=0), round(M6.carb.table$sem.pCO2, digits=0), sep=' ± ') #add combined mean and sem with ± separating them
M6.chem.table$CO2 <- paste(round(M6.carb.table$mean.CO2, digits=0), round(M6.carb.table$sem.CO2, digits=0), sep=' ± ') #add combined mean and sem with ± separating them
M6.chem.table$HCO3 <- paste(round(M6.carb.table$mean.HCO3, digits=0), round(M6.carb.table$sem.HCO3, digits=0), sep=' ± ') #add combined mean and sem with ± separating them
M6.chem.table$CO3 <- paste(round(M6.carb.table$mean.CO3, digits=0), round(M6.carb.table$sem.CO3, digits=0), sep=' ± ') #add combined mean and sem with ± separating them
M6.chem.table$DIC <- paste(round(M6.carb.table$mean.DIC, digits=0), round(M6.carb.table$sem.DIC, digits=0), sep=' ± ') #add combined mean and sem with ± separating them
M6.chem.table$Arag.Sat <- paste(round(M6.carb.table$mean.Aragonite.Sat, digits=1), round(M6.carb.table$sem.Aragonite.Sat, digits=1), sep=' ± ') #add combined mean and sem with ± separating them

##### CONTINUOUS pH MEASUREMENTS#####
#Tank NBS pH Data for Adult Exposure Experimental Period (06May14 - 17Aug14)
# read in NBS pH data from Aquacontrollers, frequency 15min
setwd(file.path(mainDir, 'Data')) #set working directory
pHs <- read.csv("Adult_Tank_NBS_pH.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
mydate.pHs <- strptime(pHs$Date.Time, format="%m/%d/%y %H:%M") #Identify date format
pHs$Tank5[is.na(pHs$Tank5)] <- pHs$Tank3[is.na(pHs$Tank5)] #merge tanks 3 and 5 into one column so the data from first 3 days when corals were in tank 3 is now showing in tank 5
tank.pHdata <-data.frame(mydate.pHs, pHs$Tank4, pHs$Tank5) #make a dataframe of temperature and time
colnames(tank.pHdata) <- c("Date.Time", "Ambient.NBS", "High.NBS") #Rename columns to describe contents
pH.data <- merge(tank.pHdata, tank.tempdata, by="Date.Time") #merge data sets by time and date
pH.data$Time <- format(pH.data$Date.Time, format = "%H:%M:%S") #separate date and time

#Plot total pH for both treatments for duration of M1
Fig16 <- ggplot(pH.data) + #plot pH total scale
  geom_line(aes(x = Date.Time, y = High.NBS, col="High")) + #plot as a line
  geom_line(aes(x = Date.Time, y = Ambient.NBS, col="Ambient")) + #plot as a line
  xlab("Date") + #Label the X Axis
  ylab("pH (NBS Scale)") + #Label the Y Axis
  ggtitle("Adult NBS pH") + #Label the graph title
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background =element_blank(), #Set the plot background
        legend.key = element_blank()) #Set plot legend key
Fig16 #View figure

quarterly.tank.pH.amb.mean <- aggregate(Ambient.NBS ~ Time, data=pH.data, mean, na.rm=TRUE) #calculate mean of pH for every 15 min interval
quarterly.tank.pH.amb.se <- aggregate(Ambient.NBS ~ Time, data=pH.data, std.error, na.rm=TRUE)  #calculate standard error of the mean of pH for every 15 min interval
quarterly.tank.pH.high.mean <- aggregate(High.NBS ~ Time, data=pH.data, mean, na.rm=TRUE) #calculate mean of pH for every 15 min interval
quarterly.tank.pH.high.se <- aggregate(High.NBS ~ Time, data=pH.data, std.error, na.rm=TRUE)  #calculate standard error of the mean of pH for every 15 min interval
adult.pH.amb.rng <- range(quarterly.tank.pH.amb.mean$Ambient.NBS) #range average ambient NBS pH
adult.pH.high.rng <-range(quarterly.tank.pH.high.mean$High.NBS) #range average high NBS pH
tank.pH.means <- data.frame(quarterly.tank.pH.amb.mean, quarterly.tank.pH.amb.se$Ambient.NBS, quarterly.tank.pH.high.mean$High.NBS, quarterly.tank.pH.high.se$High.NBS) #combine mean and standard error results
colnames(tank.pH.means) <- c("Time", "pH.Amb.NBS.mean", "pH.Amb.NBS.se", "pH.High.NBS.mean", "pH.High.NBS.se")  #Rename columns to describe contents
quarts <- as.data.frame(seq(ISOdatetime(2001,2,3,0,0,0), ISOdatetime(2001,2,4,0,0,0), by=(60*15))) #add sequence of time by 15 min
tank.pH.means$Time <- quarts[1:96,] #remove last line

Fig17 <- ggplot(tank.pH.means, aes(Time)) + # plot mean pH by tank
  geom_point(aes(x =Time, y = pH.Amb.NBS.mean, colour="Ambient"), cex=0.8) + #plot points
  geom_errorbar(aes(x=Time, ymax=pH.Amb.NBS.mean+pH.Amb.NBS.se, ymin=pH.Amb.NBS.mean-pH.Amb.NBS.se), position=position_dodge(0.9), data=tank.pH.means, col="darkgrey", width=0) + #set values for standard error bars and offset on the X axis for clarity
  geom_point(aes(x = Time, y = pH.High.NBS.mean, colour="High"), cex=0.8) + #plot points
  geom_errorbar(aes(x=Time, ymax=pH.High.NBS.mean+pH.High.NBS.se, ymin=pH.High.NBS.mean-pH.High.NBS.se), position=position_dodge(0.9), data=tank.pH.means, col="black", width=0) + #set values for standard error bars and offset on the X axis for clarity
  scale_colour_manual("Treatment", values = c("grey","black")) +
  scale_x_datetime(date_labels="%H %M") + #label hours and min
  scale_y_continuous(name="pH (NBS Scale)", breaks=c( 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8.0, 8.1), limits=c(7.35, 8.1)) + #set Y axis ticks
  ggtitle("(j)") + #Label graph
  xlab("Time") + #Label the X Axis
  ylab("pH (NBS Scale)") + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), angle = 90, vjust = 0.5, hjust=1), #set margins on labels
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        plot.title=element_text(hjust=0), #Justify the title to the top left
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        legend.position='none') #remove legend background
Fig17 #View figure


#Tank pH Data for Larval Month 1 Exposure Period (12July - 25Aug14)
# read in NBS pH data from Aquacontrollers frequency 15min
pHs.M1 <- read.csv("Month1_Tank_NBS_pH.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
mydate.pHs.M1 <- strptime(pHs.M1$Date.Time, format="%m/%d/%y %H:%M") #Identify date format
tank.pHdata.M1 <-data.frame(mydate.pHs.M1, pHs.M1$Tank4, pHs.M1$Tank5) #make a dataframe of temperature and time
colnames(tank.pHdata.M1) <- c("Date.Time", "Ambient.NBS", "High.NBS") #Rename columns to describe contents
pH.data.M1 <- merge(tank.pHdata.M1, M1.tank.tempdata, by="Date.Time") #merge data sets by time and date
pH.data.M1$Time <- format(pH.data.M1$Date.Time, format = "%H:%M:%S") #separate date and time

#Plot pH for both treatments for duration of M1
Fig18 <- ggplot(pH.data.M1) + #plot pH total scale
  geom_line(aes(x = Date.Time, y = High.NBS, col="High")) + #plot as a line
  geom_line(aes(x = Date.Time, y = Ambient.NBS, col="Ambient")) + #plot as a line
  xlab("Date") + #Label the X Axis
  ylab("pH (NBS Scale)") + #Label the Y Axis
  ggtitle("Month1 NBS pH") + #Label the graph title
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background =element_blank(), #Set the plot background
        legend.key = element_blank()) #Set plot legend key
Fig18 #View figure

quarterly.tank.pH.amb.mean <- aggregate(Ambient.NBS ~ Time, data=pH.data.M1, mean, na.rm=TRUE) #calculate mean of pH for every 15 min interval
quarterly.tank.pH.amb.se <- aggregate(Ambient.NBS ~ Time, data=pH.data.M1, std.error, na.rm=TRUE)  #calculate standard error of the mean of pH for every 15 min interval
quarterly.tank.pH.high.mean <- aggregate(High.NBS ~ Time, data=pH.data.M1, mean, na.rm=TRUE) #calculate mean of pH for every 15 min interval
quarterly.tank.pH.high.se <- aggregate(High.NBS ~ Time, data=pH.data.M1, std.error, na.rm=TRUE)  #calculate standard error of the mean of pH for every 15 min interval
M1.pH.amb.rng <- range(quarterly.tank.pH.amb.mean$Ambient.NBS) #range average ambient NBS pH
M1.pH.high.rng <-range(quarterly.tank.pH.high.mean$High.NBS) #range average high NBS pH
tank.pH.means.M1 <- data.frame(quarterly.tank.pH.amb.mean, quarterly.tank.pH.amb.se$Ambient.NBS, quarterly.tank.pH.high.mean$High.NBS, quarterly.tank.pH.high.se$High.NBS) #combine mean and standard error results
colnames(tank.pH.means.M1) <- c("Time", "pH.Amb.NBS.mean", "pH.Amb.NBS.se", "pH.High.NBS.mean", "pH.High.NBS.se")  #Rename columns to describe contents
tank.pH.means.M1$Time <- quarts[1:96,]

Fig19 <- ggplot(tank.pH.means.M1, aes(Time)) + # plot mean pH by tank
  geom_point(aes(x =Time, y = pH.Amb.NBS.mean, colour="Ambient"), cex=0.8) + #plot points
  geom_errorbar(aes(x=Time, ymax=pH.Amb.NBS.mean+pH.Amb.NBS.se, ymin=pH.Amb.NBS.mean-pH.Amb.NBS.se), position=position_dodge(0.9), data=tank.pH.means.M1, col="darkgrey", width=0) + #set values for standard error bars and offset on the X axis for clarity
  geom_point(aes(x = Time, y = pH.High.NBS.mean, colour="High"), cex=0.8) + #plot points
  geom_errorbar(aes(x=Time, ymax=pH.High.NBS.mean+pH.High.NBS.se, ymin=pH.High.NBS.mean-pH.High.NBS.se), position=position_dodge(0.9), data=tank.pH.means.M1, col="black", width=0) + #set values for standard error bars and offset on the X axis for clarity
  scale_colour_manual("Treatment", values = c("grey","black")) +
  scale_x_datetime(date_labels="%H %M") + #label hours and min
  scale_y_continuous(name="pH (NBS Scale)", breaks=c( 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8.0, 8.1), limits=c(7.35, 8.1)) + #set Y axis ticks
  ggtitle("(k)") + #Label graph
  xlab("Time") + #Label the X Axis
  ylab("pH (NBS Scale)") + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), angle = 90, vjust = 0.5, hjust=1), #set margins on labels
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        plot.title=element_text(hjust=0), #Justify the title to the top left
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        legend.position='none') #remove legend background
Fig19 #View figure


#Tank pH Data for Larval Month 6 Exposure Period (12July - 28January15)
# read in NBS pH data from Aquacontrollers frequency 15min
pHs.M6 <- read.csv("Month6_Tank_NBS_pH.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
mydate.pHs.M6 <- strptime(pHs.M6$Date.Time, format="%m/%d/%y %H:%M") #Identify date format
tank.pHdata.M6 <-data.frame(mydate.pHs.M6, pHs.M6$Tank4, pHs.M6$Tank5) #make a dataframe of temperature and time
colnames(tank.pHdata.M6) <- c("Date.Time", "Ambient.NBS", "High.NBS") #Rename columns to describe contents
pH.data.M6 <- merge(tank.pHdata.M6, M6.tank.tempdata, by="Date.Time") #merge data sets by time and date
pH.data.M6$Time <- format(pH.data.M6$Date.Time, format = "%H:%M:%S") #separate date and time

#Plot pH for both treatments for duration of m6
Fig20 <- ggplot(pH.data.M6) + #plot pH total scale
  geom_line(aes(x = Date.Time, y = High.NBS, col="High")) + #plot as a line
  geom_line(aes(x = Date.Time, y = Ambient.NBS, col="Ambient")) + #plot as a line
  xlab("Date") + #Label the X Axis
  ylab("pH (NBS Scale)") + #Label the Y Axis
  ggtitle("Month 6 Total pH") + #Label the graph title
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background =element_blank(), #Set the plot background
        legend.key = element_blank()) #Set plot legend key
Fig20 #View figure

quarterly.tank.pH.amb.mean <- aggregate(Ambient.NBS ~ Time, data=pH.data.M6, mean, na.rm=TRUE) #calculate mean of pH for every 15 min interval
quarterly.tank.pH.amb.se <- aggregate(Ambient.NBS ~ Time, data=pH.data.M6, std.error, na.rm=TRUE)  #calculate standard error of the mean of pH for every 15 min interval
quarterly.tank.pH.high.mean <- aggregate(High.NBS ~ Time, data=pH.data.M6, mean, na.rm=TRUE) #calculate mean of pH for every 15 min interval
quarterly.tank.pH.high.se <- aggregate(High.NBS ~ Time, data=pH.data.M6, std.error, na.rm=TRUE)  #calculate standard error of the mean of pH for every 15 min interval
M6.pH.amb.rng <- range(quarterly.tank.pH.amb.mean$Ambient.NBS) #range average ambient NBS pH
M6.pH.high.rng <-range(quarterly.tank.pH.high.mean$High.NBS) #range average high NBS pH
tank.pH.means.M6 <- data.frame(quarterly.tank.pH.amb.mean, quarterly.tank.pH.amb.se$Ambient.NBS, quarterly.tank.pH.high.mean$High.NBS, quarterly.tank.pH.high.se$High.NBS) #combine mean and standard error results
colnames(tank.pH.means.M6) <- c("Time", "pH.Amb.NBS.mean", "pH.Amb.NBS.se", "pH.High.NBS.mean", "pH.High.NBS.se")  #Rename columns to describe contents
tank.pH.means.M6$Time <- quarts[1:96,]

Fig21 <- ggplot(tank.pH.means.M6, aes(Time)) + # plot mean pH by tank
  geom_point(aes(x =Time, y = pH.Amb.NBS.mean, colour="Ambient"), cex=0.8) + #plot points
  geom_errorbar(aes(x=Time, ymax=pH.Amb.NBS.mean+pH.Amb.NBS.se, ymin=pH.Amb.NBS.mean-pH.Amb.NBS.se), position=position_dodge(0.9), data=tank.pH.means.M6, col="darkgrey", width=0) + #set values for standard error bars and offset on the X axis for clarity
  geom_point(aes(x = Time, y = pH.High.NBS.mean, colour="High"), cex=0.8) + #plot points
  geom_errorbar(aes(x=Time, ymax=pH.High.NBS.mean+pH.High.NBS.se, ymin=pH.High.NBS.mean-pH.High.NBS.se), position=position_dodge(0.9), data=tank.pH.means.M6, col="black", width=0) + #set values for standard error bars and offset on the X axis for clarity
  scale_colour_manual("Treatment", values = c("grey","black")) +
  scale_x_datetime(date_labels="%H %M") + #label hours and min on the X axis
  scale_y_continuous(name="pH (NBS Scale)", breaks=c( 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8.0, 8.1), limits=c(7.35, 8.1)) + #set Y axis ticks 
  ggtitle("(l)") + #Label graph
  xlab("Time") + #Label the X Axis
  ylab("pH (NBS Scale)") + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), angle = 90, vjust = 0.5, hjust=1), #set margins on labels
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        plot.title=element_text(hjust=0), #Justify the title to the top left
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        legend.position='none') #remove legend background
Fig21 #View figure

FigMT <- ggplot(tank.pH.means.M6, aes(Time)) + # plot mean pH by tank
  geom_point(aes(x =Time, y = pH.Amb.NBS.mean, colour="Ambient")) + #plot points
  geom_point(aes(x = Time, y = pH.High.NBS.mean, colour="High")) + #plot points
  scale_colour_manual("Treatment", values = c("white","white")) +
  scale_x_datetime(date_labels="%H %M") + #label hours and min on the X axis
  scale_y_continuous(name="pH (NBS Scale)", breaks=c( 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8.0, 8.1), limits=c(7.35, 8.1)) + #set Y axis ticks 
  ggtitle("(i) ") + #Label graph
  xlab("Time") + #Label the X Axis
  ylab("pH (NBS Scale)") + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), angle = 90, vjust = 0.5, hjust=1), #set margins on labels
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        plot.title=element_text(hjust=0), #Justify the title to the top left
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        legend.position='none') #remove legend background
FigMT

##### BIOLOGICAL RESPONSES #####
##### LARVAL RELEASE #####
#June
june.release.data <- read.csv("june.larval.release.data.csv", header=T, sep=",", na.string="NA", as.is=T) #load data
june.release.data <- na.omit(june.release.data) #remove NA
june.mean_larvae <- aggregate(numb.larvae ~ Lunar.Day * Treatment, data=june.release.data, FUN=mean) #calculate mean of Day * treatment
june.se_larvae <- aggregate(numb.larvae ~ Lunar.Day * Treatment, data=june.release.data, FUN=std.error) #calculate se of Day * treatment
june.larvae <- cbind(june.mean_larvae,june.se_larvae$numb.larvae) #make dataframe
colnames(june.larvae) <- c("Lunar.Day", "Treatment", "mean", "se") #rename columns
june.larvae #view data

Fig22 <- ggplot(june.larvae, aes(x=Lunar.Day, y=mean, fill=Treatment)) + #plot mean as a function of day
  geom_bar(position=position_dodge(), stat="identity") + #assign bar id and position
  scale_fill_manual(values=c("gray", "black")) + #bar fill color
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), #plot error bars
                width=0, size = 0.4,                   # Width of the error bars
                position=position_dodge(.9)) + #set bar position
  ylim(0,400) + #set y limits
  ggtitle("(a) June") + #plot title
  xlab("Lunar Day") + #x axis title
  ylab("Number of planulae released") + #y axis title
  theme_bw() + #theme black and white 
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5))+ #legend guides
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), angle = 90, vjust = 0.5, hjust=1), #set margins on labels
        axis.text=element_text(size=16), #set text size
        axis.title=element_text(size=18,face="bold"), #set axis title text size
        strip.text.x = element_text(size = 16, colour = "black", face="bold"),
        axis.line.x = element_line(color = 'black'), #Set the axes color
        axis.line.y = element_line(color = 'black'), #Set the axes color
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_blank(),  #remove legend background
        legend.position=c(0.2,0.85),  #set legend location
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0)) #set title attributes
Fig22 #view plot

june.amb <- subset(june.larvae, Treatment=="Ambient") #subset data
june.high <- subset(june.larvae, Treatment=="High") #subset data
june.ks <-ks.test(june.amb$mean, june.high$mean) #Kolmogorov-Smirnov Test Ho: differences in distribution
june.ks #view results

#July
july.release.data <- read.csv("july.larval.release.data.csv", header=T, sep=",", na.string="NA", as.is=T) #load data
july.release.data <- na.omit(july.release.data) #remove NA
july.mean_larvae <- aggregate(numb.larvae ~ Lunar.Day * Treatment, data=july.release.data, FUN=mean) #calculate mean of Day * treatment
july.se_larvae <- aggregate(numb.larvae ~ Lunar.Day * Treatment, data=july.release.data, FUN=std.error) #calculate se of Day * treatment
july.larvae <- cbind(july.mean_larvae,july.se_larvae$numb.larvae) #make dataframe
colnames(july.larvae) <- c("Lunar.Day", "Treatment", "mean", "se") #rename columns
july.larvae #view data

Fig23 <- ggplot(july.larvae, aes(x=Lunar.Day, y=mean, fill=Treatment)) + #plot mean as a function of day
  geom_bar(position=position_dodge(), stat="identity", show.legend=FALSE) + #assign bar id and position
  scale_fill_manual(values=c("gray", "black")) + #bar fill color
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), #plot error bars
                width=0, size = 0.4,                   # Width of the error bars
                position=position_dodge(.9)) + #set bar position
  ylim(0,400) + #set y limits
  ggtitle("(b) July") + #plot title
  xlab("Lunar Day") + #x axis title
 # ylab("Number of planulae released") + #y axis title
  theme_bw() + #theme black and white 
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5))+ #legend guides
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), angle = 90, vjust = 0.5, hjust=1), #set margins on labels
        axis.text=element_text(size=16), #set text size
        axis.title=element_text(size=18,face="bold"), #set axis title text size
        axis.title.y=element_blank(), #remove Y axis label
        strip.text.x = element_text(size = 16, colour = "black", face="bold"),
        axis.line.x = element_line(color = 'black'), #Set the axes color
        axis.line.y = element_line(color = 'black'), #Set the axes color
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_blank(),  #remove legend background
        legend.position="none",  #set legend location
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        plot.title = element_text(face = 'bold', 
                                size = 12, 
                                hjust = 0)) #set title attributes
Fig23

july.amb <- subset(july.larvae, Treatment=="Ambient") #subset data
july.high <- subset(july.larvae, Treatment=="High") #subset data
july.ks <-ks.test(july.amb$mean, july.high$mean) #Kolmogorov-Smirnov Test Ho: differences in distribution
july.ks #view results

#august
august.release.data <- read.csv("august.larval.release.data.csv", header=T, sep=",", na.string="NA", as.is=T) #load data
august.release.data <- na.omit(august.release.data) #remove NA
august.mean_larvae <- aggregate(numb.larvae ~ Lunar.Day * Treatment, data=august.release.data, FUN=mean) #calculate mean of Day * treatment
august.se_larvae <- aggregate(numb.larvae ~ Lunar.Day * Treatment, data=august.release.data, FUN=std.error) #calculate se of Day * treatment
august.larvae <- cbind(august.mean_larvae,august.se_larvae$numb.larvae) #make dataframe
colnames(august.larvae) <- c("Lunar.Day", "Treatment", "mean", "se") #rename columns
august.larvae #view data

Fig24 <- ggplot(august.larvae, aes(x=Lunar.Day, y=mean, fill=Treatment)) + #plot mean as a function of day
  geom_bar(position=position_dodge(), stat="identity", show.legend=FALSE) + #assign bar id and position
  scale_fill_manual(values=c("gray", "black")) + #bar fill color
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), #plot error bars
                width=0, size = 0.4,                   # Width of the error bars
                position=position_dodge(.9)) + #set bar position
  ylim(0,400) + #set y limits
  ggtitle("(c) August") + #plot title
  xlab("Lunar Day") + #x axis title
  ylab("Number of planulae released") + #y axis title
  theme_bw() + #theme black and white 
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5))+ #legend guides
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), angle = 90, vjust = 0.5, hjust=1), #set margins on labels
        axis.text=element_text(size=16), #set text size
        axis.title=element_text(size=18,face="bold"), #set axis title text size
        axis.title.y=element_blank(), #remove Y axis label
        strip.text.x = element_text(size = 16, colour = "black", face="bold"),
        axis.line.x = element_line(color = 'black'), #Set the axes color
        axis.line.y = element_line(color = 'black'), #Set the axes color
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_blank(),  #remove legend background
        legend.position="none",  #set legend location
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0)) #set title attributes
Fig24

august.amb <- subset(august.larvae, Treatment=="Ambient") #subset data
august.high <- subset(august.larvae, Treatment=="High") #subset data
august.ks <-ks.test(august.amb$mean, august.high$mean) #Kolmogorov-Smirnov Test Ho: differences in distribution
august.ks #view results

# ## Total release as a function of both treatment and time
# RM.release.data <- read.csv("RM_Release_Data.csv", header=T, sep=",", na.string="NA", as.is=T) #read in data in long format
# all.release.mean <- aggregate(Total.Release ~ Treatment + Time, data=RM.release.data, mean) #calculate mean by treatment and time
# all.release.se <- aggregate(Total.Release ~ Treatment + Time, data=RM.release.data, std.error)  #calculate se by treatment and time
# all.release <- cbind(all.release.mean, all.release.se$Total.Release) #make dataframe
# colnames(all.release) <- c("Treatment", "Time", "mean", "se") #rename columns
# 
# Fig25 <- ggplot(all.release, aes(x=Time, y=mean, colour=Treatment, group=Treatment), position=position_dodge(width=0.5)) +  #plot mean as a function of Time
#   geom_errorbar(aes(ymin=all.release$mean - all.release$se, ymax=all.release$mean + all.release$se), #plot error bars
#                 colour="black", width=0, size = 0.4, # Width of the error bars
#                 position=position_dodge(width=0.5)) + #set bar position
#   geom_point(position=position_dodge(width=0.5), size=2, shape=15) +
#   scale_colour_manual(values = c("gray","black")) + #set point fill color
#   scale_x_discrete(limits=c("June","July","August")) + #label x axis in order
#   ylab(" Total Release") + #y axis label
#   ylim(0,1700) + #y axis limits
#   ggtitle("D) Total") + #plot title
#   theme_bw() + #theme black and white 
#   theme(axis.line = element_line(color = 'black'), #Set the axes color
#         axis.text=element_text(size=10), #set text size
#         axis.title=element_text(size=12,face="bold"), #set axis title text size
#         strip.text.x = element_text(size = 12, colour = "black", face="bold"),
#         panel.border = element_blank(), #Set the border
#         axis.line.x = element_line(color = 'black'), #Set the axes color
#         axis.line.y = element_line(color = 'black'), #Set the axes color
#         axis.text.x=element_text(angle=90), #Set text angle
#         panel.grid.major = element_blank(), #Set the major gridlines
#         panel.grid.minor = element_blank(), #Set the minor gridlines
#         plot.background=element_blank(),  #Set the plot background
#         legend.key = element_blank(),  #remove legend background
#         legend.position="none",  #set legend location
#         plot.title = element_text(face = 'bold', 
#                                   size = 12, 
#                                   hjust = 0)) #set title attributes
# Fig25
# 
# #LM release by treatment and time
# RM.release.data #view data
# RM.release.data <- na.omit(RM.release.data) #remove NA
# release.lm <- lm(log10(Total.Release+1) ~ Treatment * Time, data=RM.release.data) #run generalized linear model 
# summary(release.lm) #view summary
# release.resid <-resid(release.lm)
# release.shapiro <- shapiro.test(release.resid) #runs a normality test using shapiro-wilk test on the residuals
# release.shapiro #view results
# release.qqnorm <- qqnorm(release.resid) # normal quantile plot
# release.qqline <- qqline(release.resid) # adding a qline of comparison
# hist(release.resid) #plot histogram of residuals
# plot(release.lm$fitted.values, release.lm$residuals) #plot residuals as a function of fitted data
# 
# release.posthoc <- lsmeans(release.lm, specs=c("Time")) #calculate MS means
# release.posthoc #view results
# release.posthoc.p <- contrast(release.posthoc, method="pairwise") #contrast treatment groups within a species at each time point
# release.posthoc.p #view results
# release.posthoc.lett <- cld(release.posthoc , alpha=.05, Letters=letters) #identify posthoc letter differences
# release.posthoc.lett #view results

##### SURVIVORSHIP #####
larval.data.M0 <- read.csv("Larval_Data_M0.csv", header=T, sep=",", na.string="NA", as.is=T) #load data
proportion.alive.M0 <- (larval.data.M0$Plastic + larval.data.M0$Top.Tile + larval.data.M0$Bottom.Tile +  larval.data.M0$Edge +	larval.data.M0$Swimming)/larval.data.M0$larvae.added #calculate survivorship
proportion.dead.M0 <- 1-proportion.alive.M0 #calculate mortality
larval.data.M0$Alive <- (larval.data.M0$Plastic + larval.data.M0$Top.Tile + larval.data.M0$Bottom.Tile +  larval.data.M0$Edge +	larval.data.M0$Swimming) #count alive
larval.data.M0$Dead <- larval.data.M0$larvae.added-(larval.data.M0$Plastic + larval.data.M0$Top.Tile + larval.data.M0$Bottom.Tile +  larval.data.M0$Edge +	larval.data.M0$Swimming) #calculate dead
survive.M0 <- data.frame (larval.data.M0$Chamber.num, larval.data.M0$Timepoint, larval.data.M0$Origin, larval.data.M0$Secondary, proportion.alive.M0, proportion.dead.M0, larval.data.M0$Alive, larval.data.M0$Dead) #make dataframe
colnames(survive.M0) <- c("Chamber", "Timepoint", "Origin", "Secondary", "Prop.Alive","Prop.Dead", "Alive","Dead") #rename columns
mean.survive.M0 <- aggregate(Prop.Alive ~ Origin * Secondary, data = survive.M0, FUN= "mean") #calculate mean by origin and secondary treatments
se.survive.M0 <- aggregate(Prop.Alive ~ Origin * Secondary, data = survive.M0, FUN= "std.error")  #calculate se by origin and secondary treatments
n.survive.M0 <- aggregate(Prop.Alive ~ Origin * Secondary, data = survive.M0, FUN= "length")  #calculate sample size by origin and secondary treatments
survivorship.M0 <- cbind(mean.survive.M0,se.survive.M0$Prop.Alive) #combine data
colnames(survivorship.M0) <- c("Origin", "Secondary", "mean", "se") #rename columns

#descriptive stats
(survivorship.M0[2,3]-survivorship.M0[1,3])/survivorship.M0[1,3] #percent change between treatments
(survivorship.M0[4,3]-survivorship.M0[3,3])/survivorship.M0[3,3] #percent change between treatments
mean.surs <- aggregate(Prop.Alive ~ Secondary, data = survive.M0, FUN= "mean") #calculate mean by secondary treatment
(mean.surs[1,2]-mean.surs[2,2])/mean.surs[1,2] #percent change between secondary treatments

Fig26 <- ggplot(data=survivorship.M0, aes(x=Secondary, y=mean, group=Origin, colour=Origin, shape=Origin)) + #plot data
  geom_line(size=0.7, position=position_dodge(.1)) + #plot lines
  scale_colour_manual(values=c("gray", "black")) + #set line color
  geom_point(size=4, position=position_dodge(.1), colour="black") + #set point characteristics
  scale_shape_manual(values=c(1,18)) + #set shapes
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), #plot error bars
                width=0, position=position_dodge(.1), colour="black") + #set error bar characteristics 
  ggtitle("(a) CHAMBER") + #plot title
  annotate("text", x = 0.87, y = 0.82, label = "a") + #add posthoc letters
  annotate("text", x = 0.85, y = 0.74, label = "ab") + #add posthoc letters
  annotate("text", x = 2.2, y = 0.59, label = "bc") + #add posthoc letters
  annotate("text", x = 2.15, y = 0.51, label = "cd") + #add posthoc letters
  xlab("Treatment of Offspring") + #plot x axis label
  ylab("Survivorship") + #plot y axis label
  ylim(0,1) + #Y axis limits
  theme_bw() + #theme black and white 
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.text=element_text(size=16), #set text size
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), vjust = 0.5), #set margins on labels
        axis.title=element_text(size=18,face="bold"), #set axis title text size
        axis.title.x=element_blank(), #remove x axis label
        strip.text.x = element_text(size = 16, colour = "black", face="bold"),
        axis.line.x = element_line(color = 'black'), #Set the axes color
        axis.line.y = element_line(color = 'black'), #Set the axes color
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_blank(),  #remove legend background
        legend.position="none",  #set legend location
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0)) #set title attributes

Fig26 #view plot

#Month1
larval.data.M1 <- read.csv("Larval_Data_M1.csv", header=T, sep=",", na.string="NA", as.is=T) #load data
proportion.alive.M1 <- larval.data.M1$month1/larval.data.M1$larvae.added #calculate survivorship
proportion.dead.M1 <- 1-proportion.alive.M1 #calculate mortality
larval.data.M1$Alive <- larval.data.M1$month1 #count alive
larval.data.M1$Dead <- larval.data.M1$larvae.added-larval.data.M1$month1 #calculate dead
survive.M1 <- data.frame (larval.data.M1$Chamber.num, larval.data.M1$Timepoint, larval.data.M1$Origin, larval.data.M1$Secondary, proportion.alive.M1, proportion.dead.M1, larval.data.M1$Alive, larval.data.M1$Dead) #make dataframe
colnames(survive.M1) <- c("Chamber", "Timepoint", "Origin", "Secondary", "Prop.Alive","Prop.Dead", "Alive","Dead") #rename columns
survive.M1$Timepoint <- "Time2"
mean.survive.M1 <- aggregate(Prop.Alive ~ Origin + Secondary, data = survive.M1, FUN= "mean") #calculate mean by origin and secondary treatments
se.survive.M1 <- aggregate(Prop.Alive ~ Origin + Secondary, data = survive.M1, FUN= "std.error") #calculate se by origin and secondary treatments
n.survive.M1 <- aggregate(Prop.Alive ~ Origin * Secondary, data = survive.M1, FUN= "length")  #calculate sample size by origin and secondary treatments
survivorship.M1 <- cbind(mean.survive.M1,se.survive.M1$Prop.Alive) #combine data
colnames(survivorship.M1) <- c("Origin", "Secondary", "mean", "se") #rename columns

#descriptive stats
(survivorship.M1[2,3]-survivorship.M1[1,3])/survivorship.M1[1,3] #percent change between treatments
(survivorship.M1[4,3]-survivorship.M1[3,3])/survivorship.M1[3,3] #percent change between treatments
mean.surs.m1 <- aggregate(Prop.Alive ~ Secondary, data = survive.M1, FUN= "mean") #calculate mean by secondary treatment
(mean.surs.m1[1,2]-mean.surs.m1[2,2])/mean.surs.m1[1,2] #percent change between secondary treatments

Fig27 <- ggplot(data=survivorship.M1, aes(x=Secondary, y=mean, group=Origin, colour=Origin, shape=Origin)) + #plot data
  geom_line(size=0.7, position=position_dodge(.1)) + #plot lines
  scale_colour_manual(values=c("gray", "black")) + #set line color
  geom_point(size=4, position=position_dodge(.1), colour="black") + #set point characteristics
  scale_shape_manual(values=c(1,18)) + #set shapes
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), #plot error bars
                width=0, position=position_dodge(.1), colour="black") + #set error bar characteristics 
  ggtitle("(b) MONTH 1") + #plot title
  annotate("text", x = 0.87, y = 0.52, label = "cd") + #add posthoc letters
  annotate("text", x = 0.83, y = 0.42, label = "de") + #add posthoc letters
  annotate("text", x = 2.2, y = 0.37, label = "de") + #add posthoc letters
  annotate("text", x = 2.15, y = 0.30, label = "e") + #add posthoc letters
  xlab("Treatment of Offspring") + #plot x axis label
  ylab("Survivorship") + #plot y axis label
  ylim(0,1) + #Y axis limits
  theme_bw() + #theme black and white
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.text=element_text(size=16), #set text size
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), vjust = 0.5), #set margins on labels
        axis.title=element_text(size=18,face="bold"), #set axis title text size
        axis.title.x=element_blank(), #remove x axis label
        strip.text.x = element_text(size = 16, colour = "black", face="bold"),
        axis.line.x = element_line(color = 'black'), #Set the axes color
        axis.line.y = element_line(color = 'black'), #Set the axes color
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_blank(),  #remove legend background
        legend.position="none",  #set legend location
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0)) #set title attributes

Fig27

#Month6
proportion.alive.M6 <- larval.data.M1$month6/larval.data.M1$larvae.added #claculate survival
proportion.dead.M6 <- 1-proportion.alive.M6 # calculate mortality
survive.M6 <- data.frame (larval.data.M1$Chamber.num, larval.data.M1$Timepoint, larval.data.M1$Origin, larval.data.M1$Secondary, proportion.alive.M6, proportion.dead.M6) #make dataframe
colnames(survive.M6) <- c("Chamber", "Timepoint", "Origin", "Secondary", "Prop.Alive","Prop.Dead") #rename columns
survive.M6$Timepoint <- "Time3" #identify timepoint
survive.M6$Alive <- larval.data.M1$month6 #count alive
survive.M6$Dead <- larval.data.M1$larvae.added-larval.data.M1$month6 #calculate dead
mean.survive.M6 <- aggregate(Prop.Alive ~ Origin + Secondary, data = survive.M6, FUN= "mean") #calculate mean
se.survive.M6 <- aggregate(Prop.Alive ~ Origin + Secondary, data = survive.M6, FUN= "std.error") #calculate SEM
n.survive.M6 <- aggregate(Prop.Alive ~ Origin * Secondary, data = survive.M6, FUN= "length")  #calculate sample size by origin and secondary treatments
survivorship.M6 <- cbind(mean.survive.M6,se.survive.M6$Prop.Alive) #combine descriptive statistics
colnames(survivorship.M6) <- c("Origin", "Secondary", "mean", "se") #rename columns

#descriptive stats
(survivorship.M6[2,3]-survivorship.M6[1,3])/survivorship.M6[1,3] #percent change between treatments
(survivorship.M6[4,3]-survivorship.M6[3,3])/survivorship.M6[3,3] #percent change between treatments
mean.surs.m6 <- aggregate(Prop.Alive ~ Secondary, data = survive.M6, FUN= "mean") #calculate mean by secondary treatment
(mean.surs.m6[1,2]-mean.surs.m6[2,2])/mean.surs.m6[1,2] #percent change between secondary treatments

Fig28 <- ggplot(data=survivorship.M6, aes(x=Secondary, y=mean, group=Origin, colour=Origin, shape=Origin)) + #plot data
  geom_line(size=0.7, position=position_dodge(.1)) + #plot lines
  scale_colour_manual(values=c("gray", "black"), labels=c("Ambient Parental Envt.", "High Parental Envt.")) + #set line color
  geom_point(size=4, position=position_dodge(.1), colour="black") + #set point characteristics
  scale_shape_manual(values=c(1,18), labels=c("Ambient Parental Envt.", "High Parental Envt.")) + #set shapes
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), #plot error bars
                width=0, position=position_dodge(.1), colour="black") + #set error bar characteristics 
  ggtitle("(c) MONTH 6") + #plot title
  annotate("text", x = 0.90, y = 0.19, label = "f") + #add posthoc letters
  annotate("text", x = 0.82, y = 0.12, label = "f") + #add posthoc letters
  annotate("text", x = 2.13, y = 0.06, label = "g") + #add posthoc letters
  annotate("text", x = 1.87, y = 0.0005, label = "fg") + #add posthoc letters
  xlab("Treatment of Offspring") + #plot x axis label
  ylab("Survivorship") + #plot y axis label
  ylim(0,1) + #Y axis limits
  theme_bw() + #theme black and white
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.text=element_text(size=16), #set text size
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), vjust = 0.5), #set margins on labels
        axis.title=element_text(size=18,face="bold"), #set axis title text size
        axis.title.x=element_blank(), #remove x axis label
        strip.text.x = element_text(size = 16, colour = "black", face="bold"),
        axis.line.x = element_line(color = 'black'), #Set the axes color
        axis.line.y = element_line(color = 'black'), #Set the axes color
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_blank(),  #remove legend background
        legend.position=c(0.6,0.7),  #set legend location
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0)) #set title attributes

Fig28

# #Repeated Measures Survivorship Origin = Fixed, Secondary = Fixed, Timepoint = Fixed, Chamber = Random repeated measure
All.Survivorship <- rbind(survive.M0, survive.M1, survive.M6) #combine data

#Binomial GLM
# Wald-test with H0 = 0
sur.GLM <-  glmer(cbind(Alive, Dead) ~ Origin*Secondary*Timepoint +(1|Chamber/Timepoint), data=All.Survivorship, family="binomial", na.action = "na.fail") #repeated measures ANOVA
summary(sur.GLM) #view summary
sur.mods <- dredge(sur.GLM) #describe model selection
sur.GLM <-  glmer(cbind(Alive, Dead) ~ Origin+Secondary+Timepoint + Secondary*Timepoint +(1|Chamber/Timepoint), data=All.Survivorship, family="binomial", na.action = "na.fail") #select best model for repeated measures ANOVA
Anova(sur.GLM) #view ANOVA table
Sur.Results <- Anova(sur.GLM) #view summary
dispersion_glmer(sur.GLM) #check for over dispersion
sur.resid <-resid(sur.GLM) #extract residuals
sur.shapiro <- shapiro.test(sur.resid) #runs a normality test using shapiro-wilk test on the residuals
sur.shapiro #view results
sur.qqnorm <- qqnorm(sur.resid) # normal quantile plot
sur.qqline <- qqline(sur.resid) # adding a qline of comparison
hist(sur.resid) #plot histogram of residuals
boxplot(sur.resid~ All.Survivorship$Origin * All.Survivorship$Secondary* All.Survivorship$Timepoint, ylab = "residuals", las = 2, par(mar = c(12, 5, 4, 2)+ 0.1)) #view Origin variability

#posthoc results
sur.GLM.posthoc <- summary(glht(sur.GLM, lsm(pairwise~Origin+Secondary+Timepoint)))
sur.GLM.posthoc 


##### SETTLEMENT #####
#Timepoint 1 only         
settlement.data <- larval.data.M0
settle <- (settlement.data$Plastic + settlement.data$Top.Tile + settlement.data$Bottom.Tile +  settlement.data$Edge)/(settlement.data$larvae.added)
settlement <- data.frame(settlement.data$Chamber.num, settlement.data$Origin, settlement.data$Secondary, settle)
colnames(settlement) <- c("Chamber", "Origin", "Secondary", "Prop.Settled") #rename columns
settlement$Settle <- (settlement.data$Plastic + settlement.data$Top.Tile + settlement.data$Bottom.Tile +  settlement.data$Edge) #count settlers
settlement$Not.Settle <- settlement.data$larvae.added-(settlement.data$Plastic + settlement.data$Top.Tile + settlement.data$Bottom.Tile +  settlement.data$Edge) #calculate not settled
mean.settled <- aggregate(Prop.Settled ~ Origin + Secondary, data = settlement, FUN= "mean") #calculate mean
se.settled <- aggregate(Prop.Settled ~ Origin + Secondary, data = settlement, FUN= "std.error") #calculate se
n.settled <- aggregate(Prop.Settled ~ Origin + Secondary, data = settlement, FUN= "length") #calculate se
settlement.data <- cbind(mean.settled, se.settled$Prop.Settled) #make dataframe
colnames(settlement.data) <- c("Origin", "Secondary", "mean", "se") #rename columns

#descriptive stats
(settlement.data[1,3]-settlement.data[2,3])/settlement.data[1,3] #percent change between treatments
(settlement.data[3,3]-settlement.data[4,3])/settlement.data[3,3] #percent change between treatments
mean.sets <- aggregate(Prop.Settled ~ Secondary, data = settlement, FUN= "mean") #calculate mean by secondary treatment
(mean.surs[1,2]-mean.surs[2,2])/mean.surs[1,2] #percent change between secondary treatments

Fig29 <- ggplot(data=settlement.data, aes(x=Secondary, y=mean, group=Origin, colour=Origin, shape=Origin)) + #plot data
  geom_line(size=0.7, position=position_dodge(.1)) + #plot lines
  scale_colour_manual(values=c("gray", "black")) + #set line color
  geom_point(size=4, position=position_dodge(.1), colour="black") + #set point characteristics
  scale_shape_manual(values=c(1,18)) + #set shapes
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), #plot error bars
                width=0, position=position_dodge(.1), colour="black") + #set error bar characteristics 
  annotate("text", x = 0.85, y = 0.80, label = "a") +
  annotate("text", x = 0.8, y = 0.68, label = "ab") +
  annotate("text", x = 2.25, y = 0.56, label = "bc") +
  annotate("text", x = 2.2, y = 0.49, label = "c") +
  ggtitle("(d) CHAMBER") + #plot title
  xlab("Treatment of Offspring") + #plot x axis label
  ylab("Settlement") + #plot y axis label
  ylim(0,1) + #Y axis limits
  theme_bw() + #theme black and white
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.text=element_text(size=16), #set text size
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), vjust = 0.5), #set margins on labels
        axis.title=element_text(size=18,face="bold"), #set axis title text size
        strip.text.x = element_text(size = 16, colour = "black", face="bold"),
        axis.line.x = element_line(color = 'black'), #Set the axes color
        axis.line.y = element_line(color = 'black'), #Set the axes color
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_blank(),  #remove legend background
        legend.position="none",  #set legend location
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0)) #set title attributes
Fig29

#Binomial GLM
# Wald-test with H0 = 0
set.GLM <-  glmer(cbind(Settle, Not.Settle) ~ Origin*Secondary +(1|Chamber), data=settlement, family="binomial", na.action = "na.fail") #repeated measures ANOVA with random intercept but not slope 
summary(set.GLM) #view summary
set.mods <- dredge(set.GLM) #describe model selection
set.GLM <-  glmer(cbind(Settle, Not.Settle) ~ Origin+Secondary +(1|Chamber), data=settlement, family="binomial", na.action = "na.fail") #repeated measures ANOVA with random intercept but not slope 
Set.Results <- Anova(set.GLM) #view summary
Set.Results#view ANOVA table
dispersion_glmer(set.GLM) #check for over dispersion
set.resid <-resid(set.GLM) #extract residuals
set.shapiro <- shapiro.test(set.resid) #runs a normality test using shapiro-wilk test on the residuals
set.shapiro #view results
sur.qqnorm <- qqnorm(set.resid) # normal quantile plot
sur.qqline <- qqline(set.resid) # adding a qline of comparison
hist(set.resid) #plot histogram of residuals
boxplot(set.resid~ settlement$Origin * settlement$Secondary, ylab = "residuals", las = 2, par(mar = c(12, 5, 4, 2)+ 0.1)) #view Origin variability

#posthoc results
set.GLM.posthoc <- summary(glht(set.GLM, lsm(pairwise~Origin+Secondary)))
set.GLM.posthoc 

##### GROWTH #####
data.M1 <- read.csv("Month1_Larval_Size.csv", header=T, sep=",", na.string="NA", as.is=T) #load data
data.M6 <- read.csv("Month6_Larval_Size.csv", header=T, sep=",", na.string="NA", as.is=T) #load data
growth.M1 <- aggregate(Polyp.Num.M1 ~ Date.M1 + Chamber.num, data = data.M1, FUN= "mean") #calculate size and survivorship per tile
growth.M6 <- aggregate(Polyp.Num.M6 ~ Date.M6 + Chamber.num, data = data.M6, FUN= "mean") #calculate size and survivorship per tile
growth.M6[growth.M6 == 0] <- NA #set zeros equal to NA
growth <- (cbind(growth.M1,growth.M6$Polyp.Num.M6, larval.data.M0$Origin, larval.data.M0$Secondary,larval.data.M0$Date, growth.M6$Date.M6)) #make dataframe
colnames(growth) <- c("Date.M1",  "Chamber.num",	"Polyp.Num.M1", "Polyp.Num.M6",	"Origin",	"Secondary",	"Date.M0", "Date.M6") #rename columns
growth$Date.M0<- as.Date(growth$Date.M0,format="%m/%d/%y") #set as date
growth$Date.M1<- as.Date(growth$Date.M1,format="%m/%d/%y") #set as date
growth$Date.M6<- as.Date(growth$Date.M6,format="%m/%d/%y") #set as date
growth$Days.M1 <- difftime(growth$Date.M1, growth$Date.M0, units = c("days")) #calculate the time difference in days
growth$Days.M6 <- difftime(growth$Date.M6, growth$Date.M1, units = c("days")) #calculate the time difference in days
growth$growth.rate.M1 <- (growth$Polyp.Num.M1-1)/(as.numeric(growth$Days.M1)) #calculate growth rate per day
growth$growth.rate.M6 <- (growth$Polyp.Num.M6-growth$Polyp.Num.M1)/(as.numeric(growth$Days.M6)) #calculate growth rate per day

m1.mean.growth <- aggregate(growth.rate.M1 ~ Origin + Secondary, data = growth, FUN= "mean") #calculate mean by origin and secondary treatments
m1.se.growth <- aggregate(growth.rate.M1 ~ Origin + Secondary, data = growth, FUN= "std.error") #calculate se by origin and secondary treatments
m1.n.growth <- aggregate(growth.rate.M1 ~ Origin + Secondary, data = growth, FUN= "length") #calculate se by origin and secondary treatments
m1.growth <- cbind(m1.mean.growth,m1.se.growth$growth.rate.M1) #combine data
colnames(m1.growth) <- c("Origin", "Secondary", "mean", "se") #rename columns
m6.mean.growth <- aggregate(growth.rate.M6 ~ Origin + Secondary, data = growth, FUN= "mean") #calculate mean by origin and secondary treatments
m6.se.growth <- aggregate(growth.rate.M6 ~ Origin + Secondary, data = growth, FUN= "std.error") #calculate se by origin and secondary treatments
m6.n.growth <- aggregate(growth.rate.M6 ~ Origin + Secondary, data = growth, FUN= "length") #calculate se by origin and secondary treatments
m6.growth <- cbind(m6.mean.growth,m6.se.growth$growth.rate.M6) #combine data
colnames(m6.growth) <- c("Origin", "Secondary", "mean", "se") #rename columns

Fig30 <- ggplot(data=m1.growth, aes(x=factor(Secondary), y=mean, group=Origin, colour=Origin, shape=Origin)) + #plot data
  geom_line(size=0.7, position=position_dodge(.1)) + #plot lines
  scale_colour_manual(values=c("gray", "black")) + #set line color
  geom_point(size=4, position=position_dodge(.1), colour="black") + #set point characteristics
  scale_shape_manual(values=c(1,18)) + #set shapes
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), #plot error bars
                width=0, position=position_dodge(.1), colour="black") + #set error bar characteristics 
  ggtitle("(e) MONTH 1") + #plot title
  xlab("Treatment of Offspring") + #plot x axis label
  ylab(expression(bold(~Growth~~"(polyps "*d^"1"*")"))) + #plot y axis label
  ylim(0,0.1) + #Y axis limits
  theme_bw() + #theme black and white
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.text=element_text(size=16), #set text size
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), vjust = 0.5), #set margins on labels
        axis.title=element_text(size=18,face="bold"), #set axis title text size
        strip.text.x = element_text(size = 16, colour = "black", face="bold"),
        axis.line.x = element_line(color = 'black'), #Set the axes color
        axis.line.y = element_line(color = 'black'), #Set the axes color
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_blank(),  #remove legend background
        legend.position="none",  #set legend location
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0)) #set title attributes

Fig30

Fig31 <- ggplot(data=m6.growth, aes(x=factor(Secondary), y=mean, group=Origin, colour=Origin, shape=Origin)) + #plot data
  geom_line(size=0.7, position=position_dodge(.1)) + #plot lines
  scale_colour_manual(values=c("gray", "black")) + #set line color
  geom_point(size=4, position=position_dodge(.1), colour="black") + #set point characteristics
  scale_shape_manual(values=c(1,18)) + #set shapes
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), #plot error bars
                width=0, position=position_dodge(.1), colour="black") +  #set error bar characteristics 
  ggtitle("(f) MONTH 6") +  #plot title
  xlab("Treatment of Offspring") + #plot x axis label
  ylab(expression(bold(~Growth~~"(polyps "*d^"1"*")"))) + #plot y axis label
  ylim(0,0.025) + #Y axis limits
  theme_bw() + #theme black and white
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.text=element_text(size=16), #set text size
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), vjust = 0.5), #set margins on labels
        axis.title=element_text(size=18,face="bold"), #set axis title text size
        strip.text.x = element_text(size = 16, colour = "black", face="bold"),
        axis.line.x = element_line(color = 'black'), #Set the axes color
        axis.line.y = element_line(color = 'black'), #Set the axes color
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_blank(),  #remove legend background
        legend.position="none",  #set legend location
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0)) #set title attributes

Fig31

#Repeated Measures growth
grow.M1<- cbind.data.frame(growth$Chamber.num, growth$Origin, growth$Secondary, growth$growth.rate.M1) #combine data
grow.M1$Timepoint <- "Time1" #identify time points
colnames(grow.M1) <- c( "Chamber.num", "Origin", "Secondary", "growth.rate", "Timepoint") #rename columns
grow.M6<- cbind.data.frame(growth$Chamber.num, growth$Origin, growth$Secondary, growth$growth.rate.M6) #combine data
grow.M6$Timepoint <- "Time6" #identify time points
colnames(grow.M6) <- c( "Chamber.num", "Origin", "Secondary", "growth.rate", "Timepoint") #rename columns
All.Growth <- rbind(grow.M1, grow.M6) #combine data
All.Growth <- na.omit(All.Growth) #remove NA rows

Growth.RM <- lme(log10(growth.rate+1) ~ Origin*Secondary*Timepoint, random = ~ 1|Chamber.num/Timepoint, data=All.Growth, na.action = "na.fail") #repeated measures ANOVA
summary(Growth.RM) #view results
Grow.Results <- anova(Growth.RM) #view results
anova(Growth.RM) #view results
grow.mods <- dredge(Growth.RM) #describe model selection
Growth.RM <- lme(log10(growth.rate+1) ~ Origin+Timepoint, random = ~ 1|Chamber.num/Timepoint, data=All.Growth, na.action = "na.fail") #repeated measures ANOVA
Grow.Results <- anova(Growth.RM) #view results
anova(Growth.RM) #view results
gro.resid <-resid(Growth.RM) #extract residuals
gro.shapiro <- shapiro.test(gro.resid) #runs a normality test using shapiro-wilk test on the residuals
gro.shapiro #view results
gro.qqnorm <- qqnorm(gro.resid) # normal quantile plot
gro.qqline <- qqline(gro.resid) # adding a qline of comparison
hist(gro.resid) #plot histogram of residuals
boxplot(gro.resid~ All.Growth$Origin * All.Growth$Secondary* All.Growth$Timepoint, ylab = "residuals", las = 2, par(mar = c(12, 5, 4, 2)+ 0.1)) #view Origin variability

#posthoc results
Growth.RM.posthoc <- summary(glht(Growth.RM, lsm(pairwise~Origin+Timepoint)))
Growth.RM.posthoc 

#transform and calculate descriptive stats
All.Growth$logged <- log10(All.Growth$growth.rate +1)
mean.growth <- aggregate(logged ~ Origin + Secondary + Timepoint, data = All.Growth, FUN= "mean") #calculate mean by origin and secondary treatments
se.growth <- aggregate(logged ~ Origin + Secondary + Timepoint, data = All.Growth, FUN= "std.error") #calculate se by origin and secondary treatments
growth <- cbind(mean.growth,se.growth$logged) #combine data
m1.growth <- subset(growth, Timepoint=="Time1")
colnames(m1.growth) <- c("Origin", "Secondary","Timepoint", "mean", "se") #rename columns
m6.growth <- subset(growth, Timepoint=="Time6") 
colnames(m6.growth) <- c("Origin", "Secondary","Timepoint", "mean", "se") #rename columns

#backtransform means and asymetrical error
m1.growth.bt <- m1.growth #assign data
m1.growth.bt$mean <- 10^(m1.growth.bt$mean)-1 #backtransform
m1.growth.bt$upper <- m1.growth$mean + m1.growth$se #upper sem value
m1.growth.bt$lower <- m1.growth$mean - m1.growth$se #lower sem value
m1.growth.bt$upper.bt <- 10^(m1.growth.bt$upper)-1 #backtransform
m1.growth.bt$lower.bt <- 10^(m1.growth.bt$lower)-1 #backtransform

#descriptive stats
g1 <-m1.growth.bt[2,4]/m1.growth.bt[1,4] #fold change between treatments
g2 <-m1.growth.bt[4,4]/m1.growth.bt[3,4] #fold change between treatments
mean(g1,g2)

Fig32 <- ggplot(data=m1.growth.bt, aes(x=factor(Secondary), y=mean, group=Origin, colour=Origin, shape=Origin)) + #plot data
  geom_line(size=0.7, position=position_dodge(.1)) + #plot lines
  scale_colour_manual(values=c("gray", "black")) + #set line color
  geom_point(size=4, position=position_dodge(.1), colour="black") + #set point characteristics
  scale_shape_manual(values=c(1,18)) + #set shapes
  geom_errorbar(aes(ymin=lower.bt, ymax=upper.bt), #plot error bars
                width=0, position=position_dodge(.1), colour="black") + #set error bar characteristics 
  ggtitle("(e) MONTH 1") + #plot title
  xlab("Treatment of Offspring") + #plot x axis label
  ylab(expression(bold(~Growth~~"(polyps "*d^"1"*")"))) + #plot y axis label
  ylim(0,0.1) + #Y axis limits
  theme_bw() + #theme black and white
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.text=element_text(size=16), #set text size
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), vjust = 0.5), #set margins on labels
        axis.title=element_text(size=18,face="bold"), #set axis title text size
        strip.text.x = element_text(size = 16, colour = "black", face="bold"),
        axis.line.x = element_line(color = 'black'), #Set the axes color
        axis.line.y = element_line(color = 'black'), #Set the axes color
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_blank(),  #remove legend background
        legend.position="none",  #set legend location
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0)) #set title attributes

Fig32

#backtransform means and asymetrical error
m6.growth.bt <- m6.growth #assign data
m6.growth.bt$mean <- 10^(m6.growth.bt$mean)-1 #backtransform
m6.growth.bt$upper <- m6.growth$mean + m6.growth$se #upper sem value
m6.growth.bt$lower <- m6.growth$mean - m6.growth$se #lower sem value
m6.growth.bt$upper.bt <- 10^(m6.growth.bt$upper)-1 #backtransform
m6.growth.bt$lower.bt <- 10^(m6.growth.bt$lower)-1 #backtransform

Fig33 <- ggplot(data=m6.growth.bt, aes(x=factor(Secondary), y=mean, group=Origin, colour=Origin, shape=Origin)) + #plot data
  geom_line(size=0.7, position=position_dodge(.1)) + #plot lines
  scale_colour_manual(values=c("gray", "black")) + #set line color
  geom_point(size=4, position=position_dodge(.1), colour="black") + #set point characteristics
  scale_shape_manual(values=c(1,18)) + #set shapes
  geom_errorbar(aes(ymin=lower.bt, ymax=upper.bt), #plot error bars
                width=0, position=position_dodge(.1), colour="black") +  #set error bar characteristics 
  ggtitle("(f) MONTH 6") +  #plot title
  xlab("Treatment of Offspring") + #plot x axis label
  ylab(expression(bold(~Growth~~"(polyps "*d^"1"*")"))) + #plot y axis label
  ylim(0,0.1) + #Y axis limits
  theme_bw() + #theme black and white
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.text=element_text(size=16), #set text size
        axis.ticks.length=unit(-0.2, "cm"), #turn ticks inward
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), #set margins on labels
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), vjust = 0.5), #set margins on labels
        axis.title=element_text(size=18,face="bold"), #set axis title text size
        strip.text.x = element_text(size = 16, colour = "black", face="bold"),
        axis.line.x = element_line(color = 'black'), #Set the axes color
        axis.line.y = element_line(color = 'black'), #Set the axes color
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_blank(),  #remove legend background
        legend.position="none",  #set legend location
        panel.border=element_rect(size=1.25, fill = NA), #set outer border
        plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0)) #set title attributes

Fig33

#####ALL FIGURES, TABLES, AND STATISTICAL RESULTS#####
setwd(file.path(mainDir, 'Output'))

#Capture Figures to File
pdf("Fig2.Larval.release.pdf", width = 11, height = 6)
inset <- viewport(width = 0.22, height = 0.5, x = 0.86, y = 0.65)  # the inset in upper right
grid.newpage()
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
pushViewport(viewport(layout = grid.layout(1, 3)))
print(Fig22, vp = vplayout(1, 1))
print(Fig23, vp = vplayout(1, 2))
print(Fig24, vp = vplayout(1, 3))
dev.off()

Larval.Perform <- arrangeGrob(Fig26, Fig27, Fig28,
                              Fig29, Fig32, Fig33, ncol=3)
ggsave(file="Fig3.Larval.Performance.pdf", Larval.Perform, width =12, height = 6, units = c("in"))

FigureS1.Physical <- grid.arrange(arrangeGrob(Fig2, Fig5, Fig9, Fig13, left="TEMPERATURE", ncol=4),
                                  arrangeGrob(Fig3, Fig7, Fig11, Fig15, left="IRRADIANCE", ncol=4),
                                  arrangeGrob(FigMT, Fig17, Fig19, Fig21, left="pH", ncol=4), ncol=1)
ggsave(file="FigS1.Physical_Experimental_Conditions.pdf", FigureS1.Physical, width =11, height = 8.5, units = c("in"))

# write.table(adult.chem.table, "Seawater_chemistry_table_Output_Adult.csv", sep=",", row.names = FALSE)
# write.table(M1.chem.table, "Seawater_chemistry_table_Output_M1.csv", sep=",", row.names = FALSE)
# write.table(M6.chem.table, "Seawater_chemistry_table_Output_M6.csv", sep=",", row.names = FALSE)

tt2 <- ttheme_minimal()
title1 <- "A) Adult Exposure"
title2 <- "B) Month 1 Exposure"
title3 <- "C) Month 6 Exposure"
t1 <- grid.text(title1, just="left")
t2 <- grid.text(title2, just="left")
t3 <- grid.text(title3, just="left")
SW.Chem.Tables <- grid.arrange(
  t1,
  tableGrob(adult.chem.table, theme=tt2),
  t2,
  tableGrob(M1.chem.table, theme=tt2),
  t3,
  tableGrob(M6.chem.table, theme=tt2),
  nrow=6)
ggsave(file="SW.Chemistry.Table.pdf", SW.Chem.Tables, width = 11, height = 6)

### Generate Stats Table
survivorship <- as.data.frame(Sur.Results)
survivorship <-round(survivorship[,],3)
survivorship
settlement <- as.data.frame(Set.Results)
settlement <-round(settlement[,],3)
settlement
growth <- anova(Growth.RM)
growth <-round(growth[,],3)
growth

# pdf("table2.pdf", width = 11, height = 6, nrow=3)
# sur.table <- grid.table(survivorship)
# set.table <- grid.table(settlement)
# gro.table <- grid.table(growth)
# dev.off()

#Capture statistical results to file
capture.output(june.ks, july.ks, august.ks, sur.mods, survivorship, set.mods, settlement, grow.mods, growth,  file="HI_Pdam_Parental_Stat_Results.txt")

setwd(file.path(mainDir, 'Data'))









