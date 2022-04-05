##############################################################
# ENA of Western Baffin Bay model of 2016  ###################
##############################################################

library(devtools)

#To install R tools I followed the instructions in here: https://cran.r-project.org/bin/windows/Rtools/

writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
Sys.which("make")
install.packages("jsonlite", type = "source") #this worked!

url <- "https://cran.r-project.org/src/contrib/Archive/enaR/enaR_3.0.0.tar.gz" #download the package enaR from cran
pkgFile <- "enaR_3.0.0.tar.gz"
download.file(url = url, destfile = pkgFile)
install.packages(pkgs=pkgFile, type="source", repos=NULL)
#Rtools may not have been successfully installed

#Libraries needed for this analysis
library(sna) #sna and gdata are dependencies of enaR - the package cannot be installed without these other packages
library(gdata)
library(enaR)
library(limSolve)
library(network)
library(GGally)
library(grid)
library(gridExtra)
library(gtable)
library(dplyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(extrafont)
font_import()

# This analysis was based work by Bentley et al 2019 (Ecological Indicators)
# How to load in model - need to change Ecopath model into the correct format to be read in R - we will use the SCOR format
# The Ecopath output needs to be reformatted into:
# (1) a flow matrix (Ecopath consumption matrix) 
# (2) network inputs (gross primary production) 
# (3) network exports (detritus and fisheries exports) 
# (4) respiration (calculated in EwE outputs)  
# (5) node storage (Ecopath biomass)

# As Ecopath does not provide estimates of gross primary productivity or respiration for primary producers,
# we used methods for calculating respiration and gross production for aquatic plants from a model equation proposed by Aoki (2006). 
# Gross production a() is consumed by respiration r(), production p(), and flow to detritus d( ): 
#a = r + p + d
#In accordance with Aoki (2006), they assumed that:
#r = a.rho 
#where rho is the annual respiration-gross production ratio. According to the equations:
#r = (p + d) . rho /(1-rho)
#r can therefore be estimated from values of p, d and rho, where p and d are available from Ecopath. 

# In the SCOR format file, I added diet import to the imports values (along with gross primary production - a)
 
# For gross primary production and respiration of primary producers: according to Saint-Béat et al 2019 (Elementa) the gross primary production in Baffin Bay
# is 20.4 gC/m2/month, but this does not distinguish between sea ice algae and phytoplankton. 
# Respiration was 3.15 gC/m2/month (again not distinguishing sea ice algae from phytoplankton).
# In Ecopath in ton/km2/year net production of phytoplankton is 530.4 (from GreenEdge data) and flow to detritus is 276.2 (from consumption flow matrix); 
# for sea ice algae production is 432.06 and flow to detritus is 144.1. 
# Saint-Béat et al considered respiration between 0.05*a and 0.3*a, where 0.05 and 0.3 is rho according to Aoki (2006). Based on overall gross primary production
# and respiration for western Baffin Bay I can find the rho used in this study and then from there calculate gross primary production for sea ice algae and phytoplakton:

rho <- 3.15/20.4

phy_r <- (530.4 + 276.2)*(rho/(1-rho))
SIA_r <- (432.06 + 144.1)*(rho/(1-rho))

phy_a <- phy_r + 530.4 + 276.2
SIA_a <- SIA_r + 432.06 + 144.1

#These values were added to the SCOR format text file.

# Load data

WBB <- read.scor("WBB_model12_SCOR.txt", from.file = TRUE, warn = TRUE)

#In the enaR package, a complete ecosystem network model description includes:
# F - flow matrix, oriented row-to-column
# z - vector of inputs
# r - vector of respirations
# e - vector of exports
# y - vector of outputs (respirations + exports)
# X - a vector of biomass or storage values
# Living - logical vector indicating if the node is living (TRUE) or non-living (FALSE)

#Pulling out the "vertex" (i.e. node) attributes 
WBB%v%'output'
WBB%v%'input'
WBB%v%'living'
WBB%v%'respiration' 
WBB%v%'storage'
WBB%v%'vertex.names'
WBB%v%'export'
WBB%v%'export'

#or
unpack(WBB)

#This allows me to visualize the mattrix of flows so I can check if it is correct
flows_matrix <- as.matrix(WBB,attrname="flow")
cons_flows_WBB <- as.data.frame(flows_matrix) #looks good

#Quick visualization of the data
plot(WBB) 

## Check to see if the model is balanced
ssCheck(WBB)  # not balanced...
WBB <- balance(WBB, method="AVG2") #can use different methods for balancing the model


## Ecological network analysis 

# Flow analysis - TST, APL FCI and IFI
Flow <- enaFlow(WBB)
attributes(Flow)
Flow$ns

# Control analysis - CD and system control
C <- enaControl(WBB)
attributes(C)

# Ascendency
enaAscendency(WBB)

# Mixed Trophic Impacts
mti <- enaMTI(WBB)
attributes$mti
mti <- enaMTI(WBB,eigen.check=FALSE)
attributes(mti)
mti$Relations.Table
mti$M


# Trophic aggregations (Lindeman's trophic concepts)
trop <- enaTroAgg(WBB)
attributes(trop)



##############################################################################################
## Uncertainty analysis of ENA results according to Hines et al 2018 & Bentley et al 2019 ####
##############################################################################################

# Once a network model is constructed, uncertainty data for each edge weight can be combined with the original model through a
# LIM uncertainty analysis. This uncertainty analysis can be performed by the enaUncertainty() function. This function first uses
# the original network model, along with uncertainty data tables, to generate the inputs needed for LIM. It then applies a LIM 
# analysis to generate plausible model sets using the limSolve package for R. Next, LIM results are converted back into a format
# readable by enaR and are returned to the user as a set of plausible model parameterizations for use in ENA (Hines et al).

iter = 10 # the number of plausible models to return - if testing the script use a smaller number of iterations as it takes a while to run

###############################################################################################
#percent method with similar uncertainty levels for all edges - ###############################
###############################################################################################

plausible.per.25 <- enaUncertainty(x=WBB, type="percent", p.err=25, iter=iter) # 25% error applied to all parameters
plausible.per.50 <- enaUncertainty(x=WBB, type="percent", p.err=50, iter=iter) # 50% error applied to all parameters

###############################################################################################
#percent method with different symmetric uncertainties  #######################################
###############################################################################################

# Load the uncertainty values of the allowable amount above and below each edge

f.df <- read.csv("flows.f.respCVs.zoo.csv") # based on Ecopath CVs derived from the pedigree table, except for zoo which is based in Saint-Béat et al
z.df <- read.csv("inputs.z.respCVs.csv") # based on ...
r.df <- read.csv("respiration.r.zoo.csv") # based on CVs from literature on respiration rates for some species (since this is not an input in Ecopath and no pedigree is available)
e.df <- read.csv("exports.e.csv") # based temporal variation in fisheries catches 

# Testing different levels of uncertainty applied to respiration rates to see how different the results look

r.25.df <- read.csv("respiration.r.25.csv") #applying 25 % CV to all respiration values
r.50.df <- read.csv("respiration.r.50.csv") #applying 50 % CV ...
r.10.df <- read.csv("respiration.r.10.csv") #applying 10 % CV ...



#Comparing different uncertainty levels for respiration rates (not in the Ecopath pedigree)
plausible.sym.10r <- enaUncertainty(x=WBB, F.sym=f.df, z.sym=z.df, e.sym=e.df, r.sym=r.10.df, type="sym", iter=iter)
plausible.sym.25r <- enaUncertainty(x=WBB, F.sym=f.df, z.sym=z.df, e.sym=e.df, r.sym=r.25.df, type="sym", iter=iter)
plausible.sym.50r <- enaUncertainty(x=WBB, F.sym=f.df, z.sym=z.df, e.sym=e.df, r.sym=r.50.df, type="sym", iter=iter)

plausible.sym <- enaUncertainty(x=WBB, F.sym=f.df, z.sym=z.df, e.sym=e.df, r.sym=r.df, type="sym", iter=iter)


#Plotting to compare the values using boxplots
# original model
TST.WBB <- enaFlow(WBB)$ns[2]
FCI.WBB <- enaFlow(WBB)$ns[5]

# different uncertainty levels
TST.per.25 <- unlist(lapply(plausible.per.25, function(x) enaFlow(x)$ns[2]))
TST.per.50 <- unlist(lapply(plausible.per.50, function(x) enaFlow(x)$ns[2]))
TST.per.pedigree <- unlist(lapply(plausible.sym, function(x) enaFlow(x)$ns[2]))
TST.per.pedigree.10.r <- unlist(lapply(plausible.sym.10r, function(x) enaFlow(x)$ns[2]))
TST.per.pedigree.25.r <- unlist(lapply(plausible.sym.25r, function(x) enaFlow(x)$ns[2]))
TST.per.pedigree.50.r <- unlist(lapply(plausible.sym.50r, function(x) enaFlow(x)$ns[2]))

FCI.per.25 <- unlist(lapply(plausible.per.25, function(x) enaFlow(x)$ns[5]))
FCI.per.50 <- unlist(lapply(plausible.per.50, function(x) enaFlow(x)$ns[5]))
FCI.per.pedigree <- unlist(lapply(plausible.sym, function(x) enaFlow(x)$ns[5]))
FCI.per.pedigree.10.r <- unlist(lapply(plausible.sym.10r, function(x) enaFlow(x)$ns[5]))
FCI.per.pedigree.25.r <- unlist(lapply(plausible.sym.25r, function(x) enaFlow(x)$ns[5]))
FCI.per.pedigree.50.r <- unlist(lapply(plausible.sym.50r, function(x) enaFlow(x)$ns[5]))


ENA.ind.test.5 <- as.data.frame(cbind(TST.per.pedigree, FCI.per.pedigree))

write.csv(ENA.ind.test.5, file = "ENA_indicators_10000per_comparisons_5.csv") #save to a csv file to compare results in a plot

TST.comparisons <- read.csv("ENA_indicators_10000per_comparisons_TST.csv")
FCI.comparisons <- read.csv("ENA_indicators_10000per_comparisons_FCI.csv")

# compare using different uncertainty values with boxplots

tiff("Different_uncertainty_levels_TST.tiff", width = 8, height = 6, units = 'in', res = 500)

colnames(TST.comparisons) <- c("25% all flows", "50% all flows", "Resp flow & pedigree")
par(las=1)
bp=boxplot(TST.comparisons, main="TST")
points(x=c(1,2,3,4,5), y=rep(TST.Ecopath, 5), col="red", pch=16, cex=1)
legend("bottomleft", legend=c("original model"), pch=16, col="red", bty="n", cex=0.75)

dev.off()


tiff("Different_uncertainty_levels_FCI.tiff", width = 8, height = 6, units = 'in', res = 500)

colnames(FCI.comparisons) <- c("25% all flows", "50% all flows", "Resp flow & pedigree")
par(las=1)
bp=boxplot(FCI.comparisons, main="FCI")
points(x=c(1,2,3,4,5), y=rep(FCI.Ecopath, 5), col="red", pch=16, cex=1)
legend("bottomleft", legend=c("original model"), pch=16, col="red", bty="n", cex=0.75)

dev.off()

######################

# comparison of several methods for respiration flows

tiff("Respiration_uncertainty_TST.tiff", width = 8, height = 6, units = 'in', res = 500)

TST.plot.data=cbind(TST.per.pedigree, TST.per.pedigree.25.r, TST.per.pedigree.50.r, TST.per.pedigree.10.r)
colnames(TST.plot.data) <- c("pedigree", "pedigree.r25%", "pedigree.r50%", "pedigree.r10%")
par(las=1)
bp=boxplot(TST.plot.data, main="TST")
points(x=c(1,2,3,4), y=rep(TST.WBB, 4), col="red", pch=16, cex=1)
legend("bottomleft", legend=c("original model"), pch=16, col="red", bty="n", cex=0.75)

dev.off()

tiff("Respiration_uncertainty_FCI.tiff", width = 8, height = 6, units = 'in', res = 500)

FCI.plot.data=cbind(FCI.per.pedigree, FCI.per.pedigree.25.r, FCI.per.pedigree.50.r, FCI.per.pedigree.10.r)
colnames(FCI.plot.data) <- c("pedigree", "pedigree.r25%", "pedigree.r50%", "pedigree.r10%")
par(las=1)
bp2=boxplot(FCI.plot.data, main="FCI")
points(x=c(1,2,3,4), y=rep(FCI.WBB, 4), col="red", pch=16, cex=1)
legend("bottomleft", legend=c("original model"), pch=16, col="red", bty="n", cex=0.75)

dev.off()


# extract indicators from ena results with uncertainty from pedigree

APL.per.pedigree <- unlist(lapply(plausible.sym, function(x) enaFlow(x)$ns[4])) #average path length
FCI.per.pedigree <- unlist(lapply(plausible.sym, function(x) enaFlow(x)$ns[5])) #Fin's Cycling Index
IFI.per.pedigree <- unlist(lapply(plausible.sym, function(x) enaFlow(x)$ns[8])) #Indirect flow intensity
asc.per.pedigree <- unlist(lapply(plausible.sym, function(x) enaAscendency(x)[[5]])) #ascendency
rel.asc.per.pedigree <- unlist(lapply(plausible.sym, function(x) enaAscendency(x)[[7]])) #relative ascendency
AMI.per.pedigree <- unlist(lapply(plausible.sym, function(x) enaAscendency(x)[[2]])) #average mutual information
TST.per.pedigree <- unlist(lapply(plausible.sym, function(x) enaFlow(WBB)$ns[2])) #total system control


# average values from ENAr model (based on Ecopath converted model and without uncertainty distributions)
APL.WBB <- enaFlow(WBB)$ns[4]
FCI.WBB <- enaFlow(WBB)$ns[5]
IFI.WBB <- enaFlow(WBB)$ns[8]
LD.WBB <- enaStructure(WBB)$ns[4]
asc.WBB <- enaAscendency(WBB)[[5]]
rel.asc.WBB <- enaAscendency(WBB)[[7]]
AMI.WBB <- enaAscendency(WBB)[[2]]
TST.WBB <- enaFlow(WBB)$ns[2]
robustness.WBB <- enaAscendency(WBB)[[9]]



# save a table with the permutations for all these indices
ENA.indicators <- as.data.frame(cbind(AMI.per.pedigree, APL.per.pedigree, TST.per.pedigree, FCI.per.pedigree, asc.per.pedigree, rel.asc.per.pedigree,
                                      IFI.per.pedigree))

write.csv(ENA.indicators, file = "ENA_indicators_10000per.csv")


# density plots
AMI.plot <- ggplot(ENA.indicators, aes(x=AMI.per.pedigree)) + 
            geom_density(alpha = 0.4, fill = "lightgreen", color = "darkgreen") + 
            geom_vline(aes(xintercept=mean(AMI.per.pedigree)), color="darkgreen", linetype="dashed", size=1) +
            geom_vline(aes(xintercept=AMI.WBB), color="black", linetype="dotted", size=1) +
            xlab("Average Mutual Information (AMI)") + ylab("Density") +
            theme_bw() + #xlim(4200, 4650) +
            scale_y_continuous(expand = c(0,0)) +
            theme(axis.text.y = element_text(color="black"),
               axis.text.x = element_text(color = "black"),
                 plot.margin =unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                  panel.background = element_blank())


APL.plot <- ggplot(ENA.indicators.resp.zoo, aes(x=APL.per.pedigree)) + 
  geom_density(alpha = 0.4, fill = "coral", color = "coral3") + 
  geom_vline(aes(xintercept=mean(APL.per.pedigree)), color="coral3", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=APL.WBB), color="black", linetype="dotted", size=1) +
  xlab("Average Path Length (APL)") + ylab("") +
  theme_bw() + #xlim(4200, 4650) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color = "black"),
        plot.margin =unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank())


FCI.plot <- ggplot(ENA.indicators.resp.zoo, aes(x=FCI.per.pedigree)) + 
  geom_density(alpha = 0.4, fill = "coral", color = "coral3") + 
  geom_vline(aes(xintercept=mean(FCI.per.pedigree)), color="coral3", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=FCI.WBB), color="black", linetype="dotted", size=1) +
  xlab("Finn's Cycling Index (FCI)") + ylab("") +
  theme_bw() + #xlim(4200, 4650) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color = "black"),
        plot.margin =unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank())

TST.plot <- ggplot(ENA.indicators.resp.zoo, aes(x=TST.per.pedigree)) + 
  geom_density(alpha = 0.4, fill = "lightgreen", color = "darkgreen") + 
  geom_vline(aes(xintercept=mean(TST.per.pedigree)), color="darkgreen", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=TST.WBB), color="black", linetype="dotted", size=1) +
  xlab("Total System Throughput (TST)") + ylab("Density") +
  theme_bw() + #xlim(4200, 4650) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color = "black"),
        plot.margin =unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank()) 

asc.plot <- ggplot(ENA.indicators.resp.zoo, aes(x=asc.per.pedigree)) + 
  geom_density(alpha = 0.4, fill = "lightgreen", color="darkgreen") + 
  geom_vline(aes(xintercept=mean(asc.per.pedigree)), color="darkgreen", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=asc.WBB), color="black", linetype="dotted", size=1) +
  xlab("Ascendency (A)") + ylab("Density") +
  theme_bw() + #xlim(4200, 4650) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color = "black"),
        plot.margin =unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank())

rel.asc.plot <- ggplot(ENA.indicators.resp.zoo, aes(x=rel.asc.per.pedigree)) + 
  geom_density(alpha = 0.4, fill = "lightgreen", color="darkgreen") + 
  geom_vline(aes(xintercept=mean(rel.asc.per.pedigree)), color="darkgreen", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=rel.asc.WBB), color="black", linetype="dotted", size=1) +
  xlab("Relative Ascendency (A/C)") + ylab("Density") +
  theme_bw() + #xlim(4200, 4650) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color = "black"),
        plot.margin =unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank())

IFI.plot <- ggplot(ENA.indicators.resp.zoo, aes(x=IFI.per.pedigree)) + 
  geom_density(alpha = 0.4, fill = "coral", color = "coral3") + 
  geom_vline(aes(xintercept=mean(IFI.per.pedigree)), color="coral3", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=IFI.WBB), color="black", linetype="dotted", size=1) +
  xlab("Indirect Flow Intensity (IFI)") + ylab("") +
  theme_bw() + #xlim(4200, 4650) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color = "black"),
        plot.margin =unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank())

rob.plot <- ggplot(ENA.indicators.resp.zoo, aes(x=robustness.per.pedigree)) + 
  geom_density(alpha = 0.4, fill = "lightgreen", color="darkgreen") + 
  geom_vline(aes(xintercept=mean(robustness.per.pedigree)), color="darkgreen", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=robustness.WBB), color="black", linetype="dotted", size=1) +
  xlab("Robustness") + ylab("") +
  theme_bw() + #xlim(4200, 4650) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color = "black"),
        plot.margin =unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank())

#save plots

tiff("Density_plots_1.tiff", width = 4, height = 8, units = 'in', res = 500)

panel.plots1 <- rbind(ggplotGrob(rob.plot), ggplotGrob(FCI.plot), ggplotGrob(IFI.plot), ggplotGrob(APL.plot), size = "first")
panel.plots1$widths <- unit.pmax(ggplotGrob(FCI.plot)$widths, ggplotGrob(IFI.plot)$widths, ggplotGrob(APL.plot)$widths)
grid.newpage()
grid.draw(panel.plots1)

dev.off()

tiff("Density_plots_2.tiff", width = 4, height = 8, units = 'in', res = 500)

panel.plots2 <- rbind(ggplotGrob(TST.plot), ggplotGrob(AMI.plot), ggplotGrob(asc.plot), ggplotGrob(rel.asc.plot), size = "first")
panel.plots2$widths <- unit.pmax(ggplotGrob(TST.plot)$widths, ggplotGrob(AMI.plot)$widths, ggplotGrob(asc.plot)$widths, 
                                ggplotGrob(rel.asc.plot)$widths)
grid.newpage()
grid.draw(panel.plots2)

dev.off()



# System control

CD.per.pedigree <- unlist(lapply(plausible.sym, function(x) enaControl(x)$CD)) #Schramski Control Difference Matrix
sc.per.pedigree <- unlist(lapply(plausible.sym, function(x) enaControl(x)$sc)) #System control


control.difference.mattrix <- as.data.frame(CD.per.pedigree)
write.csv(control.difference.mattrix, file = "control.difference.10000per.csv") 


f.group <- rep(c("KILLER WHALE", "POLAR BEAR", "NARWHAL", "BOWHEAD WHALE", "RINGED SEAL", "OTHER SEALS", "WALRUS", "SEABIRDS", "GREENLAND SHARK", "GREENLAND HALIBUT", "ARCTIC CHAR", "ARCTIC/POLAR COD",
             "SMALL PELAGIC FISH", "SCULPINS/EELPOUTS", "SMALL DEMERSAL FISH", "LARGE DEMERSAL FISH", "LARGE CRUSTACEANS", "CEPHALOPODS", "CARNIVOROUS ZOO", "OMNIVOROUS ZOO", "CALANUS ZOO", "MICROZOOPLANKTON",
             "POLYCHAETES", "ECHINODERMATA", "BIVALVES", "OTHER BENTHOS", "BACTERIA", "ICE ALGAE", "PHYTOPLANKTON", "DETRITUS"), times = 10000)

system.control <- data.frame(f.group, sc.per.pedigree)
write.csv(system.control, file = "system.control.10000per.csv") #not sure how the matrix is organized when doing the permutations 


# plotting the system control
C <- enaControl(WBB)
attributes(C)

# Calculating 95 % confidence intervals for system control
#plot with confidence intervals as error

sys.control <- as.data.frame(enaControl(WBB)$sc)
rownames(sys.control)[rownames(sys.control) == "MARINE WORMS"] <- "POLYCHAETES"
sys.control$f.group <- rownames(sys.control)
sys.control <- rename(sys.control, "sc" = "enaControl(WBB)$sc" )

system.control.resp.zoo.open$groupN <- rep(c(1:30),times=10000)

system.control.t <- system.control.resp.zoo.open %>%
  group_by(f.group) %>%
  summarize(sd.sc = sd(sc.per.resp.zoo), error = qnorm(0.975)*sd.sc/sqrt(iter)) # lower and upper error for each group

sys.control <- sys.control %>% left_join(system.control.t)

# respiration uncertainty levels


system.control.resp.zoo$f.group[system.control.resp.zoo$f.group=="MARINE WORMS"]<-"POLYCHAETES"

system.control.t.resp <- system.control.resp.zoo %>%
  group_by(f.group) %>%
  summarize(sd.sc = sd(sc.per.resp.zoo), error = qnorm(0.975)*sd.sc/sqrt(iter)) # lower and upper error for each group

sys.control.resp <- sys.control %>% left_join(system.control.t.resp)


# Plot

sys.cont.plot <- ggplot(sys.control.resp, aes(x = f.group, y = sc)) +
  geom_bar(color = "purple3", size = 0.3, stat = "identity", fill = "purple", alpha = 0.4) +
  geom_errorbar(aes(ymin=sc-sd.sc, ymax=sc+sd.sc), width=0.3, size = 0.3) +
  xlab("") + ylab("System Control") +
  theme_bw() +
  scale_x_discrete(limits = c("KILLER WHALE", "POLAR BEAR", "NARWHAL", "BOWHEAD WHALE", "RINGED SEAL", "OTHER SEALS", "WALRUS", "SEABIRDS", "GREENLAND SHARK", "GREENLAND HALIBUT", "ARCTIC CHAR", "ARCTIC/POLAR COD",
                            "SMALL PELAGIC FISH", "SCULPINS/EELPOUTS", "SMALL DEMERSAL FISH", "LARGE DEMERSAL FISH", "LARGE CRUSTACEANS", "CEPHALOPODS", "CARNIVOROUS ZOO", "OMNIVOROUS ZOO", "CALANUS ZOO", "MICROZOOPLANKTON",
                            "POLYCHAETES", "ECHINODERMATA", "BIVALVES", "OTHER BENTHOS", "BACTERIA", "ICE ALGAE", "PHYTOPLANKTON", "DETRITUS")) +
  theme(axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color = "black", angle = 90, vjust = 0.5, hjust=1),
        plot.margin =unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank())

tiff("system.control.resp.uncertainty.zoo.tiff", width = 6, height = 5, units = 'in', res = 500)
sys.cont.plot
dev.off()


# plotting the control difference matrix

CD.matrix <- as.data.frame(enaControl(WBB)$CD)
CD.matrix$group <- rownames(CD.matrix)

rownames(CD.matrix)[rownames(CD.matrix) == "MARINE WORMS"] <- "POLYCHAETES"
colnames(CD.matrix)[colnames(CD.matrix) == "MARINE WORMS"] <- "POLYCHAETES"

longData<-melt(CD.matrix)
longData<-longData[longData$value!=0,]

order <- rev(CD.matrix$group)


CD.plot <- ggplot(longData, aes(x = variable, y = group)) + 
  geom_tile(aes(fill=value), color = "black") + 
  scale_fill_distiller(palette = "Purples") +
  scale_y_discrete(limits=order) + 
  labs(x="", y="", title="") + 
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=90, hjust = 1, color = "black"),
                     axis.text.y=element_text(size=9, color = "black"),
                     plot.title=element_text(size=11), panel.grid = element_blank())


tiff("control.diff.colors.zoo.tiff", width = 8, height = 6, units = 'in', res = 500)
CD.plot
dev.off()



# Plotting the network
## Set the random seed to control plot output
set.seed(2)

network1 <- ggnet2(WBB, node.size = 6, node.color = "tomato3", edge.size = 0.5, edge.color = "grey", label = TRUE, label.size = 3)

ggnet2(WBB, size = nodes$biomass, label = TRUE, label.size = 3)


tiff("network1.tiff", width = 8, height = 6, units = 'in', res = 500)
network1
dev.off()

#there are a few more packages out there that I can explore for this (can also use it for the Arctic cod paper)
install.packages("igraph")
install.packages("network")
install.packages("sna")
install.packages("ndtv")

WBB%v%'output'
WBB%v%'input'
WBB%v%'living'
WBB%v%'respiration' 
WBB%v%'storage'
WBB%v%'vertex.names'
WBB%v%'export'

flows_matrix <- as.matrix(WBB,attrname="flow")
consoT<-colSums(flows_matrix)
diet_perc_WBB<-flows_matrix %*% diag(1/consoT) #convert to a proportion
cons_flows_WBB <- as.data.frame(diet_perc_WBB)
rownames(cons_flows_WBB)[rownames(cons_flows_WBB) == "MARINE WORMS"] <- "POLYCHAETES"
Species<-replace(WBB%v%'vertex.names', 23,'POLYCHAETES')
from <- rep(rownames(cons_flows_WBB),30)
to_to <- rep(rownames(cons_flows_WBB),each = 30)

TL<- c("#66C2A5", "#66C2A5", "#FC8D62", "#8DA0CB", "#FC8D62", "#FC8D62", "#8DA0CB", "#8DA0CB", "#FC8D62", "#FC8D62", "#8DA0CB", "#8DA0CB", "#8DA0CB", "#8DA0CB", "#8DA0CB", "#8DA0CB", "#E78AC3", "#8DA0CB", "#8DA0CB", "#E78AC3", "#E78AC3", "#E78AC3", "#E78AC3",
        "#E78AC3", "#E78AC3", "#E78AC3", "#E78AC3", "#FFD92F", "#FFD92F", "#FFD92F") #based on color pallet Set2, which is color friendly
brewer.pal(n = 8, name = "Set2") # "#66C2A5" "#FC8D62" "#8DA0CB" "#E78AC3" "#A6D854" "#FFD92F" "#E5C494" "#B3B3B3"

links <- cons_flows_WBB %>% gather() %>% rename(to = key, flow = value) %>% mutate(from = from, to = to_to) %>% relocate(from, to, flow)
nodes <- as.data.frame(cbind(Species = Species, biomass = WBB%v%'storage', TL = TL)) 

links <- links %>% filter(flow != 0)

net <- graph.data.frame(links, nodes, directed=T)

deg <- degree(net, mode="all")
V(net)$size <- deg/2
E(net)$width <- E(net)$flow*5

tiff("network1.tiff", width = 8, height = 8, units = 'in', res = 500)

plot(net, vertex.label.cex=.5,vertex.label.color="black", vertex.color = nodes$TL, 
     edge.color = "gray", edge.arrow.size=.2, layout = layout.auto, vertex.label.family="Arial")

legend("right", legend=c("TL 1", "TL 2","TL 3", "TL 4", "TL 5"), 
       fill=c("#FFD92F", "#E78AC3", "#8DA0CB", "#FC8D62", "#66C2A5"), cex=0.6, box.lwd = 0,box.col = "white",bg = "white")
dev.off()



