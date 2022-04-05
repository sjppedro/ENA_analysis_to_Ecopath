##############################################################
# ENA of Western Baffin Bay model of 2016  ###################
##############################################################

library(devtools)

#To install R tools I followed the instructions in here: https://cran.r-project.org/bin/windows/Rtools/

writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
Sys.which("make")
install.packages("jsonlite", type = "source") #this worked!

url <- "https://cran.r-project.org/src/contrib/Archive/enaR/enaR_3.0.0.tar.gz"
pkgFile <- "enaR_3.0.0.tar.gz"
download.file(url = url, destfile = pkgFile)
install.packages(pkgs=pkgFile, type="source", repos=NULL)
#Rtools may not have been successfully installed

library(sna) #sna and gdata are dependencies of enaR - the package cannot be installed without these other packages
library(gdata)
library(enaR)
library(limSolve)
library(network)
#install.packages("GGally")
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

#How to load in model - need to change Ecopath model into the correct format to be read in R - Bentley et al 2019 used SCOR format
#Bentley reformatted Ecopath into:
# (1) a flow matrix (Ecopath consumption matrix) 
# (2) network inputs (gross primary production in the Irish Sea) 
# (3) network exports (detritus and fisheries exports) 
# (4) respiration (calculated in EwE outputs)  
# (5) node storage (Ecopath biomass)

#As Ecopath does not provide estimates of gross primary productivity or respiration for primary producers (Heymans and Baird, 2000),
# Bentley used methods for calculating respiration and gross production for aquatic plants from a model equation proposed by Aoki (2006). 
# Gross production a() is consumed by respiration r(), production p(), and flow to detritus d( ): 
#a = r + p + d
#In accordance with Aoki (2006), they assumed that:
#r = a.rho 
#where rho is the annual respiration-gross production ratio. According to the equations:
#r = (p + d) . rho /(1-rho)
#r can therefore be estimated from values of p, d and rho, where p and d are available from Ecopath. 
#Values of rho for phytoplankton have been estimated to be 0.42 y-1 for Georges Bank (Riley, 1946), Narragansett
#Bay, Delaware Bay and Chesapeake Bay (Monaco and Ulanowicz, 1997) and 0.44 y-1 for Lake Biwa, Lake Yunoko, 
#Lake Suwa and Lake Kojima (Mori and Yamamoto, 1975). For the purpose of this model, as previously adopted by Aoki (2006), 
# Bentley used the average value of 0.4 for phytoplankton and 0.65 was for seaweed Aoki (2006).

# For my models (in the SCOR format file) I added diet import to the imports values (along with gross primary production - a)
# For gross primary production and respiration of primary producers: according to Saint-B?at et al 2019 the gross primary production is 20.4 gC/m2/month, but 
# this does not distinguish between sea ice algae and phytoplankton. Respiration was 3.15 gC/m2/month (again not distinguishing sea ice algae from phytoplankton).
# In Ecopath in ton/km2/year net production of phytoplankton is 530.4 (from GreenEdge data) and flow to detritus is 276.2 (from consumption flow matrix); 
# for sea ice algae production is 432.06 and flow to detritus is 144.1. 
# Saint-B?at et al considered respiration between 0.05*a and 0.3*a, where 0.05 and 0.3 is rho according to Aoki (2006). Based on overall gross primary production
# and respiration for western Baffin Bay I can find the rho used in this study and then from there calculate gross primary production for sea ice algae and phytoplakton
# in my models:

rho <- 3.15/20.4

phy_r <- (530.4 + 276.2)*(rho/(1-rho))
SIA_r <- (432.06 + 144.1)*(rho/(1-rho))

phy_a <- phy_r + 530.4 + 276.2
SIA_a <- SIA_r + 432.06 + 144.1


# Load data

WBB <- read.scor("C:/Users/Sara/Google Drive/PostDoc-QC/Data gathering/Models/Remarks_Model_construction/Model12/ENA/WBB_model12_SCOR.txt", from.file = TRUE, warn = TRUE)

#the file (model) below considers gross primary production (input in the system) as net primary production (biomasses of primary producers) and respiration of primary producers as zero
WBB.zero.r.pp <- read.scor("C:/Users/Sara/Google Drive/PostDoc-QC/Data gathering/Models/Remarks_Model_construction/Model12/ENA/WBB_model12_SCOR.zero.r.pp.txt", from.file = TRUE, warn = TRUE)

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
plot(WBB) #need to figure out a better package for data visualization

## Check to see if the model is balanced
ssCheck(WBB)  # not balanced...
WBB <- balance(WBB, method="AVG2")
WBB <- force.balance(WBB) #this forces the model to balance - just using here quickly to be able to generate some ENA results

## for the alternative model
ssCheck(WBB.zero.r.pp)  # not balanced...
WBB.zero.r.pp <-  balance(WBB.zero.r.pp, method="AVG2")
WBB.zero.r.pp <- force.balance(WBB.zero.r.pp) #this forces the model to balance - just using here quickly to be able to generate some ENA results

# Bentley et al used the following flow-based network statistics: total system through flow (TST), Finn's cycling index (FCI),
# indirect flow intensity (IFI) and average path length (APL). As control metrics, they used control difference (CD) and system 
# control; I could also look into keystoneness index and mixed trophic effects, ascendency, overhead or relative ascendency, and 
# Lindeman spine.

## ENA results

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


#It is possible to instantly return all of the main ENA output with enaAll:
WBB.ena <- enaAll(WBB)
names(WBB)


##############################################################################################
## Uncertainty analysis of ENA results according to Hines et al 2018 & Bentley et al 2019 ####
##############################################################################################

# Once a network model is constructed, uncertainty data for each edge weight can be combined with the original model through a
# LIM uncertainty analysis. This uncertainty analysis can be performed by the enaUncertainty() function. This function first uses
# the original network model, along with uncertainty data tables, to generate the inputs needed for LIM. It then applies a LIM 
# analysis to generate plausible model sets using the limSolve package for R. Next, LIM results are converted back into a format
# readable by enaR and are returned to the user as a set of plausible model parameterizations for use in ENA (Hines et al).

iter = 10000 # the number of plausible models to return 

###############################################################################################
#percent method with different uncertainty levels for all edges - #############################
###############################################################################################

plausible.per.25 <- enaUncertainty(x=WBB, type="percent", p.err=25, iter=iter)
plausible.per.50 <- enaUncertainty(x=WBB, type="percent", p.err=50, iter=iter)

##############################################################################################################################
#percent method with different symmetric uncertainties based on the pedigree in Ecopath and reported respiration flows #######
##############################################################################################################################

# Load the uncertainty values of the allowable amount above and below each edge

f.df <- read.csv("C:/Users/Sara/Google Drive/PostDoc-QC/Data gathering/Models/Remarks_Model_construction/Model12/ENA/flows.f.csv")
z.df <- read.csv("C:/Users/Sara/Google Drive/PostDoc-QC/Data gathering/Models/Remarks_Model_construction/Model12/ENA/inputs.z.csv")
r.df <- read.csv("C:/Users/Sara/Google Drive/PostDoc-QC/Data gathering/Models/Remarks_Model_construction/Model12/ENA/respiration.r.2.csv") # CV from literature on respiration rates for some species (see Excel)
e.df <- read.csv("C:/Users/Sara/Google Drive/PostDoc-QC/Data gathering/Models/Remarks_Model_construction/Model12/ENA/exports.e.csv")
r.25.df <- read.csv("C:/Users/Sara/Google Drive/PostDoc-QC/Data gathering/Models/Remarks_Model_construction/Model12/ENA/respiration.r.25.csv") #applying 25 % CV to all values
r.50.df <- read.csv("C:/Users/Sara/Google Drive/PostDoc-QC/Data gathering/Models/Remarks_Model_construction/Model12/ENA/respiration.r.50.csv")
r.0.df <- read.csv("C:/Users/Sara/Google Drive/PostDoc-QC/Data gathering/Models/Remarks_Model_construction/Model12/ENA/respiration.r.0.csv")
r.10.df <- read.csv("C:/Users/Sara/Google Drive/PostDoc-QC/Data gathering/Models/Remarks_Model_construction/Model12/ENA/respiration.r.10.csv") #applying 25 % CV to all values

r.df.zero.pp <- read.csv("C:/Users/Sara/Google Drive/PostDoc-QC/Data gathering/Models/Remarks_Model_construction/Model12/ENA/respiration.r.zero.pp.csv") # CV from literature on respiration rates for some species (see Excel)
f.df.respCVs <- read.csv("C:/Users/Sara/Google Drive/PostDoc-QC/Data gathering/Models/Remarks_Model_construction/Model12/ENA/flows.f.respCVs.csv")
z.df.respCVs <- read.csv("C:/Users/Sara/Google Drive/PostDoc-QC/Data gathering/Models/Remarks_Model_construction/Model12/ENA/inputs.z.respCVs.csv")

r.df.zoo <- read.csv("C:/Users/Sara/Google Drive/PostDoc-QC/Data gathering/Models/Remarks_Model_construction/Model12/ENA/respiration.r.zoo.csv") # CV from literature on respiration rates for some species (see Excel)
f.df.respCVs.zoo <- read.csv("C:/Users/Sara/Google Drive/PostDoc-QC/Data gathering/Models/Remarks_Model_construction/Model12/ENA/flows.f.respCVs.zoo.csv")


#Comparing different uncertainty levels for respiration rates (not in the Ecopath pedigree)
plausible.sym.10r <- enaUncertainty(x=WBB, F.sym=f.df, z.sym=z.df, e.sym=e.df, r.sym=r.10.df, type="sym", iter=iter)
plausible.sym.25r <- enaUncertainty(x=WBB, F.sym=f.df, z.sym=z.df, e.sym=e.df, r.sym=r.25.df, type="sym", iter=iter)
plausible.sym.50r <- enaUncertainty(x=WBB, F.sym=f.df, z.sym=z.df, e.sym=e.df, r.sym=r.50.df, type="sym", iter=iter)
plausible.sym.nor <- enaUncertainty(x=WBB, F.sym=f.df, z.sym=z.df, e.sym=e.df, r.sym=r.0.df, type="sym", iter=iter)

plausible.sym <- enaUncertainty(x=WBB, F.sym=f.df, z.sym=z.df, e.sym=e.df, r.sym=r.df, type="sym", iter=iter)


#Two other alternatives

plausible.sym.zero.pp <- enaUncertainty(x=WBB.zero.r.pp, F.sym=f.df, z.sym=z.df, e.sym=e.df, r.sym=r.df.zero.pp, type="sym", iter=iter)
# primary producers with respiration equal to zero

plausible.sym.all.r.uncertainty <- enaUncertainty(x=WBB, F.sym=f.df.respCVs, z.sym=z.df.respCVs, e.sym=e.df, r.sym=r.df, type="sym", iter=iter)
# previously defined respiration flows and respective uncertainty levels applied to all the flows 

plausible.sym.all.r.uncertainty.zoo <- enaUncertainty(x=WBB, F.sym=f.df.respCVs.zoo, z.sym=z.df.respCVs, e.sym=e.df, r.sym=r.df.zoo, type="sym", iter=iter)
# previously defined respiration flows and respective uncertainty levels applied to all the flows, but with zoo uncertainty levels based on Saint-B?at work 


#Plotting to compare the values using boxplots
# original model
TST.WBB <- enaFlow(WBB)$ns[2]
FCI.WBB <- enaFlow(WBB)$ns[5]

# different uncertainty methods
TST.per.25 <- unlist(lapply(plausible.per.25, function(x) enaFlow(x)$ns[2]))
TST.per.50 <- unlist(lapply(plausible.per.50, function(x) enaFlow(x)$ns[2]))
TST.per.pedigree <- unlist(lapply(plausible.sym, function(x) enaFlow(x)$ns[2]))
TST.per.pedigree.10.r <- unlist(lapply(plausible.sym.10r, function(x) enaFlow(x)$ns[2]))
TST.per.pedigree.25.r <- unlist(lapply(plausible.sym.25r, function(x) enaFlow(x)$ns[2]))
TST.per.pedigree.50.r <- unlist(lapply(plausible.sym.50r, function(x) enaFlow(x)$ns[2]))
TST.per.pedigree.no.r <- unlist(lapply(plausible.sym.nor, function(x) enaFlow(x)$ns[2]))

FCI.per.25 <- unlist(lapply(plausible.per.25, function(x) enaFlow(x)$ns[5]))
FCI.per.50 <- unlist(lapply(plausible.per.50, function(x) enaFlow(x)$ns[5]))
FCI.per.pedigree <- unlist(lapply(plausible.sym, function(x) enaFlow(x)$ns[5]))
FCI.per.pedigree.10.r <- unlist(lapply(plausible.sym.10r, function(x) enaFlow(x)$ns[5]))
FCI.per.pedigree.25.r <- unlist(lapply(plausible.sym.25r, function(x) enaFlow(x)$ns[5]))
FCI.per.pedigree.50.r <- unlist(lapply(plausible.sym.50r, function(x) enaFlow(x)$ns[5]))


TST.per.pedigree.zero.resp.pp <- unlist(lapply(plausible.sym.zero.pp, function(x) enaFlow(x)$ns[2]))
FCI.per.pedigree.zero.resp.pp <- unlist(lapply(plausible.sym.zero.pp, function(x) enaFlow(x)$ns[5]))

TST.per.pedigree.all.resp <- unlist(lapply(plausible.sym.all.r.uncertainty, function(x) enaFlow(x)$ns[2]))
FCI.per.pedigree.all.resp <- unlist(lapply(plausible.sym.all.r.uncertainty, function(x) enaFlow(x)$ns[5]))


ENA.ind.test.5 <- as.data.frame(cbind(TST.per.pedigree.all.resp, FCI.per.pedigree.all.resp))

write.csv(ENA.ind.test.5, file = "ENA_indicators_10000per_comparisons_5.csv")

TST.comparisons <- read.csv("C:/Users/Sara/Google Drive/PostDoc-QC/Data gathering/Models/Remarks_Model_construction/Model12/ENA/ENA_indicators_10000per_comparisons_TST.csv")
FCI.comparisons <- read.csv("C:/Users/Sara/Google Drive/PostDoc-QC/Data gathering/Models/Remarks_Model_construction/Model12/ENA/ENA_indicators_10000per_comparisons_FCI.csv")

# compare the methods with boxplots

# comparison of several methods for all flows

tiff("Different_uncertainty_levels_TST.tiff", width = 8, height = 6, units = 'in', res = 500)

colnames(TST.comparisons) <- c("25%", "50%", "Resp pp = 0", "Resp flow unc" ,"Pedigree")
par(las=1)
bp=boxplot(TST.comparisons, main="TST")
points(x=c(1,2,3,4,5), y=rep(TST.Ecopath, 5), col="red", pch=16, cex=1)
legend("bottomleft", legend=c("original model"), pch=16, col="red", bty="n", cex=0.75)

dev.off()


tiff("Different_uncertainty_levels_FCI.tiff", width = 8, height = 6, units = 'in', res = 500)

colnames(FCI.comparisons) <- c("25%", "50%", "Resp pp = 0", "Resp flow unc" ,"Pedigree")
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

# compare the methods with density plots

a <- as.data.frame(TST.plot.data)

ggplot(a, aes(x=TST.per.pedigree)) + 
  geom_density(alpha = 0.4, fill = "tomato3") + 
  geom_vline(aes(xintercept=TST.WBB), color="black", linetype="dashed", size=1) +
  theme_bw() + xlim(4200, 4650) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color = "black"),
        plot.margin =unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank())

ggplot(a, aes(x=TST.per.pedigree.25.r)) + 
  geom_density(alpha = 0.4, fill = "tomato3") + 
  geom_vline(aes(xintercept=TST.WBB), color="black", linetype="dashed", size=1) +
  theme_bw() + xlim(4200, 4650) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color = "black"),
        plot.margin =unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank())

# extract other indicators from ena results with uncertainty from pedigree

APL.per.pedigree <- unlist(lapply(plausible.sym, function(x) enaFlow(x)$ns[4])) #average path length
FCI.per.pedigree <- unlist(lapply(plausible.sym, function(x) enaFlow(x)$ns[5])) #Fin's Cycling Index
IFI.per.pedigree <- unlist(lapply(plausible.sym, function(x) enaFlow(x)$ns[8])) #Indirect flow intensity
asc.per.pedigree <- unlist(lapply(plausible.sym, function(x) enaAscendency(x)[[5]])) #ascendency
rel.asc.per.pedigree <- unlist(lapply(plausible.sym, function(x) enaAscendency(x)[[7]])) #relative ascendency
AMI.per.pedigree <- unlist(lapply(plausible.sym, function(x) enaAscendency(x)[[2]])) #average mutual information
TST.per.pedigree <- unlist(lapply(plausible.sym, function(x) enaFlow(WBB)$ns[2])) #total system control

# extract other indicators from ena results with uncertainty from respiration flows from literature

APL.per.resp <- unlist(lapply(plausible.sym.all.r.uncertainty, function(x) enaFlow(x)$ns[4])) #average path length
FCI.per.resp <- unlist(lapply(plausible.sym.all.r.uncertainty, function(x) enaFlow(x)$ns[5])) #Fin's Cycling Index
IFI.per.resp <- unlist(lapply(plausible.sym.all.r.uncertainty, function(x) enaFlow(x)$ns[8])) #Indirect flow intensity
asc.per.resp <- unlist(lapply(plausible.sym.all.r.uncertainty, function(x) enaAscendency(x)[[5]])) #ascendency
rel.asc.per.resp <- unlist(lapply(plausible.sym.all.r.uncertainty, function(x) enaAscendency(x)[[7]])) #relative ascendency
AMI.per.resp <- unlist(lapply(plausible.sym.all.r.uncertainty, function(x) enaAscendency(x)[[2]])) #average mutual information
TST.per.resp <- unlist(lapply(plausible.sym.all.r.uncertainty, function(x) enaFlow(x)$ns[2])) #total system control

# extract other indicators from ena results with uncertainty from respiration flows from literature and Saint-B?at for zoo

APL.per.resp.zoo <- unlist(lapply(plausible.sym.all.r.uncertainty.zoo, function(x) enaFlow(x)$ns[4])) #average path length
FCI.per.resp.zoo <- unlist(lapply(plausible.sym.all.r.uncertainty.zoo, function(x) enaFlow(x)$ns[5])) #Fin's Cycling Index
IFI.per.resp.zoo <- unlist(lapply(plausible.sym.all.r.uncertainty.zoo, function(x) enaFlow(x)$ns[8])) #Indirect flow intensity
asc.per.resp.zoo <- unlist(lapply(plausible.sym.all.r.uncertainty.zoo, function(x) enaAscendency(x)[[5]])) #ascendency
rel.asc.per.resp.zoo <- unlist(lapply(plausible.sym.all.r.uncertainty.zoo, function(x) enaAscendency(x)[[7]])) #relative ascendency
AMI.per.resp.zoo <- unlist(lapply(plausible.sym.all.r.uncertainty.zoo, function(x) enaAscendency(x)[[2]])) #average mutual information
TST.per.resp.zoo <- unlist(lapply(plausible.sym.all.r.uncertainty.zoo, function(x) enaFlow(x)$ns[2])) #total system control
robustness.per.resp.zoo <- unlist(lapply(plausible.sym.all.r.uncertainty.zoo, function(x) enaAscendency(x)[[9]])) #robustness


# average values from ENAr model (based on Ecopath converted model and without uncertainty)
APL.WBB <- enaFlow(WBB)$ns[4]
FCI.WBB <- enaFlow(WBB)$ns[5]
IFI.WBB <- enaFlow(WBB)$ns[8]
LD.WBB <- enaStructure(WBB)$ns[4]
asc.WBB <- enaAscendency(WBB)[[5]]
rel.asc.WBB <- enaAscendency(WBB)[[7]]
AMI.WBB <- enaAscendency(WBB)[[2]]
TST.WBB <- enaFlow(WBB)$ns[2]
robustness.WBB <- enaAscendency(WBB)[[9]]

# values from Ecopath
TST.Ecopath <- 3536
FCI.Ecopath <- 0.1277
APL.Ecopath <- 3.67
asc.Ecopath <- 4105
rel.asc.Ecopath <- 0.226
AMI.Ecopath <- 1.161 #called information (in bits) in Ecopath
#no idea where to get IFI in Ecopath
  
# density plot
# save a table with the permutations for all these indices
ENA.indicators <- as.data.frame(cbind(AMI.per.pedigree, APL.per.pedigree, TST.per.pedigree, FCI.per.pedigree, asc.per.pedigree, rel.asc.per.pedigree,
                                      IFI.per.pedigree))

write.csv(ENA.indicators, file = "ENA_indicators_10000per_2.csv")


ENA.indicators.resp <- as.data.frame(cbind(AMI.per.resp, APL.per.resp, TST.per.resp, FCI.per.resp, asc.per.resp, rel.asc.per.resp,
                                      IFI.per.resp))

write.csv(ENA.indicators.resp, file = "ENA_indicators_10000per_resp_uncertainty.csv")


ENA.indicators.resp.zoo <- as.data.frame(cbind(AMI.per.resp.zoo, APL.per.resp.zoo, TST.per.resp.zoo, FCI.per.resp.zoo, asc.per.resp.zoo, rel.asc.per.resp.zoo,
                                           IFI.per.resp.zoo, robustness.per.resp.zoo))

write.csv(ENA.indicators.resp.zoo, file = "ENA_indicators_10000per_resp_uncertainty.zoo.csv")


AMI.plot <- ggplot(ENA.indicators.resp.zoo, aes(x=AMI.per.resp.zoo)) + 
            geom_density(alpha = 0.4, fill = "lightgreen", color = "darkgreen") + 
            geom_vline(aes(xintercept=mean(AMI.per.resp.zoo)), color="darkgreen", linetype="dashed", size=1) +
            geom_vline(aes(xintercept=AMI.WBB), color="black", linetype="dotted", size=1) +
            xlab("Average Mutual Information (AMI)") + ylab("Density") +
            theme_bw() + #xlim(4200, 4650) +
            scale_y_continuous(expand = c(0,0)) +
            theme(axis.text.y = element_text(color="black"),
               axis.text.x = element_text(color = "black"),
                 plot.margin =unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                  panel.background = element_blank())


APL.plot <- ggplot(ENA.indicators.resp.zoo, aes(x=APL.per.resp.zoo)) + 
  geom_density(alpha = 0.4, fill = "coral", color = "coral3") + 
  geom_vline(aes(xintercept=mean(APL.per.resp.zoo)), color="coral3", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=APL.WBB), color="black", linetype="dotted", size=1) +
  xlab("Average Path Length (APL)") + ylab("") +
  theme_bw() + #xlim(4200, 4650) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color = "black"),
        plot.margin =unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank())


FCI.plot <- ggplot(ENA.indicators.resp.zoo, aes(x=FCI.per.resp.zoo)) + 
  geom_density(alpha = 0.4, fill = "coral", color = "coral3") + 
  geom_vline(aes(xintercept=mean(FCI.per.resp.zoo)), color="coral3", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=FCI.WBB), color="black", linetype="dotted", size=1) +
  xlab("Finn's Cycling Index (FCI)") + ylab("") +
  theme_bw() + #xlim(4200, 4650) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color = "black"),
        plot.margin =unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank())

TST.plot <- ggplot(ENA.indicators.resp.zoo, aes(x=TST.per.resp.zoo)) + 
  geom_density(alpha = 0.4, fill = "lightgreen", color = "darkgreen") + 
  geom_vline(aes(xintercept=mean(TST.per.resp.zoo)), color="darkgreen", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=TST.WBB), color="black", linetype="dotted", size=1) +
  xlab("Total System Throughput (TST)") + ylab("Density") +
  theme_bw() + #xlim(4200, 4650) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color = "black"),
        plot.margin =unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank()) 

asc.plot <- ggplot(ENA.indicators.resp.zoo, aes(x=asc.per.resp.zoo)) + 
  geom_density(alpha = 0.4, fill = "lightgreen", color="darkgreen") + 
  geom_vline(aes(xintercept=mean(asc.per.resp.zoo)), color="darkgreen", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=asc.WBB), color="black", linetype="dotted", size=1) +
  xlab("Ascendency (A)") + ylab("Density") +
  theme_bw() + #xlim(4200, 4650) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color = "black"),
        plot.margin =unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank())

rel.asc.plot <- ggplot(ENA.indicators.resp.zoo, aes(x=rel.asc.per.resp.zoo)) + 
  geom_density(alpha = 0.4, fill = "lightgreen", color="darkgreen") + 
  geom_vline(aes(xintercept=mean(rel.asc.per.resp.zoo)), color="darkgreen", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=rel.asc.WBB), color="black", linetype="dotted", size=1) +
  xlab("Relative Ascendency (A/C)") + ylab("Density") +
  theme_bw() + #xlim(4200, 4650) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color = "black"),
        plot.margin =unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank())

IFI.plot <- ggplot(ENA.indicators.resp.zoo, aes(x=IFI.per.resp.zoo)) + 
  geom_density(alpha = 0.4, fill = "coral", color = "coral3") + 
  geom_vline(aes(xintercept=mean(IFI.per.resp.zoo)), color="coral3", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=IFI.WBB), color="black", linetype="dotted", size=1) +
  xlab("Indirect Flow Intensity (IFI)") + ylab("") +
  theme_bw() + #xlim(4200, 4650) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color = "black"),
        plot.margin =unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank())

rob.plot <- ggplot(ENA.indicators.resp.zoo, aes(x=robustness.per.resp.zoo)) + 
  geom_density(alpha = 0.4, fill = "lightgreen", color="darkgreen") + 
  geom_vline(aes(xintercept=mean(robustness.per.resp.zoo)), color="darkgreen", linetype="dashed", size=1) +
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

tiff("Density_plots_based_on_resp_uncertainty_1_zoo.tiff", width = 4, height = 8, units = 'in', res = 500)

panel.plots1 <- rbind(ggplotGrob(rob.plot), ggplotGrob(FCI.plot), ggplotGrob(IFI.plot), ggplotGrob(APL.plot), size = "first")
panel.plots1$widths <- unit.pmax(ggplotGrob(FCI.plot)$widths, ggplotGrob(IFI.plot)$widths, ggplotGrob(APL.plot)$widths)
grid.newpage()
grid.draw(panel.plots1)

dev.off()

tiff("Density_plots_based_on_resp_uncertainty_2_zoo.tiff", width = 4, height = 8, units = 'in', res = 500)

panel.plots2 <- rbind(ggplotGrob(TST.plot), ggplotGrob(AMI.plot), ggplotGrob(asc.plot), ggplotGrob(rel.asc.plot), size = "first")
panel.plots2$widths <- unit.pmax(ggplotGrob(TST.plot)$widths, ggplotGrob(AMI.plot)$widths, ggplotGrob(asc.plot)$widths, 
                                ggplotGrob(rel.asc.plot)$widths)
grid.newpage()
grid.draw(panel.plots2)

dev.off()




# max and min values

min(ENA.indicators$AMI.per.pedigree)
max(ENA.indicators$AMI.per.pedigree)

min(ENA.indicators$APL.per.pedigree)
max(ENA.indicators$APL.per.pedigree)

min(ENA.indicators$FCI.per.pedigree)
max(ENA.indicators$FCI.per.pedigree)

min(ENA.indicators$IFI.per.pedigree)
max(ENA.indicators$IFI.per.pedigree)

min(ENA.indicators$asc.per.pedigree)
max(ENA.indicators$asc.per.pedigree)

min(ENA.indicators$rel.asc.per.pedigree)
max(ENA.indicators$rel.asc.per.pedigree)

min(ENA.indicators$TST.per.pedigree)
max(ENA.indicators$TST.per.pedigree)

# from respiration uncertainty

min(ENA.indicators.resp$AMI.per.resp)
max(ENA.indicators.resp$AMI.per.resp)

min(ENA.indicators.resp$APL.per.resp)
max(ENA.indicators.resp$APL.per.resp)

min(ENA.indicators.resp$FCI.per.resp)
max(ENA.indicators.resp$FCI.per.resp)

min(ENA.indicators.resp$IFI.per.resp)
max(ENA.indicators.resp$IFI.per.resp)

min(ENA.indicators.resp$TST.per.resp)
max(ENA.indicators.resp$TST.per.resp)

min(ENA.indicators.resp$asc.per.resp)
max(ENA.indicators.resp$asc.per.resp)

min(ENA.indicators.resp$rel.asc.per.resp)
max(ENA.indicators.resp$rel.asc.per.resp)


# from respiration uncertainty and zoo from Saint-B?at

min(ENA.indicators.resp.zoo$AMI.per.resp.zoo)
max(ENA.indicators.resp.zoo$AMI.per.resp.zoo)
mean.AMI<-mean(ENA.indicators.resp.zoo$AMI.per.resp.zoo)

min(ENA.indicators.resp.zoo$APL.per.resp.zoo)
max(ENA.indicators.resp.zoo$APL.per.resp.zoo)
mean.APL<-mean(ENA.indicators.resp.zoo$APL.per.resp.zoo)

min(ENA.indicators.resp.zoo$FCI.per.resp.zoo)
max(ENA.indicators.resp.zoo$FCI.per.resp.zoo)
mean.FCI<-mean(ENA.indicators.resp.zoo$FCI.per.resp.zoo)


min(ENA.indicators.resp.zoo$IFI.per.resp.zoo)
max(ENA.indicators.resp.zoo$IFI.per.resp.zoo)
mean.IFI<-mean(ENA.indicators.resp.zoo$IFI.per.resp.zoo)


min(ENA.indicators.resp.zoo$TST.per.resp.zoo)
max(ENA.indicators.resp.zoo$TST.per.resp.zoo)
mean.TST<-mean(ENA.indicators.resp.zoo$TST.per.resp.zoo)

min(ENA.indicators.resp.zoo$asc.per.resp.zoo)
max(ENA.indicators.resp.zoo$asc.per.resp.zoo)
mean.asc<-mean(ENA.indicators.resp.zoo$asc.per.resp.zoo)


min(ENA.indicators.resp.zoo$rel.asc.per.resp.zoo)
max(ENA.indicators.resp.zoo$rel.asc.per.resp.zoo)
mean.rel.asc<-mean(ENA.indicators.resp.zoo$rel.asc.per.resp.zoo)


min(ENA.indicators.resp.zoo$robustness.per.resp.zoo)
max(ENA.indicators.resp.zoo$robustness.per.resp.zoo)
mean.rob<-mean(ENA.indicators.resp.zoo$robustness.per.resp.zoo)


# Variation between original model and uncertainty based ENA indices

APL.var <- (APL.WBB-mean.APL)/(mean(APL.WBB,mean.APL))*100

TST.var <- (TST.WBB-mean.TST)/(mean(TST.WBB,mean.TST))*100

FCI.var <- (FCI.WBB-mean.FCI)/(mean(FCI.WBB,mean.FCI))*100

AMI.var <- (AMI.WBB-mean.AMI)/(mean(AMI.WBB,mean.AMI))*100

IFI.var <- (IFI.WBB-mean.IFI)/(mean(IFI.WBB,mean.IFI))*100

asc.var <- (asc.WBB-mean.asc)/(mean(asc.WBB,mean.asc))*100

rel.asc.var <- (rel.asc.WBB-mean.rel.asc)/(mean(rel.asc.WBB,mean.rel.asc))*100

rob.var <- (robustness.WBB-mean.rob)/(mean(robustness.WBB,mean.rob))*100

# System control

CD.per.pedigree <- unlist(lapply(plausible.sym, function(x) enaControl(x)$CD)) #Schramski Control Difference Matrix
sc.per.pedigree <- unlist(lapply(plausible.sym, function(x) enaControl(x)$sc)) #System control

CD.per.resp <- unlist(lapply(plausible.sym.all.r.uncertainty, function(x) enaControl(x)$CD)) #Schramski Control Difference Matrix
sc.per.resp <- unlist(lapply(plausible.sym.all.r.uncertainty, function(x) enaControl(x)$sc)) #System control

CD.per.resp.zoo <- unlist(lapply(plausible.sym.all.r.uncertainty.zoo, function(x) enaControl(x)$CD)) #Schramski Control Difference Matrix
sc.per.resp.zoo <- unlist(lapply(plausible.sym.all.r.uncertainty.zoo, function(x) enaControl(x)$sc)) #System control


control.difference.mattrix <- as.data.frame(CD.per.pedigree)
write.csv(control.difference.mattrix, file = "control.difference.10000per.csv") 

control.difference.mattrix.resp <- as.data.frame(CD.per.resp)
write.csv(control.difference.mattrix.resp, file = "control.difference.10000per.resp.uncertainty.csv")  

control.difference.mattrix.resp <- as.data.frame(CD.per.resp.zoo)
write.csv(control.difference.mattrix.resp, file = "control.difference.10000per.resp.uncertainty.zoo.csv")  


f.group <- rep(c("KILLER WHALE", "POLAR BEAR", "NARWHAL", "BOWHEAD WHALE", "RINGED SEAL", "OTHER SEALS", "WALRUS", "SEABIRDS", "GREENLAND SHARK", "GREENLAND HALIBUT", "ARCTIC CHAR", "ARCTIC/POLAR COD",
             "SMALL PELAGIC FISH", "SCULPINS/EELPOUTS", "SMALL DEMERSAL FISH", "LARGE DEMERSAL FISH", "LARGE CRUSTACEANS", "CEPHALOPODS", "CARNIVOROUS ZOO", "OMNIVOROUS ZOO", "CALANUS ZOO", "MICROZOOPLANKTON",
             "POLYCHAETES", "ECHINODERMATA", "BIVALVES", "OTHER BENTHOS", "BACTERIA", "ICE ALGAE", "PHYTOPLANKTON", "DETRITUS"), times = 10000)

system.control <- data.frame(f.group, sc.per.pedigree)
write.csv(system.control, file = "system.control.10000per.csv") #not sure how the matrix is organized when doing the permutations 

system.control.resp <- data.frame(f.group, sc.per.resp)
write.csv(system.control.resp, file = "system.control.10000per.resp.uncertainty.csv") #not sure how the matrix is organized when doing the permutations 

system.control.resp.zoo <- data.frame(f.group, sc.per.resp.zoo)
write.csv(system.control.resp.zoo, file = "system.control.10000per.resp.uncertainty.zoo.csv") #not sure how the matrix is organized when doing the permutations 


system.control.resp.zoo.open <- read.csv("C:/Users/Sara/Google Drive/PostDoc-QC/Data gathering/Models/Remarks_Model_construction/Model12/ENA/system.control.10000per.resp.uncertainty.zoo.csv")


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



