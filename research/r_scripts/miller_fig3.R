##############
# Spencer Fox
# Read two strain joel miller model output
# Graph to recreate figure 3
##############

rm(list=ls())
library(ggplot2)
library(dplyr)
library(reshape2)
library(gridExtra)
vplayout<-function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
setwd("~/projects/network/EpiFire/research/r_scripts/")

#####################
# Heatmaps
names <- c("theta", "pI1", "pI2", "I1", "I2", "R1", "R2")
heat_dat = read.table("../data/miller_fig3.txt")
colnames(heat_dat) = c("rho1", "rho2", names)

fac_labels = rep(" ", length(unique(heat_dat$rho1)))
fac_labels[round(seq(1,101,length.out = 7))] = c("10^-6", "10^-5", "10^-4", "10^-3", "10^-2", "10^-1", "10^0")

disease_1 <- ggplot(heat_dat, aes(as.factor(rho1), as.factor(rho2), fill = R1)) + geom_tile() +
  scale_fill_gradient(low = "blue", high = "red", name="Epidemic Size 1") + 
  scale_x_discrete(labels = fac_labels) + 
  scale_y_discrete(labels = fac_labels) +
  labs(y = "Proportion Initially Infected Disease 2", 
       x = "Proportion Initially Infected Disease 1")

disease_2 <- ggplot(heat_dat, aes(as.factor(rho1), as.factor(rho2), fill = R2)) + geom_tile() +
  scale_fill_gradient(low = "blue", high = "red", name="Epidemic Size 2") + 
  scale_x_discrete(labels = fac_labels) + 
  scale_y_discrete(labels = fac_labels) + 
  labs(y = "Proportion Initially Infected Disease 2", 
       x = "Proportion Initially Infected Disease 1")

pdf("miller_fig3.pdf", height = 12, width = 10)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1)))
print(disease_1, vp = vplayout(1, 1))
print(disease_2, vp = vplayout(2, 1))
dev.off()

