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
heat_dat = read.table("../data/pl2_fig3.txt")
colnames(heat_dat) = c("rho1", "rho2", names)

regime <- data.frame(rho1s= unique(heat_dat$rho1))
fac_labels = rep(" ", length(regime$rho1s))
fac_labels[round(seq(1,101,length.out = 7))] = c("10^-6", "10^-5", "10^-4", "10^-3", "10^-2", "10^-1", "10^0")

# rate1 <- 1.2
# rate2 <- 0.616666
# c_val_min <- -6
# c_val_max <- 3*rate2/rate1
# 
# #Returns rho 2 values for the overlapping region based on other parameters
# overlapping <- function(rho1, rate1, rate2, C_val){
#   return(exp(log(rho1)*(rate2/rate1) - 6))
# }
# regime$rho2s_min <- overlapping(regime$rho1s, rate1, rate2, c_val_min)
# regime$rho2s_max <- overlapping(regime$rho1s, rate1, rate2, c_val_max)


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

pdf("pl2_fig3.pdf", height = 12, width = 10)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1)))
print(disease_1, vp = vplayout(1, 1))
print(disease_2, vp = vplayout(2, 1))
dev.off()
