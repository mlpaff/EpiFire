##############
# Spencer Fox
# Read two strain joel miller model output
# Use output from models to create figures for presentation
##############

rm(list=ls())
library(ggplot2)
library(dplyr)
library(reshape2)
setwd("~/projects/network/EpiFire/R scripts/")

data <- read.table("../examples/test.txt")
names <- c("theta", "pI1", "pI2", "I1", "I2", "R1", "R2")
colnames(data) <- names
data$S = 1 - data$I1 - data$I2 - data$R1 - data$R2 
data$index = seq(0, nrow(data)-1)
melt(data,  id.vars=c("index"), measure.vars = c("S", "I1", "I2", "R1", "R2")) %>%
  ggplot(., aes(x=index, y = value, color = variable)) + geom_line() + theme_bw() +
  geom_hline(y=1)
