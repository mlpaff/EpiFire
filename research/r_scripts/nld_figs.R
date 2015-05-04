##############
# Spencer Fox
# Read two strain joel miller model output
# Use output from models to create figures for presentation
##############

rm(list=ls())
library(ggplot2)
library(dplyr)
library(reshape2)
setwd("~/projects/network/EpiFire/research/data/")

data <- read.table("ex_time_course2.txt")
names <- c("theta", "pI1", "pI2", "I1", "I2", "R1", "R2")
colnames(data) <- names
data$S = 1 - data$I1 - data$I2 - data$R1 - data$R2 
data$index = seq(0, nrow(data)-1)
data = data[1:50,] 

pdf("../r_scripts/time_course2.pdf", width = 10, height = 8)
melt(data,  id.vars=c("index"), measure.vars = c("S", "I1", "I2", "R1", "R2")) %>%
  ggplot(., aes(x=index, y = value, color = variable)) + geom_line(size = 2) + 
  labs(x="Time Step", y="Proportion of Population")+
  geom_hline(y=1)+
  theme_bw() %+replace% theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                              axis.text.x = element_text(size = 24), axis.text.y = element_text(size = 24, angle=90),
                              axis.title.x = element_text(size=24), axis.title.y=element_text(size=24, angle = 90),
                              legend.position = c(0.8, 0.7), legend.text = element_text(size=20),
                              legend.title = element_blank(), legend.key.size = unit(1.5, unit="cm"))
dev.off()
