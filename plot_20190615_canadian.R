rm(list=ls())

# loading libraries
library(fda)
library(ggplot2)
library(reshape2)
library(stringr)

data(CanadianWeather)
data <- CanadianWeather

data.mlt <- melt(data)
View(data.mlt)
data.mlt <- data.mlt[, c(1, 2, 4)]
colnames(data.mlt) <- c("day", "region", "temp")

data <- dcast(data.mlt, day+region~temp)

plot(data[1, ], type="o", col=1, main="temperature of 90 regions in Korea",
     xlab="month", ylab="temperature", ylim=c(-10, 40))
for (i in 1:nrow(data)){
  lines(data[i, ], type="o", col=i)
}






