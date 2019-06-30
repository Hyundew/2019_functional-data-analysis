rm(list=ls())
getwd()
setwd("/Volumes/Hyunjoo/Korea Univ Master/1학기/이론특수통계학")
data <- read.csv("KoreanDataSheet.csv", fileEncoding = "CP949", encoding="UTF-8")
data0 <- read.csv("KoreanDataSheet_2018.csv", fileEncoding = "CP949", encoding="UTF-8")

head(data)
head(data0)

# loading libraries
library(fda)
library(ggplot2)
library(reshape2)
library(corrplot)


# 0. data preprocessing
region <- data[, 1]
regionID <- data[, 2]

data0 <- data0[, -c(1:2)]
data <- data[, -c(1:2)]
colnames(data) <- paste(2019, colnames(data)); colnames(data)
colnames(data0) <- paste(2018, colnames(data0)); colnames(data0)

data <- cbind(data0, data)
head(data)

# scale
scale0 <- function(x){(x-mean(x))/sd(x)}
range0 <- function(x){(x-min(x))/(max(x)- min(x))}
center0 <- function(x){x-mean(x)}
# data <- t(apply(data[, -c(1:2)], 1, range0)) - 차이 별로 없어 ㅜ
data <- t(apply(data, 1, center0))

corrplot(cor(data))





# 1. raw data
data.mlt <- melt(data)
View(data.mlt)
data.mlt[, 1] <- rep(region, 12*2)
str(data.mlt)
colnames(data.mlt) <- c("region", "month", "temp")
data.mlt$region <- as.factor(data.mlt$region)

## use ggplot - fd data는 못 그림ㅠ
ggplot(data.mlt, aes(x=month, y=temp, group=region, color=region)) +
  geom_line() +
  geom_point()

## use plot
plot(data[1, ], type="o", col=1, main="temperature of 90 regions in Korea",
     xlab="month", ylab="temperature", ylim=c(-2, 2))
for (i in 1:nrow(data)){
  lines(data[i, ], type="o", col=i)
}


# 2. smoothing using regression analysis
monthbasis65 <- create.fourier.basis(c(0, 24), nbasis=5, period=24)

# example
# datafd <- smooth.basis(1:12, data[1,], monthbasis65)
# plot(datafd$fd)

smooth.basis.fun <- function(x) {
  result <- smooth.basis(1:24, x, monthbasis65)
  return(result$fd)
}
data.smooth <- apply(data, 1, smooth.basis.fun)
# str(data.smooth)
plot(data.smooth[[1]], xlab="month", ylab="temperature",
     col=1, main="smooth temperature", ylim=c(-2, 2))
for (i in 2:nrow(data)) lines(data.smooth[[i]], col=i) 

# 2.1. first derivative of smoothing data

# ex
# plot(deriv.fd(data.smooth[[1]], 1))
# plot(deriv.fd(data.smooth[[1]]))

data.smooth.1 <- list()
for (i in 1:90){
  data.smooth.1[[i]] <- deriv.fd(data.smooth[[i]], 1)
}
plot(data.smooth.1[[1]], xlab="month", ylab="temperature",
     col=1, main="the first derivative curves", ylim=c(-11, 11))
for (i in 2:90) lines(data.smooth.1[[i]], col=i)

# 2.2. second derivative of smoothing data

data.smooth.2 <- list()
for (i in 1:90){
  data.smooth.2[[i]] <- deriv.fd(data.smooth[[i]], 2)
}
plot(data.smooth.2[[1]], xlab="month", ylab="temperature",
     col=1, main="the second derivative curves", ylim=c(-10, 10))
for (i in 2:90) lines(data.smooth.2[[i]], col=i)





