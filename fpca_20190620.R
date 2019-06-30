rm(list=ls())
getwd()

# loading packages
library(fda)
library(ggplot2)
library(ggmap)
library(reshape2)
library(dplyr)
library(ggrepel)
library(data.table)
library(fdapace)
library(caret)
library(dbscan)



# 0. data preprocessing
data <- read.csv("KoreanDataSheet_지점추가.csv",
                 fileEncoding = "CP949", encoding = "UTF-8") # 2018 data

region <- data$Region
rownames(data) <- region; data <- data[,-1]
head(data)

# scaling
range0 <- function(x){(x-min(x))/(max(x)- min(x))}
center0 <- function(x) x - mean(x) 

ctr_data <- t(apply(data[,-1], 1, center0))
nor_data <- t(apply(data[,-1], 1, scale))



# 1. Plotting
# 1.1 Raw data
data.mlt <- melt(data[,-1]) # remain only month(variable) and value(temp)
head(data.mlt)
data.mlt <- cbind(region = rep(region, 12), data.mlt)
colnames(data.mlt) <- c("region", "month", "temp")
data.mlt$region <- as.factor(data.mlt$region)
head(data.mlt)

mean.dat <- data.mlt %>% group_by(month) %>% summarise(m = mean(temp))
with(mean.dat, month <- as.character(month))

ggplot(data.mlt, aes(x = month, y = temp, group = region, color = region)) +
  geom_line() + geom_point(size = 2) + ggtitle("Raw Data") + theme_bw() + 
  geom_line(data = mean.dat, aes(x = month, y = m, group = 1), inherit.aes = FALSE, size = 1.3)

# 1.2 Normalized data
data.mlt <- melt(nor_data)
head(data.mlt)
colnames(data.mlt) <- c("region", "month", "temp")
data.mlt$region <- as.factor(data.mlt$region)
head(data.mlt)

mean.dat <- data.mlt %>% group_by(month) %>% summarise(m = mean(temp))
with(mean.dat, month <- as.character(month))

ggplot(data.mlt, aes(x = month, y = temp, group = region, color = region)) +
  geom_line() + geom_point(size = 2) + ggtitle("Scaled Data") + theme_bw() + 
  geom_line(data = mean.dat, aes(x = month, y = m, group = 1), inherit.aes = FALSE, size = 1.3)

# 1.3 Observatory - key 변경하기
register_google(key='AIzaSyDLbJFdxvyywERerFK2piCIDPLjNLMlbQk')

obs <- fread("관측지점.csv")
df <- obs[종료일 == ""&지점 %in% data$RegionID, .(name = 지점명, lon = 경도, lat = 위도)]
head(df)

cen <- c(mean(df$lon),mean(df$lat))
map <- get_googlemap(center=cen,
                     maptype="roadmap",
                     zoom=7)

ggmap(map) + geom_point(data = df, aes(x = lon, y = lat), size = 4.5, col = "red") 
#+ geom_text_repel(data = df,aes(label = name), col = "red", size= 4.5) 



# 2. smoothing using regression analysis
z <- NULL
for(i in 1:5){
  monthbasis <- create.fourier.basis(c(0, 12), nbasis = (2*i +1), period = 12)
  fd_obj <- smooth.basis(1:12, t(data[,-1]), monthbasis, fdnames = list("month","region","Deg C"))$fd
  est <- eval.fd(1:12, fd_obj)
  mse <- mean(abs(t(data[,-1])-est))
  z <- c(z, mse)
}
z # 2*3+1= 7개의 basis로 선택

monthbasis <- create.fourier.basis(c(0, 12), nbasis = 7, period = 12)
fd_obj <- smooth.basis(1:12, t(data[,-1]), monthbasis, fdnames = list("month","region","Deg C"))$fd

month <- colnames(data)[-1]


smooth.basis.fun <- function(x) {
  result <- smooth.basis(1:12, x, monthbasis)
  return(result$fd)
}
data.smooth <- apply(data[,-1], 1, smooth.basis.fun)

# str(data.smooth)
plot(data.smooth[[1]], xlab="month", ylab="temperature",
     col=1, main="smooth temperature", ylim=c(-10, 40), xaxt="n")
for (i in 2:nrow(data)) lines(data.smooth[[i]], col=i) 
axis(1, at=1:12, labels=month)


# 2.1. first derivative of smoothing data

# ex
# plot(deriv.fd(data.smooth[[1]], 1))
# plot(deriv.fd(data.smooth[[1]]))

data.smooth.1 <- list()
for (i in 1:90){
  data.smooth.1[[i]] <- deriv.fd(data.smooth[[i]], 1)
}
plot(data.smooth.1[[1]], xlab="month", ylab="temperature",
     col=1, main="the first derivative curves", ylim=c(-11, 11), xaxt="n")
for (i in 2:90) lines(data.smooth.1[[i]], col=i)
axis(1, at=1:12, labels=month)

# 2.2. second derivative of smoothing data
data.smooth.2 <- list()
for (i in 1:90){
  data.smooth.2[[i]] <- deriv.fd(data.smooth[[i]], 2)
}
plot(data.smooth.2[[1]], xlab="month", ylab="temperature",
     col=1, main="the second derivative curves", ylim=c(-10, 10), xaxt="n")
for (i in 2:90) lines(data.smooth.2[[i]], col=i)
axis(1, at=1:12, labels=month)

# 3. Principal Component Analysis

monthbasis <- create.fourier.basis(c(0,12), nbasis = 7, period = 12)
plot(monthbasis)
fd_obj <- smooth.basis(1:12, t(data[,-1]), monthbasis, fdnames = list("month","region","Deg C"))$fd

# PCA - mean +/- eigenfunctions
pca_obj <- pca.fd(fd_obj, nharm = 7, centerfns = T )
plot.pca.fd(pca_obj)

# eigenfunctions
fdmat <- eval.fd(1:12, pca_obj[[1]])
plot(fdmat[,1], type = "b", ylim = c(-0.6, 0.6)); abline(h = 0, lty = 2)

# eigenvlaues
# number of eigenvalues
eigvals <- pca_obj[[2]]

plot(1:7, eigvals[1:7], type="b",
     xlab="Eigenvalue Number", ylab="Eigenvalues")

# score functions
harmscr <- pca_obj[[3]]

plot(harmscr[,1], harmscr[,2], xlab="Harmonic 1", ylab="Harmonic 2")
text(harmscr[,1], harmscr[,2], rownames(data), col=4)


# 4. clustering - region
set.seed(123)
str(data)
data.clustering <- data[, -1]

data.kmeans <- kmeans(data.clustering, centers=3)
data.kmeans$centers

data.clustering$cluster <- as.factor(data.kmeans$cluster)

qplot(data.clustering$Jan, data.clustering$Aug, colour=cluster, 
      data=data.clustering) +
  geom_text(label=rownames(data.clustering), nudge_x = 0.25, nudge_y = 0.25, check_overlap = T)




