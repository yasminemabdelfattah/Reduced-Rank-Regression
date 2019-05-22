load("~/Google Drive/PhD Thesis & Analysis/My PhD R work/R workfile/Spatiotemporal modeling/tropicalpcdf.RData")
View(tropicalpcdf)
# canonical correlation analysis--- Specification-----
library(candisc)
library(data.table)
setDT(tropicalpcdf)[, paste0('PC1', 1:3) := shift(PC1, 1:3)][]
setDT(tropicalpcdf)[, paste0('PC2', 1:3) := shift(PC2, 1:3)][]
setDT(tropicalpcdf)[, paste0('PC3', 1:3) := shift(PC3, 1:3)][]

## reduced rank regression 
library(rrr)
library(dplyr)
galaxy=cbind.data.frame(tropicalpcdf, rainpcdf)
galaxy <- na.omit(galaxy)
galaxy$date=NULL
galaxy=as.data.frame(galaxy)
galaxy_x <- select(galaxy, -PC1rain,-PC2rain, -PC3rain, -date)
galaxy_y <- select(galaxy, PC1rain:PC3rain)
pairwise_plot(galaxy_x, galaxy_y, type = "cva")
rank_trace(galaxy_x, galaxy_y, type = "cva")
rank_trace(galaxy_x, galaxy_y, type = "cva", plot = FALSE)
residuals(galaxy_x, galaxy_y, type = "cva", rank = 2, k = 0.001)
rrr(galaxy_x, galaxy_y, type = "cva", rank =2, k = 0.001)

scores=scores(galaxy_x, galaxy_y, type = "cva", rank = 2)
cca.y=as.data.frame(scores[c(-1,-3)])
cca.x=as.data.frame(scores[c(-2,-3)])
cca=cbind.data.frame(cca.x, cca.y)
library(zoo)
cca.y$date<- seq(from = as.Date("1998-04-01", tz = 'UTC'), to = as.Date("2016-12-01", tz = 'UTC'), by = "month")
cca.x$date<- seq(from = as.Date("1998-04-01", tz = 'UTC'), to = as.Date("2016-12-01", tz = 'UTC'), by = "month")
cca=merge(cca.y, cca.x, by = "date")
str(cca)
cca1<- select(cca, -omega.omega2,-xi.xi2)
cca2<- select(cca, omega.omega2,xi.xi2, date)

library(reshape2)
library(scales)
require(ggplot2)
ccaMelted1<- reshape2::melt(cca1, id.var='date')
cca1=ggplot(ccaMelted1, aes(x=date, y=value, col=variable)) + geom_line()
cca1=cca1+ scale_x_date(date_breaks = "1 year", date_labels =  "%Y")  
cca1
ccaMelted2<- reshape2::melt(cca2, id.var='date')
cca2=ggplot(ccaMelted2, aes(x=date, y=value, col=variable)) + geom_line()
cca2=cca2+ scale_x_date(date_breaks = "1 year", date_labels =  "%Y")  
cca2

sqrt_matrix <- function(matr){
  eigens <- eigen(matr)
  vecs <- eigens[["vectors"]]
  vals <- eigens[["values"]]
  vecs %*% diag(sqrt(vals)) %*% t(vecs)
}
organize <- function(vars, type = "cov"){
  matr <- as.matrix(vars)
  if(type == "cov"){
    matr %>%
      scale(scale = FALSE) %>%
      t()
  } else if(type == "cor"){
    matr %>%
      scale() %>%
      t()
  } else {
    stop("type not recognized")
  }
}

cov_matrix <- function(x, y, type = "cov"){
  n <- dim(x)[1]
  x_org <- organize(x, type)
  y_org <- organize(y, type)
  x %*% t(y) / n
}


reduce_rank_regression <- function(x, y, gamma_matrix, rank = "full", k = 0){
  full_rank <- min(dim(x)[2], dim(y)[2])
  if(rank == "full"){
    reduce_rank <- full_rank
  } else if(rank <= full_rank){
    reduce_rank <- rank
  } else {
    stop("rank out of bounds")
  }
  cov_x <- cov(x) + k * diag(1, dim(x)[2])
  cov_yx <- cov(y, x)
  cov_y <- cov(y) + k * diag(1, dim(y)[2])
  cov_xy <- t(cov_yx)
  sqrtm <- sqrt_matrix(gamma_matrix)
  weighted_matrix <- sqrtm %*%
    cov_yx %*%
    solve(cov_x) %*%
    cov_xy %*%
    sqrtm
  eigens <- eigen(weighted_matrix)
  eigen_values <- eigens[["values"]]
  V_t <- eigens[["vectors"]][,1:reduce_rank] %>%
    as.matrix(ncol = reduce_rank)
  A_t <- solve(sqrtm) %*% V_t
  rownames(A_t) <- names(y)
  B_t <- t(V_t) %*%
    sqrtm %*%
    cov_yx %*%
    solve(cov_x)
  C_t <- A_t %*% B_t
  mu_y <- colMeans(y)
  mu_x <- colMeans(x)
  mu_t <- mu_y - C_t %*% mu_x
  list(mean = mu_t, A = A_t, B = B_t, C = C_t, eigenvalues = eigen_values)
}


cva <- function(x, y, rank = "full", k = 0) {
  gamma <- solve(cov(y) + k * diag(1, dim(y)[2]))
  rrr_object <- reduce_rank_regression(x, y, gamma, rank, k)
  H <- ginv(rrr_object[["A"]])
  colnames(H) <- names(y)
  list(mean = rrr_object[["mean"]],
       G = rrr_object[["B"]],
       H = H,
       canonical_corr = rrr_object[["eigenvalues"]])
}

cva_error <- function(x, y, x_new, y_new, rank = "full", k = 0){
  cva_object <- cva(x, y, rank, k)
  index <- data_frame(index = 1:dim(y_new)[1])
  error <- as_data_frame(t(cva_object[["H"]] %*%
                             organize(y_new) - cva_object[["G"]] %*% organize(x_new)))
  names(error) <- paste("CV", 1:dim(error)[2], sep = "")
  dplyr::bind_cols(index, error)
}


# galaxy_train <-
series <- ts(galaxy, frequency = 12, start = c(1998, 4))
# Training set
# Use data from 1949 to 1955 for forecasting
train= window(series, start= c(1998, 4), end=c(2010,12))
train=as.data.frame(train)
# Test set
# Use remaining data from 1956 to 1960 to test accuracy
test= window(series, start= c(2011, 1), end=c(2016,12))
test=as.data.frame(test)

train_x <- select(train, -PC1rain,-PC2rain, -PC3rain, -date)
train_y <- select(train, PC1rain:PC3rain)
test_x <- select(test, -PC1rain,-PC2rain, -PC3rain, -date)
test_y <- select(test, PC1rain:PC3rain)
cvaerror=cva_error(train_x, train_y, test_x, test_y, rank = 2, k = 0.001)
cvaerror=as.data.frame(cvaerror)
RMSE1=sqrt(mean((cvaerror$CV1)^2))
RMSE1
RMSE2=sqrt(mean((cvaerror$CV2)^2))
RMSE2

# canonical correlation analysis--- lead 1-----
load("~/Google Drive/PhD Thesis & Analysis/My PhD R work/R workfile/Spatiotemporal modeling/tropicalpcdf.RData")
View(tropicalpcdf)
library(candisc)
library(data.table)
setDT(tropicalpcdf)[, paste0('PC1', 1:4) := shift(PC1, 1:4)][]
setDT(tropicalpcdf)[, paste0('PC2', 1:4) := shift(PC2, 1:4)][]
setDT(tropicalpcdf)[, paste0('PC3', 1:4) := shift(PC3, 1:4)][]

## reduced rank regression 
library(rrr)
library(dplyr)
galaxy=cbind.data.frame(tropicalpcdf, rainpcdf)
galaxy <- na.omit(galaxy)
galaxy$date=NULL
galaxy=as.data.frame(galaxy)
galaxy_x <- select(galaxy, -PC1rain,-PC2rain, -PC3rain, -date)
galaxy_y <- select(galaxy, PC1rain:PC3rain)
pairwise_plot(galaxy_x, galaxy_y, type = "cva")
rank_trace(galaxy_x, galaxy_y, type = "cva")
rank_trace(galaxy_x, galaxy_y, type = "cva", plot = FALSE)
residuals(galaxy_x, galaxy_y, type = "cva", rank = 2, k = 0.001)
rrr(galaxy_x, galaxy_y, type = "cva", rank =2, k = 0.001)

scores=scores(galaxy_x, galaxy_y, type = "cva", rank = 2)
cca.y=as.data.frame(scores[c(-1,-3)])
cca.x=as.data.frame(scores[c(-2,-3)])
library(zoo)
cca.y$date<- seq(from = as.Date("1998-05-01", tz = 'UTC'), to = as.Date("2016-12-01", tz = 'UTC'), by = "month")
cca.x$date<- seq(from = as.Date("1998-05-01", tz = 'UTC'), to = as.Date("2016-12-01", tz = 'UTC'), by = "month")
cca=merge(cca.y, cca.x, by = "date")
str(cca)
cca1<- select(cca, -omega.omega2,-xi.xi2)
cca2<- select(cca, omega.omega2,xi.xi2, date)

library(reshape2)
library(scales)
require(ggplot2)
ccaMelted1<- reshape2::melt(cca1, id.var='date')
cca1=ggplot(ccaMelted1, aes(x=date, y=value, col=variable)) + geom_line()
cca1=cca1+ scale_x_date(date_breaks = "1 year", date_labels =  "%Y")  
cca1
ccaMelted2<- reshape2::melt(cca2, id.var='date')
cca2=ggplot(ccaMelted2, aes(x=date, y=value, col=variable)) + geom_line()
cca2=cca2+ scale_x_date(date_breaks = "1 year", date_labels =  "%Y")  
cca2

# galaxy_train <-
series <- ts(galaxy, frequency = 12, start = c(1998, 5))
# Training set
# Use data from 1949 to 1955 for forecasting
train= window(series, start= c(1998, 5), end=c(2010,12))
train=as.data.frame(train)
# Test set
# Use remaining data from 1956 to 1960 to test accuracy
test= window(series, start=c(2011, 1), end=c(2016,12))
test=as.data.frame(test)

train_x <- select(train, -PC1rain,-PC2rain, -PC3rain, -date)
train_y <- select(train, PC1rain:PC3rain)
test_x <- select(test, -PC1rain,-PC2rain, -PC3rain, -date)
test_y <- select(test, PC1rain:PC3rain)
cvaerror=cva_error(train_x, train_y, test_x, test_y, rank = 2, k = 0.001)
cvaerror=as.data.frame(cvaerror)
RMSE1=sqrt(mean((cvaerror$CV1)^2))
RMSE1
RMSE2=sqrt(mean((cvaerror$CV2)^2))
RMSE2
# canonical correlation analysis--- lead 2-----
load("~/Google Drive/PhD Thesis & Analysis/My PhD R work/R workfile/Spatiotemporal modeling/tropicalpcdf.RData")
View(tropicalpcdf)
library(candisc)
library(data.table)
setDT(tropicalpcdf)[, paste0('PC1', 2:5) := shift(PC1, 2:5)][]
setDT(tropicalpcdf)[, paste0('PC2', 2:5) := shift(PC2, 2:5)][]
setDT(tropicalpcdf)[, paste0('PC3', 2:5) := shift(PC3, 2:5)][]

## reduced rank regression 
library(rrr)
library(dplyr)
galaxy=cbind.data.frame(tropicalpcdf, rainpcdf)
galaxy <- na.omit(galaxy)
galaxy$date=NULL
galaxy=as.data.frame(galaxy)
galaxy_x <- select(galaxy, -PC1rain,-PC2rain, -PC3rain, -date)
galaxy_y <- select(galaxy, PC1rain:PC3rain)
pairwise_plot(galaxy_x, galaxy_y, type = "cva")
rank_trace(galaxy_x, galaxy_y, type = "cva")
rank_trace(galaxy_x, galaxy_y, type = "cva", plot = FALSE)
residuals(galaxy_x, galaxy_y, type = "cva", rank = 2, k = 0.001)
rrr(galaxy_x, galaxy_y, type = "cva", rank =2, k = 0.001)

scores=scores(galaxy_x, galaxy_y, type = "cva", rank = 2)
cca.y=as.data.frame(scores[c(-1,-3)])
cca.x=as.data.frame(scores[c(-2,-3)])
library(zoo)
cca.y$date<- seq(from = as.Date("1998-05-01", tz = 'UTC'), to = as.Date("2016-12-01", tz = 'UTC'), by = "month")
cca.x$date<- seq(from = as.Date("1998-05-01", tz = 'UTC'), to = as.Date("2016-12-01", tz = 'UTC'), by = "month")
cca=merge(cca.y, cca.x, by = "date")
str(cca)
cca1<- select(cca, -omega.omega2,-xi.xi2)
cca2<- select(cca, omega.omega2,xi.xi2, date)

library(reshape2)
library(scales)
require(ggplot2)
ccaMelted1<- reshape2::melt(cca1, id.var='date')
cca1=ggplot(ccaMelted1, aes(x=date, y=value, col=variable)) + geom_line()
cca1=cca1+ scale_x_date(date_breaks = "1 year", date_labels =  "%Y")  
cca1
ccaMelted2<- reshape2::melt(cca2, id.var='date')
cca2=ggplot(ccaMelted2, aes(x=date, y=value, col=variable)) + geom_line()
cca2=cca2+ scale_x_date(date_breaks = "1 year", date_labels =  "%Y")  
cca2

# galaxy_train <-
series <- ts(galaxy, frequency = 12, start = c(1998, 6))
# Training set
# Use data from 1949 to 1955 for forecasting
train= window(series, start= c(1998, 6), end=c(2010,12))
train=as.data.frame(train)
# Test set
# Use remaining data from 1956 to 1960 to test accuracy
test= window(series, start= c(2011, 1), end=c(2016,12))
test=as.data.frame(test)

train_x <- select(train, -PC1rain,-PC2rain, -PC3rain, -date)
train_y <- select(train, PC1rain:PC3rain)
test_x <- select(test, -PC1rain,-PC2rain, -PC3rain, -date)
test_y <- select(test, PC1rain:PC3rain)
cvaerror=cva_error(train_x, train_y, test_x, test_y, rank = 2, k = 0.001)
cvaerror=as.data.frame(cvaerror)
RMSE1=sqrt(mean((cvaerror$CV1)^2))
RMSE1
RMSE2=sqrt(mean((cvaerror$CV2)^2))
RMSE2
# canonical correlation analysis--- lead 3-----
load("~/Google Drive/PhD Thesis & Analysis/My PhD R work/R workfile/Spatiotemporal modeling/tropicalpcdf.RData")
View(tropicalpcdf)
library(candisc)
library(data.table)
setDT(tropicalpcdf)[, paste0('PC1', 3:6) := shift(PC1, 3:6)][]
setDT(tropicalpcdf)[, paste0('PC2', 3:6) := shift(PC2, 3:6)][]
setDT(tropicalpcdf)[, paste0('PC3', 3:6) := shift(PC3, 3:6)][]

## reduced rank regression 
library(rrr)
library(dplyr)
galaxy=cbind.data.frame(tropicalpcdf, rainpcdf)
galaxy <- na.omit(galaxy)
galaxy$date=NULL
galaxy=as.data.frame(galaxy)
galaxy_x <- select(galaxy, -PC1rain,-PC2rain, -PC3rain, -date)
galaxy_y <- select(galaxy, PC1rain:PC3rain)
pairwise_plot(galaxy_x, galaxy_y, type = "cva")
rank_trace(galaxy_x, galaxy_y, type = "cva")
rank_trace(galaxy_x, galaxy_y, type = "cva", plot = FALSE)
residuals(galaxy_x, galaxy_y, type = "cva", rank = 2, k = 0.001)
rrr(galaxy_x, galaxy_y, type = "cva", rank =2, k = 0.001)

scores=scores(galaxy_x, galaxy_y, type = "cva", rank = 2)
cca.y=as.data.frame(scores[c(-1,-3)])
cca.x=as.data.frame(scores[c(-2,-3)])
library(zoo)
cca.y$date<- seq(from = as.Date("1998-07-01", tz = 'UTC'), to = as.Date("2016-12-01", tz = 'UTC'), by = "month")
cca.x$date<- seq(from = as.Date("1998-07-01", tz = 'UTC'), to = as.Date("2016-12-01", tz = 'UTC'), by = "month")
cca=merge(cca.y, cca.x, by = "date")
str(cca)
cca1<- select(cca, -omega.omega2,-xi.xi2)
cca2<- select(cca, omega.omega2,xi.xi2, date)

library(reshape2)
library(scales)
require(ggplot2)
ccaMelted1<- reshape2::melt(cca1, id.var='date')
cca1=ggplot(ccaMelted1, aes(x=date, y=value, col=variable)) + geom_line()
cca1=cca1+ scale_x_date(date_breaks = "1 year", date_labels =  "%Y")  
cca1
ccaMelted2<- reshape2::melt(cca2, id.var='date')
cca2=ggplot(ccaMelted2, aes(x=date, y=value, col=variable)) + geom_line()
cca2=cca2+ scale_x_date(date_breaks = "1 year", date_labels =  "%Y")  
cca2

# galaxy_train <-
series <- ts(galaxy, frequency = 12, start = c(1998, 7))
# Training set
# Use data from 1949 to 1955 for forecasting
train= window(series, start= c(1998, 07), end=c(2010,12))
train=as.data.frame(train)
# Test set
# Use remaining data from 1956 to 1960 to test accuracy
test= window(series, start= c(2011, 1), end=c(2016,12))
test=as.data.frame(test)

train_x <- select(train, -PC1rain,-PC2rain, -PC3rain, -date)
train_y <- select(train, PC1rain:PC3rain)
test_x <- select(test, -PC1rain,-PC2rain, -PC3rain, -date)
test_y <- select(test, PC1rain:PC3rain)
cvaerror=cva_error(train_x, train_y, test_x, test_y, rank = 2, k = 0.001)
cvaerror=as.data.frame(cvaerror)
RMSE1=sqrt(mean((cvaerror$CV1)^2))
RMSE1
RMSE2=sqrt(mean((cvaerror$CV2)^2))
RMSE2
