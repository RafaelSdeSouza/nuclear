# Rafael de Souza, UNC
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License version 3 as published by
#the Free Software Foundation.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#
#' @title  Estimate reaction rates over a temperature grid
#' @description Provides a reaction rate table
#' @aliases NumRateHe3dpLoop
#' @usage NumRateHe3dpLoop (mat = mat, vars=vars, N = N, T9  =T9)
#' @format \describe{
#' \item{x}{
#' The function has 4 arguments: mat, vars, N, T9}
#' }
#' @param mat mat
#' @param vars  vars
#' @param N    N
#' @param T9  T9
#' @return NumRateHe3dpLoop
#' @import gsl
#' @author Rafael de Souza, UNC,  and Christian Illiadis, UNC
#'
#' @keywords NumRateHe3dpLoop
#' @export
#'
#'
NumRateHe3dpLoop <- function(mat, vars=c("e0","gi","gf","ri","rf"), N = 1000,T9){
#  mcdat_I <- as.data.frame(do.call(rbind, as.mcmc(mat)[,vars]))
  mcdat_I <- mat[,vars]
  index <- sample(1:nrow(mcdat_I ),size=N,replace=FALSE)
  mcdat_I  <- mcdat_I [index,]

  gdat <-  sapply(Tgrid,function(Tgrid){NumRate3Hedp(mcdat_I,T9=Tgrid)})

  gg <-  as.data.frame(gdat)

  gg2 <- apply(gg, 1, quantile, probs=c(0.16, 0.5, 0.84), na.rm=TRUE)

  fu <- function(x){exp(sqrt(log(1+var(x)/mean(x)^2)))}

  fu_I<-apply(gg, 1, fu)

  gg2data <- data.frame(T9 =T9, lower = gg2["16%",], mean = gg2["50%",], upper = gg2["84%",] )
  gg2data$fu <- fu_I
  return(gg2data)
}
