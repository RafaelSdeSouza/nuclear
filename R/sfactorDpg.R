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
#' @title  Estimate Astrophysical S-factor
#' @description Provides a confusion matrix of classification statistics following logistic regression.
#' @aliases sfactorDpg
#' @usage sfactorDpg(ecm = ecm, a.scale = a.scale)
#' @format \describe{
#' \item{x}{
#' The function has two  arguments:  ecm, a.scale}
#' }
#' @param ecm ecm
#' @param a.scale a.scale
#' @return S-factor
#' @examples
#' library(nuclear)
#'
#' N <- 300
#' obsx1 <- exp(seq(log(1e-3), log(1),length.out=N))
#' plot(obsx1,sfactorDpg(obsx1),
#' col="red",cex=1.25,type="l",ylab="S-factor",xlab="E",log="x")
#'
#' @author Rafael de Souza, UNC,  and Christian Illiadis, UNC
#'
#' @keywords S-factor
#' @export
#'
#'
sfactorDpg <- function(ecm,a.scale=1){
  data(dpg)
  th <- approxfun(dpg[,1], dpg[,2])
  SF <- a.scale*th(ecm)*1e-6
  return(SF = SF)
}

