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
#' @aliases sfactorDdp
#' @usage sfactorDdp(ecm = ecm)
#' @format \describe{
#' \item{x}{
#' The function has one arguments:  ecm}
#' }
#' @param ecm ecm
#' @return S-factor
#' @examples
#' library(nuclear)
#'
#' N <- 300
#' obsx1 <- exp(seq(log(1e-3), log(1),length.out=N))
#' plot(obsx1,sfactorDdp(obsx1),
#' col="red",cex=1.25,type="l",ylab="S-factor",xlab="E",log="x")
#'
#' @author Rafael de Souza, UNC,  and Christian Illiadis, UNC
#'
#' @keywords S-factor
#' @export
#'
#'
sfactorDdp <- function(ecm){
  data(ddp_th)
  th <- approxfun(ddp_th[,1], ddp_th[,2])
  SF <- th(ecm)
  return(SF = SF)
}

