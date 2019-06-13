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
#' @aliases sfactorHe3He3
#' @usage sfactorHe3He3(ecm = ecm, S0 = S0, S1 = S1,S2 = S2, i.screen = i.screen)
#' @format \describe{
#' \item{x}{
#' The function has 5  arguments: ecm,S0, S1, S2,i.screen}
#' }
#' @param ecm ecm
#' @param S0  S0
#' @param S1  S1
#' @param S2  S2
#' @param i.screen i.screen
#' @return S-factor
#' @examples
#' library(nuclear)
#'
#' N <- 300
#' obsx1 <- exp(seq(log(1e-2), log(1),length.out=N))
#' plot(obsx1,sfactorHe3He3(obsx1,5.14,-2.69,2.14,325*1e-6),
#' col="red",cex=1.25,type="l",ylab="S-factor",xlab="E",log="x")
#'
#' @author Rafael de Souza, UNC,  and Christian Illiadis, UNC
#'
#' @keywords S-factor
#' @export
#'
#'
#'
#'
sfactorHe3He3 <- function(ecm,S0, S1, S2,i.screen){
  M0 = 3.0160293; M1 = 3.0160293;		# masses (amu) of p and d
  Z0 = 2; Z1 = 2 ;			# charges of t and d
  mue <- (M0*M1)/(M0 + M1) # Reduced Center of Mass

SF <- (exp(0.98951013*Z0*Z1*sqrt(mue)/(2*ecm^(3/2))*i.screen))*(S0 + S1*ecm + 0.5*S2*ecm^2)

  return(SF = SF)
}
