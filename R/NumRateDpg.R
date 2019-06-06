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
#' @title  Estimate reaction rates
#' @description Calculate numerical reaction rates
#' @aliases NumRateDpg
#' @usage NumRateDpg(a.scale = a.scale,T9 = T9)
#' @format \describe{
#' \item{x}{
#' The function has two  arguments: a.scale, T9}
#' }
#' @param a.scale  a.scale
#' @param T9  T9
#' @return nuclear_rate
#' @import gsl
#' @examples
#' library(nuclear)
#'
#'  NumRateDpg(1,10)
#'
#' @author Rafael de Souza, UNC,  and Christian Illiadis, UNC
#'
#' @keywords Nuclear rate
#' @export
#'
#'

NumRateDpg  <-  Vectorize(function(a.scale,T9){
  # Constants
  M0 = 1.007276466812; M1 = 2.01355318262;		# masses (amu) of p and d
  Z0 = 1; Z1 = 1 ;			# charges of t and d

  #   DEFINITIONS
  mue <- (M0*M1)/(M0 + M1)
  dpieta <- function(E){0.98951013*Z0*Z1*sqrt(mue/E)}

  #     ----------------------------------------------------
  #     Integrand
  #     ----------------------------------------------------

  integrand <- function(E,T9) {exp(-dpieta(E))*sfactorDpg(E,a.scale)*exp(-E/(0.086173324*T9))}

  # CALCULATE Nuclear rate

  Nasv <- function(Temp){(3.7318e10/Temp^{3/2})*sqrt(1/mue)*integrate(integrand, lower = 1e-5, upper = 2,
                                                                      abs.tol = 0L,T9 = Temp)$value}

  out <- Nasv(T9)
  return(Nasv=out)
}
)

