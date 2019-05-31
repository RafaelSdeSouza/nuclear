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
#' @aliases NumRate3Hedp
#' @usage NumRate3Hedp(e0 = e0, gi = gi, gf = gf, ri = ri, rf = rf, T9 = T9)
#' @format \describe{
#' \item{x}{
#' The function has six arguments: e0, gi, gf, ri, rf, T9}
#' }
#' @param e0  e0
#' @param gi  gi
#' @param gf  gf
#' @param ri  ri
#' @param rf  rf
#' @param T9  T9
#' @return nuclear_rate
#' @import gsl
#' @examples
#' library(nuclear)
#'
#'  NumRate3Hedp(0.35,1.0085,0.025425,5,5,10)
#'
#' @author Rafael de Souza, UNC,  and Christian Illiadis, UNC
#'
#' @keywords Nuclear rate
#' @export
#'
#'

NumRate3Hedp  <-  Vectorize(function(e0,gi,gf,ri,rf,T9){
  er = e0
  # Constants
  M0 = 3.01493216; M1 = 2.01355332;		# masses (amu) of projectile and target
  Z0 = 2; Z1 = 1 ;			# charge projectile and target

  #   DEFINITIONS
  mue <- (M0*M1)/(M0 + M1)
  dpieta <- function(E){0.98951013*Z0*Z1*sqrt(mue/E)}

  #     ----------------------------------------------------
  #     Integrand
  #     ----------------------------------------------------

  integrand <- function(E,T9) {exp(-dpieta(E))*sfactor3Hedp(E,e0,gi,gf,ri,rf)*exp(-E/(0.086173324*T9))}

  # CALCULATE Nuclear rate

  Nasv <- function(Temp){(3.7318e10/Temp^{3/2})*sqrt(1/mue)*integrate(integrand, lower = 1e-5, upper = Inf,
          abs.tol = 0L,T9 = Temp)$value}

  out <- Nasv(T9)
  return(Nasv=out)
}
)
