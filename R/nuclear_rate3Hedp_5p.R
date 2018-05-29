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
#' @aliases nuclear_rate3Hedp_5p
#' @usage nuclear_rate3Hedp_5p(ECM = ECM, gd = gd, gp = gp, rd = rd, rp = rp, T9 = T9)
#' @format \describe{
#' \item{x}{
#' The function has four arguments: ECM, gd, gp, rd, rp, T9}
#' }
#' @param ER  ECM
#' @param gd  gd
#' @param gp  gp
#' @param rd  rd
#' @param rp  rp
#' @param T9  T9
#' @return nuclear_rate
#' @import gsl
#' @examples
#' library(nuclear)
#'
#' N <- 300
#' obsx1 <- exp(seq(log(1e-3), log(1),length.out=N))
#' plot(obsx1,sfactor3Hedp(obsx1,0.35,1.0085,0.025425),
#' col="red",cex=1.25,type="l",ylab="S-factor",xlab="E",log="x")
#'
#' @author Rafael de Souza, UNC,  and Christian Illiadis, UNC
#'
#' @keywords Nuclear_rate
#' @export
#'
#'

nuclear_rate3Hedp_5p <- function(ECM, gd, gp, rd, rp, T9){
  # Constants
  M0 = 3.01493216; M1 = 2.01355332;		# masses (amu) of t and d
  Z0 = 2; Z1 = 1 ;			# charges of t and d

  #   DEFINITIONS
  mue <- (M0*M1)/(M0 + M1)
  dpieta <- function(E){0.98951013*Z0*Z1*sqrt(mue/E)}

  #     ----------------------------------------------------
  #     Integrand
  #     ----------------------------------------------------

  integrand <- function(E,T9) {exp(-dpieta(E))*sfactor3Hedp_5p(E,ECM,gd,gp,rd,rp)*exp(-E/(0.086173324*T9))}

  # CALCULATE Nuclear rate

  Nasv <- function(Temp){(3.7318e10/Temp^{3/2})*sqrt(1/mue)*integrate(integrand, lower = 1e-5, upper = Inf,
          abs.tol = 0L,T9 = Temp)$value}

  out <- Nasv(T9)
  return(Nasv=out)
}
