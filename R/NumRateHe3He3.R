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
#' @aliases NumRateHe3He3
#' @usage NumRateHe3He3 (S0 = S0,S1 = S1,S2 = S2,i.screen = i.screen,T9 = T9)
#' @format \describe{
#' \item{x}{
#' The function has two  arguments: S0,S1,S2,i.screen, T9}
#' }
#' @param S0  S0
#' @param S1  S1
#' @param S2  S2
#' @param i.screen  i.screen
#' @param T9  T9
#' @return nuclear_rate
#' @examples
#' library(nuclear)
#'
#'  NumRateHe3He3 (1,1)
#'
#' @author Rafael de Souza, UNC,  and Christian Illiadis, UNC, and HK
#'
#' @keywords Nuclear rate
#' @export
#'
#'
NumRateHe3He3 <- Vectorize(function(S0,S1,S2,i.screen,T9){

  # Constants
  M0 = 3.0160293; M1 = 3.0160293;		# masses (amu) of p and d
  Z0 = 2; Z1 = 2;			# charges of t and d

  #   DEFINITIONS
  mue <- (M0*M1)/(M0 + M1) # Reduced Center of Mass
  dpieta <- function(E){0.98951013*Z0*Z1*sqrt(mue/E)} # 2 pi eta

  #     ----------------------------------------------------
  #     Integrand
  #     ----------------------------------------------------

  integrand <- function(E,T9) {exp(-dpieta(E))*sfactorHe3He3(E,S0,S1,S2,i.screen)*exp(-E/(0.086173324*T9))}
  # Integrand, this means that we are looking at Gamow*S-factor*Boltzmann Factor

  # CALCULATE Nuclear rate

  Nasv <- function(Temp){(3.7318e10/Temp^{3/2})*sqrt(1/mue)*integrate(integrand, lower = 1e-5, upper = Inf,
                                                                      T9 = Temp)$value}

  # Note to self, the limits of integration, in some sense, the scale should be appropriate.
  # From HELP, the first argument MUST BE integrated. The optional argument T9 is used to be substituted in
  # Nasv <-> N A <sigma v >

  out <- Nasv(T9)
  return(Nasv = out)
}
)




NumRateHe3He3(5.190775,-2.768626,1.12796,0,1)





