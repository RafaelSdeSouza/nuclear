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
#' @aliases NumRate7Benp
#' @usage NumRate7Benp(e0 = e0, gi = gi, gf = gf, ri = ri, rf = rf, T9 = T9)
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
#' NumRate7Benp(0.35,1.0085,0.025425,5,5,10)
#'
#' @author Rafael de Souza, UNC,  and Christian Illiadis, UNC
#'
#' @keywords Nuclear rate
#' @export
#'
#'

NumRate7Benp   <- function(x, T9){


  #     ----------------------------------------------------
  #     Integrand
  #     ----------------------------------------------------

  integrand <- function(E,T9) {(sqrt(E) * sigma7Benp7mod(E,x) + x["hbg"])*
      exp(-E/(0.086173324*T9))}

  # CALCULATE Nuclear rate

  m1 = 7.01473482886
  m2 = 1.00866491582   # masses (amu) of 7Be and n
  mue = (m1*m2)/(m1+m2)

  Nasv <- function(Temp){(3.7318e10/Temp^{3/2})*sqrt(1/mue)*integrate(integrand, lower = 1e-10, upper = 2,
                                                                      abs.tol = 0L,
                                                                      T9 = Temp)$value}

  # Note to self, the limits of integration, in some sense, the scale should be appropriate.
  # From HELP, the first argument MUST BE integrated. The optional argument T9 is used to be substituted in
  # Nasv <-> N A <sigma v >

  out <- Nasv(T9)
  return(Nasv=out)
}

