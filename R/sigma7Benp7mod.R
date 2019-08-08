

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
#' @description Provides a sqrt(E)*sigma for Be7(n,p).
#' @aliases sigma7Benp7mod
#' @usage sigma7Benp7mod(ecm = ecm, e_01 = e_01, ga_1 = ga_1, gb_1 = gb_1,
#' ra_1 = ra_1, rb_1 = rb_1, jr_1 = jr_1,la_1 = la_1, lb_1 = lb_1)
#' @format \describe{
#' \item{x}{
#' The function has 57 arguments: ecm, e0_1, ga_1, gb_1, ra_1, rb_1, jr_1, la_1, lb_1,
#'                                     e0_2, ga_2, gb_2, ra_2, rb_2, jr_2, la_2, lb_2,
#'                                     e0_3, ga_3, gb_3, ra_3, rb_3, jr_3, la_3, lb_3,
#'                                     e0_4, ga_4, gb_4, ra_4, rb_4, jr_4, la_4, lb_4,
#'                                     e0_5, ga_5, gb_5, ra_5, rb_5, jr_5, la_5, lb_5,
#'                                     e0_6, ga_6, gb_6, ra_6, rb_6, jr_6, la_6, lb_6,
#'                                     e0_7, ga_7, gb_7, ra_7, rb_7, jr_7, la_7, lb_7}
#' }
#' @param ecm ecm
#' @param e0_1  e0_1
#' @param ga_1  ga_1
#' @param g_b  g_b
#' @param ra_1 ra_1
#' @param rb_1 rb_1
#' @param jr_1 jr_1
#' @param la_1 la_1
#' @param lb_1 lb_1
#' @return Sigma
#' @import gsl
#' @author Rafael de Souza, UNC,  and Christian Illiadis, UNC
#'
#' @keywords sqrt(E) * SIGMA
#' @export
#'
#'
sigma7Benp7mod <- function(ecm,M){

  SF1 <-  sigma7Benp(ecm, e0 = M[["e0_1"]], ga = M[["ga_1"]], gb = M[["gb_1"]], ra = M[["ra"]], rb = M[["rb"]], jr = 2, la = 0, lb = 0)
  SF2 <-  sigma7Benp(ecm, e0 = M[["e0_2"]], ga = M[["ga_2"]], gb = M[["gb_2"]], ra = M[["ra"]], rb = M[["rb"]], jr = 3, la = 1, lb = 1)
  SF3 <-  sigma7Benp(ecm, e0 = M[["e0_3"]], ga = M[["ga_3"]], gb = M[["gb_3"]], ra = M[["ra"]], rb = M[["rb"]], jr = 3, la = 1, lb = 1)
  SF4 <-  sigma7Benp(ecm, e0 = M[["e0_4"]], ga = M[["ga_4"]], gb = M[["gb_4"]], ra = M[["ra"]], rb = M[["rb"]], jr = 1, la = 0, lb = 0)
  SF5 <-  sigma7Benp(ecm, e0 = M[["e0_5"]], ga = M[["ga_5"]], gb = M[["gb_5"]], ra = M[["ra"]], rb = M[["rb"]], jr = 4, la = 3, lb = 3)
  SF6 <-  sigma7Benp(ecm, e0 = M[["e0_6"]], ga = M[["ga_6"]], gb = M[["gb_6"]], ra = M[["ra"]], rb = M[["rb"]], jr = 2, la = 1, lb = 1)
  SF7 <-  sigma7Benp(ecm, e0 = M[["e0_7"]], ga = M[["ga_7"]], gb = M[["gb_7"]], ra = M[["ra"]], rb = M[["rb"]], jr = 0, la = 1, lb = 1)
  SF <- SF1 + SF2 + SF3 + SF4 + SF5 + SF6 + SF7
  return(SF = SF)
}


