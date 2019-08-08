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
#' @aliases sigma7Benp
#' @usage sigma7Benp(ecm = ecm, e0 = e0, ga = ga, gb = gb,
#' ra = ra, rb = rb, jr = jr,la = la, lb = lb)
#' @format \describe{
#' \item{x}{
#' The function has 9  arguments: ecm, e0, ga, gb, ra, rb, jr,la,lb}
#' }
#' @param ecm ecm
#' @param e0  e0
#' @param ga  ga
#' @param gb  gb
#' @param ra ra
#' @param rb rb
#' @param jr jr
#' @param la la
#' @param lb lb
#' @return Sigma
#' @import gsl
#' @examples
#' library(nuclear)
#'
#' N <- 300
#' obsx1 <- exp(seq(log(1e-5), log(5),length.out=N))
#' plot(obsx1,sqrt(obsx1)*sigma7Benp(obsx1,0.021,1.6,0.4,5,5,2,0,0),
#' col="red",cex=1.25,type="l",ylab="Sigma",xlab="E",log="x")
#'
#' @author Rafael de Souza, UNC,  and Christian Illiadis, UNC
#'
#' @keywords sqrt(E) * SIGMA
#' @export
#'
#'
sigma7Benp  <-  function(ecm, e0, ga, gb, ra, rb, jr, la, lb){
  # input masses, charges, angular momenta
  m1_i = 7.01473482886
  m2_i = 1.00866491582  # masses (amu) of 7Be and n
  m1_f = 7.01435791572
  m2_f = 1.00727646658  # masses (amu) of 7Li and p
  z1_i = 4
  z2_i = 0              # charges of 7Be and n
  z1_f = 3
  z2_f = 1              # charges of 7Li and p
  jt = 1.5              # spins of target, projectile
  jp = 0.5
  Q = 1.644425          # reaction Q-value (MeV) [from nuclear masses]

  # reduced masses
  mue_i <- (m1_i*m2_i)/(m1_i+m2_i)
  mue_f <- (m1_f*m2_f)/(m1_f+m2_f)

  # constants
  pek <- 6.56618216e-1/mue_i
  omega <- (2*jr+1)/((2*jt+1)*(2*jp+1))

  ### CALCULATE S-FACTOR
  ## incoming channel
  eta_a=0.1574854*z2_i*z1_i*sqrt(mue_i)
  rho_a=0.218735*ra*sqrt(mue_i)
  eta_i=eta_a/(sqrt(ecm))
  rho_i=rho_a*(sqrt(ecm))
  P3 <- coulomb_wave_FG(eta_i, rho_i, la, k=0)
  # penetration and shift factor
  p_i <- rho_i/(P3$val_F^2 + P3$val_G^2)
  s_i <- rho_i*(P3$val_F*P3$val_Fp + P3$val_G*P3$val_Gp)/(P3$val_F^2 + P3$val_G^2)
  # shift factor at energy E0 [eigenvalue]
  xeta_i=eta_a/(sqrt(e0))
  xrho_i=rho_a*(sqrt(e0))
  PX1 <- coulomb_wave_FG(xeta_i, xrho_i, la, k=0)
  b_i <- xrho_i*(PX1$val_F*PX1$val_Fp + PX1$val_G*PX1$val_Gp)/(PX1$val_F^2 + PX1$val_G^2)
  # partial width
  Ga <- 2*ga*p_i

  ## outgoing channel
  eta_b=0.1574854*z2_f*z1_f*sqrt(mue_f)
  rho_b=0.218735*rb*sqrt(mue_f)
  eta_f=eta_b/(sqrt(ecm+Q))
  rho_f=rho_b*(sqrt(ecm+Q))
  P4 <- coulomb_wave_FG(eta_f, rho_f, lb, k=0)
  # penetration and shift factor
  p_f <- rho_f/(P4$val_F^2 + P4$val_G^2)
  s_f <- rho_f*(P4$val_F*P4$val_Fp + P4$val_G*P4$val_Gp)/(P4$val_F^2 + P4$val_G^2)
  # shift factor at energy E0+Q
  xeta_f=eta_b/(sqrt(e0+Q))
  xrho_f=rho_b*(sqrt(e0+Q))
  PX2 <- coulomb_wave_FG(xeta_f, xrho_f, lb, k=0)
  b_f <- xrho_f*(PX2$val_F*PX2$val_Fp + PX2$val_G*PX2$val_Gp)/(PX2$val_F^2 + PX2$val_G^2)
  # partial width
  Gb <- 2*gb*p_f

  tapp <- (s_i-b_i)*ga + (s_f-b_f)*gb

  s1=pek*omega*Ga*Gb
  s2=((e0-ecm-tapp)^2)+0.25*((Ga+Gb)^2)
  SG <- (s1/s2)*(1/ecm)

  return(SG = SG)
}



