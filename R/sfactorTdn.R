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
#' @aliases sfactorTdn
#' @usage sfactorTdn(ECM = ECM, ER = ER, gi = gi, gf = gf)
#' @format \describe{
#' \item{x}{
#' The function has four arguments: ECM, ER, gi, gf}
#' }
#' @param ECM ECM
#' @param ER  ER
#' @param gi  gi
#' @param gf  gf
#' @return S-factor
#' @import gsl
#' @examples
#' library(nuclear)
#'
#' N <- 300
#' obsx1 <- exp(seq(log(1e-3), log(1),length.out=N))
#' plot(obsx1,sfactorTdn(obsx1,0.0912,2.93,0.0794),ylim=c(0,32.5),
#'      xlim=c(1e-3,1),col="red",cex=1.25,type="l",ylab="S-factor",xlab="E",log="x")
#'
#' @author Rafael de Souza, UNC,  and Christian Illiadis, UNC
#'
#' @keywords S-factor
#' @export
#'
#'
sfactorTdn <- function(ECM,ER,gi,gf){

# Constants
m1_i = 3.016; m2_i = 2.014;		# masses (amu) of t and d
m1_f = 4.0026; m2_f = 1.0087;	# masses (amu) of n and 4He
z1_i = 1; z2_i = 1;			# charges of t and d
z1_f = 2; z2_f = 0;				#charges of n and 4He
r_i = 6.0; r_f = 5.0;			# channel radii (fm)
jt = 0.5; jp = 1.0; jr = 1.5;			#spins of target, projectile, resonance
Q = 17.589;						#reaction Q-value (MeV)
la = 0; lb = 2;					#orbital angular momenta of d and n

#   DEFINITIONS

mue_i <- (m1_i*m2_i)/(m1_i + m2_i);
mue_f <- (m1_f*m2_f)/(m1_f + m2_f);
pek <- 6.56618216e-1/mue_i;
omega <- (2*jr + 1)/((2*jt + 1)*(2*jp + 1));

#     ----------------------------------------------------
#     PENETRABILITY AND SHIFT FUNCTION AT ER
#     ----------------------------------------------------


  eta_a = 0.15748927*z2_i*z1_i*sqrt(mue_i)
  rho_a = 0.218735097*r_i*sqrt(mue_i)
  reta_i = eta_a/(sqrt(ER))
  rrho_i = rho_a*(sqrt(ER))

  P1 <- coulomb_wave_FG(reta_i, rrho_i, la, k = 0)
  pr_i <- rrho_i/(P1$val_F^2 + P1$val_G^2)
  sr_i <- rrho_i*(P1$val_F*P1$val_Fp + P1$val_G*P1$val_Gp)/(P1$val_F^2 + P1$val_G^2)
  ga <- 2*gi*pr_i


  eta_b = 0.15748927*z2_f*z1_f*sqrt(mue_f)
  rho_b = 0.218735097*r_f*sqrt(mue_f)
  reta_f = eta_b/(sqrt(ER + Q))
  rrho_f = rho_b*(sqrt(ER + Q))

  P2 <- coulomb_wave_FG(reta_f, rrho_f, lb, k = 0)
  pr_f <- rrho_f/(P2$val_F^2 + P2$val_G^2)
  sr_f <- rrho_f*(P2$val_F*P2$val_Fp + P2$val_G*P2$val_Gp)/(P2$val_F^2 + P2$val_G^2)
  gb <- 2*gf*pr_f

# CALCULATE S-FACTOR

  etpe_i = exp(0.989534267*z1_i*z2_i*sqrt(mue_i/ECM))
  eta_i = eta_a/(sqrt(ECM))
  rho_i = rho_a*(sqrt(ECM))
  P3 <- coulomb_wave_FG(eta_i, rho_i, la, k = 0)
  p_i <- rho_i/(P3$val_F^2 + P3$val_G^2)
  s_i <- rho_i*(P3$val_F*P3$val_Fp + P3$val_G*P3$val_Gp)/(P3$val_F^2 + P3$val_G^2)
  prat_i = p_i/pr_i

  eta_f = eta_b/(sqrt(ECM + Q))
  rho_f = rho_b*(sqrt(ECM + Q))
  P4 <- coulomb_wave_FG(eta_f, rho_f, lb, k = 0)
  p_f <- rho_f/(P4$val_F^2 + P4$val_G^2)
  s_f <- rho_f*(P4$val_F*P4$val_Fp + P4$val_G*P4$val_Gp)/(P4$val_F^2 + P4$val_G^2)
  prat_f = p_f/pr_f

  tapp <- (s_i - sr_i)*gi + (s_f - sr_f)*gf

  s1 = pek*etpe_i*omega*prat_i*prat_f*ga*gb
  s2 = ((ER - ECM - tapp)^2) + 0.25*((ga*prat_i + gb*prat_f)^2)
  SF <- s1/s2
  return(SF = SF)
}

