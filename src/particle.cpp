/*
 * File: particle.cpp
 * Description: Implementation of the particle classes
 * -------------------------------------------------------
 *
 * John R. Dowdle
 * The University of Texas at Austin
 * June 2008
 *
 */

#include "particle.h"

/********************************************************
 *  Particle                                            *
 *******************************************************/

/************** Public *********************************/

Particle::Particle(){
  // assume 3D LJ
  //zero_coords(3);
  assign_params(0, HS_SIGMA, LJ_SIGMA, LJ_EPSILON, JG_LAMBDA0, 
		JG_LAMBDA1, JG_LAMBDA2, JG_EPSILON1, JG_EPSILON2);
}

Particle::Particle(int p, double hs_sigma, double lj_sigma,
		   double lj_epsilon, double jg_lambda0,
		   double jg_lambda1, double jg_lambda2,
		   double jg_epsilon1, double jg_epsilon2){
  // ptype = p;
  init_params();
  assign_params(p, hs_sigma, lj_sigma, lj_epsilon, jg_lambda0, 
		jg_lambda1, jg_lambda2, jg_epsilon1, jg_epsilon2);
}


Particle::~Particle(){

}

void Particle::init_params(){
  // initialize all parameters
  lambda_0 = -1.;
  lambda_1 = -1.;
  lambda_2 = -1.;
  sigma = -1.;
  epsilon = -1.;
  epsilon_1 = -1.;
  epsilon_2 = -1.;
  mass = -1.;
  fixed = false;
}

void Particle::assign_params(int p, double hs_sigma, double lj_sigma, double lj_epsilon,
			     double jg_lambda0, double jg_lambda1, double jg_lambda2,
			     double jg_epsilon1, double jg_epsilon2){
  ptype = p;

  switch (p){
  case 0:
    // Argon params by default
    epsilon = lj_epsilon;   // [kJ/mol]
    sigma = lj_sigma;     // [nm]
    mass = 39.94800;         // [u]
    break;
  case 1:
    //Hard Sphere
    sigma = hs_sigma;         // [nm]
    mass = 18.0153;         // [u]
    break;
  case 2:
    // Jagla
    // (3D params from Jagla's 2001 LL paper)
    epsilon_2 = jg_epsilon2;         // [kJ/mol]
    epsilon_1 = jg_epsilon1;  // [kJ/mol]
    lambda_0 = jg_lambda0;           // [nm]  first solvation shell of water
    lambda_1 = jg_lambda1;    // [nm]  second ""
    lambda_2 = jg_lambda2;     // [nm]
    mass = 18.0153;              // [u]
    sigma = lambda_0;            // for mixtures
    epsilon = fabs(epsilon_2);   // "
    break;
  case 3:
    // Square-Well
    sigma = hs_sigma;
    lambda = 1.5*sigma;
    epsilon = 1.0;
    mass = 18.0153;         // [u]
    break;
  default:
    //Lennard-Jones OPLS-AA params for Argon
    epsilon = 9.78638e-01;   // [kJ/mol]
    sigma = 3.40100e-01;     // [nm]
    mass = 39.94800;         // [u]
  }
}


/************** Private ********************************/

