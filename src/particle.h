/*
 * File: particle.h
 * Description: Description of the particle classes
 * -------------------------------------------------------
 *
 * John R. Dowdle
 * The University of Texas at Austin
 * June 2008
 *
 */
#ifndef _particle_h
#define _particle_h

#include <iostream>
#include <string>
#include <vector>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "utils.h"

using namespace std;


// default particle potential parameters
const double HS_SIGMA = 0.2800;
const double LJ_SIGMA = 0.3415;
const double LJ_EPSILON = 1.03931;
const double SW_SIGMA = 0.2800;
const double SW_EPSILON = 1.00;
const double SW_LAMBDA = 0.42;   // 1.5*sigma
const double JG_LAMBDA0 = 0.2800;
const double JG_LAMBDA1 = 0.4816;  // 1.72*lambda0
const double JG_LAMBDA2 = 0.8400;  // 3.*lambda0
const double JG_EPSILON1 = 3.5;
const double JG_EPSILON2 = -1.0;

class Particle {

 public:

  //int index;
  int ptype;
  //vector<double> coords;
  double mass;
  double sigma, epsilon;  // LJ, HS, SW
  double lambda;          // SW
  double lambda_0, lambda_1, lambda_2, epsilon_2, epsilon_1; // JG
  bool fixed;
  vector<int> bonds;

  Particle();
  Particle(int p, double hs_sigma=HS_SIGMA, double lj_sigma=LJ_SIGMA, 
	   double lj_epsilon=LJ_EPSILON, double jg_lambda0=JG_LAMBDA0,
	   double jg_lambda1=JG_LAMBDA1, double jg_lambda2=JG_LAMBDA2,
	   double jg_epsilon1=JG_EPSILON1, double jg_epsilon2=JG_EPSILON2);

  ~Particle();

  void init_params();
  void assign_params(int p, double hs_sigma=HS_SIGMA, double lj_sigma=LJ_SIGMA, 
	   double lj_epsilon=LJ_EPSILON, double jg_lambda0=JG_LAMBDA0,
	   double jg_lambda1=JG_LAMBDA1, double jg_lambda2=JG_LAMBDA2,
	   double jg_epsilon1=JG_EPSILON1, double jg_epsilon2=JG_EPSILON2);
};


#endif
