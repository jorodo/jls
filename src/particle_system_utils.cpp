/*
 * File: particle_system_utils.cpp
 * Description: Implementation of particle system base class
 * for MC/MD sims.  Utility functions for the PS class.
 * ---------------------------------------------------------
 *
 * John R. Dowdle
 * The University of Texas at Austin
 * June 2008
 * GPLv3.0+
 *
 */

#include "particle_system.h"

/******************* Class utils *******************************************/

void ParticleSystem::Run(){
  if (calculation=="mc") {
    if (ensemble=="nvt"){
      MC_NVT();
    } else if (ensemble=="npt"){
      MC_NPT();
    }
  } else if (calculation=="test_particle"){
    TestParticleInsertion();
  } else if (calculation=="cavity" || calculation=="cavity_pmf"){
    CavityAnalysis();
  //} else if (calculation=="cavity_pmf"){
  //  CavityPMF();
  } else if (calculation=="cavity_dist"){
    CavitySizeDistributions();
  } else if (calculation=="bennett"){
    BennettMethodTrajectory();
  }
}


double ParticleSystem::get_system_volume(){
  double V = 1.;
  for (int k=0; k<ndims; k++){
    V *= side_lengths[k];
  }
  return V;
}


double ParticleSystem::pair_distance(int i, int j){
  double dq_k, dist_sqrd=0.;
  for (int k=0; k<ndims; k++){
    dq_k = gsl_matrix_get(coords,i,k) - gsl_matrix_get(coords,j,k);
    dq_k -= side_lengths[k] * (double)round_(dq_k / side_lengths[k]);
    // if (dq_k > side_lengths[k]/2.0) dq_k -= side_lengths[k];        //minimum image
    // else if (dq_k < -side_lengths[k]/2.0) dq_k += side_lengths[k];
    dist_sqrd += dq_k*dq_k;
  }
  
  return sqrt(dist_sqrd);
}


double ParticleSystem::pair_distance_sqrd(int i, int j){
  double dq_k, dist_sqrd=0.;
  for (int k=0; k<ndims; k++){
    dq_k = gsl_matrix_get(coords,i,k) - gsl_matrix_get(coords,j,k);
    dq_k -= side_lengths[k] * (double)round_(dq_k / side_lengths[k]);
    // if (dq_k > side_lengths[k]/2.0) dq_k -= side_lengths[k];        //minimum image
    // else if (dq_k < -side_lengths[k]/2.0) dq_k += side_lengths[k];
    dist_sqrd += dq_k*dq_k;
  }
  
  return dist_sqrd;
}


void ParticleSystem::add_particles(){
  int count, index;
  index = 0;
  count = 0;

  if (polymer){
    assign_bonds();
  }

  for (map<int,int>::iterator ii=ptypes.begin(); ii!=ptypes.end(); ++ii){
    for (int j=0; j<(*ii).second; j++){
      particles.push_back(Particle((*ii).first,hs_sigma,lj_sigma,lj_epsilon,
				   jg_lambda0, jg_lambda1, jg_lambda2, jg_epsilon1,
				   jg_epsilon2));
      index = count + j;
    }
    count += (*ii).second;
  }

  if (calculation=="test_particle" || calculation=="cavity" || calculation=="cavity_dist"){
    particles.push_back(Particle(test_ptype,hs_sigma,test_sigma,test_epsilon));
  } else if (calculation=="cavity_pmf" ){
    particles.push_back(Particle(test_ptype,hs_sigma,test_sigma,test_epsilon));
    particles.push_back(Particle(test_ptype,hs_sigma,test_sigma,test_epsilon));
  }
}


double ParticleSystem::get_max_sigma(){
  double sigma_max = 0.;
  for (int i=0; i<natoms; i++){
    if (particles[i].sigma>sigma_max){
      sigma_max = particles[i].sigma;
    }
  }

  return sigma_max;
}


void ParticleSystem::dir_manage(){
  string new_name;
  if (!dir_exists(out_dir)){
    mkdir(out_dir.c_str(),0777);
  } else {
    bkup_file(coordfile_name);
    bkup_file(xdrcoordfile_name);
    bkup_file(corcoordfile_name);
    bkup_file(engfile_name);
    bkup_file(seqEngFileName);
    bkup_file(tpi_engfile_name);
    bkup_file(bm_engfile_name);
    bkup_file(growth_engfile_name);
    bkup_file(logfile_name);
  }

  if (!dir_exists(rst_dir))
    mkdir(rst_dir.c_str(),0777);

}


void ParticleSystem::lj_mixing_rule(int i, int j){
  double sigma_i, sigma_j, epsilon_i, epsilon_j;
  sigma_i = particles[i].sigma;
  epsilon_i = particles[i].epsilon;
  sigma_j = particles[j].sigma;
  epsilon_j = particles[j].epsilon;
  epsilon_ij = sqrt(epsilon_i*epsilon_j);
  if (use_geom_sigma){
    sigma_ij_sqrd = sigma_i*sigma_j;
    sigma_ij = sqrt(sigma_ij_sqrd);
  } else {
    sigma_ij = (sigma_i + sigma_j)/2.;
    sigma_ij_sqrd = sigma_ij*sigma_ij;
  }
}


void ParticleSystem::scale_coords(vector<double> &new_boxsize){
  double old_pos, scale_factor;

  if (hard_wall && hw_adjust){
    hw_thickness *= new_boxsize[ndims-1]/side_lengths[ndims-1];
    hw_high_limit = new_boxsize[ndims-1]/2. + hw_thickness/2.;
    hw_low_limit = new_boxsize[ndims-1]/2. - hw_thickness/2.;
  }

  for (int i=0; i<natoms; i++){
    for (int j=0; j<ndims; j++){
      scale_factor = new_boxsize[j]/side_lengths[j];
      old_pos = gsl_matrix_get(coords,i,j);
      gsl_matrix_set(coords,i,j,old_pos*scale_factor);
    }
  }
}


void ParticleSystem::update_boxsize(vector<double> &new_boxsize){
  for (int k=0; k<ndims; k++){
    side_lengths[k] = new_boxsize[k];
  }
}


void ParticleSystem::save_boxsize(){
  for (int k=0; k<ndims; k++){
    old_side_lengths[k] = side_lengths[k];
  }
}


void ParticleSystem::rescale_center_of_mass(){
  double com_diff_k, old_pos;

  compute_center_of_mass();

  for (int k=0; k<ndims; k++){
    com_diff_k = center_of_mass[k] - center_of_mass_ref[k];
    for (int i=0; i<natoms; i++){
      old_pos = gsl_matrix_get(coords,i,k);
      gsl_matrix_set(coords,i,k,old_pos-com_diff_k);
    }
  }

  compute_center_of_mass();
  applyPBC();
}


void ParticleSystem::compute_center_of_mass(){
  // this is actually the geometric center

  // zero COM
  for (int k=0; k<ndims; k++){
    center_of_mass[k] = 0.;
  }

  // compute numerator
  for (int i=0; i<natoms; i++){
    for (int k=0; k<ndims; k++){
      center_of_mass[k] += gsl_matrix_get(coords,i,k);
    }
  }

  // finalize
  for (int k=0; k<ndims; k++){
    center_of_mass[k] /= natoms;
  }
}




void ParticleSystem::applyPBC(){
  double old_pos;
  for (int i=0; i<natoms; i++){
    for (int k=0; k<ndims; k++){
      old_pos = gsl_matrix_get(coords, i, k);
      gsl_matrix_set(coords, i, k, modulo(old_pos, side_lengths[k]));
      // if (gsl_matrix_get(coords,i,k) < 0.){
      // 	old_pos = gsl_matrix_get(coords,i,k);
      // 	gsl_matrix_set(coords,i,k,old_pos+side_lengths[k]);
      // }
      // if (gsl_matrix_get(coords,i,k) >= side_lengths[k]){
      // 	old_pos = gsl_matrix_get(coords,i,k);
      // 	gsl_matrix_set(coords,i,k,old_pos-side_lengths[k]);
      // }
    }
  }
}


void ParticleSystem::sidelengths2xdrbox(matrix box){
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++){
      box[i][j] = 0.;
    }
  }

  for (int i=0; i<ndims; i++){
    box[i][i] = side_lengths[i];
  }
}


void ParticleSystem::xdrbox2sidelengths(matrix box){
  for (int i=0; i<ndims; i++){
    side_lengths[i] = box[i][i];
  }
}


void ParticleSystem::coords2xdr(rvec *xdrcoords){

  for (int i=0; i<natoms; i++){
    for (int j=0; j<ndims; j++){
      xdrcoords[i][j] = (float)gsl_matrix_get(coords,i,j);
    }
  }
}


void ParticleSystem::xdr2coords(rvec *xdrcoords){

  for (int i=0; i<natoms; i++){
    for (int j=0; j<ndims; j++){
      gsl_matrix_set(coords,i,j,(double)xdrcoords[i][j]);
    }
  }
}
