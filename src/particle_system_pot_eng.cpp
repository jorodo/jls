/*
 * File: particle_system_eng.cpp
 * Description: Implementation of particle system base class
 * for MC/MD sims.  Potential energy and related functions.
 * ---------------------------------------------------------
 *
 * John R. Dowdle
 * The University of Texas at Austin
 * June 2008
 * GPLv3.0+
 *
 */

#include "particle_system.h"


/***************** Potential Energy ******************************************/

double ParticleSystem::get_pair_energy(int i, int j){
  double eng;
  int iptype, jptype;

  iptype = particles[i].ptype;
  jptype = particles[j].ptype;

  // either we have a lj-lj interaction, a jg-jg interaction,
  // or a hs-hs interaction.  no others.
  if (iptype == 0 && jptype == 0){
    eng = LJ_LJ_Energy(i,j);
  } else if (iptype == 2 && jptype == 2){
    eng = JG_JG_Energy(i,j);
  } else if (iptype == 3 && jptype == 3){
    eng = SW_SW_Energy(i,j);
  } else {
    eng = HS_HS_Energy(i,j);
  }
  
  return eng;;
}


double ParticleSystem::get_particle_energy(int i, int jmin){
  // double zpos;
  double eng = 0.0;

  if (use_cell_list){  
    eng = get_cell_list_particle_energy(i, jmin);
  } else {   
    eng = get_full_particle_energy(i);
  }

  return eng;
}

// simple N^2 algorithm
double ParticleSystem::get_full_particle_energy(int i){
  double eng = 0.;

  for (int j=0; j<natoms; j++){
    if (j==i){
      continue;
    }
    else {
      eng += get_pair_energy(i,j);
    }
  }

  return eng;
}

// use the cell list method for energy evaluation
double ParticleSystem::get_cell_list_particle_energy(int i, int jmin){
  double eng = 0.0;
  std::vector<double> r;
  int j;

  r.resize(ndims);

  // create a double vector of the particle position to pass to the cell list
  for (int k=0; k<ndims; k++){
    r[k] = gsl_matrix_get(coords, i, k);
  }

  // get the neighbors of the particle (includes particle itself)
  cell_list.getNeighbors(r, neighbor_list, nneighbors);

  // sanity check
  assert(nneighbors <= cell_list.maxNeighbors);

  // now evaluate the energy from the neighbor list
  for (int n = 0; n < nneighbors; n++){
    j = neighbor_list[n];
    if (i != j && j >= jmin){
      eng += get_pair_energy(i, j);
    }
  }

  return eng;
}


double ParticleSystem::get_system_energy(){
  double eng = 0.0;
  int jmin = 0;

  // if (compute_pressure){  
  //   p_flag = true;
  //   virial = 0.;
  //   if ((particles[0].ptype==1 || particles[0].ptype==2) && is_homogeneous) {
  //     gsl_vector_set_zero(ndelta_r_cyc);
  //   } else if (particles[0].ptype==3 && is_homogeneous){
  //     gsl_vector_set_zero(ndelta_r_cyc);
  //     gsl_vector_set_zero(ndelta_rlm_cyc);
  //     gsl_vector_set_zero(ndelta_rlp_cyc);
  //   }
  // }

  if (use_cell_list){
    for (int i=0; i<natoms; i++){
      jmin = i + 1;
      eng += get_particle_energy(i, jmin);
    }
  } else {
    for (int i=0; i<natoms; i++){
      for (int j=i+1; j<natoms; j++){
	eng += get_pair_energy(i,j);
      }
    }
  }
  eng += (tail_cor?(natoms*eng_cor):0.0);
  // if (compute_pressure)
  //   pressure = get_system_pressure();

  return eng;
}


double ParticleSystem::get_sequence_energy(){
  double eng = 0.;
  int jmin = 0;

  if (use_cell_list){
    for (int i=0; i < seqLength; i++){
      jmin = i + 1;
      eng += get_particle_energy(i, jmin);
    }
  } else {
    for (int i=0; i < seqLength; i++){
      // interaction with other sequence members
      for (int j=i+1; j < seqLength; j++){
	eng += get_pair_energy(i,j);
      }

      // interaction with solvent
      for (int j=seqLength; j < natoms; j++){
	eng += get_pair_energy(i,j);
      }
    }
  }

  return eng;
}


double ParticleSystem::wall_overlap_test(){
  double zpos, eng;

  for (int i=0; i<natoms; i++){
    zpos = gsl_matrix_get(coords,i,ndims-1);
    if (zpos > hw_low_limit && zpos < hw_high_limit){
      overlap = true;
      eng = numeric_limits<double>::infinity();
      return eng;
    }
  }

  return 0.;
}


double ParticleSystem::get_system_pressure(){
  double pres, c1, cov11, sumsq; 
  //double ndrk, ndrk1, ndrcyck;

  pres = 0.;

  // discontinous corrections to virial
  if ((particles[0].ptype==1 || particles[0].ptype==2) && is_homogeneous){ // hard-sphere or jagla
      //for (int k=0; k<int(delta_r_size); k++){
      //	ndrcyck = gsl_vector_get(ndelta_r_cyc,k);
      //	ndrk = gsl_vector_get(ndelta_r,k);
      //	ndrk1 = ((ncycles-start_cycle-1)*ndrk + ndrcyck)/(ncycles-start_cycle);
      //	gsl_vector_set(ndelta_r,k,ndrk1);
      //}
      //gsl_fit_mul(delta_r->data,1,ndelta_r->data,1,delta_r_size,&c1,&cov11,&sumsq);
    gsl_fit_mul(delta_r->data,1,ndelta_r_cyc->data,1,delta_r_size,&c1,&cov11,&sumsq);
    virial += sigma_ij*c1/beta;
  } else if (particles[0].ptype==3 && is_homogeneous){  // square-well
    gsl_fit_mul(delta_r->data,1,ndelta_r_cyc->data,1,delta_r_size,&c1,&cov11,&sumsq);
    virial += sigma_ij*c1/beta;
    gsl_fit_mul(delta_rl->data,1,ndelta_rlm_cyc->data,1,delta_r_size,&c1,&cov11,&sumsq);
    virial -= lambda_ij*c1/beta;
    gsl_fit_mul(delta_rl->data,1,ndelta_rlp_cyc->data,1,delta_r_size,&c1,&cov11,&sumsq);
    virial += lambda_ij*c1/beta;
  }

  pres = density/beta
    + 1./(ndims*volume)*virial + (tail_cor?p_cor:0.);

  return pres;
}


void ParticleSystem::set_tail_corrections(int ptype){
    fprintf(logfile,"applying tail correction ...\n");
    if (ptype==0){
      LJ_LJ_tail_corrections();
    } else {
      fprintf(logfile,"warning:  tail corrections undefined for ptype %d.\n",ptype);
    }
}


double ParticleSystem::LJ_LJ_Energy(int i, int j){
  double r_ij_sqrd, eng;
  double core_term, decay_term;
  vector< int >::iterator location;

  eng = 0.;

  // if (!is_homogeneous){  # FIXME:  allow for mixtures
  //   lj_mixing_rule(i,j);
  // }

  r_ij_sqrd = pair_distance_sqrd(i,j);

  // if (polymer) {  // if there is a polymer in the system
  //   if ((i < polymer_size) && (j < polymer_size)){  // if both i and j are monomers
  //     // look for a bond between i and j
  //     location = find( particles[i].bonds.begin(), particles[i].bonds.end(), j );
  //     if (location != particles[i].bonds.end()){  // if a bond was found
  // 	if (r_ij_sqrd > max_bond_length*max_bond_length){  // test for overextension
  // 	  overlap = true;
  // 	  eng = numeric_limits<double>::infinity();
  // 	  return eng;
  // 	}
  //     }
  //   }
  // }

  if (r_ij_sqrd<cutoff_sqrd){
    core_term = gsl_pow_int((sigma_ij_sqrd/r_ij_sqrd),6); 
    decay_term = gsl_pow_int((sigma_ij_sqrd/r_ij_sqrd),3);
    eng = 4*epsilon_ij*(core_term - decay_term) - (shift?cutoff_eng:0.0);
    // update the virial if required
    if (p_flag){
      virial += 24.*epsilon_ij*(2.*core_term - decay_term);
    }
  } else {
    eng = 0.;
  }
  
  return eng;
}


double ParticleSystem::JG_JG_Energy(int i, int j){
  double dr, r_ij, drk, ndrcyck;
  double eng;
  vector< int >::iterator location;

  eng = 0.;

  r_ij = pair_distance(i,j);
  dr = 0.01*lambda_0_ij;  // FIXME: must be a better way ...

  // if (polymer) {  // if there is a polymer in the system
  //   if ((i < polymer_size) && (j < polymer_size)){  // if both i and j are monomers
  //     // look for a bond between i and j
  //     location = find( particles[i].bonds.begin(), particles[i].bonds.end(), j );
  //     if (location != particles[i].bonds.end()){  // if a bond was found
  // 	if (r_ij > max_bond_length){  // test for overextension
  // 	  overlap = true;
  // 	  eng = numeric_limits<double>::infinity();
  // 	  return eng;
  // 	}
  //     }
  //   }
  // }

  if (r_ij < lambda_0_ij){
    // --------------------- debug --------------------------
    // cout << "------ rank " << rank << " ------- ncycles: " << ncycles << endl;
    // cout << "i,j: " << i << "," << j << endl;
    // cout << "ptypei,ptypej: " << particles[i].ptype << "," << particles[j].ptype << endl;
    // cout << "xi: " << gsl_matrix_get(coords, i, 0) << endl;
    // cout << "yi: " << gsl_matrix_get(coords, i, 1) << endl;
    // cout << "zi: " << gsl_matrix_get(coords, i, 2) << endl;
    // cout << "xj: " << gsl_matrix_get(coords, j, 0) << endl;
    // cout << "yj: " << gsl_matrix_get(coords, j, 1) << endl;
    // cout << "zj: " << gsl_matrix_get(coords, j, 2) << endl;
    // cout << "r_ij: " << r_ij << endl;
    // cout << "lambda_0_ij: " << lambda_0_ij << endl;
    // cout << "--------------------" << endl;
    // cout.flush();
    // -------------------------------------------------------
    overlap = true;
    eng = numeric_limits<double>::infinity();
  } else if (r_ij < lambda_1_ij){
    eng = m1*r_ij + b1;
    if (p_flag) {
      virial += -m1*r_ij;
      for (int k=int(delta_r_size)-1; k>=0; k--){
	drk = gsl_vector_get(delta_r,k);
	ndrcyck = gsl_vector_get(ndelta_r_cyc,k);
	if (r_ij < lambda_0_ij + drk)
	  gsl_vector_set(ndelta_r_cyc,k,ndrcyck+1);
	else 
	  break;
      }
    }
  } else if (r_ij < lambda_2_ij){
    eng = m2*r_ij + b2;
    if (p_flag) {
      virial += -m2*r_ij;
    }
  } else {
    eng = 0.;
  }

  return eng;
}


// only homogenous SW systems right now
double ParticleSystem::SW_SW_Energy(int i, int j){
  double r_ij, eng, drk, ndrk, ndrlk;
  double drlk;
  vector< int >::iterator location;
  
  //double dr1;

  eng = 0.;



  r_ij = pair_distance(i,j);

  // if (polymer) {  // if there is a polymer in the system
  //   if ((i < polymer_size) && (j < polymer_size)){  // if both i and j are monomers
  //     // look for a bond between i and j
  //     location = find( particles[i].bonds.begin(), particles[i].bonds.end(), j );
  //     if (location != particles[i].bonds.end()){  // if a bond was found
  // 	if (r_ij > max_bond_length){  // test for overextension
  // 	  overlap = true;
  // 	  eng = numeric_limits<double>::infinity();
  // 	  return eng;
  // 	}
  //     }
  //   }
  // }

  if (r_ij < sigma_ij){
    overlap = true;
    eng = numeric_limits<double>::infinity();
  } else if (r_ij < lambda_ij){
    eng = -epsilon_ij;
    if (p_flag) {
      for (int k=int(delta_r_size)-1; k>=0; k--){
	drk = gsl_vector_get(delta_r,k);
	drlk = gsl_vector_get(delta_rl,k);
	ndrk = gsl_vector_get(ndelta_r_cyc,k);
	ndrlk = gsl_vector_get(ndelta_rlm_cyc,k);
	if (r_ij < sigma_ij + drk)
	  gsl_vector_set(ndelta_r_cyc,k,ndrk+1);
	else if (r_ij > lambda_ij - drlk)
	  gsl_vector_set(ndelta_rlm_cyc,k,ndrlk+1);
	else
	  break;
      }
    }
  } else {
    if (p_flag) {
      for (int k=int(delta_r_size)-1; k>=0; k--){
	drlk = gsl_vector_get(delta_rl,k);
	ndrlk = gsl_vector_get(ndelta_rlp_cyc,k);
	if (r_ij < lambda_ij + drlk)
	  gsl_vector_set(ndelta_rlp_cyc,k,ndrlk+1);
      }
    }
    eng = 0.;
  }

  return eng;
}


double ParticleSystem::HS_HS_Energy(int i, int j) {
  double r_ij, sigma_i, sigma_j, eng, drk, dr;
  double ndrcyck, exclusion_dist;
  vector< int >::iterator location;

  eng = 0.;

  sigma_i = particles[i].sigma;
  sigma_j = particles[j].sigma;

  exclusion_dist = (sigma_i + sigma_j)/2.;
  dr = 0.01*exclusion_dist;
  r_ij = pair_distance(i,j);

  // if (polymer) {  // if there is a polymer in the system
  //   if ((i < polymer_size) && (j < polymer_size)){  // if both i and j are monomers
  //     // look for a bond between i and j
  //     location = find( particles[i].bonds.begin(), particles[i].bonds.end(), j );
  //     if (location != particles[i].bonds.end()){  // if a bond was found
  // 	if (r_ij > max_bond_length){  // test for overextension
  // 	  overlap = true;
  // 	  eng = numeric_limits<double>::infinity();
  // 	  return eng;
  // 	}
  //     }
  //   }
  // }

  if (r_ij<exclusion_dist){
    // --------------------- debug --------------------------
    // cout << "------ rank " << rank << " ------- ncycles: " << ncycles << endl;
    // cout << "i,j: " << i << "," << j << endl;
    // cout << "ptypei,ptypej: " << particles[i].ptype << "," << particles[j].ptype << endl;
    // cout << "xi: " << gsl_matrix_get(coords, i, 0) << endl;
    // cout << "yi: " << gsl_matrix_get(coords, i, 1) << endl;
    // cout << "zi: " << gsl_matrix_get(coords, i, 2) << endl;
    // cout << "xj: " << gsl_matrix_get(coords, j, 0) << endl;
    // cout << "yj: " << gsl_matrix_get(coords, j, 1) << endl;
    // cout << "zj: " << gsl_matrix_get(coords, j, 2) << endl;
    // cout << "r_ij: " << r_ij << endl;
    // cout << "exclusion_dist: " << exclusion_dist << endl;
    // cout << "--------------------" << endl;
    // cout.flush();
    // -------------------------------------------------------
    overlap = true;
    eng = numeric_limits<double>::infinity();
  } else {
    eng = 0.;
    if (p_flag) {
      for (int k=int(delta_r_size)-1; k>=0; k--){
	drk = gsl_vector_get(delta_r,k);
	ndrcyck = gsl_vector_get(ndelta_r_cyc,k);
	if (r_ij < exclusion_dist + drk)
	  gsl_vector_set(ndelta_r_cyc,k,ndrcyck+1);
	else 
	  break;
      }
    }
  }

  return eng;
}


void ParticleSystem::LJ_LJ_tail_corrections(){
  double Cn, decay_term, core_term;

  Cn = nsphere_SA_constant(ndims);

  if (ndims>6){
    fprintf(stderr,"error:  cannot use tail correction for LJ with ndims > 6");
    fprintf(logfile,"error:  cannot use tail correction for LJ with ndims > 6");
    exit(USAGE);
  }

  decay_term = gsl_pow_int(cutoff,ndims-6)/(ndims-6);
  core_term = gsl_pow_int(cutoff,ndims-12)/(ndims-12);
  eng_cor = 2.*epsilon_ij*density*ndims*Cn*(gsl_pow_int(sigma_ij,6)*decay_term
					- gsl_pow_int(sigma_ij,12)*core_term);
  p_cor = 12.*epsilon_ij*gsl_pow_int(density,2)*Cn
    *(gsl_pow_int(sigma_ij,6)*decay_term - 2.*gsl_pow_int(sigma_ij,12)*core_term);
}



void ParticleSystem::update_cutoff_energies(){
  double min_sl, core_term, decay_term;
  if (box_length_cutoff){ // only if cutoff = L/2
    min_sl = *min_element(side_lengths.begin(), side_lengths.end());
    cutoff = min_sl/2.;
    cutoff_sqrd = cutoff*cutoff;
    if (shift){
      core_term = gsl_pow_int((sigma_ij_sqrd/cutoff_sqrd),6); 
      decay_term = gsl_pow_int((sigma_ij_sqrd/cutoff_sqrd),3);
      cutoff_eng = 4.*epsilon_ij*(core_term - decay_term);
    }
  }
  if (tail_cor){
    set_tail_corrections();
  }
}
