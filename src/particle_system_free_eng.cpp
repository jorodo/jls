/*
 * File: particle_system_free_eng.cpp
 * Description: Implementation of particle system base class
 * for free energy calculations
 * ---------------------------------------------------------
 *
 * John R. Dowdle
 * The University of Texas at Austin
 * June 2008
 * GPLv3.0+
 *
 */

#include "particle_system.h"


void ParticleSystem::TestParticleInsertion(){
  xdrcoordfile = xdrfile_open(traj_analysis_file.c_str(),"r");
  int result;

  v_sum = 0.;
  start_cycle = 0;
  ncycles = 0;
  tpi_eng_sum = 0.;
  tpi_eng = 0.;
  logfile = fopen(logfile_name.c_str(),"a");
  fprintf(logfile,"Running Test Particle insertion simulation on trajectory file %s",
	  traj_analysis_file.c_str());

  if ( xdrcoordfile != NULL ){
    write_tpi_eng();
    while (true){
      result = ReadXTCFrame(xdrcoordfile);
      if (result != 0) break;  // eof
      v_sum += volume;
      ncycles++;

      if (frame_num%insertion_freq==0){
	// attempt 'ninsertions' randomly located insertions
	tpi_eng_sum = 0.;
	for (int i=0; i<ninsertions; i++){
	  InsertParticle();
	}
	write_tpi_eng();
      }
    }
  } else {
    fileio_error(traj_analysis_file.c_str());
  }

  FinalInsertionStats();
  fflush(NULL);
  fclose(logfile);
}


// calculate \mu_{ex} for inserting a cavity volume defined by two
// spherical cavities separated by a cavity pair distance.
// void ParticleSystem::CavityPMF(){
//   xdrcoordfile = xdrfile_open(traj_analysis_file.c_str(),"r");
//   int result;
// 
//   start_cycle = 0;
//   ncycles = 0;
//   logfile = fopen(logfile_name.c_str(),"a");
//   fprintf(logfile,"Running cavity pmf analysis on trajectory file %s",
// 	  traj_analysis_file.c_str());
// 
//   if ( xdrcoordfile != NULL ){
//     write_tpi_eng();
//     while (true){
//       result = ReadXTCFrame(xdrcoordfile);
//       if (result != 0) break;  // eof
//       ncycles++;
// 
//       if (frame_num%insertion_freq==0){
// 	// attempt 'ninsertions' randomly located insertions
// 	nwsucc_cfg = 0;
// 	for (int i=0; i<ninsertions; i++){
// 	  CavityPairTest();
// 	}
// 	write_tpi_eng();
//       }
//     }
//   } else {
//     fileio_error(traj_analysis_file.c_str());
//   }
// 
//   fflush(NULL);
//   fclose(logfile);
// }


void ParticleSystem::CavityAnalysis(){
  xdrcoordfile = xdrfile_open(traj_analysis_file.c_str(),"r");
  int result;
  start_cycle = 0;
  ncycles = 0;
  logfile = fopen(logfile_name.c_str(),"a");
  fprintf(logfile,"Running Cavity Analysis simulation on trajectory file %s",
	  traj_analysis_file.c_str());

  if ( xdrcoordfile != NULL ){
    write_tpi_eng();

    while (true){
      result = ReadXTCFrame(xdrcoordfile);
      if (result != 0) break;  // eof
      // zero the insertion params for frame
      ncenters_cfg_sum = 0;
      ncenters_sqrd_cfg_sum = 0;
      for (int i=0; i<ncenters_max; i++){
	pndist_cfg[i] = 0.;
      }
	
      ncycles++;

      if (frame_num%insertion_freq==0){
	// attempt 'ninsertions' randomly located cavity searches
	for (int i=0; i<ninsertions; i++){
	  if (calculation=="cavity_pmf"){
	    CavityPairTest();
	  }
	  else
	    CavityTest();
	  ncenters_cfg_sum += ncenters;
	  ncenters_sqrd_cfg_sum += gsl_pow_int(ncenters,2);
	  ncenters_sum += ncenters;
	  ncenters_sqrd_sum += gsl_pow_int(ncenters,2);
	}
	// get insertion means for current frame
	ncenters_mean = ncenters_cfg_sum/(double)ninsertions;
	ncenters_sqrd_mean = ncenters_sqrd_cfg_sum/(double)ninsertions;
	for (unsigned int i=0; i<pndist.size(); i++){
	  pndist_cfg[i] /= (double)ninsertions;
	}
      }

      write_tpi_eng();

    }
  } else {
    fileio_error(traj_analysis_file.c_str());
  }

  FinalCavityStats();
  fflush(NULL);
  fclose(logfile);
}


void ParticleSystem::InsertParticle(){
  overlap = false;
  for (int k=0; k<ndims; k++){
    gsl_matrix_set(coords,natoms,k,gsl_rng_uniform(rng)*side_lengths[k]);
  }
  test_eng = get_particle_energy(natoms);
  if (particles[natoms].ptype==0 && tail_cor)  // FIXME:  gerry-rigged
    test_eng += 2.*eng_cor;

  nwtrials++;
  if (overlap){
    //if (verbose && ncycles%stdout_freq==0)
    //fprintf(logfile,"insertion rejected");
    tpi_eng = 0.;
  }
  else {
    //if (verbose && ncycles%stdout_freq==0)
    //fprintf(logfile,"insertion accepted");
    nwsucc++;
    tpi_eng = ensemble=="npt"?volume*exp(-beta*test_eng):exp(-beta*test_eng);
    tpi_eng_sum += tpi_eng;
  }      
}


void ParticleSystem::CavityTest(){
  double r_ij;
  ncenters = 0;
  overlap = false;

  for (int k=0; k<ndims; k++){
    gsl_matrix_set(coords,natoms,k,gsl_rng_uniform(rng)*side_lengths[k]);
  }
  for (int j=0; j<natoms; j++){
    r_ij = pair_distance(natoms,j);
    if (r_ij < cavity_radius){
      ncenters++;
    }
  }
  if (ncenters < (int)pndist.size()){
    pndist[ncenters]++;
    pndist_cfg[ncenters]++;
  } else {
    fprintf(logfile,"ncenters %d out of histogram range",ncenters);
  }
  if (ncenters > 0)
    overlap = true;
  nwtrials++;
  if (overlap){
    //if (verbose && ncycles%stdout_freq==0)
    //fprintf(logfile,"insertion rejected");
    tpi_eng = 0.;
  }
  else {
    //if (verbose && ncycles%stdout_freq==0)
    //fprintf(logfile,"insertion accepted");
    nwsucc++;
    //tpi_eng = ensemble=="npt"?volume*exp(-beta*test_eng):exp(-beta*test_eng);
    //tpi_eng_sum += tpi_eng;
  }
}


void ParticleSystem::CavityPairTest(){
  double r_ij0, r_ij1, r, xk, ran_gauss;
  gsl_vector *x;
  double sum=0.;

  ncenters = 0;
  overlap = false;

  x = gsl_vector_alloc(ndims);
  gsl_vector_set_zero(x);

  // center of the first cavity search
  for (int k=0; k<ndims; k++){
    gsl_matrix_set(coords,natoms,k,gsl_rng_uniform(rng)*side_lengths[k]);
  }  
  
  // center of the second cavity search
  // uniformly randomly generate a point on the surface of the n-ball
  for (int k=0; k<ndims; k++){
    ran_gauss = gsl_ran_gaussian(rng,1.);
    gsl_vector_set(x,k,ran_gauss);
    sum += ran_gauss*ran_gauss;
  }
  r = pow(sum,1./2.);
  gsl_vector_scale(x,1./r);
  gsl_vector_scale(x,cavity_pair_dist);

  // center the point on first cavity
  for (int k=0; k<ndims; k++){
    xk = gsl_matrix_get(coords,natoms,k) + gsl_vector_get(x,k);
    xk = modulo(xk, side_lengths[k]);
    // if (xk < 0.0) xk += side_lengths[k];     // PBC
    // if (xk >= side_lengths[k]) xk -= side_lengths[k];
    gsl_matrix_set(coords, natoms+1, k, xk);
  }

  // test for occupancy
  for (int j=0; j<natoms; j++){
    r_ij0 = pair_distance(natoms,j);
    r_ij1 = pair_distance(natoms+1,j);
    if (r_ij0 < cavity_radius)
      ncenters++;
    if (r_ij1 < cavity_radius)
      ncenters++;
  }

  // record stats
  if (ncenters < (int)pndist.size()){
    pndist[ncenters]++;
    pndist_cfg[ncenters]++;
  } else {
    fprintf(logfile,"ncenters %d out of histogram range\n",ncenters);
  }
  if (ncenters > 0)
    overlap = true;
  nwtrials++;
  if (overlap){
    tpi_eng = 0.;
  }
  else {
    nwsucc++;
  }

  gsl_vector_free(x);
}


void ParticleSystem::FinalInsertionStats(){
  engfile = fopen(tpi_engfile_name.c_str(),"a");
  v_mean = v_sum/ncycles;
  mu_ex = ensemble=="npt"?-1./beta*log(1./v_mean*tpi_eng_sum/nwtrials):
    -1./beta*log(tpi_eng_sum/nwtrials);
  fprintf(engfile,"# excess_chemical_potential:  %0.6e [kJ mol^-1]\n",mu_ex);
  fprintf(logfile,"# done with tpi insertion calculation");
  fclose(engfile);
}


void ParticleSystem::FinalCavityStats(){
  engfile = fopen(tpi_engfile_name.c_str(),"a");
  ncenters_mean = ncenters_sum/nwtrials;
  ncenters_sqrd_mean = ncenters_sqrd_sum/nwtrials;
  fprintf(engfile,"# ---------------------------------------\n");
  fprintf(engfile,"# pn distribution:\n");
  for (unsigned int i=0; i<pndist.size(); i++){
    pndist[i] /= (double)nwtrials;
    fprintf(engfile,"# \tp%d = %0.6f\n",i,pndist[i]);
  }
  fprintf(engfile,"# ---------------------------------------\n");
  fclose(engfile);
  //fprintf(engfile,"# <n> = %0.6e\n",ncenters_mean);
  //fprintf(engfile,"# <n^2> = %0.6e\n",ncenters_sqrd_mean);
}


//void ParticleSystem::StagedParticleGrowth(){
//  double eng;
//  ifstream traj_file(traj_analysis_file.c_str());
//
//  ng_att = 0;
//  ng_acc = 0;
//  ncycles = 0;
//  logfile = fopen(logfile_name.c_str(),"a");
//  fprintf(logfile,"nattempts\t naccepted\n");
//  fprintf(logfile,"-----------------------\n");
//
//  if (traj_file.is_open()){
//    while (true){
//      ReadPDBFrame(traj_file);
//      if (traj_file.eof()) break;
//      ncycles++;
//      overlap = false;
//      eng = get_particle_energy(growth_idx);
//      if (!overlap){
//	ng_att++;
//	ng_acc++;
//      } else {
//	ng_att++;
//      }
//      fprintf(logfile,"%d\t%d\n",ng_att,ng_acc);
//
//    }
//  } else {
//    fileio_error(traj_analysis_file.c_str());
//  }
//
//  fprintf(logfile,"-----------------------\n");
//  fprintf(logfile,"nframes:  %d\n",ncycles);
//  fprintf(logfile,"ng_acc:  %d\n",ng_acc);
//  fprintf(logfile,"ng_att:  %d\n",ng_att);
//  fprintf(logfile,"p_0(d_n|d_{n-1}):  %0.6f\n",(ng_acc/(double)ng_att));
//  fclose(logfile);
//}


void ParticleSystem::TrialGrowth(){
  double sigma0, part_eng;
  
  ng_att++;
  ng_att_cycle++;
  sigma0 = particles[growth_idx].sigma;
  particles[growth_idx].sigma = growth_sigma;
  overlap = false;
  part_eng = get_particle_energy(growth_idx);
  if (!overlap){
    ng_acc++;
    ng_acc_cycle++;
    //cout << "no ovrlp" << endl;
    //cout << "ncycles " << ncycles << endl;
    //cout << "rst_freq " << rst_freq << endl;
    //cout << "growth_sigma " << growth_sigma << endl;
    //cout << "growth_sigma_avg " << growth_sigma_avg << endl;
    write_growth_rst_coords();
    write_rst();
    //if (ncycles%rst_freq==0 && growth_sigma>=growth_sigma_avg){
    //  cout << "writing gr rst coords" << endl;
    //  write_growth_rst_coords();
    //}
  }
  particles[growth_idx].sigma = sigma0;
  p_growth = ng_acc/(double)ng_att;
  //p_growth = ng_acc_cycle/(double)ng_att_cycle;
  growth_sigma_sum += growth_sigma;
}


void ParticleSystem::BennettMethodTrajectory(){
  xdrcoordfile = xdrfile_open(traj_analysis_file.c_str(),"r");
  double bm_eng0, bm_eng1, hw_eng;
  double old_volume;
  int result;

  if (bennett_transform==""){
    fprintf(logfile,"error: must specify bennett_transform\n");
    exit(USAGE);
  } else if (bennett_transform=="scale"){
    if (new_side_lengths.size()==0){
      fprintf(logfile,"error: must specify new_side_lengths\n");
      exit(USAGE);
    }
  }

  bm_eng0 = 0.;
  bm_eng1 = 0.;
  start_cycle = 0;
  ncycles = 0;
  logfile = fopen(logfile_name.c_str(),"a");
  fprintf(logfile,"Running Bennett's Method simulation on trajectory file %s",
	  traj_analysis_file.c_str());

  if ( xdrcoordfile != NULL ){
    write_bm_eng(bm_eng0,bm_eng1);
    while (true){
      result = ReadXTCFrame(xdrcoordfile);
      if (result != 0) break;  // eof

      old_volume = volume;
      ncycles++;
      overlap = false;
      bm_eng0 = get_system_energy();
      //if (overlap) continue;  // FIXME: possible overlap due to round errors
      if (bennett_transform=="scale"){
	save_boxsize();
	scale_coords(new_side_lengths);
	update_boxsize(new_side_lengths);
	update_cutoff_energies();
	//volume = get_system_volume();
      }
      overlap = false;
      bm_eng1 = get_system_energy();
      if (hard_wall)
	hw_eng = wall_overlap_test();
      //if (!overlap)          // don't write inf energies to file
      write_bm_eng(bm_eng0,bm_eng1);
      scale_coords(old_side_lengths);
      update_boxsize(old_side_lengths);
      update_cutoff_energies();
    }
  } else {
    fileio_error(traj_analysis_file.c_str());
  }

  fflush(NULL);
  fclose(logfile);
}


// BAR on a given frame
void ParticleSystem::BennettMethodFrame(){
  double bm_eng, hw_eng;

  bm_eng = 0.;
  frame_num = ncycles;

  if (bennett_transform=="scale"){
    save_boxsize();
    scale_coords(new_side_lengths);
    update_boxsize(new_side_lengths);
    update_cutoff_energies();
  }
  overlap = false;
  bm_eng = get_system_energy();
  if (hard_wall)
    hw_eng = wall_overlap_test();
  if (!overlap)          // don't write inf energies to file
    write_bm_eng(pot_eng,bm_eng);
  scale_coords(old_side_lengths);
  update_boxsize(old_side_lengths);
  update_cutoff_energies();

}


// generate a distribution of nearest particles distances
// from a randomly located point
void ParticleSystem::CavitySizeDistributions(){
  xdrcoordfile = xdrfile_open(traj_analysis_file.c_str(),"r");
  int result;
  double r_ij, min_distance;

  gsl_hist = gsl_histogram_alloc(nbins);
  double bins_min = 0.;   // [nm]
  double bins_max = 1.;  // [nm]

  gsl_histogram_set_ranges_uniform(gsl_hist, bins_min, bins_max);

  logfile = fopen(logfile_name.c_str(),"a");
  fprintf(logfile,"Running cavity cavity size distribution analysis on trajectory file %s",
	  traj_analysis_file.c_str());

  if ( xdrcoordfile != NULL ){
    while (true){
      // load the frame
      result = ReadXTCFrame(xdrcoordfile);
      if (result != 0) break;  // eof
      
      if (frame_num%insertion_freq==0){
	// attempt 'ninsertions' randomly located points
	// and bin each min_distance
	for (int i=0; i<ninsertions; i++){
	  // center of random point
	  for (int k=0; k<ndims; k++){
	    gsl_matrix_set(coords,natoms,k,gsl_rng_uniform(rng)*side_lengths[k]);
	  }

	  min_distance = side_lengths[0];
	  
	  // get the minimum distance
	  for (int j=0; j<natoms; j++){
	    r_ij = pair_distance(natoms,j);
	    if (r_ij < min_distance){
	      min_distance = r_ij;
	    }
	  }
	  gsl_histogram_increment(gsl_hist, min_distance);  // bin the min_distance
	}
	write_cav_dist();  // write the updated cavity dist to file
      }
    }
  } else {
    fileio_error(traj_analysis_file.c_str());
  }

  fflush(NULL);
  fclose(logfile);
}
