/*
 * File: particle_system_mc.cpp
 * Description: Implementation of particle system base class
 * for molec. sims.  MC calculations.
 * ---------------------------------------------------------
 *
 * John R. Dowdle
 * The University of Texas at Austin
 * June 2008
 * GPLv3.0+
 *
 */

#include "particle_system.h"

/******************* Monte Carlo calculations ********************************/


void ParticleSystem::MC_NVT(void){
  double minVisits, histMean;
  int imon, isolv, allSolvFreq, nmonSamples=10;

  logfile = fopen(logfile_name.c_str(),"a");
  
  if (ncycles==0){
    write_coords();
    write_eng();
  }

  if (seqMC){
    nHS = 0;  // only need to update nHS after mutations
    for (int iseq = 0; iseq < seqLength; iseq++){
      if (particles[iseq].ptype == 1) nHS++;
    }
    nmonSamples = seqLength;
    allSolvFreq = 10;
  }

  fprintf(logfile,"Running canonical ensemble MC sampling for %d cycles.\n",
	  ncycles_target);

  for (ncycles=start_cycle+1; ncycles<end_cycle+1; ncycles++){
    if (seqMC){
      assert((seqEng >= engHistMin) && (seqEng < engHistMax));
      // ----------------- sequence energy landscape sampling --------
      if (ncycles%allSolvFreq == 0) { 
	// -- MC moves for all solvent --
	for (int i=0; i < natoms - seqLength; i++){
	  isolv = seqLength + (int)gsl_rng_uniform_int(rng, natoms - seqLength);
	  trial_move_part(isolv);
	}
      } else {
	// -- MC moves for solvation shell solvent only --
	for (int i=0; i < seqLength; i++){
	  imon = (int)gsl_rng_uniform_int(rng, seqLength);
	  sampleNeighbors(imon);
	}
      }

      if (ncycles%seqMoveFreq == 0){
	// -- try some sequence moves --
	for (int i=0; i<seqLength; i++){
	  if (gsl_rng_uniform(rng) < 0.5)
	    trial_move_sequence_swap();
	  else
	    trial_move_sequence_mutate();
	}
      }

      if (seqWL || seqWLTM){  // Wang-Landau bookkeepping
	write_seq_eng();
	write_seq_entropy_hist();
	write_seq_eng_hist();
	write_seq_comp_hist();
	minVisits = gsl_histogram_min_val(seqEngHist);
	histMean = 0.;
	for (int bin=0; bin < nbins; bin++){
	  histMean += gsl_histogram_get(seqEngHist, bin);
	}
	histMean /= nbins;
	// if (nVisits > 20){  // Shell & Deb. criterion
	if (minVisits / histMean > 0.8){   // W&L criterion
	  gsl_histogram_reset(seqEngHist);
	  if (gModFactor*0.5 > numeric_limits<double>::epsilon())
	    gModFactor *= 0.5;
	}
	// if (gModFactor < gModTol) break; # <-- FIXME: only master process breaks now
      } else 
	  write_seq_eng();
      // --------------------------------------------------------------
    } else {
      // -------------------- usual MC sampling --------------------------------
      for (int ntrials=0; ntrials<nmoves; ntrials++){
	trial_move_part();
      }
      // -----------------------------------------------------------------------
    }

    pot_eng = get_system_energy();

    if (ncycles%cor_out_freq == 0) write_coords();
    if (ncycles%rst_freq == 0) write_rst();
    if (ncycles%eng_out_freq == 0) write_eng();
    if (ncycles%npart_adjust == 0 && adjust) {
      adjust_delta_max();
      if (!cubic_system) adjust_delta_zmax();
    }
    if (verbose && ncycles%stdout_freq==0) talk();

    // if (ncycles%rescaleCOM_freq == 0 && fixCOM) rescale_center_of_mass();
    // if (bennett_method) BennettMethodFrame();

    fflush(NULL);

  }

  fprintf(logfile,"\n");
  talk();
  fclose(logfile);
}


void ParticleSystem::MC_NPT(){

  logfile = fopen(logfile_name.c_str(),"a");

  if (ncycles==0){
    write_coords();
    write_eng();
  }

  fprintf(logfile,
	  "Running isothermal-isobaric ensemble MC sampling for %d cycles.\n"
	  ,ncycles_target);
  fprintf(logfile,"-------------------------------\n");

  for (ncycles=start_cycle+1; ncycles<end_cycle+1; ncycles++){    
    for (int ntrials = 0; ntrials<nmoves; ntrials++){
      if ((int)gsl_rng_uniform_int(rng,natoms+1)<natoms){
	trial_move_part();
      } else {
	pot_eng = get_system_energy();
	(cubic_system?trial_move_vol():trial_move_vol_rect());
      }
    }
    
    pot_eng = get_system_energy();

    if (ncycles%cor_out_freq == 0) write_coords();
    if (ncycles%rst_freq == 0) write_rst();
    if (ncycles%eng_out_freq == 0) write_eng();
    if (growth) {
      TrialGrowth();
      if (growth_equil && !overlap) break;
      if (ncycles%eng_out_freq == 0) write_growth_eng();
      if (adjust_growth && ng_att%ng_adjust==0) 
	adjust_growth_sigma();
    }
    if (adjust){
      if (npart_attempted%npart_adjust == 0)
	adjust_delta_max();
      if (nvol_attempted%nvol_adjust == 0 && nvol_att_cycle>0){
	if (cubic_system)
	  adjust_delta_lnV_max();
	else
	  adjust_delta_L_max();
      }
    }
    if (verbose && ncycles%stdout_freq==0)
      talk();
    fflush(NULL);
  }

  fclose(logfile);
}


void ParticleSystem::trial_move_part(int part){
  int atom_idx;
  double eng_old=0., eng_new=0., new_pos, delta;
  bool move_accepted, isSeqNeighborBefore, isSeqNeighborAfter, isSeqNeighbor;
  double seqEngOld, seqEngNew, S0, S1, tol=1.e-10;
  size_t binNum0, binNum1;
  int old_cell;
  int new_cell;
  std::vector<double> r;
  double distance;

  // avoid warnings
  isSeqNeighborBefore = false;
  isSeqNeighborAfter = false;
  move_accepted = false;
  seqEngOld = 0.;
  seqEngNew = 0.;

  // User submits the index of a specific particle
  // for trial displacement of specific particle,
  // or we randomly select a particle to displace.
  if (part >=0)
    atom_idx = part;
  else
    atom_idx = gsl_rng_uniform_int(rng,natoms);

  // Check for the "fixed" attribute on a particle,
  // do not attempt to displace it if set true.
  if (particles[atom_idx].fixed){
    return;
  }

  // determine if the particle is a sequence neighbor before the trial move,
  // because if it is, it will affect the seqEng
  if (seqWL || seqWLTM){
    isSeqNeighborBefore = false;
    // create a double vector of the particle position to pass to the cell list
    r.resize(ndims);
    for (int k=0; k<ndims; k++){
      r[k] = gsl_matrix_get(coords, atom_idx, k);
    }
    // get the neighbors of the particle (includes particle itself)
    cell_list.getNeighbors(r, neighbor_list, nneighbors);
    // check if any neighbors are part of the sequence
    for (int n = 0; n < nneighbors; n++){
      int j = neighbor_list[n];
      if (atom_idx != j && j < seqLength){
	distance = pair_distance(atom_idx,j);
	if (distance <= cutoff){
	  isSeqNeighborBefore = true;
	  break;
	}
      }
    }
  }
  
  // now get the particle energy
  eng_old = get_particle_energy(atom_idx);
  overlap = false;

  for (int k=0; k<ndims; k++){
    old_coords[k] = gsl_matrix_get(coords,atom_idx,k);
    new_coords[k] = gsl_matrix_get(coords,atom_idx,k); // init. new_coords to old
  }
  
  // allow large displacements in z-direction for LVcoex simulations
  // to facilitate equilibration between liquid and vapor phases
  if (!cubic_system){
    if (ncycles%10 == 0){  // try large displacements in z-dir
      delta = side_lengths[2] * (0.5 - gsl_rng_uniform(rng));
      new_pos = old_coords[2] + delta;
      new_pos = modulo(new_pos, side_lengths[2]);
      gsl_matrix_set(coords,atom_idx,2,new_pos);
      new_coords[2] = new_pos;
    } else if (ncycles%3 == 0) {  // try small displacements in z-dir
      npart_zattempted++;
      npart_zatt_cycle++;
      delta = delta_zmax * (0.5 - gsl_rng_uniform(rng));
      new_pos = old_coords[2] + delta;
      new_pos = modulo(new_pos, side_lengths[2]);
      gsl_matrix_set(coords,atom_idx,2,new_pos);
      new_coords[2] = new_pos;
    } else {
      npart_attempted++;
      npart_att_cycle++;
      for (int k=0; k<2; k++){   // small displacements in x-y
	delta = delta_max*(0.5-gsl_rng_uniform(rng));  // random displacement
	new_pos = old_coords[k] + delta;
	new_pos = modulo(new_pos, side_lengths[k]);
	gsl_matrix_set(coords,atom_idx,k,new_pos);
	new_coords[k] = new_pos;
      }
    }
  } else {  // cubic system, displace symmetrically
    npart_attempted++;
    npart_att_cycle++;

    for (int k=0; k<ndims; k++){
      delta = delta_max*(0.5-gsl_rng_uniform(rng));  // random displacement
      new_pos = old_coords[k] + delta;
      new_pos = modulo(new_pos, side_lengths[k]);
      gsl_matrix_set(coords,atom_idx,k,new_pos);
      new_coords[k] = new_pos;
    }
  }

  // determine if the particle is a sequence neighbor after the trial move,
  // because if it is, it will affect the seqEng
  if (seqWL || seqWLTM){
    isSeqNeighborAfter = false;
    // create a double vector of the particle position to pass to the cell list
    r.resize(ndims);
    for (int k=0; k<ndims; k++){
      r[k] = gsl_matrix_get(coords, atom_idx, k);
    }
    // get the neighbors of the particle (includes particle itself)
    cell_list.getNeighbors(r, neighbor_list, nneighbors);
    // check if any neighbors are part of the sequence
    for (int n = 0; n < nneighbors; n++){
      int j = neighbor_list[n];
      if (atom_idx != j && j < seqLength){
	distance = pair_distance(atom_idx,j);
	if (distance <= cutoff){
	  isSeqNeighborAfter = true;
	  break;
	}
      }
    }
  }

  isSeqNeighbor = (isSeqNeighborBefore || isSeqNeighborAfter);

  eng_new = get_particle_energy(atom_idx);
  if (large_hs) get_pair_energy(atom_idx, large_hs_idx);

  if (overlap){
    move_accepted = false;
  } else {
    if ((seqWL || seqWLTM) && isSeqNeighbor){
      seqEngOld = seqEng;
      if (use_cell_list){   // --- update cell list before getting seqEng --
    	old_cell = cell_list.cellIndexFromPosition(old_coords);
    	new_cell = cell_list.cellIndexFromPosition(new_coords);
    	if (old_cell != new_cell)
    	  cell_list.build(side_lengths, cutoff, coords, natoms);
      }
      seqEngNew = get_sequence_energy();
      if (seqEngNew < engHistMin || seqEngNew >= engHistMax){
	// -- reject if not in range --
	for (int k=0; k<ndims; k++){ // --- put particle back ---
	  gsl_matrix_set(coords,atom_idx,k,old_coords[k]);
	}
	if (use_cell_list){// -- rebuild cell list --
	  old_cell = cell_list.cellIndexFromPosition(old_coords);
	  new_cell = cell_list.cellIndexFromPosition(new_coords);
	  if (old_cell != new_cell) {
	    cell_list.build(side_lengths, cutoff, coords, natoms);
	  }
	}
      } else { // ----- accept based on DOS -----------------
	gsl_histogram_find(seqEntropyHist, seqEngOld, &binNum0);
	gsl_histogram_find(seqEntropyHist, seqEngNew, &binNum1);
	S0 = gsl_histogram_get(seqEntropyHist, binNum0);
	S1 = gsl_histogram_get(seqEntropyHist, binNum1);
	if (fabs(S1) < tol){  // don't divide by (nearly) zero
	  move_accepted = true;
	} else {
	  move_accepted = accept_dos(rng, S0, S1);
	}
      }
    } else {
      move_accepted = accept_metropolis(rng, beta, eng_old, eng_new);
    }
  }

  if (move_accepted){
    // -- stuff for LVcoex ---
    if (!cubic_system) {
      if (ncycles%10==0) {
      } else if (ncycles%3==0){
	npart_zaccepted++;
	npart_zacc_cycle++;
      } else {
	npart_accepted++;
	npart_acc_cycle++;
      }
    } else {
      npart_accepted++;
      npart_acc_cycle++;
    }
    // ---------------------

    // ---------------------
    if ((seqWL || seqWLTM) && isSeqNeighbor){  // --- seq. sampling ---
      // already updated cell_list
      seqEng = seqEngNew;  
    } else if (use_cell_list){ // ---- normal stuff ------
      old_cell = cell_list.cellIndexFromPosition(old_coords);
      new_cell = cell_list.cellIndexFromPosition(new_coords);
      if (old_cell != new_cell)
	cell_list.build(side_lengths, cutoff, coords, natoms);
    }
    // ---------------------
  } else {  // ---- move was rejected ---
    for (int k=0; k<ndims; k++){
      gsl_matrix_set(coords,atom_idx,k,old_coords[k]);
    }
  }

  // // update all bookkeeping stuff for WL/WLTM
  // if ((seqWL || seqWLTM) && (isSeqNeighbor)) {
  //   gsl_histogram_accumulate(seqEntropyHist, seqEng, gModFactor);
  //   gsl_histogram_increment(seqEngHist, seqEng);
  //   gsl_histogram_increment(seqCompEngHists[nHS].hist, seqEng);
  // }  
}


void ParticleSystem::trial_move_vol(){
  double sl, lnV0, lnV;
  double V, V0, rho0, eng0, vir0, arg;

  nvol_attempted++;
  nvol_att_cycle++;
  overlap = false;

  // old params
  V0 = volume;
  save_boxsize();
  rho0 = density;
  eng0 = pot_eng;
  vir0 = virial;

  // random volume change
  lnV0 = log(V0);
  lnV = lnV0 + delta_lnV_max*(gsl_rng_uniform(rng)-0.5);

  // assign new parameters
  volume = exp(lnV);
  V = volume;
  density = natoms/volume;

  // cubic system
  sl = pow(volume,1./ndims);
  for (int k=0; k<ndims; k++){
    new_side_lengths[k] = sl;
  }
  
  scale_coords(new_side_lengths);
  update_boxsize(new_side_lengths);
  update_cutoff_energies();  // update the cutoff, tail corrs, shift eng, etc.

  if (use_cell_list){
    cell_list.build(side_lengths, cutoff, coords, natoms);
  }
  pot_eng = get_system_energy();    // get new eng
  if (large_hs) get_full_particle_energy(large_hs_idx);

  // apply acc criterion
  arg = -beta*(target_pressure*(V-V0)-(natoms+1)*(lnV-lnV0)\
		      /beta+pot_eng-eng0);
  if (gsl_rng_uniform(rng) < exp(arg) && !overlap) {
    // ACCEPTED
    nvol_accepted++;
    nvol_acc_cycle++;
  } else {
    // REJECTED
    // reassign old params
    density = rho0;
    volume = V0;
    pot_eng = eng0;
    virial = vir0;
    scale_coords(old_side_lengths);
    update_boxsize(old_side_lengths);
    update_cutoff_energies();
    if (use_cell_list){
      cell_list.build(side_lengths, cutoff, coords, natoms);
    }
  }
}

// use this for a rectangular cell
void ParticleSystem::trial_move_vol_rect(){
  double L0, L1, delta;
  double V, V0, lnV, lnV0, rho0, eng0, vir0, arg, sl_ratio;

  nvol_attempted++;
  nvol_att_cycle++;
  overlap = false;

  // old params
  L0 = side_lengths[0];
  V0 = volume;
  lnV0 = log(V0);
  save_boxsize();
  rho0 = density;
  eng0 = pot_eng;
  vir0 = virial;

  // random side length change
  delta = delta_L_max*(0.5-gsl_rng_uniform(rng));  // random displacement
  L1 = L0 + delta;
  sl_ratio = L1/L0;
  
  for (int k=0; k<ndims; k++){
    new_side_lengths[k] = old_side_lengths[k]*sl_ratio;
  }

  // update coords, side_lengths, and cutoff
  scale_coords(new_side_lengths);
  update_boxsize(new_side_lengths);
  update_cutoff_energies();  // update the cutoff, tail corrs, shift eng, etc.

  // assign new parameters
  volume = get_system_volume();
  V = volume;
  lnV = log(V);
  density = natoms/volume;

  if (use_cell_list){
    cell_list.build(side_lengths, cutoff, coords, natoms);
  }
  pot_eng = get_system_energy();    // get new eng

  if (!hw_adjust){
    pot_eng += wall_overlap_test();
  }

  // apply acc criterion
  arg = -beta*(target_pressure*(V-V0)-(natoms+1)*(lnV-lnV0)\
	       /beta+pot_eng-eng0);
  if (gsl_rng_uniform(rng) < exp(arg) && !overlap) {
    nvol_accepted++;
    nvol_acc_cycle++;
  } else {
    // REJECTED
    // reassign old params
    density = rho0;
    volume = V0;
    pot_eng = eng0;
    virial = vir0;
    scale_coords(old_side_lengths);
    update_boxsize(old_side_lengths);
    update_cutoff_energies();
    if (use_cell_list){
      cell_list.build(side_lengths, cutoff, coords, natoms);
    }
  }
}


void ParticleSystem::sampleNeighbors(int i){
  std::vector<double> r;
  int j;

  // create a double vector of the particle position to pass to the cell list
  r.resize(ndims);
  for (int k=0; k<ndims; k++){
    r[k] = gsl_matrix_get(coords, i, k);
  }

  // get the neighbors of the particle (includes particle itself)
  cell_list.getNeighbors(r, neighbor_list, nneighbors);
  // sanity check
  assert(nneighbors <= cell_list.maxNeighbors);

  // now sample the neighbors
  for (int n = 0; n < nneighbors; n++){
    j = neighbor_list[n];
    if (i != j ){
      trial_move_part(j);
    }
  }
}


void ParticleSystem::putParticleBack(int atom_idx){
  int old_cell;
  int new_cell;
  for (int k=0; k<ndims; k++){ // --- put particle back ---
    gsl_matrix_set(coords,atom_idx,k,old_coords[k]);
  }
  if (use_cell_list){// -- rebuild cell list --
    old_cell = cell_list.cellIndexFromPosition(old_coords);
    new_cell = cell_list.cellIndexFromPosition(new_coords);
    if (old_cell != new_cell) {
      cell_list.build(side_lengths, cutoff, coords, natoms);
    }
  }
}


void ParticleSystem::trial_move_sequence_swap(){
  double seqEngOld=0., seqEngNew=0.;
  int atom_idx0, atom_idx1;
  int maxits = 1000;
  int its=0, ptypeTmp;
  size_t binNum0, binNum1;
  bool move_accepted=false;
  double tol=1.e-10, S0, S1;

  // select two different particle types in the sequence
  // FIXME:  dirty hack.  right now I just iterate guessing
  // particles until I get two different partticles of different
  // types
  while (true){
    if (its > maxits)
      break;
    atom_idx0 = gsl_rng_uniform_int(rng, seqLength);
    atom_idx1 = gsl_rng_uniform_int(rng, seqLength);
    if (atom_idx0 != atom_idx1)
      if (particles[atom_idx0].ptype != particles[atom_idx1].ptype)
	break;
    its++;
  }

  seqEngOld = seqEng;

  ptypeTmp = particles[atom_idx0].ptype;
  particles[atom_idx0].assign_params(particles[atom_idx1].ptype, hs_sigma);
  particles[atom_idx1].assign_params(ptypeTmp, hs_sigma);

  overlap = false;
  seqEngNew = get_sequence_energy();

  if (overlap){ // -- reject anything that leads to an overlap ---
    move_accepted = false;
  } else if (seqWL || seqWLTM){ // -- binning for seqWL and seqWLTM --
    if (seqEngNew < engHistMin || seqEngNew >= engHistMax) {
      move_accepted = false; // -- reject if not in eng range --
    } else { // -- new eng is in range --
      // -- accept based on DOS estimates ---
      gsl_histogram_find(seqEntropyHist, seqEngOld, &binNum0);  
      gsl_histogram_find(seqEntropyHist, seqEngNew, &binNum1);
      S0 = gsl_histogram_get(seqEntropyHist, binNum0);
      S1 = gsl_histogram_get(seqEntropyHist, binNum1);
      if (fabs(S1) < tol){  // don't divide by (nearly) zero
	move_accepted = true;
      } else {
	move_accepted = accept_dos(rng, S0, S1);
      }
    }
  } else {
    move_accepted = accept_metropolis(rng, seqBeta, seqEngOld, seqEngNew);
  }

  if (move_accepted){
    // keep the swap move
    // cout << "rank " << rank << " swap move accepted between " << atom_idx0 << " and " << atom_idx1 << endl;
    // cout << "new ptypes:  " << particles[atom_idx0].ptype << " " << particles[atom_idx1].ptype << endl;
    // cout.flush();
    seqEng = seqEngNew;
    // no need to rebuild cell lists since swap move
    // does not change a particle's cell
  } else {
    // unswap
    ptypeTmp = particles[atom_idx0].ptype;
    particles[atom_idx0].assign_params(particles[atom_idx1].ptype, hs_sigma);
    particles[atom_idx1].assign_params(ptypeTmp, hs_sigma);
  }

  // // update all bookkeeping stuff for WL/WLTM
  // if (seqWL || seqWLTM){
  //   gsl_histogram_accumulate(seqEntropyHist, seqEng, gModFactor);
  //   gsl_histogram_increment(seqEngHist, seqEng);
  //   gsl_histogram_increment(seqCompEngHists[nHS].hist, seqEng);
  // }
}

void ParticleSystem::trial_move_sequence_mutate(){
  double seqEngOld=0., seqEngNew=0.;
  int atom_idx;
  size_t binNum0, binNum1;
  bool move_accepted=false;
  double tol=1.e-10, S0, S1;

  seqEngOld = seqEng;

  atom_idx = gsl_rng_uniform_int(rng, seqLength);
  
  // for now assume only JG and HS ptypes, so just switch
  // from one to the other for mutations.  For larger alpha-
  // bets, will have to use random(0,len(alphabet)) to select
  // a type to which to mutate
  if (particles[atom_idx].ptype == 1){
    particles[atom_idx].assign_params(2, hs_sigma);
  } else {
    particles[atom_idx].assign_params(1, hs_sigma);
  }
  
  overlap = false;
  seqEngNew = get_sequence_energy();

  if (overlap){ // -- reject anything that leads to overlap ---
    move_accepted = false;
  } else if (seqWL || seqWLTM){ // -- binning for seqWL and seqWLTM --
    if (seqEngNew < engHistMin || seqEngNew >= engHistMax) {
      move_accepted = false; // -- reject move if not in eng range --
    } else { // -- new eng is in range --------------------
      // -- accept based on DOS estimate ----
      gsl_histogram_find(seqEntropyHist, seqEngOld, &binNum0);
      gsl_histogram_find(seqEntropyHist, seqEngNew, &binNum1);
      S0 = gsl_histogram_get(seqEntropyHist, binNum0);
      S1 = gsl_histogram_get(seqEntropyHist, binNum1);
      if (fabs(S1) < tol){  // don't divide by (nearly) zero
	move_accepted = true;
      } else {
	move_accepted = accept_dos(rng, S0, S1);
      }
    }
  } else {
    move_accepted = accept_metropolis(rng, seqBeta, seqEngOld, seqEngNew);
  }

  if (move_accepted){
    // keep the mutation move
    // no need to rebuild cell lists since mutation move
    // does not change a particle's cell cell
    // cout << "rank " << rank << " mutation move accepted for particle " << atom_idx << endl;
    // cout << "new ptype:  " << particles[atom_idx].ptype << endl;
    // cout.flush();
    seqEng = seqEngNew;
    nHS = 0;
    for (int iseq = 0; iseq < seqLength; iseq++){
      if (particles[iseq].ptype == 1) nHS++;
    }
  } else {
    // unmutate
    // for now assume only JG and HS ptypes, so just switch
    // from one to the other for mutations.  For larger alpha-
    // bets, will have to use random(0,len(alphabet)) to select
    // a type to which to mutate
    if (particles[atom_idx].ptype == 1){
      particles[atom_idx].assign_params(2, hs_sigma);
    } else {
      particles[atom_idx].assign_params(1, hs_sigma);
    }
  }
  
  // // update all bookkeeping stuff for WL/WLTM
  // if (seqWL || seqWLTM){
  //   gsl_histogram_accumulate(seqEntropyHist, seqEng, gModFactor);
  //   gsl_histogram_increment(seqEngHist, seqEng);
  //   gsl_histogram_increment(seqCompEngHists[nHS].hist, seqEng);
  // }
}


// FIXME:  consolidate adjust_* into one function, pass pointers,
// and put in utils.cpp
void ParticleSystem::adjust_delta_max(){
  double ratio;
  ratio = npart_acc_cycle/((double)npart_att_cycle);
  //if (verbose && ncycles%stdout_freq==0) {
  //  fprintf(logfile,"adjusting delta_max ... \n");
  //  fprintf(logfile,"npart_acc_cycle: %d\n",npart_acc_cycle);
  //  fprintf(logfile,"npart_adjust: %d\n",npart_adjust);
  //  fprintf(logfile,"ratio:  %f\n",ratio);
  //}
  if (ratio > 0.5) delta_max *= 1.05;
  else delta_max *= 0.95;
  if (delta_max > 10.0*delta_max0) delta_max = 10.0*delta_max0;
  if (delta_max < 0.1*delta_max0) delta_max = 0.1*delta_max0;
  npart_acc_cycle = 0;
  npart_att_cycle = 0;
  //if (verbose && ncycles%stdout_freq==0) 
  //  fprintf(logfile,"delta_max:  %f\n",delta_max);
}


void ParticleSystem::adjust_delta_zmax(){
  double ratio;
  ratio = npart_zacc_cycle/((double)npart_zatt_cycle);
  //if (verbose && ncycles%stdout_freq==0) {
  //  fprintf(logfile,"adjusting delta_max ... \n");
  //  fprintf(logfile,"npart_acc_cycle: %d\n",npart_acc_cycle);
  //  fprintf(logfile,"npart_adjust: %d\n",npart_adjust);
  //  fprintf(logfile,"ratio:  %f\n",ratio);
  //}
  if (ratio > 0.5) delta_zmax *= 1.05;
  else delta_zmax *= 0.95;
  if (delta_zmax > 10.0*delta_zmax0) delta_zmax = 10.0*delta_zmax0;
  if (delta_zmax < 0.1*delta_zmax0) delta_zmax = 0.1*delta_zmax0;
  npart_zacc_cycle = 0;
  npart_zatt_cycle = 0;
  //if (verbose && ncycles%stdout_freq==0) 
  //  fprintf(logfile,"delta_max:  %f\n",delta_max);
}


void ParticleSystem::adjust_delta_lnV_max(){
  double ratio;
  ratio = nvol_acc_cycle/(double)nvol_att_cycle;    
  //if (verbose && ncycles%stdout_freq==0) {
  //  fprintf(logfile,"adjusting delta_lnV_max ... \n");
  //  fprintf(logfile,"nvol_acc_cycle:  %d\n",nvol_acc_cycle);
  //  fprintf(logfile,"nvol_att_cycle:  %d\n",nvol_att_cycle);
  //  fprintf(logfile,"ratio:  %f\n",ratio);
  //}
  if (ratio > 0.5) delta_lnV_max *= 1.05;
  else delta_lnV_max *= 0.95;
  if (delta_lnV_max > 10.0*delta_lnV_max0) delta_lnV_max = 10.0*delta_lnV_max0;
  if (delta_lnV_max < 0.1*delta_lnV_max0) delta_lnV_max = 0.1*delta_lnV_max0;
  nvol_acc_cycle = 0;
  nvol_att_cycle = 0;
  //if (verbose && ncycles%stdout_freq==0)   
  //  fprintf(logfile,"delta_lnV_max:  %f\n",delta_lnV_max);
}


void ParticleSystem::adjust_delta_L_max(){
  double ratio;
  ratio = nvol_acc_cycle/(double)nvol_att_cycle;    
  //if (verbose && ncycles%stdout_freq==0) {
  //  fprintf(logfile,"adjusting delta_L_max ... \n");
  //  fprintf(logfile,"nvol_acc_cycle:  %d\n",nvol_acc_cycle);
  //  fprintf(logfile,"nvol_att_cycle:  %d\n",nvol_att_cycle);
  //  fprintf(logfile,"ratio:  %f\n",ratio);
  //}
  if (ratio > 0.5) delta_L_max *= 1.05;
  else delta_L_max *= 0.95;
  if (delta_L_max > 10.0*delta_L_max0) delta_L_max = 10.0*delta_L_max0;
  if (delta_L_max < 0.1*delta_L_max0) delta_L_max = 0.1*delta_L_max0;
  nvol_acc_cycle = 0;
  nvol_att_cycle = 0;
  //if (verbose && ncycles%stdout_freq==0)
  //  fprintf(logfile,"delta_L_max:  %f\n",delta_L_max);
}


void ParticleSystem::adjust_growth_sigma(){
  double ratio, eps=0.05*p_growth_target;
  bool changed_dr;

  ratio = ng_acc_cycle/(double)ng_att_cycle;
  if (ratio>p_growth_target+eps){
    growth_sigma *= 1 + (ratio - p_growth_target);
    changed_dr = true;
  }
  else if (ratio<p_growth_target-eps){
    growth_sigma *= 1 + (ratio - p_growth_target);
    changed_dr = true;
  }
  if (growth_sigma > 10.*growth_sigma0) growth_sigma = 10.*growth_sigma0;
  if (growth_sigma < 0.1*growth_sigma0) growth_sigma = 0.1*growth_sigma0;
  if (changed_dr){
    ng_acc_cycle = 0;
    ng_att_cycle = 0;
  }
}


// dirty hack
void ParticleSystem::adjust_config(bool openlog){
  int cnt, maxits = 100*natoms;
  double pair_eng;
  // double dr0;
  std::vector<double> r;
  
  // dr0 = delta_max;
  // delta_max = 0.01;
  r.resize(ndims);

  if (openlog)
    logfile = fopen(logfile_name.c_str(),"w");

  fprintf(logfile,"adjusting config to remove overlap ...\n");
  fflush(logfile);

  for (int i=0; i<natoms; i++){
    for (int j=i+1; j<natoms; j++){
      overlap = false;
      pair_eng = get_pair_energy(i,j);
      cnt = 0;
      while (overlap) {
	if (cnt==0){
	  fprintf(logfile,"overlap detected on pair %d, %d\n",i,j);
	}
	if (cnt>maxits) {
	  fprintf(logfile,"error:  could not remove overlap from configuration\n");
	  exit(OVERLAP);
	}
	// overlap = false;
	trial_move_part(i);
	trial_move_part(j);
	// cout << "j x coords:  " << gsl_matrix_get(coords,j,0) << endl;
	// cout.flush();
	// pair_eng = get_pair_energy(i,j);
	// exit(1);
	// if (ensemble=="npt")
	//   trial_move_vol();
	cnt++;
	fflush(logfile);
      }
      if (cnt>0){
	fprintf(logfile,"removed overlap on particle %d after %d iterations.\n",
		i,cnt);
	if (use_cell_list){
	  cell_list.build(side_lengths, cutoff, coords, natoms);
	}
      }
    }
  }
  fflush(logfile);
  if (openlog)
    fclose(logfile);

  // delta_max = dr0;
}
