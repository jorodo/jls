/*
 * File: particle_system_init.cpp
 * Description: Implementation of particle system base class
 * for MC/MD sims.  Class initialization functions.
 * ---------------------------------------------------------
 *
 * John R. Dowdle
 * The University of Texas at Austin
 * June 2008
 * GPLv3.0+
 *
 */

#include "particle_system.h"

/*
 * Function: init_coords()
 * Usage: ParticleSystem::init_coords()
 * -----------------------
 * Initialize the positions of the particles
 */
void ParticleSystem::init_coords(){

  if (ic_polymer) {
    init_polymer();
  } else if (ic_fcc == true){   // check to see if fcc lattice requested
    if (ndims==3)
      init_fcc();
    else {
      fprintf(logfile,"warning:  currently can only initialize fcc in d=3.\n");
      fprintf(logfile,"\t defaulting to square lattice.\n");
      init_square();
    }
  } else {
    init_square();
  }
}


/*
 * Function: init_square()
 * Usage: ParticleSystem::init_square()
 * -----------------------
 * Initialize the positions of the particles to a square
 * lattice.
 */
void ParticleSystem::init_square(){
  int nside=0;
  vector<int> n;
  double delta, zpos;

  delta = 0.001;  // [nm]  FIXME: this is a hack to prevent particle placement
                  // on boundary
  // square lattice
  for (int i=0; i<ndims; i++){
    n.push_back(0);
  }

  while (gsl_pow_int(nside,ndims)<natoms) nside++;

  // note that I am using only Lx here
  for (int i=0; i<natoms; i++) {
    for (int j=0; j<ndims; j++) {
      gsl_matrix_set(coords,i,j,((double)n[j]+0.5)*(side_lengths[0]-delta)/nside);
    }
    n[0]++;
    RecursiveIncrement(0,nside,ndims,n);
  }

  // move COM to center of cell
  if (!cubic_system){
    for (int i=0; i < natoms; i++){
      zpos = gsl_matrix_get(coords,i,2);
      gsl_matrix_set(coords, i, 2, zpos + (side_lengths[2]/2. - side_lengths[0]/2.));
    }
  }
}


void ParticleSystem::init_fcc(){
  int ncells, i, j, k, m, n, nplaced;
  double cell_length;
  double delta, zpos;
  double rFCC[4][3] = {{0.0, 0.0, 0.0}, {0.0, 0.5, 0.5},
		       {0.5, 0.0, 0.5}, {0.5, 0.5, 0.0}};
  double rCell[3];

  delta = 0.001;  // [nm]  FIXME: this is a hack to prevent particle placement
                  // on boundary

  // find number ncells needed to place all particles
  for (ncells = 1; ; ncells++){
    if (4*ncells*ncells*ncells >= natoms)
      break;
  }
  cell_length = (side_lengths[0]-delta) / ncells;	// FIXME: assumes cubic system
  nplaced = 0;			        
  for (i = 0; i < ncells; i++) {
    rCell[0] = i * cell_length;
    for (j = 0; j < ncells; j++) {
      rCell[1] = j * cell_length;
      for (k = 0; k < ncells; k++) {
        rCell[2] = k * cell_length;
        for (m = 0; m < 4; m++) 
          if (nplaced < natoms) {
            for (n = 0; n < 3; n++)
              gsl_matrix_set(coords,nplaced,n,rCell[n] + cell_length * rFCC[m][n]);
            ++nplaced;
          }
      }
    }
  }

  // move COM to center of cell
  if (!cubic_system){
    for (int i=0; i < natoms; i++){
      zpos = gsl_matrix_get(coords,i,2);
      gsl_matrix_set(coords, i, 2, zpos + (side_lengths[2]/2. - side_lengths[0]/2.));
    }
  }
}


/*
 * Function: init_random()
 * Usage: ParticleSystem::init_random()
 * --------------------------------------
 * Initialize the coords randomly, except
 * always place the polymer at the center
 * of the box
 *
 */
void ParticleSystem::init_random(){

}


/*
 * Function: init_polymer()
 * Usage: ParticleSystem::init_polymer()
 * --------------------------------------
 * Initialize the coords linearly and equil-
 * ibrate.  Always place the polymer at the 
 * center of the box
 *
 */
void ParticleSystem::init_polymer(){
  double dx;
  vector<double> center_of_mass, com_diff;

  dx = 0.001 * hs_sigma;
  for (int i=0; i<polymer_size; i++){
    gsl_matrix_set(coords, i, 0, i * (hs_sigma + dx));
    for (int k=1; k<ndims; k++){
      gsl_matrix_set(coords, i, k, side_lengths[k]/2.);
    }
  }
    
  // // polymer has been successfully placed
  // // so move center of mass to cell center
  // center_of_mass.resize(ndims);  // COM of polymer
  // com_diff.resize(ndims);
  // 
  // // loop through and calc com polymer
  // for (int k = 0; k < ndims; k++){
  //   for (int i = 0; i < polymer_size; i++){
  //     center_of_mass[k] += gsl_matrix_get(coords, i, k);
  //   }
  //   center_of_mass[k] /= polymer_size;
  //   // difference in com with cell center
  //   com_diff[k] = center_of_mass[k] - side_lengths[k] / 2.;  
  // }
  //   
  // // move to cell center
  // for (int i = 0; i < polymer_size; i++){
  //   for (int k = 0; k < ndims; k++){
  //     old_pos = gsl_matrix_get(coords, i, k);
  //     gsl_matrix_set(coords, i, k, old_pos - com_diff[k]);
  //   }
  // }
}

// // now randomly place solvent particles
// for (int i=polymer_size; i<natoms; i++){
//   while (true){
//     for (int k = 0; k < ndims; k++){
//       gsl_matrix_set(coords, i, k, side_lengths[k] * gsl_rng_uniform(rng));
//     }
//     overlap = false;
//     get_particle_energy(i);
//     if (!overlap) break;
//   }
//  }



/*
 * Function: init_mc_vars()
 * Usage: ParticleSystem::init_mc_vars()
 * --------------------------------------
 * Initialize the MC specific vars
 *
 */
void ParticleSystem::init_mc_vars(){
  old_coords.resize(ndims);
  new_coords.resize(ndims);
  old_side_lengths.resize(ndims);
  if (new_side_lengths.size() == 0)
    new_side_lengths.resize(ndims);
  if (npart_adjust==0){
    npart_adjust = 2;
  }
  if (ensemble=="npt"){
    if (nvol_adjust==0){
      nvol_adjust = 10*npart_adjust;
    }
  }

  if (!cubic_system && delta_zmax == 0.){
    delta_zmax = delta_max;
    delta_zmax0 = delta_max;
  }
}


/* Function:  init_vars
 * Usage:  init_vars()
 * --------------------------------------------------------------
 * Initialize all the particle system members to their default
 * values.
 */
void ParticleSystem::init_vars(){
  ensemble = "";
  calculation = "";
  natoms = 0;
  temp = 0.;
  density = 0.;
  nreplex = 0;
  nmoves = 0;
  ncycles = 0;
  ncycles_target = 0;
  npart_accepted = 0;
  npart_acc_cycle = 0;
  npart_att_cycle = 0;
  npart_adjust = 0;
  npart_attempted = 0;
  nvol_accepted = 0;
  nvol_attempted = 0;
  nvol_adjust = 0;
  nvol_att_cycle = 0;
  nvol_acc_cycle = 0;
  ndims = 3;
  rst_freq = 0;
  rst_name = "";
  rst_coords = "";
  proj_name = "";
  tpifile_name = "";
  tpi_engfile_name = "";
  bm_engfile_name = "";
  growth_engfile_name = "";
  traj_analysis_file = "";
  restart_coords = "";
  coordfile_name = "";
  xdrcoordfile_name = "";
  logfile_name = "";
  engfile_name = "";
  cor_rst_name = "";
  gr_cor_rst_name = "";
  out_dir = "";
  rst_dir = "";
  rst_path = "";
  log_dir = "";
  bennett_transform = "";
  xdrfile = true;
  corfile = false;
  bennett_method = false;
  seed_random = false;
  test_ptype = 1;     // HS
  adjust = true;
  adjust_growth = false;
  shift = false;
  tail_cor = false;
  verbose = false;
  ic_fcc = false;
  growth = false;
  cor_out_freq = 0;
  eng_out_freq = 0;
  stdout_freq = 0;
  rescaleCOM_freq = 10;
  target_pressure = 0.;
  delta_max = 0.;
  delta_max0 = 0.;
  delta_lnV_max = 0.;
  delta_lnV_max0 = 0.;
  delta_L_max = 0.;
  delta_L_max0 = 0.;
  lj_sigma = 0.3415;
  lj_epsilon = 1.03931;
  test_sigma = 0.3415;
  test_epsilon = 1.039313750;
  hs_sigma = 0.2800;
  jg_lambda0 = 0.2800;
  jg_lambda1 = 0.4816;
  jg_lambda2 = 0.8400;
  jg_epsilon1 = 3.5;
  jg_epsilon2 = -1.0;
  lambda_0_ij = 0.2800;
  lambda_1_ij = 0.4816;
  lambda_2_ij = 0.8400;
  epsilon_1_ij = 3.5;
  epsilon_2_ij = -1.0;
  cutoff = -1.;
  cutoff_sqrd = -1.;
  epsilon_ij = 0.;
  sigma_ij = 0.;
  sigma_ij_sqrd = 0.;
  growth_sigma = 0.;
  //growth_sigma_avg = 0.;
  growth_sigma_sum = 0.;
  ng_att = 0;
  ng_acc = 0;
  ng_att_cycle = 0;
  ng_acc_cycle = 0;
  ng_adjust = 100;
  p_growth = 0.;
  p_growth_target = 0.1;
  growth_equil = false;
  has_cutoff = false;
  box_length_cutoff = false;
  fixCOM = false;
  compute_pressure = true;
  cavity_radius = 0.;
  cavity_pair_dist = 0.;
  insertion_freq = 1;
  frame_num = -1;
  ninsertions = 0;
  eng_cor = 0.;
  p_cor = 0.;
  virial = 0.;
  pot_eng = 0.;
  test_eng = 0.;
  nwtrials = 0;
  nwsucc = 0;
  nwsucc_cfg = 0;
  nwtrials_cfg = 0;
  tpi_eng = 0.;
  tpi_eng_sum = 0.;
  nbins = -1; // 1000;
  ncenters = 0;
  ncenters_max = 0;
  ncenters_mean = 0.;
  ncenters_sqrd_mean = 0.;
  ncenters_sum = 0.;
  ncenters_sqrd_sum = 0.;
  ncenters_cfg_sum = 0.;
  ncenters_sqrd_cfg_sum = 0.;
  cutoff_eng = 0.;
  overlap = false;
  is_homogeneous = true;
  use_geom_sigma = false;
  p_flag = false;
  start_cycle = 0;
  end_cycle = 0;
  beta = 0.;
  cubic_system = false;
  mu_ex = 0.;
  v_sum = 0.;
  v_mean = 0.;
  delta_r_size = 50;
  xdrcoordfile = NULL;
  coords = NULL;
  delta_r = NULL;
  ndelta_r = NULL;
  ndelta_r_cyc = NULL;
  delta_rl = NULL;
  ndelta_rl = NULL;
  ndelta_rl_cyc = NULL;
  ndelta_rlm_cyc = NULL;
  ndelta_rlp_cyc = NULL;
  reduced_dr_bin_size = 0.0001;
  hard_wall = false;
  hw_adjust = true;
  hw_thickness = 0.;
  ic_polymer = false;
  polymer = false;
  max_bond_length = -1.;
  polymer_size = -1;
  gsl_hist = NULL;
  use_cell_list = false;
  nneighbors = -1;
  large_hs_idx = 0;
  large_hs = false;
  rng_seed = 1361;
  delta_zmax = 0.;
  delta_zmax0 = 0.;
  npart_zatt_cycle = 0.;
  npart_zacc_cycle = 0.;
  npart_zattempted = 0.;
  npart_zaccepted = 0.;
  // ------ additions for WL-TM sequence sampling --------------------------------
  // only include "H" and "P" monomers (JG and HS, resp.)
  seqMC = false;                 // perform sequence swap and mutation moves
  seqWL = false;                 // flag for performing WL sequence sampling
  seqWLTM = false;               // " " " WL-TM sequence sampling
  seqLength = 0;                 // number of monomers in sequence
  seqTemp = 1.;                  // the sequence temperature
  seqEng = 0.;                   // energy of the sequence
  engHistMin = 0.;               // minimum of energy range
  engHistMax = 0.;               // maximum of energy range
  seqMoveFreq = 1;               // Number of canonical MC sweeps btwn seq. moves
  seqEntropyHist = NULL;         // the dimensionless sequence entropy
  seqEngHist = NULL;             // visited energies histogram
  // seqCompEngHists = NULL;        // " " " for each possible seq. composition
  gModFactor = 1.;               // modification factor for entropy
  gModTol = 1.e-10;              // stop the simulation when g drops below this
  transitionMatrix = NULL;       // transition matrix for WL-TM
  nSeqTrials = 0;                // number of times seq moves have been tried
  nFixed = 0;                    // number of particles with fixed coords
  nHS = 0;                       // number of hard spheres in the sequence
  seqMutations = false;          // allow mutation moves
  seqEngFileName = "";           // name of the file to which to write seq. eng.
  seqEntropyHistFileName = "";   // name of density of states output file
  seqCompHistFileName = "";      // name of x_i(E) output file
  restartDOS = false;            // if set to true, read DOS from file
  seqEntRstName = "";            // the estimate of the DOS from which to restart
  // -----------------------------------------------------------------------------
}

