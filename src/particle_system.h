/*
 * File: particle_system.h
 * Description: The system of particles for the MC/MD sims
 * -------------------------------------------------------
 *
 * John R. Dowdle
 * The University of Texas at Austin
 * June 2008
 * GPLv3.0+
 *
 *
 * The class is organized in the following way:
 *
 * - Construction/Destruction
 * - MC calculations
 * - Free Energy calculations
 * - Object utils
 * - File I/O
 * - Energy functions
 * - Initialization
 * 
 * Each section is broken into its own .cpp file to 
 * increase modularity.
 *
 *
 */

#ifndef __particle_system_h__
#define __particle_system_h__

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_histogram.h>
#include <sys/stat.h>
#include <unistd.h>
#include <algorithm>
#include <assert.h>
#include <limits>
unsigned int sleep(unsigned int seconds);

extern "C" {
#include <xdrfile/xdrfile_xtc.h>
}

#include "utils.h"
#include "particle.h"
#include "constants.h"
#include "cell_list.h"

using namespace std;

/*
 * Class:  ParticleSystem
 * Usage:  ParticleSystem* System;
 *         System = new ParticleSystemNVT(settings)
 *         System->Run();
 * -----------------------------------------------------------------------
 * All of the parameters for the particle system are defined in settings.
 * 
 */
class ParticleSystem {

 public:
  
  /*************************** Parameters ************************************/
  // These are all of the possible parameters of the ParticleSystem
  double temp, beta, density, volume, sigma_ij, lambda_0_ij, lambda_1_ij;
  double pot_eng, pressure, epsilon_ij, lambda_2_ij;
  double lambda_ij;
  double virial, hs_sigma, cutoff_sqrd, sigma_ij_sqrd, lj_sigma, lj_epsilon;
  double jg_lambda0, jg_lambda1, jg_lambda2, jg_epsilon1, jg_epsilon2;
  double m1, b1, m2, b2, epsilon_1_ij, epsilon_2_ij, test_sigma, test_epsilon;
  int cor_out_freq, eng_out_freq, rst_freq, stdout_freq, ncycles_target;
  int natoms, ndims, nmoves, ncycles, start_cycle, end_cycle, rank, nreplex;
  int npart_accepted,npart_attempted, npart_adjust, npart_acc_cycle;
  string proj_name, ensemble, restart_coords, out_dir, rst_dir, rst_coords;
  string coordfile_name, engfile_name, tpifile_name, rst_name, logfile_name;
  string xdrcoordfile_name, corcoordfile_name;
  string tpi_engfile_name, cor_rst_name, log_dir, growth_engfile_name;
  string rst_path, calculation, bm_engfile_name, distfile_name;
  bool seed_random, p_flag, tail_cor, overlap, use_geom_sigma, ic_fcc, adjust;
  bool shift, verbose, is_homogeneous, has_cutoff, box_length_cutoff, growth;
  bool cubic_system, fixCOM, compute_pressure, hard_wall, adjust_growth;
  bool growth_equil, use_cell_list, large_hs;
  double delta_max, delta_max0, target_pressure, delta_lnV_max, delta_lnV_max0;
  double delta_L_max, delta_L_max0;
  int nvol_accepted,nvol_attempted, ng_acc, ng_att, npart_att_cycle;
  int nvol_adjust,nvol_acc_cycle, growth_idx, nvol_att_cycle, large_hs_idx;
  int ng_att_cycle, ng_acc_cycle, ng_adjust, frame_num, rescaleCOM_freq;
  double p_growth_target, mu_ex, v_sum, v_mean, reduced_dr_bin_size;
  double test_eng, tpi_eng, tpi_eng_sum, ncenters_sum; 
  double ncenters_sqrd_sum, cavity_radius, ncenters_cfg_sum, cavity_pair_dist;
  double ncenters_sqrd_cfg_sum, growth_sigma, growth_sigma0, growth_sigma_avg;
  double growth_sigma_sum, p_growth;
  string traj_analysis_file, gr_cor_rst_name, bennett_transform;
  int nwtrials, nwsucc, insertion_freq, ninsertions, test_ptype, ncenters, ncenters_max;
  int nwsucc_cfg, nwtrials_cfg, nbins, polymer_size, nneighbors;
  int npart_zatt_cycle, npart_zacc_cycle, npart_zattempted, npart_zaccepted;
  double delta_zmax, delta_zmax0;
  unsigned long int rng_seed;
  double cutoff, eng_cor, cutoff_eng, p_cor, ncenters_mean, ncenters_sqrd_mean;
  double hw_thickness, hw_high_limit, hw_low_limit, max_bond_length;
  bool hw_adjust, bennett_method, xdrfile, polymer, ic_square, ic_polymer, corfile;
  int nSeqTrials, seqMoveFreq, seqLength, nFixed, nHS;
  double gModFactor, seqTemp, seqBeta, seqEng, engHistMin, engHistMax, gModTol;
  bool seqMC, seqWL, seqWLTM, seqMutations, restartDOS;
  string seqEntropyHistFileName, seqEngHistFileName, seqEngFileName, seqCompHistFileName;
  string seqEntRstName;
  gsl_histogram *gsl_hist, *seqEngHist, *seqEntropyHist;
  FILE *tpifile;
  FILE *coordfile;
  XDRFILE *xdrcoordfile;
  FILE *corcoordfile;
  FILE *engfile;
  FILE *distfile;
  FILE *logfile;
  FILE *rst_file;
  FILE *rst_coordfile;
  FILE *seqEngFile;
  FILE *seqEngHistFile;
  FILE *seqEntropyHistFile;
  FILE *seqCompHistFile;
  gsl_rng *rng;
  gsl_vector *ndelta_r, *ndelta_rl, *ndelta_rlm_cyc, *ndelta_rlp_cyc;
  gsl_vector *delta_r, *delta_rl, *ndelta_r_cyc, *ndelta_rl_cyc;
  gsl_matrix *coords, *transitionMatrix;
  size_t delta_r_size;
  vector<double> pndist, pndist_cfg;
  vector<double> center_of_mass, center_of_mass_ref;
  vector<double> old_coords, new_coords, side_lengths, old_side_lengths, new_side_lengths;
  vector<int> fixed_idx;
  vector<int> neighbor_list;
  vector<HistogramT> seqCompEngHists;
  vector<Particle> particles;
  CellList cell_list;
  map<int, int> ptypes;
  /***************************************************************************/

  /****************** Construction/Destruction *******************************/

  // Default constructor
  ParticleSystem();
  // Constructor with an input file
  ParticleSystem(string input_file_name, int myrank=0);
  // destructor
  ~ParticleSystem();
  
  /* Function:  setup
   * Usage:  setup(fname)
   * -----------------------------------------------------------
   * Setup the ParticleSystem object with all of the parameters
   * from the input file.  This determines everything about the
   * system, including the type of calculation (MC, MD, Tpi,
   * etc.).  See the example input file in files/.
   */
  void setup(string input_file_name);

  /* Function:  InitSysParams
   * Usage:  InitSysParams()
   * -----------------------------------------------------------
   * Initialize system parameters, read user input file
   */
  void InitSysParams(string input_file_name);

  /* Function:  SeedRNG
   * Usage:  SeedRNG()
   * -----------------------------------------------------------
   * Seed the RNG
   */
  void SeedRNG();

  /* Function:  AssignParticleParams
   * Usage:  AssignParticleParams()
   * -----------------------------------------------------------
   * Add particles to the system, assign particle-specific 
   * params
   */
  void AssignParticleParams();

  /* Function:  ConstructSystem
   * Usage:  ConstructSystem()
   * -----------------------------------------------------------
   * Build the system geometry and initial configuration from
   * either restart file(s) or initial lattice.
   */
  void ConstructSystem();

  /* Function:  InitEnergyParams
   * Usage:  InitEnergyParams()
   * -----------------------------------------------------------
   * Initialize system-specific parameters for the evaluation of
   * potential energies
   */
  void InitEnergyParams();

  /* Function:  InitCavityAnalysisParams
   * Usage:  InitCavityAnalysisParams()
   * -----------------------------------------------------------
   * Initialize parameters for cavity analysis calculations
   */
  void InitCavityAnalysisParams();

  /* Function:  assign_geometry
   * Usage:  assign_geometry()
   * -----------------------------------------------------------
   * Build the system geometry from user inputs
   */
  void assign_geometry();

  /* Function:  LJ_setup
   * Usage:  LJ_setup()
   * -----------------------------------------------------------
   * Assign LJ-specific params
   */
  void LJ_setup();

  /* Function:  HS_setup
   * Usage:  HS_setup()
   * -----------------------------------------------------------
   * Assign HS-specific params
   */
  void HS_setup();

  /* Function:  JG_setup
   * Usage:  JG_setup()
   * -----------------------------------------------------------
   * Assign JG-specific params
   */
  void JG_setup();

  /* Function:  SW_setup
   * Usage:  SW_setup()
   * -----------------------------------------------------------
   * Assign SW-specific params
   */
  void SW_setup();

  /* Function:  parse_input_file
   * Usage:  parse_input_file(fname)
   * -----------------------------------------------------------
   * Parse the given input file for parameters to be assigned to
   * the system object.
   */
  void parse_input_file(string fname);

  /* Function:  assign_param
   * Usage:  assign_param(keyword,value)
   * -----------------------------------------------------------
   * Assign value to the object parameter keyword.
   * 
   */
  void assign_param(string keyword, vector<string> &values);

  /*************************** MC calculations *******************************/

  /* Function:  MC_NVT
   * Usage:  MC_NVT()
   * -----------------------------------------------------------
   * Run canonical ensemble Monte Carlo
   * 
   */
  void MC_NVT(void);


  /* Function:  trial_move_part
   * Usage:  trial_move_part()
   * -----------------------------------------------------------
   * Attempt a particle move
   * 
   */
  void trial_move_part(int part=-1);

  /* Function:  trial_move_sequence_swap
   * Usage:  trial_move_sequence_swap()
   * -----------------------------------------------------------
   * Attempt a particle swap along backbone sequence
   * 
   */
  void trial_move_sequence_swap();

  /* Function:  sampleNeighbors
   * Usage:  sampleNeighbors(i)
   * -----------------------------------------------------------
   * Make an MC sweep over all of the neighbors of particle i
   * 
   */
  void sampleNeighbors(int i);

  /* Function:  putParticleBack
   * Usage:  putParticleBack(i)
   * -----------------------------------------------------------
   * Put particle i back in its previous position (old_coords)
   * and rebuild the cell list if necessary
   * 
   */
  void putParticleBack(int atom_idx);

  /* Function:  trial_move_sequence_mutate
   * Usage:  trial_move_sequence_mutate()
   * -----------------------------------------------------------
   * Attempt a particle mutation along backbone sequence
   * 
   */
  void trial_move_sequence_mutate();

  /* Function:  adjust_delta_max
   * Usage:  adjust_delta_max()
   * -----------------------------------------------------------
   * Adjust the maximum allowed displacement for particles moves
   * 
   */
  void adjust_delta_max();

  /* Function:  adjust_delta_zmax
   * Usage:  adjust_delta_zmax()
   * -----------------------------------------------------------
   * Adjust the maximum allowed displacement for particles moves
   * in the z-direction
   */
  void adjust_delta_zmax();

  /* Function:  MC_NPT
   * Usage:  MC_NPT()
   * -----------------------------------------------------------
   * Run isothermal-isobaric ensemble Monte Carlo
   * 
   */
  void MC_NPT(void);

  /* Function:  trial_move_vol
   * Usage:  trial_move_vol()
   * -----------------------------------------------------------
   * Attempt a volume change
   * 
   */
  void trial_move_vol();

  /* Function:  trial_move_vol_rect
   * Usage:  trial_move_vol_rect()
   * -----------------------------------------------------------
   * Attempt a volume change
   * 
   */
  void trial_move_vol_rect();

  /* Function:  adjust_delta_lnV_max
   * Usage:  adjust_delta_lnV_max();
   * -----------------------------------------------------------
   * Adjust the maximum allowed change in volume.
   * 
   */
  void adjust_delta_lnV_max();

  /* Function:  adjust_delta_L_max
   * Usage:  adjust_delta_L_max();
   * -----------------------------------------------------------
   * Adjust the maximum allowed change in volume.
   * 
   */
  void adjust_delta_L_max();

  /* Function:  adjust_growth_sigma
   * Usage:  adjust_growth_sigma();
   * -----------------------------------------------------------
   * Adjust the dr for trial growths.
   * 
   */
  void adjust_growth_sigma();

  /* Function:  adjust_config
   * Usage:  adjust_config();
   * -----------------------------------------------------------
   * Increase the volume to remove overlapping hard cores.
   * 
   */
  void adjust_config(bool openlog=false);

  /***************************************************************************/

  /************************* TestParticle calculations ******************************/

  /* Function:  TestParticleInsertion
   * Usage:  TestParticleInsertion()
   * -----------------------------------------------------------
   * Run a test particle insertion analysis on the 
   * trajectory.
   * 
   */
  void TestParticleInsertion();


  /* Function:  CavityPMF
   * Usage:  CavityPMF()
   * -----------------------------------------------------------
   * Run a test particle insertion analysis on the 
   * trajectory.
   * 
   */
  //void CavityPMF();


  /* Function:  CavityAnalysis
   * Usage:  CavityAnalysis()
   * -----------------------------------------------------------
   * Run a cavity analysis on the trajectory.
   * 
   */
  void CavityAnalysis();


  /* Function:  InsertParticle
   * Usage:  InsertParticle()
   * -----------------------------------------------------------
   * Attempt to insert a particle.
   * 
   * 
   */
  void InsertParticle();


  /* Function:  CavityTest
   * Usage:  CavityTest()
   * -----------------------------------------------------------
   * Test for the existence of a cavity of size 'cavity_radius' 
   * at a random location in the fluid.  If no cavity is found, 
   * count the number of solvent centers inside the volume.
   */
  void CavityTest();


  /* Function:  CavitySizeDistributions
   * Usage:  CavitySizeDistributions()
   * -----------------------------------------------------------
   * Generate a histogram of minimum distances from a randomly
   * located point to a particle.
   */
  void CavitySizeDistributions();


  /* Function:  CavityPairTest
   * Usage:  CavityPairTest()
   * -----------------------------------------------------------
   * Test for the existence of two cavities of size 'cavity_radius' 
   * separated by a distance 'cavity_pair_dist' at a random location 
   * in the fluid.  
   */
  void CavityPairTest();


  /* Function:  FinalInsertionStats
   * Usage:  FinalInsertionStats()
   * -----------------------------------------------------------
   * Calculate and write the final statistics for the test
   * particle insertion or cavity analysis calculation
   * 
   */
  void FinalInsertionStats();


  /* Function:  FinalCavityStats
   * Usage:  FinalCavityStats()
   * -----------------------------------------------------------
   * Calculate and write the final statistics for the test
   * particle insertion or cavity analysis calculation
   * 
   */
  void FinalCavityStats();


  /* Function:  StagedParticleGrowth
   * Usage:  StagedParticleGrowth()
   * -----------------------------------------------------------
   * Attempt to grow the particle of index growth_idx by dr.
   * To be used on an existing trajectory.
   * 
   */
  //void StagedParticleGrowth();


  /* Function:  TrialGrowth
   * Usage:  TrialGrowth()
   * -----------------------------------------------------------
   * Attempt to grow the particle of index growth_idx by dr.
   * To be used during the generation of a trajectory.
   * 
   */
  void TrialGrowth();


  /* Function:  BennettMethodTrajectory
   * Usage:  BennettMethodTrajectory()
   * -----------------------------------------------------------
   * For each frame in traj_analysis_file (sys0):
   *   - compute the energy of sys0
   *   - transform sys0 --> sys1 according to bennett_transform
   *   - compute the energy of sys1
   *   - output both energies to file proj_name.beng
   * To complete the Bennet Method calculation, a second 
   * simulation must be performed at the sys1 conditions
   * and transformed to sys0 according to the inverse of 
   * bennett_transform.
   */
  void BennettMethodTrajectory();


  /* Function:  BennettMethodFrame
   * Usage:  BennettMethodFrame()
   * -----------------------------------------------------------
   * For each frame in traj_analysis_file (sys0):
   *   - compute the energy of sys0
   *   - transform sys0 --> sys1 according to bennett_transform
   *   - compute the energy of sys1
   *   - output both energies to file proj_name.beng
   * To complete the Bennet Method calculation, a second 
   * simulation must be performed at the sys1 conditions
   * and transformed to sys0 according to the inverse of 
   * bennett_transform.
   */
  void BennettMethodFrame();


  /***************************************************************************/

  /********************************* Object utils ****************************/

  /* Function:  Run
   * Usage:  Run()
   * -----------------------------------------------------------
   * Run the calculation specified for the system in the input
   * file.
   * 
   * 
   */
  void Run();


  /* Function:  get_system_volume
   * Usage:  get_system_volume()
   * -----------------------------------------------------------
   * Returns the system volume as the product of side lengths
   * 
   */
  double get_system_volume();


  /* Function:  pair_distance
   * Usage:  pair_distance(i,j)
   * -----------------------------------------------------------
   * Calculate the Euclidean distance between particles i and j.
   * 
   */
  double pair_distance(int i, int j);


  /* Function:  pair_distance_sqrd
   * Usage:  pair_distance_sqrd(i,j)
   * -----------------------------------------------------------
   * Calculate the square of the Euclidean distance between particles 
   * i and j.
   */
  double pair_distance_sqrd(int i, int j);


  /* Function:  lj_mixing_rule
   * Usage:  lj_mixing_rule(i,j)
   * -----------------------------------------------------------
   * Apply the mixing rule (geom or berthelot) to lj particles.
   * 
   */
  void lj_mixing_rule(int i, int j);


  /* Function:  scale_coords
   * Usage:  scale_coords(new_boxsize)
   * -----------------------------------------------------------
   * Scale all coordinates in the system by the new side lengths
   * given in new_boxsize
   * 
   */
  void scale_coords(vector<double> &new_boxsize);


  /* Function:  rescale_center_of_mass
   * Usage:  rescale_center_of_mass()
   * -----------------------------------------------------------
   * Move all coordinates a distance -(rcom - rcom_ref) so that
   * the center of mass of the system stays at the reference
   * location.  Useful for inhomogenous systems (e.g., liquid-
   * vapor equilibrium)
   * 
   */
  void rescale_center_of_mass();


  /* Function:  compute_center_of_mass
   * Usage:  compute_center_of_mass()
   * -----------------------------------------------------------
   * Compute the center of mass of the current frame and assign
   * it to center_of_mass
   * 
   */
  void compute_center_of_mass();


  /* Function:  applyPBC
   * Usage:  applyPBC()
   * ------------------------------------------------------------
   * Apply the minimum image convention to all particle positions
   * 
   */
  void applyPBC();


  /* Function:  update_boxsize
   * Usage:  update_boxsize(new_boxsize)
   * -----------------------------------------------------------
   * Update the side lengths to those in the vector new_boxsize
   * 
   */
  void update_boxsize(vector<double> &new_boxsize);


  /* Function:  save_boxsize
   * Usage:  save_boxsize(new_boxsize)
   * ---------------------------------------------------------------
   * Archive the side lengths to the vector old_side_lengths
   * 
   */
  void save_boxsize();


  /* Function:  update_cutoff_energies
   * Usage:  update_cutoff_energies()
   * -----------------------------------------------------------
   * If the cutoff is set to L/2, then update the cutoff, tail
   * corrections, and shift energy if required.
   * 
   */
  void update_cutoff_energies();


  /* Function:  get_temp
   * Usage:  get_temp()
   * -----------------------------------------------------------
   * Get the system temperature.
   * 
   */
  double get_temp();

  /* Function:  dir_manage
   * Usage:  dir_manage()
   * -----------------------------------------------------------
   * Backup existing trajectory files before writing new ones.
   * 
   */
  void dir_manage();

  /* Function:  get_max_sigma
   * Usage:  get_max_sigma()
   * -----------------------------------------------------------
   * Get the largest hard core diameter in the system.
   * 
   */
  double get_max_sigma();

  /* Function:  name_files
   * Usage:  name_files()
   * -----------------------------------------------------------
   * Name the output files for the simulation.
   * 
   */
  void name_files();

  /* Function:  get_natoms
   * Usage:  get_natoms()
   * -----------------------------------------------------------
   * Obtain the number of atoms in the system.
   * 
   */
  void get_natoms();

  /* Function:  add_particles
   * Usage:  add_particles()
   * -----------------------------------------------------------
   * Add all of the particles to the system's coordinates vector.
   * 
   */
  void add_particles();

  /* Function:  sidelengths2xdrbox
   * Usage:  sidelengths2xdrbox(box)
   * -----------------------------------------------------------
   * Assign the matrix box from the system side_lengths so that
   * write_xtc (from xdrfile) may have its desired matrix type
   * input.
   * 
   */
  void sidelengths2xdrbox(matrix box);

  /* Function:  xdrbox2sidelengths
   * Usage:  xdrbox2sidelengths(box)
   * -----------------------------------------------------------
   * Assign the matrix box from the system side_lengths so that
   * write_xtc (from xdrfile) may have its desired matrix type
   * input.
   * 
   */
  void xdrbox2sidelengths(matrix box);

  /* Function:  coords2xdr
   * Usage:  coords2xdr(xdrcoords);
   * -----------------------------------------------------------
   * Assign the coordinates of the system to the rvec* type
   * xdrcoords, to be passed to write_xtc from xdrfile.
   * 
   */
  void coords2xdr(rvec *xdrcoords);

  /* Function:  xdr2coords
   * Usage:  xdr2coords(xdrcoords);
   * -----------------------------------------------------------
   * Assign the coordinates of the system to the rvec* type
   * xdrcoords, to be passed to write_xtc from xdrfile.
   * 
   */
  void xdr2coords(rvec *xdrcoords);

  /***************************************************************************/

  /****************************** File I/O ***********************************/

  /* Function:  ReadXTCFrame
   * Usage:  ReadXTCFrame(traj_file)
   * --------------------------------------------------------------
   * Read in the next single frame from the .xtc file
   * traj_file should already be opened
   */
  int ReadXTCFrame(XDRFILE *xdrtraj, string xdrtraj_path="");
  
  /* Function:  ReadPDBFrame
   * Usage:  ReadPDBFrame(traj_file)
   * --------------------------------------------------------------
   * Read in the next single frame from the .pdb file traj_file.
   * 
   */
  void ReadPDBFrame(ifstream &traj_file);

  /* Function:  loadRstDOS
   * Usage:  loadRstDOS(seqEntRstFile)
   * --------------------------------------------------------------
   * Load the estimate for the DOS
   * 
   */
  void loadRstDOS();

  /* Function:  write_eng
   * Usage:  write_eng()
   * --------------------------------------------------------------
   * Write the current energy information to the system's .eng file
   * 
   */
  void write_eng();

  /* Function:  write_seq_eng
   * Usage:  write_seq_eng()
   * --------------------------------------------------------------
   * Write the current sequence energy to the system's .seqEng file
   * 
   */
  void write_seq_eng();

  /* Function:  talk
   * Usage:  talk(fsnum)
   * --------------------------------------------------------------
   * Print system information to the given file stream (default is
   * stdout).
   * 
   */
  void talk(int fsnum=0);

  /* Function:  write_coords
   * Usage:  write_coords()
   * --------------------------------------------------------------
   * Write the current coordinates to the trajectory file in PDB
   * format.
   * 
   */
  void write_coords();

  /* Function:  write_tpi
   * Usage:  write_tpi()
   * --------------------------------------------------------------
   * Write the current coordinates plus the test particle position
   * to a trajectory file.
   * 
   */
  void write_tpi();

  /* Function:  write_tpi_eng
   * Usage:  write_tpi_eng()
   * --------------------------------------------------------------
   * Write the test particle insertion information to the .weng 
   * file.
   *
   */
  void write_tpi_eng();  

  /* Function:  write_cav_dist
   * Usage:  write_cav_dist()
   * --------------------------------------------------------------
   * Write the distribution of cavity sizes to file.
   *
   */
  void write_cav_dist();

  /* Function:  write_seq_entropy_hist
   * Usage:  write_seq_entropy_hist()
   * --------------------------------------------------------------
   * Write the entropy histogram to file
   *
   */
  void write_seq_entropy_hist();

  /* Function:  write_seq_eng_hist
   * Usage:  write_seq_eng_hist()
   * --------------------------------------------------------------
   * Write the current energy histogram to file
   *
   */
  void write_seq_eng_hist();

  /* Function:  write_seq_comp_hist
   * Usage:  write_seq_comp_hist()
   * --------------------------------------------------------------
   * Write the composition energy histograms to file
   *
   */
  void write_seq_comp_hist();

  /* Function:  write_bm_eng
   * Usage:  write_bm_eng()
   * --------------------------------------------------------------
   * Write the energies of system 0 and system 1 in the Bennett's
   * Acceptance Ratio method calculation to file.
   *
   */
  void write_bm_eng(double bm_eng0, double bm_eng1);  

  /* Function:  write_growth_eng
   * Usage:  write_growth_eng()
   * --------------------------------------------------------------
   * Write the trial growth information to the .weng 
   * file.
   *
   */
  void write_growth_eng();  

  /* Function:  write_rst
   * Usage:  write_rst()
   * --------------------------------------------------------------
   * Write restart info.
   *
   */
  void write_rst();

  /* Function:  write_rst_coords
   * Usage:  write_rst_coords()
   * --------------------------------------------------------------
   * Write restart coordinates.
   *
   */
  void write_rst_coords();


  /* Function:  write_growth_rst_coords
   * Usage:  write_growth_rst_coords()
   * --------------------------------------------------------------
   * Write restart coordinates.
   *
   */
  void write_growth_rst_coords();

  /***************************************************************************/

  /********************************** Energy functions ***********************/

  /* Function:  LJ_LJ_Energy
   * Usage:  LJ_LJ_Energy()
   * --------------------------------------------------------------
   * Return Lennard-Jones/Lennard-Jones interaction potential energy 
   * between particles i and j.
   *
   */
  double LJ_LJ_Energy(int i, int j);

  /* Function:  HS_HS_Energy
   * Usage:  HS_HS_Energy()
   * --------------------------------------------------------------
   * Determine whether or not hard spheres i and j overlap.
   *
   */
  double HS_HS_Energy(int i, int j);

  /* Function:  JG_JG_Energy
   * Usage:  JG_JG_Energy()
   * --------------------------------------------------------------
   * Return Jagla/Jagla interaction energy.  
   *
   */
  double JG_JG_Energy(int i, int j);

  /* Function:  SW_SW_Energy
   * Usage:  SW_SW_Energy()
   * --------------------------------------------------------------
   * Return Square-Well interaction energy.  
   *
   */
  double SW_SW_Energy(int i, int j);

  /* Function:  LJ_LJ_tail_corrections
   * Usage:  LJ_LJ_tail_corrections(sigma_ij,epsilon_ij)
   * --------------------------------------------------------------
   * Calculate the tail corrections for a homogeneous Lennard-Jones
   * system of params sigma_ij and epsilon_ij.
   *
   */
  void LJ_LJ_tail_corrections();

  /* Function:  get_pair_energy
   * Usage:  get_pair_energy(i,j)
   * --------------------------------------------------------------
   * Calculate interaction energy between particle i and j.  
   * Select the appropriate ineraction function from below.
   *
   */
  double get_pair_energy(int i, int j);

  /* Function:  get_particle_energy
   * Usage:  get_particle_energy(i)
   * --------------------------------------------------------------
   * Calculate the total interaction energy for particle i.
   *
   */
  double get_particle_energy(int i, int jmin=0);

  /* Function:  get_particle_energy
   * Usage:  get_particle_energy(i)
   * --------------------------------------------------------------
   * Calculate the total interaction energy for particle i.
   *
   */
  double get_full_particle_energy(int i);

  /* Function:  get_particle_energy
   * Usage:  get_particle_energy(i)
   * --------------------------------------------------------------
   * Calculate the total interaction energy for particle i.
   *
   */
  double get_cell_list_particle_energy(int i, int jmin);

  /* Function:  get_system_energy
   * Usage:  get_system_energy()
   * --------------------------------------------------------------
   * Calculate the total system energy.
   *
   */
  double get_system_energy();

  /* Function:  get_sequence_energy
   * Usage:  get_sequence_energy()
   * --------------------------------------------------------------
   * Calculate the energy of the heteropolymer sequence in its
   * current configuration, including any solvent-polymer
   * interactions.
   *
   */
  double get_sequence_energy();

  /* Function:  wall_overlap_test
   * Usage:  wall_overlap_test()
   * --------------------------------------------------------------
   * Test for hard core overlaps on the wall
   *
   */
  double wall_overlap_test();

  /* Function:  get_system_pressure
   * Usage:  get_system_pressure()
   * --------------------------------------------------------------
   * Calculate the total system pressure.
   *
   */
  double get_system_pressure();

  /* Function:  set_tail_corrections
   * Usage:  set_tail_corrections()
   * --------------------------------------------------------------
   * Set the tail correction for the system energy evaluation.
   *
   */
  void set_tail_corrections(int ptype=0);

  /***************************************************************************/


  /********************************** Initialization *************************/

  /* Function:  init_vars
   * Usage:  init_vars()
   * --------------------------------------------------------------
   * Initialize all the particle system members.
   *
   */
  void init_vars();

  /* Function:  init_mc_vars
   * Usage:  init_mc_vars()
   * --------------------------------------------------------------
   * Initialize all the MC-specific vars.
   *
   */
  void init_mc_vars();

  /* Function:  init_coords
   * Usage:  init_coords()
   * --------------------------------------------------------------
   * Initialize the coordinates of the system to a square lattice
   * in d-dimensions unless fcc is specified.
   *
   */
  void init_coords();
 
 /* Function:  init_fcc
   * Usage:  init_fcc()
   * --------------------------------------------------------------
   * Initialize the coordinates of the system to a fcc lattice.
   * Only works in three dimensions.
   *
   */
  void init_fcc();

 /* Function:  init_square
   * Usage:  init_square()
   * --------------------------------------------------------------
   * Initialize the coordinates of the system to a square lattice.
   *
   */
  void init_square();

 /* Function:  init_polymer
   * Usage:  init_polymer()
   * --------------------------------------------------------------
   * Initialize the coordinates of a polymer
   *
   */
  void init_polymer();

 /* Function:  init_random
   * Usage:  init_random()
   * --------------------------------------------------------------
   * Initialize the coordinates of the system randomly
   *
   */
  void init_random();

  /* Function:  assign_restart_coords
   * Usage:  assign_restart_coords()
   * --------------------------------------------------------------
   * Assign the coordinates from the system's restart file.
   * 
   */
  void assign_restart_coords();

  /* Function:  assign_bonds
   * Usage:  assign_bonds()
   * --------------------------------------------------------------
   * Assign the polymer's bonds
   * 
   */
  void assign_bonds();

  /***************************************************************************/

};


#endif
