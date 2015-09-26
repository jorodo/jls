/*
 * File: particle_system_constructors.cpp
 * Description: Implementation of particle system base class
 * for MC/MD sims.  Class constructor/destructor functions.
 * ---------------------------------------------------------
 *
 * John R. Dowdle
 * The University of Texas at Austin
 * June 2008
 * GPLv3.0+
 *
 */

#include "particle_system.h"

/******************* Object Construction **********************************/

/*
 * Constructor: ParticleSystem
 * Usage: encoding = new ParticleSystem();
 * ---------------------------------
 * Initializes a new empty particle system.
 */
ParticleSystem::ParticleSystem(){

}


ParticleSystem::ParticleSystem(string input_file_name, int myrank){
  rank = myrank;
  setup(input_file_name);
}


/*
 * Destructor: ~Particleystem
 * Usage: delete System
 * -----------------------
 * Frees all storage associated with a particle system.
 */
ParticleSystem::~ParticleSystem(){

  gsl_rng_free(rng);

  if (coords != NULL)
    gsl_matrix_free(coords);
  if (delta_r != NULL)
    gsl_vector_free(delta_r);
  if (delta_rl != NULL)
    gsl_vector_free(delta_rl);
  if (ndelta_r != NULL)
    gsl_vector_free(ndelta_r);
  if (ndelta_rl != NULL)
    gsl_vector_free(ndelta_rl);
  if (ndelta_r_cyc != NULL)
    gsl_vector_free(ndelta_r_cyc);
  if (ndelta_rl_cyc != NULL)
    gsl_vector_free(ndelta_rl_cyc);
  if (ndelta_rlm_cyc != NULL)
    gsl_vector_free(ndelta_rlm_cyc);
  if (ndelta_rlp_cyc != NULL)
    gsl_vector_free(ndelta_rlp_cyc);
  if (gsl_hist != NULL)
    gsl_histogram_free(gsl_hist);
  if (transitionMatrix != NULL)
    gsl_matrix_free(transitionMatrix);
  if (seqEntropyHist != NULL)
    gsl_histogram_free(seqEntropyHist);
  if (seqEngHist != NULL)
    gsl_histogram_free(seqEngHist);
  if (seqCompEngHists.size() > 0){
    for (int i=0; i < seqLength+1; i++){
      gsl_histogram_free(seqCompEngHists[i].hist);
    }
  }

}

// Initialization
void ParticleSystem::setup(string input_file_name){

  InitSysParams(input_file_name);            // initial params, user input
  logfile = fopen(logfile_name.c_str(),"w");    
  SeedRNG();
  AssignParticleParams();                    // add particles to the system
  ConstructSystem();                         // build geom. and config
  InitEnergyParams();                        // parameters for potentials
  if (calculation=="cavity" || calculation=="cavity_pmf"){                // cavity analysis params
    InitCavityAnalysisParams();
  }

  if (fixCOM){
    if (center_of_mass.size() == 0){
      // find the COM
      center_of_mass.resize(ndims);
      compute_center_of_mass();
    }
    // assign COM reference if not already input by user
    if (center_of_mass_ref.size() == 0){
      for (unsigned int i=0; i<center_of_mass.size(); i++){
	center_of_mass_ref.push_back(center_of_mass[i]);
      }
    }
  }

  // get initial pressure and energy
  pot_eng = get_system_energy();
  pressure = density/beta + 1./(ndims*volume)*virial + (tail_cor?p_cor:0.);

  if (overlap) {
    adjust_config();  // adjust config for overlap due to rounding errors
    pot_eng = get_system_energy();
    pressure = density/beta + 1./(ndims*volume)*virial + (tail_cor?p_cor:0.);  
  }

  // BAR stuff
  if (bennett_transform != ""){
    bennett_method = true;
    if (bennett_transform=="scale"){
      if (new_side_lengths.size()==0){
	fprintf(logfile,"error: must specify new_side_lengths\n");
	exit(USAGE);
      }
    }
    write_bm_eng(pot_eng,pot_eng);
  }

  talk();

  fclose(logfile);

}


void ParticleSystem::InitSysParams(string input_file_name){
  init_vars();
  parse_input_file(input_file_name);
  init_mc_vars();
  name_files();
  dir_manage();

  if (rst_name.length()>0){
    // set the start cycle to the number of completed cycles 
    // in the restart file
    start_cycle = ncycles;
  }

  end_cycle = start_cycle + ncycles_target;

}


void ParticleSystem::SeedRNG(){
  if (seed_random)
    rng = InitRNG(true,rank,rng_seed);
  else
    rng = InitRNG(false,rank,rng_seed);
}


void ParticleSystem::AssignParticleParams(){
  int ncoords;
  natoms = 0;

  for (map<int,int>::iterator ii=ptypes.begin(); ii!=ptypes.end(); ++ii){
    natoms += (*ii).second;
  }

  if (calculation=="test_particle" || calculation=="cavity" || calculation=="cavity_dist")
    ncoords = natoms + 1;
  else if (calculation=="cavity_pmf")
    ncoords = natoms + 2;
  else
    ncoords = natoms;
  coords = gsl_matrix_alloc(size_t(ncoords),size_t(ndims));
  gsl_matrix_set_zero(coords);

  if (restart_coords.length() > 0) {  // restarting from previous run
    assign_restart_coords();          // assign volume, config, and geometry
    if (polymer){
      assign_bonds();
    }
  } else {                            // ^-- assumes never analyze from rst.cor file
    add_particles();
  }

  // set the 'fixed' flag on the specified particles
  // to prevent trial displacements in MC
  for (int i=0; i<int(fixed_idx.size()); i++){
    particles[fixed_idx[i]].fixed = true;
  }

  beta = 1.0/(K_BOLTZMANN*temp);
  if (nmoves==0) {
    if (ensemble=="npt")
      nmoves = natoms+1;
    else
      nmoves = natoms;
  }
}


void ParticleSystem::ConstructSystem(){
  if (restart_coords.length() > 0) {  // restarting from previous run
    density = natoms/volume;
  } else {
    assign_geometry();                // build sys geometry
    init_coords();                    // initialize configuration
  }

  if (box_length_cutoff){
    // assign the box_length_cutoff to 1/2 the minimum box length
    double min_sl = *min_element(side_lengths.begin(), side_lengths.end());
    cutoff = min_sl/2.;
  }
  
  if (use_cell_list){
    assert(cutoff>0.);
    cell_list.build_new(side_lengths, cutoff, coords, natoms);
    neighbor_list.resize(cell_list.maxNeighbors);
  }

  if (hard_wall){
    if (hw_thickness==0.)
      hw_thickness = side_lengths[ndims-1]/5.;  // vapor-liquid-hw-liquid-vapor
    hw_high_limit = side_lengths[ndims-1]/2. + hw_thickness/2.;
    hw_low_limit = side_lengths[ndims-1]/2. - hw_thickness/2.;
  }

  // if (hs_sigma > cutoff) large_hs = true;
}


void ParticleSystem::InitEnergyParams(){
  map<int,int>::iterator it_lj, it_hs, it_jg, it_sw;

  // test for homogeneity
  for (int i=1; i<(int)particles.size(); i++){
    if (particles[i].ptype != particles[i-1].ptype){
      is_homogeneous = false;
      break;
    }
  }

  // search for different particle types
  it_lj = ptypes.find(0);  // LJ particles
  it_hs = ptypes.find(1);  // HS particles
  it_jg = ptypes.find(2);  // JG particles
  it_sw = ptypes.find(3);  // SW particles

  // assign lj-specific params
  if (it_lj != ptypes.end()){
    LJ_setup();
  }

  // assign hs-specific params
  if (it_hs != ptypes.end()){
    HS_setup();
  }

  // assign jagla-specific parameters
  if (it_jg!=ptypes.end()){
    JG_setup();
  }

  // assign square-well-specific parameters
  if (it_sw!=ptypes.end()){
    SW_setup();
  }

  if (seqWL || seqWLTM){
    assert(nbins > 0);
    seqEngHist = gsl_histogram_alloc(nbins);
    gsl_histogram_set_ranges_uniform(seqEngHist, engHistMin, engHistMax);
    seqEntropyHist = gsl_histogram_alloc(nbins);
    gsl_histogram_set_ranges_uniform(seqEntropyHist, engHistMin, engHistMax);
    seqCompEngHists.resize(seqLength+1);
    for (int i=0; i < seqLength+1; i++){
      seqCompEngHists[i].hist = gsl_histogram_alloc(nbins);
      gsl_histogram_set_ranges_uniform(seqCompEngHists[i].hist, engHistMin, engHistMax);
      seqCompEngHists[i].nbins = nbins;
    }
    transitionMatrix = gsl_matrix_alloc(size_t(nbins),size_t(nbins));
    gsl_matrix_set_zero(transitionMatrix);
    if (restartDOS)
      loadRstDOS();
  }
  
  if (seqMC){
    overlap = false;
    seqEng = get_sequence_energy();
    if (overlap) {
      adjust_config();  // adjust config for overlap due to rounding errors
      seqEng = get_sequence_energy();
      // cout << "----------" << endl;
      // for (int i=0; i<seqLength; i++){
      // 	cout << "ptype " << i << ": " << particles[i].ptype << endl;
      // 	cout.flush();
      // }
      if (use_cell_list){
	cell_list.build(side_lengths, cutoff, coords, natoms);
      }
    }
    seqBeta = 1./seqTemp;
    if (seqWL || seqWLTM){
      assert(seqEng > engHistMin);
      assert(seqEng < engHistMax);
    }
  }

}


void ParticleSystem::InitCavityAnalysisParams(){
  // initialize cavity analysis parameters
  if (ncenters_max==0) ncenters_max = 100;
  for (int i=0; i<ncenters_max; i++){
    pndist.push_back(0.0);
    pndist_cfg.push_back(0.0);
  }
  is_homogeneous = false;
}


void ParticleSystem::assign_geometry(){
  // user must specify either the density or the side lengths
  if (side_lengths.size()==0){    // no side lengths given, so assume cubic
    cubic_system = true;
    double sl = pow(natoms/density,1./ndims);
    side_lengths.resize(ndims);
    for (int k=0; k<ndims; k++){
      side_lengths[k] = sl;
    }
    volume = get_system_volume();
  } else {  // side lengths are specified
    volume = get_system_volume();
    density = natoms/volume;
  }
}


void ParticleSystem::LJ_setup(){
  double core_term, decay_term;

  for (int i=0; i<natoms; i++){
    if (particles[i].ptype==0){
      epsilon_ij = particles[i].epsilon;
      sigma_ij = particles[i].sigma;
      sigma_ij_sqrd = sigma_ij*sigma_ij;
      break;
    }
  }

  if (cutoff>0.){
    has_cutoff = true;
    cutoff_sqrd = cutoff*cutoff;
    if (tail_cor) set_tail_corrections();
    if (shift){
      core_term = gsl_pow_int((sigma_ij_sqrd/cutoff_sqrd),6); 
      decay_term = gsl_pow_int((sigma_ij_sqrd/cutoff_sqrd),3);
      cutoff_eng = 4.*epsilon_ij*(core_term - decay_term);
    }
  } 
}


void ParticleSystem::HS_setup(){
  double stride;

  delta_r = gsl_vector_alloc(delta_r_size);
  ndelta_r = gsl_vector_alloc(delta_r_size);
  ndelta_r_cyc = gsl_vector_alloc(delta_r_size);

  gsl_vector_set_zero(ndelta_r);
  gsl_vector_set_zero(ndelta_r_cyc);

  stride = reduced_dr_bin_size*sigma_ij;
  // setup dr array for pressure calculation
  for (int i=0; i<int(delta_r_size); i++){
    gsl_vector_set(delta_r,i,stride + i*stride);
  }

}


void ParticleSystem::JG_setup(){
  double stride;
  // FIXME:  for now, assuming homogeneous JG system
  // and assigning homogeneous JG params.  Update if
  // for some reason mixing JG particles is desired.

  delta_r = gsl_vector_alloc(delta_r_size);
  ndelta_r = gsl_vector_alloc(delta_r_size);
  ndelta_r_cyc = gsl_vector_alloc(delta_r_size);

  gsl_vector_set_zero(ndelta_r);
  gsl_vector_set_zero(ndelta_r_cyc);

  //if (is_homogeneous){
  for (int i=0; i<natoms; i++){
    if (particles[i].ptype==2){
      sigma_ij = particles[i].lambda_0;
      lambda_0_ij = particles[i].lambda_0;
      lambda_1_ij = particles[i].lambda_1;
      lambda_2_ij = particles[i].lambda_2;
      epsilon_1_ij = particles[i].epsilon_1;
      epsilon_2_ij = particles[i].epsilon_2;
      m1 = (epsilon_2_ij-epsilon_1_ij)/(lambda_1_ij-lambda_0_ij);
      b1 = epsilon_2_ij - m1*lambda_1_ij;
      m2 = -epsilon_2_ij/(lambda_2_ij-lambda_1_ij);
      b2 = epsilon_2_ij - m2*lambda_1_ij;
      break;
    }
  }

  stride = reduced_dr_bin_size*sigma_ij;

  // setup dr array for pressure calculation
  for (int i=0; i<int(delta_r_size); i++){
    gsl_vector_set(delta_r,i,stride + i*stride);
  }

}


void ParticleSystem::SW_setup(){
  double stride;
  //double denom;

  delta_r = gsl_vector_alloc(delta_r_size);
  delta_rl = gsl_vector_alloc(delta_r_size);
  ndelta_r = gsl_vector_alloc(delta_r_size);
  ndelta_r_cyc = gsl_vector_alloc(delta_r_size);
  ndelta_rlm_cyc = gsl_vector_alloc(delta_r_size);
  ndelta_rlp_cyc = gsl_vector_alloc(delta_r_size);

  gsl_vector_set_zero(ndelta_r);
  gsl_vector_set_zero(ndelta_r_cyc);  
  gsl_vector_set_zero(ndelta_rlm_cyc);
  gsl_vector_set_zero(ndelta_rlp_cyc);

  for (int i=0; i<natoms; i++){
    if (particles[i].ptype==3){
      sigma_ij = particles[i].sigma;
      lambda_ij = particles[i].lambda;
      epsilon_ij = particles[i].epsilon;
      break;
    } 
  }

  stride = reduced_dr_bin_size*sigma_ij;

  // setup dr array for pressure calculation
  for (int i=0; i<int(delta_r_size); i++){
    gsl_vector_set(delta_r,i,stride + i*stride);
  }

  gsl_vector_memcpy(delta_rl,delta_r);
  gsl_vector_scale(delta_rl,lambda_ij/sigma_ij);

}


void ParticleSystem::parse_input_file(string fname){
  //string::size_type cmt_idx, delim_idx;
  //string::size_type ipos=0;
  string            line;
  string            keyword, params;
  const string      comment("#");
  const string      delim("=");
  vector<string>    values;

  ifstream input_file(fname.c_str());
    
  //if(!input_file) {
  //  fprintf(stderr,"error: could not open input file %s\n",fname.c_str());
  //  fprintf(logfile,"error: could not open input file %s",fname.c_str());
  //  exit(FILEIO);
  //}

  if (input_file.is_open()){
    while (!input_file.eof()){
      getline(input_file,line);
      if(line.empty());                     // Ignore empty lines
      else {
	strip_string(line);
	if (line[0]=='#');                  // Ignore commented lines
	else {
	  Tokenize(line,values);
	  keyword = values[0];
	  assign_param(keyword,values);
	  values.clear();
	}
      }
    }
  } else {
    fileio_error(fname);
  }
}

void ParticleSystem::assign_param(string keyword, vector<string> &values){
  // values[0] is the keyword, values[1] is "=" or some delimiter,
  // values[2] ... are the values to be assigned
  if (keyword=="ensemble"){
    ensemble = values[2].c_str();
  } else if (keyword=="calculation"){
    calculation = values[2].c_str();
  } else if (keyword=="lj_particles" || keyword=="LJ"){
    ptypes[0] = atoi(values[2].c_str());
  } else if (keyword=="hs_particles" || keyword=="HS"){
    ptypes[1] = atoi(values[2].c_str());
    // int nhs;
    // nhs = atoi(values[2].c_str());
    // map<int,int>::iterator it_hs;
    // it_hs = ptypes.find(1);  // HS particles
    // if (it_hs != ptypes.end()){
    //   ptypes[1] += nhs;
    // } else {
    //   //ptypes[1] = nhs;  // FIXME: polymer / hs
    // }
  } else if (keyword=="jg_particles" || keyword=="JG"){
    ptypes[2] = atoi(values[2].c_str());
  } else if (keyword=="sw_particles" || keyword=="SW"){
    ptypes[3] = atoi(values[2].c_str());
  } else if (keyword=="temp"){
    temp = atof(values[2].c_str());
  } else if (keyword=="density"){
    density = atof(values[2].c_str()); 
  } else if (keyword=="side_lengths"){
    if (ndims==0){
      fprintf(logfile,"error: ndims must be specified before side_lengths.\n");
      exit(USAGE);
    } else {
      side_lengths.resize(ndims);
      for (int k=0; k<ndims; k++){
	side_lengths[k] = atof(values[2+k].c_str()); 
      }
    }
  } else if (keyword=="new_side_lengths"){
    if (ndims==0){
      fprintf(logfile,"error: ndims must be specified before new_side_lengths.\n");
      exit(USAGE);
    } else {
      new_side_lengths.resize(ndims);
      for (int k=0; k<ndims; k++){
	new_side_lengths[k] = atof(values[2+k].c_str()); 
      }
    }
  } else if (keyword=="center_of_mass_ref"){
    if (ndims==0){
      fprintf(logfile,"error: ndims must be specified before center_of_mass_ref.\n");
      exit(USAGE);
    } else {
      center_of_mass_ref.resize(ndims);
      for (int k=0; k<ndims; k++){
	center_of_mass_ref[k] = atof(values[2+k].c_str()); 
      }
    }
  } else if (keyword=="fixed"){
    for (int k=2; k<int(values.size()); k++){
      if (values[k][0] == '#')
	break;
      nFixed++;
    }
    if (values[3]=="to"){
      int first = atoi(values[2].c_str());
      int last = atoi(values[4].c_str());
      int idx = 0;
      nFixed = last - first;
      fixed_idx.resize(last-first);
      for (int k=first; k<last; k++){
	fixed_idx[idx] = k;
	idx++;
      }
    } else {
      fixed_idx.resize(nFixed);
      for (int k=0; k<int(fixed_idx.size()); k++){
	fixed_idx[k] = atoi(values[2+k].c_str());
      }
    }
  } else if (keyword=="nreplex"){
    nreplex = atoi(values[2].c_str());
  } else if (keyword=="nmoves"){
    nmoves = atoi(values[2].c_str());
  } else if (keyword=="ncycles"){
    ncycles = atoi(values[2].c_str());
  } else if (keyword=="ncycles_target"){
    ncycles_target = atoi(values[2].c_str());
  } else if (keyword=="npart_accepted"){
    npart_accepted = atoi(values[2].c_str());
  } else if (keyword=="npart_attempted"){
    npart_attempted = atoi(values[2].c_str());
  } else if (keyword=="npart_adjust"){
    npart_adjust = atoi(values[2].c_str());
  } else if (keyword=="nvol_accepted"){
    nvol_accepted = atoi(values[2].c_str());
  } else if (keyword=="nvol_attempted"){
    nvol_attempted = atoi(values[2].c_str());
  } else if (keyword=="nvol_adjust"){
    nvol_adjust = atoi(values[2].c_str());
  } else if (keyword=="ndims"){
    ndims = atoi(values[2].c_str());
  } else if (keyword=="rst_freq"){
    rst_freq = atoi(values[2].c_str());
  } else if (keyword=="rst_coords"){
    rst_coords = values[2].c_str();
  } else if (keyword=="rst_name"){
    rst_name = values[2].c_str();
    string rst_ext = ".rst";
    string rst_filename = rst_name + rst_ext;
    rst_path = rst_dir + rst_filename;
    parse_input_file(rst_path);
  } else if (keyword=="proj_name"){
    proj_name = values[2].c_str();
  } else if (keyword=="traj_analysis_file"){
    traj_analysis_file = values[2].c_str();
  } else if (keyword=="restart_coords"){
    restart_coords = values[2].c_str();
  } else if (keyword=="output_dir"){
    out_dir = values[2].c_str();
  } else if (keyword=="rst_dir"){
    rst_dir = values[2].c_str();
  } else if (keyword=="log_dir"){
    log_dir = values[2].c_str();
  } else if (keyword=="seed_random"){
    if (values[2]=="true"){
      seed_random = true;
    } else {
      seed_random = false;
    }
  } else if (keyword=="use_geom_sigma"){
    if (values[2]=="true"){
      use_geom_sigma = true;
    } else {
      use_geom_sigma = false;
    }
  } else if (keyword=="test_ptype"){
    test_ptype = atoi(values[2].c_str());
  } else if (keyword=="adjust"){
    if (values[2]=="true"){
      adjust = true;
    } else {
      adjust = false;
    }
  } else if (keyword=="adjust_growth"){
    if (values[2]=="true"){
      adjust_growth = true;
    } else {
      adjust_growth = false;
    }
  } else if (keyword=="growth_equil"){
    if (values[2]=="true"){
      growth_equil = true;
    } else {
      growth_equil = false;
    }
  } else if (keyword=="shift"){
    if (values[2]=="true"){
      shift = true;
    } else {
      shift = false;
    }
  } else if (keyword=="tail_cor"){
    if (values[2]=="true"){
      tail_cor = true;
    } else {
      tail_cor = false;
    }
  } else if (keyword=="verbose"){
    if (values[2]=="true"){
      verbose = true;
    } else {
      verbose = false;
    }
  } else if (keyword=="ic_fcc"){
    if (values[2]=="true"){
      ic_fcc = true;
    } else {
      ic_fcc = false;
    }
  } else if (keyword=="growth"){
    if (values[2]=="true"){
      growth = true;
    } else {
      growth = false;
    }
  } else if (keyword=="is_homogeneous"){
    if (values[2]=="true"){
     is_homogeneous = true;
    } else {
      is_homogeneous = false;
    }
  } else if (keyword=="cor_out_freq"){
    cor_out_freq = atoi(values[2].c_str());
  } else if (keyword=="eng_out_freq"){
    eng_out_freq = atoi(values[2].c_str());
  } else if (keyword=="stdout_freq"){
    stdout_freq = atoi(values[2].c_str());
  } else if (keyword=="rescaleCOM_freq"){
    rescaleCOM_freq = atoi(values[2].c_str());
  } else if (keyword=="target_pressure" || keyword=="pressure"){
    target_pressure = atof(values[2].c_str());
  } else if (keyword=="delta_max"){
    delta_max = atof(values[2].c_str());
    delta_max0 = delta_max;
  } else if (keyword=="delta_zmax"){
    delta_zmax = atof(values[2].c_str());
    delta_zmax0 = delta_zmax;
  } else if (keyword=="delta_lnV_max"){
    delta_lnV_max = atof(values[2].c_str());
    delta_lnV_max0 = delta_lnV_max;
  } else if (keyword=="delta_L_max"){
    delta_L_max = atof(values[2].c_str());
    delta_L_max0 = delta_L_max;
  } else if (keyword=="hs_sigma"){
    hs_sigma = atof(values[2].c_str());
  } else if (keyword=="lj_sigma"){
    lj_sigma = atof(values[2].c_str());
  } else if (keyword=="lj_epsilon"){
    lj_epsilon = atof(values[2].c_str());
  } else if (keyword=="jg_lambda0"){
    jg_lambda0 = atof(values[2].c_str());
  } else if (keyword=="jg_lambda1"){
    jg_lambda1 = atof(values[2].c_str());
  } else if (keyword=="jg_lambda2"){
    jg_lambda2 = atof(values[2].c_str());
  } else if (keyword=="jg_epsilon1"){
    jg_epsilon1 = atof(values[2].c_str());
  } else if (keyword=="jg_epsilon2"){
    jg_epsilon2 = atof(values[2].c_str());
  } else if (keyword=="reduced_dr_bin_size"){
    reduced_dr_bin_size = atof(values[2].c_str());
  } else if (keyword=="test_sigma"){
    test_sigma = atof(values[2].c_str());
  } else if (keyword=="test_epsilon"){
    test_epsilon = atof(values[2].c_str());
  } else if (keyword=="growth_idx"){
    growth_idx = atoi(values[2].c_str());
  } else if (keyword=="cutoff"){
    cutoff = atof(values[2].c_str());
  } else if (keyword=="growth_sigma"){
    growth_sigma = atof(values[2].c_str());
    growth_sigma0 = atof(values[2].c_str());
  } else if (keyword=="p_growth_target"){
    p_growth_target = atof(values[2].c_str());
  } else if (keyword=="ng_adjust"){
    ng_adjust = atoi(values[2].c_str());
  } else if (keyword=="box_length_cutoff"){
    if (values[2]=="true"){
      box_length_cutoff = true;
    } else {
      box_length_cutoff = false;
    }
  } else if (keyword=="fixCOM"){
    if (values[2]=="true"){
      fixCOM = true;
    } else {
      fixCOM = false;
    }
  } else if (keyword=="compute_pressure"){
    if (values[2]=="true"){
      compute_pressure = true;
    } else {
      compute_pressure = false;
    }
  } else if (keyword=="hard_wall"){
    if (values[2]=="true"){
      hard_wall = true;
    } else {
      hard_wall = false;
    }
  } else if (keyword=="hw_adjust"){
    if (values[2]=="true"){
      hw_adjust = true;
    } else {
      hw_adjust = false;
    }
  } else if (keyword=="xdrfile"){
    if (values[2]=="true"){
      xdrfile = true;
    } else {
      xdrfile = false;
    }
  } else if (keyword=="corfile"){
    if (values[2]=="true"){
      corfile = true;
    } else {
      corfile = false;
    }
  } else if (keyword=="cell_list"){
    if (values[2]=="true"){
      use_cell_list = true;
    } else {
      use_cell_list = false;
    }
  } else if (keyword=="seqMC"){
    if (values[2]=="true"){
      seqMC = true;
    } else {
      seqMC = false;
    }
  } else if (keyword=="seqWL"){
    if (values[2]=="true"){
      seqWL = true;
    } else {
      seqWL = false;
    }
  } else if (keyword=="seqWLTM"){
    if (values[2]=="true"){
      seqWLTM = true;
    } else {
      seqWLTM = false;
    }
  } else if (keyword=="seqMutations"){
    if (values[2]=="true"){
      seqMutations = true;
    } else {
      seqMutations = false;
    }
  } else if (keyword=="cell_list"){
    if (values[2]=="true"){
      use_cell_list = true;
    } else {
      use_cell_list = false;
    }
  } else if (keyword=="restartDOS"){
    if (values[2]=="true"){
      restartDOS = true;
      // seqEntRstName = out_dir + rst_name + ".seqEnt.hist";
    } else {
      restartDOS = false;
    }
  } else if (keyword=="seqEntRstName"){
    seqEntRstName = values[2].c_str();
  } else if (keyword=="hw_thickness"){
    hw_thickness = atof(values[2].c_str());
  } else if (keyword=="ncenters_max"){
    ncenters_max = atoi(values[2].c_str());
  } else if (keyword=="cavity_radius"){
    cavity_radius = atof(values[2].c_str());
  } else if (keyword=="cavity_pair_dist"){
    cavity_pair_dist = atof(values[2].c_str());
  } else if (keyword=="insertion_freq"){
    insertion_freq = atoi(values[2].c_str());
  } else if (keyword=="ninsertions"){
    ninsertions = atoi(values[2].c_str());
  } else if (keyword=="nbins"){
    nbins = atoi(values[2].c_str());
  } else if (keyword=="seqLength"){
    seqLength = atoi(values[2].c_str());
  } else if (keyword=="seqMoveFreq"){
    seqMoveFreq = atoi(values[2].c_str());
  } else if (keyword=="seqTemp"){
    seqTemp = atof(values[2].c_str());
  } else if (keyword=="engHistMin"){
    engHistMin = atof(values[2].c_str());
  } else if (keyword=="engHistMax"){
    engHistMax = atof(values[2].c_str());
  } else if (keyword=="rng_seed"){
    rng_seed = strtoul(values[2].c_str(), NULL, 10);
  } else if (keyword=="bennett_transform"){
    bennett_transform = values[2].c_str();
  } else if (keyword=="polymer"){
    polymer = true;
    ic_polymer = true;
  } else if (keyword=="polymer_size"){
    polymer_size = atoi(values[2].c_str());
    if (rst_name.length() == 0){
      map<int,int>::iterator it_hs;
      it_hs = ptypes.find(1);  // HS particles
      if (it_hs != ptypes.end()){
	ptypes[1] += polymer_size;
      } else {
	ptypes[1] = polymer_size;
      }
    }
  } else if (keyword=="max_bond_length"){
    max_bond_length = atof(values[2].c_str());
  } else if (keyword=="large_hs"){
    if (values[2]=="true"){
      large_hs = true;
    } else {
      large_hs = false;
    }
  } else if (keyword=="large_hs_idx"){
    large_hs_idx = atoi(values[2].c_str());
  } else fprintf(stderr,"Could not assign parameter %s\n",keyword.c_str());

}
