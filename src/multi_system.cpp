/*
 * File: replica_exchange.cpp
 * Description: Implementation of the replica exchange
 * algorithm using MPI.
 * ---------------------------------------------------------
 *
 * John R. Dowdle
 * The University of Texas at Austin
 * June 2008
 * GPLv3.0+
 *
 */


#include "multi_system.h"

MultiSystem::MultiSystem(){
  
}

MultiSystem::MultiSystem(string input_file_name, int num_systems, 
			 int num_procs, int myrank, int mytnum,
			 string mypfunc_fname, bool write_mode){

  write_mode_temps = write_mode;
  tnum = mytnum;
  pfunc_fname = mypfunc_fname;
  rank = myrank;
  nsystems = num_systems;
  nprocs = num_procs;
  System = new ParticleSystem(input_file_name,rank);

  pfunc = new int[nsystems];
  pfunc_inv = new int[nsystems];
  pot_eng = new double[nsystems];
  vol = new double[nsystems];
  pres = new double[nsystems];
  temp = new double[nsystems];
  beta = new double[nsystems];
  dr_max = new double[nsystems];
  dlnV_max = new double[nsystems];
  dL_max = new double[nsystems];
  temp_inv = new double[nsystems];
  beta_inv = new double[nsystems];
  nvol_att = new int[nsystems];
  nvol_acc = new int[nsystems];
  npart_att = new int[nsystems];
  npart_acc = new int[nsystems];
  dr_max_inv = new double[nsystems];
  dlnV_max_inv = new double[nsystems];
  dL_max_inv = new double[nsystems];
  pair_attempts = new int[nsystems-1];
  pair_exchanges = new int[nsystems-1];
  pair_probs = new double[nsystems-1];
  if (System->seqWL || System->seqWLTM){
    gModFactors = new double[nsystems];
    engHistMins = new double[nsystems];
    engHistMaxs = new double[nsystems];
    allSeqEng = new double[nsystems];
    allCoords = new double[nsystems*System->natoms*System->ndims];
    seqPtypes = new int[System->seqLength];
    allSeqPtypes = new int[nsystems*System->seqLength];
    for (int i=0; i<System->seqLength; i++){
      seqPtypes[i] = System->particles[i].ptype;
    }

    gModTol = 1.e-10;   // stop the simulation when all g's drop below this
    // for (int i=0; i<nsystems*natoms*ndims; i++){
    //   allCoords[i] = 0.;
    // }
  }

  for (int i=0; i<nsystems-1; i++){
    pair_attempts[i] = 0;
    pair_exchanges[i] = 0;
    pair_probs[i] = 0.;
  }

  // output file renaming
  if (tnum==-1)
    tnum = rank;
  pvs_tnum_str = int2str(tnum);

  int proj_name_length = (int)System->proj_name.length();
  string pfunc_str = System->proj_name;
  pfunc_str.erase(proj_name_length-pvs_tnum_str.length(),
		  pvs_tnum_str.length());
  pfunc_fname = System->out_dir + pfunc_str + "_pfunc.dat";

  MPI_Allgather(&(System->beta),1,MPI_DOUBLE,beta_inv,1,MPI_DOUBLE,
		MPI_COMM_WORLD);
  MPI_Allgather(&(System->temp),1,MPI_DOUBLE,temp_inv,1,MPI_DOUBLE,
		MPI_COMM_WORLD);
  MPI_Allgather(&(System->delta_max),1,MPI_DOUBLE,dr_max_inv,1,MPI_DOUBLE,
		MPI_COMM_WORLD);
  if (System->ensemble=="npt"){
    if (System->cubic_system){
      MPI_Allgather(&(System->delta_lnV_max),1,MPI_DOUBLE,dlnV_max_inv,1,MPI_DOUBLE,
		    MPI_COMM_WORLD);
    } else {
      MPI_Allgather(&(System->delta_L_max),1,MPI_DOUBLE,dL_max_inv,1,MPI_DOUBLE,
		    MPI_COMM_WORLD);
    }
  }

  MPI_Allgather(&tnum,1,MPI_INT,pfunc_inv,1,MPI_INT,MPI_COMM_WORLD);

  if (System->seqWL || System->seqWLTM){
    MPI_Allgather(&(System->gModFactor), 1, MPI_DOUBLE, gModFactors, 1, 
		  MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(&(System->engHistMin), 1, MPI_DOUBLE, engHistMins, 1, 
		  MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(&(System->engHistMax), 1, MPI_DOUBLE, engHistMaxs, 1, 
		  MPI_DOUBLE, MPI_COMM_WORLD);
  }

  // update pfunc after gathering pfunc_inv
  for (int i=0; i<nsystems; i++){
    pfunc[pfunc_inv[i]] = i;
    temp[pfunc_inv[i]] = temp_inv[i];
    beta[pfunc_inv[i]] = beta_inv[i];
    dr_max[pfunc_inv[i]] = dr_max_inv[i];
    if (System->ensemble=="npt"){
      if (System->cubic_system)
	dlnV_max[pfunc_inv[i]] = dlnV_max_inv[i];
      else
	dL_max[pfunc_inv[i]] = dL_max_inv[i];
    }
  }
}


MultiSystem::~MultiSystem(){
  delete System;

  delete[] pfunc;
  delete[] pot_eng;
  delete[] vol;
  delete[] pres;
  delete[] temp;
  delete[] beta;
  delete[] temp_inv;
  delete[] beta_inv;
  delete[] pair_attempts;
  delete[] pair_exchanges;
  delete[] pair_probs;
  delete[] dr_max;
  delete[] dlnV_max;
  delete[] dL_max;
  delete[] nvol_att;
  delete[] nvol_acc;
  delete[] npart_att;
  delete[] npart_acc;
  delete[] dr_max_inv;
  delete[] dlnV_max_inv;
  delete[] dL_max_inv;
  delete[] pfunc_inv;
  if (System->seqWL || System->seqWLTM){
    delete[] gModFactors;
    delete[] engHistMins;
    delete[] engHistMaxs;
    delete[] allCoords;
  }

}


void MultiSystem::RunAll(){
  System->Run();
}


void MultiSystem::ReplicaExchange(){

  // open permutation function file
  if (rank==0){
    if (System->ncycles==0){
      pfunc_file.open(pfunc_fname.c_str());
      // WritePfunc();
      WritePfuncInv();
    } else {
      pfunc_file.open(pfunc_fname.c_str(),ios::app);
    }
  }

  // update run and total cycle counts
  remc_end_cycle = System->end_cycle;
  System->ncycles_target = System->nreplex;
  nrun = System->ncycles/System->nreplex;

  while (System->ncycles < remc_end_cycle){
    System->end_cycle = System->start_cycle + System->nreplex; 

    // sample
    System->Run();
    System->ncycles -= 1;  // dirty hack
    // bookkeep for file renaming later
    pvs_tnum_str = int2str(tnum);    
    // send new data to master process
    // double mpi_start_time = MPI_Wtime();
    MPISendToMaster();
    // fprintf(stdout,"rank%d: mpi gather time:  %f\n",rank,MPI_Wtime() - mpi_start_time);
    
    // master evaluate exchanges
    if (rank==0){
      // cout << endl;
      // cout << "=============== REMC run " << nrun << " =====================" << endl;
      AttemptREMCExchanges();
    }

    // scatter new temps
    MPIMasterSend();
    // fprintf(stdout,"rank%d: mpi total time: %f \n",rank,MPI_Wtime()-mpi_start_time);
    UpdateParams();

    // write pfunc
    if (rank==0){
      // WritePfunc();
      WritePfuncInv();
    }

    //prepare for next run
    System->start_cycle = System->ncycles;
    nrun++;
    fflush(stdout);
  }

  if (rank==0){
    // cout << "===========================================" << endl;
    // cout << endl;
    WriteExchangeStats();
    pfunc_file.close();
  }
}

void MultiSystem::WangLandau(){

  // open permutation function file
  if (rank==0){
    if (System->ncycles==0){
      pfunc_file.open(pfunc_fname.c_str());
      // WritePfunc();
      WritePfuncInv();
    } else {
      pfunc_file.open(pfunc_fname.c_str(), ios::app);
    }
  }

  // update run and total cycle counts
  remc_end_cycle = System->end_cycle;
  System->ncycles_target = System->nreplex;
  nrun = System->ncycles/System->nreplex;

  while (System->ncycles < remc_end_cycle){
    System->end_cycle = System->start_cycle + System->nreplex; 
    // logfile = fopen((System->logfile_name).c_str(),"a");   
    // fprintf(logfile,"=============== WL run %d =====================\n",nrun);
    // fclose(logfile);

    // sample
    System->Run();
    System->ncycles -= 1;  // dirty hack

    // update the sequence composition
    for (int i=0; i<System->seqLength; i++){
      seqPtypes[i] = System->particles[i].ptype;
    }

    // send new data to master process
    // double mpi_start_time = MPI_Wtime();
    MPISendToMasterWL();
    // fprintf(stdout,"rank%d: mpi gather time:  %f\n",rank,MPI_Wtime() - mpi_start_time);

    // master evaluate exchanges
    if (rank==0){
      // cout << endl;
      // cout << "=============== WL run " << nrun << " =====================" << endl;
      // check g's to see if we have converged
      CheckModFactors();
      if (gModsLTtol && !System->restartDOS) break;
      AttemptWLExchanges();
    }
    // scatter new engs
    MPIMasterSendWL();
    // fprintf(stdout,"rank%d: mpi total time: %f \n",rank,MPI_Wtime()-mpi_start_time);
    UpdateParamsWL();

    // write pfunc
    if (rank==0){
      // WritePfunc();
      WritePfuncInv();
    }

    //prepare for next run
    System->start_cycle = System->ncycles;
    nrun++;
    fflush(stdout);
  }

  if (rank==0){
    cout << "===========================================" << endl;
    cout << endl;
    WriteExchangeStats();
    pfunc_file.close();
  }
}


void MultiSystem::MPISendToMaster(){
  MPI_Gather(&(System->pot_eng),1,MPI_DOUBLE,pot_eng,1,MPI_DOUBLE,0,
	     MPI_COMM_WORLD);
  MPI_Gather(&(System->volume),1,MPI_DOUBLE,vol,1,MPI_DOUBLE,0,
	     MPI_COMM_WORLD);
}


void MultiSystem::MPIMasterSend(){
  MPI_Scatter(pfunc_inv,1,MPI_INT,&tnum,1,MPI_INT,0,MPI_COMM_WORLD);
}


void MultiSystem::UpdateParams(){
  System->temp = temp[tnum];
  System->beta = beta[tnum];
  System->delta_max = dr_max[tnum];
  if (System->ensemble=="npt"){
    if (System->cubic_system)
      System->delta_lnV_max = dlnV_max[tnum];
    else
      System->delta_L_max = dL_max[tnum];
  }
}


void MultiSystem::UpdateParamsWL(){
  for (int i=0; i < System->seqLength; i++){
    System->particles[i].assign_params(seqPtypes[i], System->hs_sigma);
  }
  if (System->use_cell_list){
    System->cell_list.build(System->side_lengths, System->cutoff, 
			    System->coords, System->natoms);
  }
}


void MultiSystem::MPISendToMasterWL(){
  MPI_Gather(&(System->seqEng),1,MPI_DOUBLE,allSeqEng,1,MPI_DOUBLE,0,
	     MPI_COMM_WORLD);
  MPI_Gather(&(System->gModFactor),1,MPI_DOUBLE,gModFactors,1,MPI_DOUBLE,0,
	     MPI_COMM_WORLD);
  MPI_Gather(System->coords->data, System->natoms*System->ndims,
	     MPI_DOUBLE, allCoords, System->natoms*System->ndims,
	     MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(seqPtypes, System->seqLength, MPI_INT, allSeqPtypes,
	     System->seqLength, MPI_INT, 0, MPI_COMM_WORLD);
}


void MultiSystem::MPIMasterSendWL(){
  MPI_Scatter(pfunc_inv,1,MPI_INT,&tnum,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Scatter(allCoords, System->natoms*System->ndims, MPI_DOUBLE, System->coords->data,
  	      System->natoms*System->ndims, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatter(allSeqPtypes, System->seqLength, MPI_INT, seqPtypes,
	      System->seqLength, MPI_INT, 0, MPI_COMM_WORLD);
}


void MultiSystem::CheckModFactors(){
  // stop the simulation if the g mod factors are less than the tolerance
  // for all windows
  gModsLTtol = true;
  for (int i=0; i<nsystems; i++){
    if (gModFactors[i] > gModTol){
      gModsLTtol = false;
      break;
    }
  }
}


void MultiSystem::AttemptREMCExchanges(){
  int i,j,m,n;
  double delta, p_acc, rand_real, targ_pres;

  targ_pres = System->target_pressure;

  if (nrun%2==0 || nsystems==2){
    m = 0;
    n = 1;
  } else {
    m = 1;
    n = 2;
  }
  
  while (n<nsystems){
    i = pfunc[m];
    j = pfunc[n];

    // cout << "-----------------------------------" << endl;
    // cout << "attempting exchange:" << endl;
    // cout << "replicas:  " << i << " and " << j << endl;
    // cout << "at temperatures:  " << temp[m] << " (" << m << ") ";
    // cout << " and " << temp[n] << " (" << n << ") " << ", resp." << endl;
    // cout << "beta[n]:  " << beta[n] << endl;  // DEBUG
    // cout << "beta[m]:  " << beta[m] << endl;
    // cout << "eng[i]:  " << pot_eng[i] << endl;
    // cout << "eng[j]:  " << pot_eng[j] << endl;
    // cout << "vol[i]:  " << vol[i] << endl;
    // cout << "vol[j]:  " << vol[j] << endl;
    
    pair_attempts[m]++;
    if (System->ensemble=="npt"){
      delta = (beta[n] - beta[m])*(pot_eng[i] - pot_eng[j])
	+ (beta[n] - beta[m])*targ_pres*(vol[i] - vol[j]);
    } else {
      delta = (beta[n]-beta[m])*(pot_eng[i]-pot_eng[j]);
    }
    p_acc = min(1.,exp(-delta));
    rand_real = gsl_rng_uniform(System->rng);

    // cout << "p_acc:  " << p_acc << endl;           // DEBUG
    // cout << "rand_real:  " << rand_real << endl;

    if (rand_real<p_acc){
      // cout << "EXCHANGE ACCEPTED" << endl;
      pfunc[m] = j;
      pfunc[n] = i;
      pfunc_inv[j] = m;
      pfunc_inv[i] = n;
      pair_exchanges[m]++;
    } 
    // cout << "-----------------------------------" << endl;

    m += 2;
    n += 2;
  }
}


void MultiSystem::AttemptWLExchanges(){
  // Assume that Emin always corresponds to replica 0,
  // and Emax always corresponds to replica M (for replicas
  // enumerated [0, ... , M]).
  int i,j,m,n;
  double tmp;
  int pTmp;

  if (nrun%2==0 || nsystems==2){
    i = 0;
    j = 1;
  } else {
    i = 1;
    j = 2;
  }
  
  while (j < nsystems){
    m = pfunc_inv[i];
    n = pfunc_inv[j];

    // cout << "-----------------------------------" << endl;
    // cout << "attempting exchange:" << endl;
    // cout << "replicas:  " << i << " and " << j << endl;
    // cout << "seqEng[i]:  " << allSeqEng[i] << endl;
    // cout << "seqEng[j]:  " << allSeqEng[j] << endl;
    
    pair_attempts[i]++;

    if ((allSeqEng[i] > engHistMins[j]) && (allSeqEng[j] < engHistMaxs[i])){
      // perform swap
      // cout << "EXCHANGE ACCEPTED" << endl;
      // cout << "... swapping " << i << " and " << j << "..." << endl;
      // swap coordinates
      for (int atomn=0; atomn<System->natoms; atomn++){
        for (int dim=0; dim<System->ndims; dim++){
          tmp = allCoords[i*System->natoms*System->ndims \
                          + atomn*System->ndims + dim];
          allCoords[i*System->natoms*System->ndims \
                    + atomn*System->ndims + dim] = allCoords[j*System->natoms*System->ndims\
                                                     + atomn*System->ndims + dim];
          allCoords[j*System->natoms*System->ndims + atomn*System->ndims + dim] = tmp;
        }
      }
      // swap compositions
      for (int atomn=0; atomn<System->seqLength; atomn++){
	pTmp = allSeqPtypes[i*System->seqLength + atomn];
	allSeqPtypes[i*System->seqLength + atomn] = allSeqPtypes[j*System->seqLength + atomn];
	allSeqPtypes[j*System->seqLength + atomn] = pTmp;
      }
      // update bookkeeping
      pfunc_inv[j] = m;
      pfunc_inv[i] = n;
      pfunc[m] = j;
      pfunc[n] = i;
      pair_exchanges[i]++;
    }

    // cout << "-----------------------------------" << endl;

    i += 2;
    j += 2;
  }
}


void MultiSystem::WriteExchangeStats(){

  cout << "-----------------------------------" << endl;
  cout << "Exchange Statistics: " << endl;
  cout << "pair:  ";
  for (int m=0; m<nsystems-1; m++){
    pair_probs[m] = pair_exchanges[m] / (double) pair_attempts[m];
    cout << m << "-" << m+1 << "  ";
  }
  cout << endl;
  cout << "attempts:  ";
  for (int m=0; m<nsystems-1; m++){
    cout << pair_attempts[m] << "  ";
  }
  cout << endl;
  cout << "exchanges:  ";
  for (int m=0; m<nsystems-1; m++){
    cout << pair_exchanges[m] << "  ";
  }
  cout << endl;
  cout << "probs:  ";
  for (int m=0; m<nsystems-1; m++){
    cout << pair_probs[m] << "  ";
  }
  cout << endl;
  cout << "-----------------------------------" << endl;
}


void MultiSystem::WritePfunc(){
  pfunc_file << System->ncycles << ":  ";
  for (int i=0; i<nsystems; i++){
    pfunc_file << pfunc[i] << " ";
  }
  pfunc_file << endl;
}

void MultiSystem::WritePfuncInv(){
  pfunc_file << System->ncycles << ":  ";
  for (int i=0; i<nsystems; i++){
    pfunc_file << pfunc_inv[i] << " ";
  }
  pfunc_file << endl;
}

