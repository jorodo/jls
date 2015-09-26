#ifndef __replica_exchange_h__
#define __replica_exchange_h__

// Libraries
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <sys/stat.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <time.h>
#include <assert.h>
#include "mpi.h"

#include "particle_system.h"
#include "constants.h"

// Namespaces
using namespace std;

class MultiSystem {

 public:
  
  int remc_end_cycle, nrun, proj_name_length, tnum;
  double *pot_eng, *vol, *pres, *beta, *temp;
  double *pair_probs, *dr_max, *dlnV_max, *dL_max, *beta_inv;
  double *temp_inv, *dr_max_inv, *dlnV_max_inv, *dL_max_inv;
  double *engHistMins, *engHistMaxs, *allCoords, *allSeqEng;
  double *gModFactors, gModTol;
  bool gModsLTtol;
  int *nvol_att, *nvol_acc, *npart_acc, *npart_att, *pfunc;
  int *nvol_att_inv, *npart_acc_inv, *npart_att_inv;
  int *pfunc_inv, *nvol_acc_inv;
  int *pair_attempts, *pair_exchanges, *seqPtypes, *allSeqPtypes;
  int nsystems, nprocs, rank;
  string pfunc_fname, tnum_str, pvs_tnum_str;
  ofstream pfunc_file;
  bool write_mode_temps;
  FILE* logfile;
  ParticleSystem* System;


  MultiSystem();
  MultiSystem(string input_file_name, int num_systems, 
	      int num_procs, int myrank, int mytnum=-1,
	      string mypfunc_fname="", bool write_mode=false);
  ~MultiSystem();
  void ReplicaExchange();
  void WangLandau();
  void RunAll();

 private:

  void MPISendToMaster();
  void MPIMasterSend();
  void UpdateParams();
  void MPISendToMasterWL();
  void MPIMasterSendWL();
  void UpdateParamsWL();
  void AttemptREMCExchanges();
  void AttemptWLExchanges();
  void CheckModFactors();
  void WriteExchangeStats();
  void WritePfunc();
  void WritePfuncInv();

};



#endif
