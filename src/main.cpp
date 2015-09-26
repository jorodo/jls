/* File: main.cpp
 * Description:  A simple atomistic Monte Carlo code
 * --------------------------------------------------
 * Author:  John R. Dowdle
 * Date:  June 2008
 * GPLv3.0+
 *
*/

// Libraries
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <sys/stat.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <time.h>
#include <assert.h>
#include "mpi.h"

#include "particle_system.h"
#include "constants.h"
#include "multi_system.h"

// Namespaces
using namespace std;

// Prototypes
void MultiSystemRun(string input_file_name, int nsystems, int nprocs, int rank,
		    string pfunc_fname, bool write_mode, bool infile_list);
void LoadPfunc(string pfunc_fname, vector<int>& pfunc,
	       vector<int>& pfunc_inv, int nsystems, int rank);
void StagedGrowth(string input_file_name, string growth_array_fname);
void LoadGrowthArray(string growth_array_fname, vector<double> &growth_array);

// Main
// int main(int argc, char *argv[]){
int main(int argc, char *argv[]){  // valgrind ??
  int nprocs, rank;
  int nsystems=1;
  string input_file_name, growth_array_fname="", pfunc_fname="";
  bool write_mode = false;
  bool infile_list = false;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  ParseCommandLine(argc,argv,input_file_name,nsystems,pfunc_fname,
		   growth_array_fname,write_mode,infile_list);

  if (nsystems>1){
    MultiSystemRun(input_file_name,nsystems,nprocs,rank,pfunc_fname,write_mode,
  		   infile_list);
  } else if (growth_array_fname.length()>0){
    StagedGrowth(input_file_name,growth_array_fname);
  } else {
    ParticleSystem* System;
    System = new ParticleSystem(input_file_name);
    System->Run();
    delete System;
  }

  // ParticleSystem* System;
  // System = new ParticleSystem(input_file_name);
  // System->Run();
  // delete System;

  MPI_Finalize();

  return 0;
}

// Function definitions
void MultiSystemRun(string input_file_name, int nsystems, int nprocs, int rank,
		    string pfunc_fname, bool write_mode, bool infile_list){
  string stnum;
  vector<int> pfunc, pfunc_inv;
  int tnum;

  if (pfunc_fname.length()>0){
    pfunc.resize(nsystems);
    pfunc_inv.resize(nsystems);
    //LoadPfunc(pfunc_fname, pfunc, pfunc_inv, nsystems, rank);
    LoadPfunc(pfunc_fname, pfunc_inv, pfunc, nsystems, rank);
    tnum = pfunc_inv[rank];
  } else {
    tnum = rank;
  }
  if (infile_list){  // the input file in the command line contains a list of files
    vector<string> lines;
    load_file(input_file_name,lines);
    input_file_name = lines[rank];
  } else {  // the input file in the command line is a basename
    stnum = int2str(tnum);
    input_file_name += stnum;
  }
  MultiSystem* MSystem;
  MSystem = new MultiSystem(input_file_name,nsystems,nprocs,rank,tnum,
			    pfunc_fname, write_mode);

  if (MSystem->System->nreplex > 0){
    if (nsystems != nprocs){
    	cerr << "error: replica exchange and Wang-Landau can currently only be run ";
    	cerr << "with processors and replicas in 1-to-1 correspondence." << endl;
    	exit(USAGE);
    }
    if (MSystem->System->seqWL || MSystem->System->seqWLTM){
	MSystem->WangLandau();
    } else {
      MSystem->ReplicaExchange();
    }
  } else {
    MSystem->RunAll();
  }

  delete MSystem;
}


void LoadPfunc(string pfunc_fname, vector<int>& pfunc,
	       vector<int>& pfunc_inv, int nsystems, int rank){
  string line, pvs_line;
  vector<string> tokens;

  if (rank==0)
    cout << "restarting REMC from " << pfunc_fname << endl;
  
  ifstream pfunc_in(pfunc_fname.c_str());
  if (pfunc_in.is_open()){
    while (!pfunc_in.eof()){
      pvs_line = line;
      getline(pfunc_in,line);
    }
    pfunc_in.close();
    Tokenize(pvs_line,tokens);
    for (int m=0; m<nsystems; m++){
      pfunc[m] = atof(tokens[m+1].c_str());
    }
    pfunc_in.close();
  } else {
    cout << "warning:  unable to open remc pfunc file." << endl;
  }

  for (int m=0; m < nsystems; m++){
    pfunc_inv[pfunc[m]] = m;
  }
}


void StagedGrowth(string input_file_name, string growth_array_fname){
  vector<double> growth_array;
  int nstages, proj_name_length;
  string stagen_str, stagen1_str, rst_name, rst_coord_name;
  ParticleSystem* System;
  double eng, hs_sigma;
  FILE* logfile;
  
  LoadGrowthArray(growth_array_fname,growth_array);
  nstages = (int)growth_array.size();

  System = new ParticleSystem(input_file_name);
  System->start_cycle = 0;
  System->ncycles = 0;
  System->end_cycle = System->ncycles_target;
  rst_coord_name = System->rst_coords;

  for (int stagen=0; stagen<nstages-1; stagen++){
    logfile = fopen(System->logfile_name.c_str(),"w");
    fprintf(logfile,"beginning insertion stage %d ...\n",stagen);
    fclose(logfile);
    System->Run();

    System->parse_input_file(System->rst_path);  // get adjusted growth_sigma
    hs_sigma = System->growth_sigma;
    System->hs_sigma = hs_sigma;
    System->particles[System->growth_idx].sigma = hs_sigma;
    
    System->ncycles = System->start_cycle;
    System->growth_sigma = growth_array[stagen+1];
    System->growth_sigma_avg = System->growth_sigma;
    System->ng_att = 0;
    System->ng_acc = 0;
    System->ng_att_cycle = 0;
    System->ng_acc_cycle = 0;
    System->p_growth = 0.;
    System->growth_sigma_sum = 0.;

    rst_coord_name = System->rst_dir + System->proj_name + "_gr.rst.cor";
    proj_name_length = (int)System->proj_name.length();
    stagen_str = int2str(stagen);
    System->proj_name.erase(proj_name_length-stagen_str.length(),
			    stagen_str.length());
    stagen1_str = int2str(stagen+1);
    System->proj_name.append(stagen1_str);
    System->name_files();
    System->restart_coords = rst_coord_name;
    System->assign_restart_coords();
    eng = System->get_system_energy();
    if(System->overlap){
      System->adjust_config(true);
    }
  }

  delete System;
}


void LoadGrowthArray(string growth_array_fname, vector<double> &growth_array){
  ifstream growth_file(growth_array_fname.c_str());
  vector<string> tokens;
  string line;

  if (growth_file.is_open()){
    while (!growth_file.eof()){
      getline(growth_file,line);
      Tokenize(line,tokens);
      for (int i=0; i<(int)tokens.size(); i++){
	growth_array.push_back(atof(tokens[i].c_str()));
      }
      tokens.clear();
    }
  } else {
    fileio_error(growth_array_fname);
  }

}


