#ifndef _utils_h
#define _utils_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <math.h>
#include <assert.h>
#include <getopt.h>
#include <time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_histogram.h>
#include <sys/stat.h>

#include "random_seed.h"
#include "constants.h"

using namespace std;


struct _HistogramT {
  gsl_histogram *hist;
  int nbins;
};

// The HistogramT type
typedef struct _HistogramT HistogramT;

/* Function:  Usage
 * Usage:  Usage()
 * -----------------------------------------------------------
 * Tell the user how to correctly use JLS
 */
void Usage(char *str);


/* Function:  ParseCommandLine
 * Usage:  ParseCommandLine(argc,argv,input_fname,nsystems,
 *                          pfunc_fname, growth_array_fname,
 *                          write_mode)
 * -----------------------------------------------------------
 * Read the command line arguments and assign to simulation
 * parameters
 */
void ParseCommandLine(int argc, char *argv[], string& input_file_name, \
		      int& nsystems, string& pfunc_fname, 
		      string &growth_array_fname, bool &write_mode, 
		      bool &infile_list);


/* Function:  ParticleLabel
 * Usage:  ParticleLabel()
 * -----------------------------------------------------------
 * Return the label string for the given particle type
 */
string ParticleLabel(int p);


/* Function:  InverseParticleLabel
 * Usage:  InverseParticleLabel(s)
 * -----------------------------------------------------------
 * Return the label string for the given particle type
 */
int InverseParticleLabel(string s);


/* Function:  InitRNG
 * Usage:  InitRNG(seed_random_true,rank)
 * -----------------------------------------------------------
 * Return the RNG
 */
gsl_rng* InitRNG(bool seed_random, int rank=0,
		 unsigned long int seed=2833311595);


/* Function:  RecursiveIncrement
 * Usage:  RecursiveIncrement(idx,lim,sup,n)
 * -----------------------------------------------------------
 * Recursively increment n idxs for lattice positions
 */
void RecursiveIncrement(int idx, int lim, int sup, vector<int> &n);



/* Function:  dir_exists
 * Usage:  dir_exists(dir_name)
 * -----------------------------------------------------------
 * True if directory exists, false otherwise
 */
bool dir_exists(string fname);



/* Function:  file_exists
 * Usage:  file_exists(file_name)
 * -----------------------------------------------------------
 * True if file exists, false otherwise
 */
bool file_exists(string fname);


/* Function:  bkup_file
 * Usage:  bkup_file(fpath)
 * -----------------------------------------------------------
 * Backup file fpath.ext to fpath.ext<time>, where <time> is 
 * the unix clock time in seconds since the beginning of time
 * (midnight Jan 1, 1970, of course)
 */
void bkup_file(string fname);


/* Function:  nsphere_SA_constant
 * Usage:  nsphere_SA_constant(ndims)
 * -----------------------------------------------------------
 * Return the prefactor for the surface area of the n-sphere
 */
double nsphere_SA_constant(int ndims);


/* Function:  strip_string
 * Usage:  strip_string(str)
 * -----------------------------------------------------------
 * Strip leading and trailing white space from str
 */
void strip_string(string& StringToModify);


/* Function:  GetDistance
 * Usage:  GetDistance(coords1,coords2,side_lengths)
 * -----------------------------------------------------------
 * Return the Euclidean distance between two points in a 
 * periodic cell
 */
double GetDistance(vector<double> &coords_i, vector<double> &coords_j,
		   vector<double> &side_lengths);


/* Function:  Tokenize
 * Usage:  Tokenize(str,tokens)
 * -----------------------------------------------------------
 * Tokenize based on whitespaces and tabs as delimiters.
 * Similar to pythons str.split() method
 */
void Tokenize(const string& str, vector<string>& tokens);


/* Function:  fileio_error
 * Usage:  fileio_error(s)
 * -----------------------------------------------------------
 * Exit with a file I/O error for file s
 */
void fileio_error(string s);


/* Function:  int2str
 * Usage:  int2str(n)
 * -----------------------------------------------------------
 * Return a string of the integer n
 */
string int2str(int n);


/* Function:  sign
 * Usage:  sign(x)
 * -----------------------------------------------------------
 * Return the sign of the double x
 */
double sign(double x);


/* Function:  modulo
 * Usage:  modulo(x)
 * -----------------------------------------------------------
 * Modulo of doubles x and y.  Works identically to Python's
 * "%" operator.  Works for y > 0.
 */
double modulo(double x, double y);

/* Function:  round_
 * Usage:  round_(x)
 * -----------------------------------------------------------
 * Round the double "x" to the nearest integer
 * 
 */
double round_(double x);

/* Function:  load_file
 * Usage:  load_file(fname,lines)
 * -----------------------------------------------------------------
 * Load all of the lines in the file fname into the vector of strings,
 * lines.  Each line in the file is an element in the vector.
 */
void load_file(string fname, vector<string> &lines);

/* Function:  accept_metropolis
 * Usage:  accept_metropolis()
 * -----------------------------------------------------------
 * Accept the move to system with eng_new with the Metropolis
 * criterion.
 * 
 */
bool accept_metropolis(const gsl_rng *rng, double beta, double eng_old, 
		       double eng_new);

/* Function:  accept_dos
 * Usage:  accept_dos()
 * -----------------------------------------------------------
 * Accept the move based on the current estimate of the 
 * density of states.  Used for flat-histogram methods.
 * 
 */
bool accept_dos(const gsl_rng *rng, double log_gE_old, double log_gE_new);

/****************** inlined functions ************************/

inline double sign(double x){
  return (x>=0)?1.0:-1.0;
}

inline double modulo(double x, double y){
  double r;

  r = fmod(x,y);
  return (sign(r) == 1.) ? r : r + y;
}

inline double round_(double x)
{
  return x < 0.0 ? ceil(x - 0.5) : floor(x + 0.5);
}

#endif
