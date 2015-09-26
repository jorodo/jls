#include "utils.h"

void Usage(char *str) {
  cout << "Usage: " << str;
  cout << "  [-r <pfunc_fname>] [-n <num_systems>] [-f <input_file>]" << endl;
  cout << "  [-l <input_file_list>] [-g <growth_array_fname>]" << endl;
  cout << "  [-w <toggle_write_mode>] [-h print help]" << endl;
  cout << endl;
  cout << "See the example input file in jls/files/." << endl;
}


void ParseCommandLine(int argc, char *argv[], string& input_file_name, \
		      int& nsystems, string& pfunc_fname, string &growth_array_fname,
		      bool &write_mode, bool &infile_list){
  int option_char;
  int num_opts;
  opterr = 0;
  num_opts = 0;
  
  while ((option_char = 
	  getopt(argc, argv, "r:n:f:w:g:l:h")) != -1){
    switch (option_char){
    case 'r':
      pfunc_fname = optarg;
      break;
    case 'n':
      nsystems = atoi(optarg);
      break;
    case 'f':
      input_file_name = optarg;
      break;
    case 'w':
      write_mode = false;
      break;
    case 'g':
      growth_array_fname = optarg;
      break;
    case 'l':
      input_file_name = optarg;
      infile_list = true;
      break;
    case 'h':
      Usage(argv[0]);
      exit(SUCCESS);
      break;
    case '?':
      if (isprint (optopt))
	fprintf (stderr, "Error:  Unknown option `-%c'.\n", optopt);
      else
	fprintf (stderr,
		 "Error:  Unknown option character `\\x%x'.\n",
		 optopt);
      Usage(argv[0]);
      abort();
    default:
      Usage(argv[0]);
      exit(USAGE);
    }
    num_opts++;
  }
  if (num_opts==0){
      Usage(argv[0]);
      exit(USAGE);
  }
}


string ParticleLabel(int p){
  switch (p) {
  case 0:
    return "LJ";  // Lennard-Jones
  case 1:
    return "HS";  // Hard Sphere
  case 2:
    return "JG";  // Jagla
  case 3:
    return "SW";  // Square-Well
  default:
    cout << endl;
    cout << "error:  unknown particle type" << endl;
    exit(USAGE);
  }
}

int InverseParticleLabel(string s){

  if (s=="LJ")
    return 0;
  else if (s=="HS")
    return 1;
  else if (s=="JG")
    return 2;
  else if (s=="SW")
    return 3;
  else {
    cout << endl;
    cout << "error:  unknown particle label" << endl;
    exit(USAGE);
  }
}

gsl_rng* InitRNG(bool seed_random, int rank, unsigned long int seed){
  const gsl_rng_type *rng_type;
  gsl_rng *rng;
  int ntest;

  ntest = 4;
  rng_type = gsl_rng_default;
  rng = gsl_rng_alloc(rng_type);
  //cout <<"rank" << rank << ":  Initializing the RNG ...\n" << endl;

  if (seed_random){
    seed = random_seed(rank);
  } 
  // else {
  //   gsl_rng_env_setup();
  // }
  gsl_rng_set(rng,seed);

  // talk to user
  // printf("generator type: %s\n", gsl_rng_name(rng));
  // printf("seed = %lu\n", seed);
  // printf("first value = %lu\n", gsl_rng_get(rng));
  // printf("\n");

  return rng;

}

void RecursiveIncrement(int idx, int lim, int sup, vector<int> &n){
  if (idx==sup)
    return;
  else {
    if (n[idx]==lim){
      n[idx] = 0;
      n[idx+1]++;
      RecursiveIncrement(idx+1,lim,sup,n);
    }
    else return;
  }
}

// Check to see if directory exists
// by attempting to open a new file
// for output within it.
bool dir_exists(string fname) {
  size_t len = fname.length();
  if(fname[len-1] != '/' && fname[len-1] != '\\')
    fname.append("/");
  fname.append("000.tmp");
  ofstream outf(fname.c_str());
  bool existFlag = outf;
  if(outf) {
    outf.close();
    remove(fname.c_str());
  }
  return existFlag;
}



// Check to see if file exists
bool file_exists(string fname) {
  struct stat stFileInfo;
  bool blnReturn;
  int intStat;

  // Attempt to get the file attributes
  intStat = stat(fname.c_str(),&stFileInfo);
  if(intStat == 0) {
    // We were able to get the file attributes
    // so the file obviously exists.
    blnReturn = true;
  } else {
    // We were not able to get the file attributes.
    // This may mean that we don't have permission to
    // access the folder which contains this file. If you
    // need to do that level of checking, lookup the
    // return values of stat which will give you
    // more details on why stat failed.
    blnReturn = false;
  }
  
  return(blnReturn);
}

void bkup_file(string fname){
  int result;
  string new_name;
  string ext;

  if (file_exists(fname)){
    ext = int2str((int)time(NULL));
    new_name = fname + ext;
    result = rename(fname.c_str(),new_name.c_str());
    if ( result == 0 )
      fprintf(stdout,"backed up file %s to %s\n",
	      fname.c_str(),new_name.c_str());
    else
      fprintf(stderr,"warning: could not move %s to %s\n",
	      fname.c_str(),new_name.c_str());
  }
}

// return the surface area constant, Cn, for ndims
double nsphere_SA_constant(int ndims){
  int k;
  double num, denom, Cn;
  if (ndims < 1){
    fprintf(stderr,"error: ndims less than 1");
    exit(USAGE);
  } else if (ndims == 1){
    Cn = 1.;
  }
  if (ndims%2 == 0){
    k = ndims/2;
    Cn = gsl_pow_int(M_PI,k)/gsl_sf_fact(k);
  } else {
    k = (ndims-1)/2;
    num = gsl_pow_int(2.,2*k+1)*gsl_sf_fact(k)*gsl_pow_int(M_PI,k);
    denom = gsl_sf_fact(2*k+1);
    Cn = num/denom;
  }

  return Cn;
}


// remove leading and trailing whitespace
void strip_string(string& string_to_modify){
  if(string_to_modify.empty()) return;

  int start_index = string_to_modify.find_first_not_of(" ");
  int end_index = string_to_modify.find_last_not_of(" ");
  string temp_string = string_to_modify;
  string_to_modify.erase();

  string_to_modify = temp_string.substr(start_index,(end_index-start_index+1));
}
    
// Euclidean distance between two points i and j in a periodic cell
double GetDistance(vector<double> &coords_i, vector<double> &coords_j,
		   vector<double> &side_lengths){
  double distance, dist_sqrd, delta;
  int ndims;
  ndims = coords_i.size();
  dist_sqrd = 0.0;

  for (int k=0; k<ndims; k++){
    delta = coords_i[k] - coords_j[k];
    delta -= side_lengths[k] * (double)round_(delta / side_lengths[k]);
    // if (delta > side_lengths[k]/2.0) delta -= side_lengths[k];        //minimum image
    // else if (delta < -side_lengths[k]/2.0) delta += side_lengths[k];
    dist_sqrd += delta*delta;
  }
  distance = sqrt(dist_sqrd);
  
  return distance;
}

//void Tokenize(const string& str, vector<string>& tokens, const string& delimiters) {
//  // Skip delimiters at beginning.
//  string::size_type lastPos = str.find_first_not_of(delimiters, 0);
//  // Find first "non-delimiter".
//  string::size_type pos     = str.find_first_of(delimiters, lastPos);
//
//  while (string::npos != pos || string::npos != lastPos)
//    {
//      // Found a token, add it to the vector.
//      tokens.push_back(str.substr(lastPos, pos - lastPos));
//      // Skip delimiters.  Note the "not_of"
//      lastPos = str.find_first_not_of(delimiters, pos);
//      // Find next "non-delimiter"
//      pos = str.find_first_of(delimiters, lastPos);
//    }
//}


// Tokenize based on either whitespaces or tabs
// Similar to pythons str.split() method
void Tokenize(const string& str, vector<string>& tokens){
  char white_space = ' ';
  char tab = '\t';
  int start_pos, stop_pos;
  int str_length = int(str.size());

  start_pos = 0;
  stop_pos = 0;

    while (true){
      // move forward until non-whitespace char found
      while (true){
	if (start_pos >= str_length) break;
	if (str[start_pos] == white_space){
	  start_pos += 1;
	} else if (str[start_pos] == tab){
	  start_pos += 1;
	} else {
	  break;
	}
      }

      if (start_pos >= str_length) break;
    
      stop_pos = start_pos + 1;  // sub-string must have lenght 1 or greater
      // move the stop position forward until whitespace/tab/end is found
      while (true){
	if (stop_pos > str_length){
	  break;
	} else if (str[stop_pos] == white_space){
	  break;
	} else if (str[stop_pos] == tab){
	  break;
	} else {
	  stop_pos += 1;
	}
      }
  
      tokens.push_back(str.substr(start_pos,stop_pos-start_pos));
      start_pos = stop_pos;
    }
}

void fileio_error(string s){
  cerr << "error: unable to open file " << s << endl;
  exit(FILEIO);
}


string int2str(int n){
  ostringstream oss;
  oss << n;
  return oss.str();
}

void load_file(string fname, vector<string> &lines){
  ifstream infile(fname.c_str());
  string line;

  if (infile.is_open()){
    while (!infile.eof()){
      getline(infile,line);
      lines.push_back(line);
    }
  } else {
    fileio_error(fname);
  }
}

// apply the Metropolis acceptance criterion
bool accept_metropolis(const gsl_rng *rng, double beta, 
		       double eng_old, double eng_new){
  double p_acc,rand_real;
  p_acc = min(1.0,exp(-beta*(eng_new-eng_old)));
  rand_real = gsl_rng_uniform(rng);
  //if (verbose && ncycles%stdout_freq==0){
  //  fprintf(logfile,"eng_old:  %f\n",eng_old);
  //  fprintf(logfile,"eng_new:  %f\n",eng_new);
  //  fprintf(logfile,"p_acc:  %f\n",p_acc);
  //  fprintf(logfile,"rand_real:  %f\n",rand_real);
  //}

  return (rand_real < p_acc);
}

// accept based on density of states (DOS) for flat 
// histogram methods
bool accept_dos(const gsl_rng *rng, double log_gE_old, 
		double log_gE_new){
  double p_acc,rand_real;
  p_acc = min(1.0, exp(log_gE_old - log_gE_new));
  rand_real = gsl_rng_uniform(rng);

  return (rand_real < p_acc);
}
