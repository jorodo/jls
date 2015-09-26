/*
 * File: random_seed.cpp
 * Description: Random seed generator
 * -------------------------------------------
 *
 * John R. Dowdle
 * June 2008
 *
 */
#include "random_seed.h"

/*
 * Function: random_seed
 * Usage: random_seed();
 * ----------------------
 * Generate a random seed from /dev/random.  If
 * /dev/random is unavailable, use the system clock.
 *
 */
unsigned long int random_seed(int rank){
 unsigned int seed;
 struct timeval tv;
 FILE *devrandom;
 size_t result;

 if ((devrandom = fopen("/dev/random","r")) == NULL) {
   gettimeofday(&tv,0);
   seed = tv.tv_sec + tv.tv_usec;
   printf("rank%d:  Got seed %u from gettimeofday()\n",rank,seed);
 } else {
   result = fread(&seed,sizeof(seed),1,devrandom);
   printf("rank%d:  Got seed %u from /dev/random\n",rank,seed);
   fclose(devrandom);
 }

 return seed;

}
