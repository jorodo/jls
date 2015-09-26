/*
 * File: random_seed.h
 * Description: Random seed generator
 * -------------------------------------------
 *
 * John R. Dowdle
 * June 2008
 *
 */
#ifndef _random_seed_h
#define _random_seed_h

#include <stdio.h>
#include <sys/time.h>


/*
 * Function: random_seed
 * Usage: random_seed();
 * ----------------------
 * Generate a random seed from /dev/random.  If
 * /dev/random is unavailable, use the system clock.
 *
 */
unsigned long int random_seed(int rank);

#endif
