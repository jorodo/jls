/*
 * File: cell.h
 * Description: Class for a single cell in a cell list
 * --------------------------------------------------------
 *
 * John R. Dowdle
 * The University of Texas at Austin
 * June 2008
 *
 */

#ifndef _cell_h
#define _cell_h

#include <vector>
#include <assert.h>
#include <iostream>
#include "constants.h"

// namespace JLS_NS
// {

  class Cell
  {

  public:

    /*
     *  Maximum number of particles allowed in a given cell.
     *  
     */
    static const int maxParticles = 40;

    /* 
     * Vector containing the particle indices of particles in the cell
     */
    std::vector<int> particles;          

    /*
     * Constructor
     */
    Cell();

    /*
     * Destructor
     */
    virtual ~Cell();

    /* Member function: nParticles()
     * -----------------------------------------------------
     * Return the number of particles in the cell.
     */ 
    int nParticles() const;


    /* Member function: addParticle(particleIndex)
     * -----------------------------------------------------
     * Add a particle to the cell
     */
    void addParticle(int i);

    /* Member function:  clearParticles()
     * -----------------------------------------------------
     * Remove all particles from the cell
     */
    void clearParticles();

  private:

    int nParticles_;  // number of particles currently in cell

  };


  inline void Cell::addParticle(int i)
  { 

    try {
      if (nParticles_ + 1 > maxParticles) {
	throw USAGE;  // FIXME: make this a proper error and exception handling
      }
    } catch(int e) {
      std::cerr << "Error: maximum number of atoms in cell has been exceeded" << std::endl;
      // exit(e);
    }

    particles[nParticles_] = i;
    nParticles_++;
  }

  inline int Cell::nParticles() const 
  {
    return nParticles_;
  }

  inline void Cell::clearParticles()
  {
    nParticles_ = 0;
  }

// }
#endif
