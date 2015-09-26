/*
 * File: cell_list.h
 * Description: Class for implementing a cell list to speed
 * up energy routines.
 * --------------------------------------------------------
 *
 * John R. Dowdle
 * The University of Texas at Austin
 * June 2008
 *
 */
#ifndef _cell_list_h
#define _cell_list_h

#include <iostream>
#include <vector>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_matrix.h>
#include "cell.h"
#include "constants.h"
#include "utils.h"

// namespace JLS_NS
// {

  class CellList 
  {

  public:

    /* Constructor
     *
     */
    CellList();

    /* Destructor
     *
     */
    virtual ~CellList();

    /* Maximum possible number of neighboring atoms   
     * Assumes d=3
     */
    static const int maxNeighbors = 27 * Cell::maxParticles;

    /* 
     * Maximum possible number of cells in list
     */
    static const int maxCells = 1000000;

    /* Member function:  build_new()
     * -----------------------------------------------------
     * Create a new cell list.
     */
    void build_new(const std::vector<double> &sideLengths, double cutoff, 
		   gsl_matrix *coords, int nparticles);
  
    /* Member function:  build()
     * -----------------------------------------------------
     * Build the cell list.
     */
    void build(const std::vector<double> &sideLengths, double cutoff, 
	       gsl_matrix *coords, int nparticles);

    /* Member function:  getNeighbors(r, neighbors)
     * ------------------------------------------------------------
     * Find the particles in cell containing the position vector
     * r and all of the particles in neighboring cells.
     */
    void getNeighbors(const std::vector<double> &r, std::vector<int> &neighbors,
		      int &nNeighbors);

    /* Member function:  cellIndexFromPosition(r)
     * ----------------------------------------------------
     * Return the index of the cell containing the position
     * "r".
     */
    int cellIndexFromPosition(const std::vector<double> &r) const;

    /* Member function:  getCellParticles(cellIndex, neighbors)
     * --------------------------------------------------------
     * Fill the vector "neighbors" with all of the particles in
     * the current cell.
     */
    void getCellParticles(int cellIndex, std::vector<int> &neighbors,
			  int &nNeighbors) const;

    /* Member function:  ntotcells()
     * --------------------------------------------------------
     * Return the total number of cells in the cell list.
     * 
     */
    int ntotcells() const;

  private:

    int ntotcells_;                     // total number of cells in list
    std::vector<int> ncells_;           // number of cells in each dimension
    std::vector<double> cellLengths_;   // size of the cells in each dimension
    std::vector<Cell> cells_;           // the cells
    std::vector<double> positionVector_;// a dummy vector used to store positions

    /* Member function: applyPBC(r)
     * ---------------------------------------------------
     * Apply periodic boundary conditions to the position
     * vector "r"
     */
    void applyPBC(std::vector<double> &vec);

  };


  /*
   * Return the index of the cell containing the position
   * "r".
   */
  inline int CellList::cellIndexFromPosition(const std::vector<double> &r) const 
  {
    // FIXME: make d-dimensional
    int n, nx, ny, nz;

    nx = int(r[0]/cellLengths_[0]);
    ny = int(r[1]/cellLengths_[1]);
    nz = int(r[2]/cellLengths_[2]);

    n = nx \
      + ny * ncells_[0] \
      + nz * ncells_[0] * ncells_[1];

    return n;
  }

  inline int CellList::ntotcells() const
  {
    return ntotcells_;
  }

  // }
#endif


