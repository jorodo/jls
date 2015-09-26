/* jrd
 * GPLv3+
 */

// FIXME: make sure the position vector is 'const' and can't get changed

#include "cell_list.h"

// using namespace JLS_NS;

/*
 * Default constructor
 */
CellList::CellList()
{

}

/*
 * Destructor (virtual)
 */
CellList::~CellList()
{

}

/*
 * Create a new cell list
 */
void CellList::build_new(const std::vector<double> &sideLengths, double cutoff, 
			 gsl_matrix *coords, int nparticles)
{

  int ndims = sideLengths.size();

  // Only use cell lists for d=3 at the moment
  assert(ndims==3);

  ntotcells_ = 0;
  positionVector_.resize(ndims);
  ncells_.resize(ndims);
  cellLengths_.resize(ndims);
  cells_.resize(maxCells);
  build(sideLengths, cutoff, coords, nparticles);

}

/*
 * build the cell list
 */
void CellList::build(const std::vector<double> &sideLengths, double cutoff, 
		     gsl_matrix *coords, int nparticles)
{
  int nk, cellIndex;
  int ndims = sideLengths.size();

  // clear the particles from the current cells
  for (int i = 0; i < ntotcells_; i++){
    cells_[i].clearParticles();
  }

  // positionVector_.resize(ndims); ////
  ntotcells_ = 1;

  // determine the number of cells in each direction and the dimensions of 
  // the cells.
  for (int k = 0; k < ndims; k++){
    nk = int(sideLengths[k] / cutoff);
    ntotcells_ *= nk;
    ncells_[k] = nk;
    cellLengths_[k] = sideLengths[k] / ncells_[k];
  }
  
  // cells_.clear();
  // cells_.resize(ntotcells_);

  // sort the particles into their cells
  for (int i = 0; i < nparticles; i++){
    for (int k = 0; k < ndims; k++){
      positionVector_[k] = gsl_matrix_get(coords, i, k);
    }
    cellIndex = cellIndexFromPosition(positionVector_);
    cells_[cellIndex].addParticle(i);
  }

}

/*
 * Fill the neighbors array with all particles from the cell cellIndex
 * and all of the particles from its neighboring cells
 */
void CellList::getNeighbors(const std::vector<double> &r, std::vector<int> &neighbors,
			    int &nNeighbors)
{
  double xDelta, yDelta, zDelta, cellIndex;
  // int ndims = (int)ncells_.size();

  // positionVector_.resize(ndims);
  nNeighbors = 0;

  // FIXME:  make d-dimensional.  right now only works in d=3
  for (int i=-1; i< 2; i++){
    xDelta = i * cellLengths_[0];
    positionVector_[0] = r[0] + xDelta;
    for (int j=-1; j< 2; j++){
      yDelta = j * cellLengths_[1];
      positionVector_[1] = r[1] + yDelta;
      for (int k=-1; k< 2; k++){
	zDelta = k * cellLengths_[2];
	positionVector_[2] = r[2] + zDelta;
	applyPBC(positionVector_);
	cellIndex = cellIndexFromPosition(positionVector_);
	getCellParticles(cellIndex, neighbors, nNeighbors);
      }
    }
  }

  assert(nNeighbors <= maxNeighbors);

}

/*
 * Fill the neighbors vector with all of the particles from the
 * cell cellIndex
 */
void CellList::getCellParticles(int cellIndex, std::vector<int> &neighbors,
				int &nNeighbors) const
{
  int idx;

  for (int i = 0; i < cells_[cellIndex].nParticles(); i++){
    idx = i + nNeighbors;
    neighbors[idx] = cells_[cellIndex].particles[i];
  }

  nNeighbors += cells_[cellIndex].nParticles();
}

/*
 * Aply periodic boundary condition if particle moves outsie of box
 */
void CellList::applyPBC(std::vector<double> &vec)
{
  double sideLength;
  int ndims = ncells_.size();
  
  for (int k = 0; k < ndims; k++){
    sideLength = ncells_[k] * cellLengths_[k];
    vec[k] = modulo(vec[k], sideLength);
    // if (vec[k] >= sideLength) vec[k] -= sideLength;
    // else if (vec[k] < 0) vec[k] += sideLength;
  }
}
