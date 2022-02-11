#include "0_staggered_grid.h"

StaggeredGrid::StaggeredGrid(std::array<int,2> nCells, std::array<double,2> meshWidth) :
  u_({nCells[0] + 1, nCells[1] + 2},   {0.0,      meshWidth[1]/2.}, meshWidth),
  v_({nCells[0] + 2, nCells[1] + 1},   {meshWidth[0]/2.,   0.},    meshWidth),
  p_({nCells[0] + 2, nCells[1] + 2},   {meshWidth[0]/2.,   meshWidth[1]/2.}, meshWidth),
  f_({nCells[0] + 1, nCells[1] + 2},   {0.,      meshWidth[1]/2.}, meshWidth),
  g_({nCells[0] + 2, nCells[1] + 1},   {meshWidth[0]/2.,   0.},    meshWidth),
  rhs_({nCells[0] + 2, nCells[1] + 2}, {meshWidth[0]/2.,   meshWidth[1]/2.}, meshWidth),
  isObstacleCell_({nCells[0] + 2, nCells[1] + 2},   {meshWidth[0]/2.,   meshWidth[1]/2.}, meshWidth),
  hasFluidNeighbourLeft_({nCells[0] + 2, nCells[1] + 2},   {meshWidth[0]/2.,   meshWidth[1]/2.}, meshWidth),
  hasFluidNeighbourRight_({nCells[0] + 2, nCells[1] + 2},   {meshWidth[0]/2.,   meshWidth[1]/2.}, meshWidth),
  hasFluidNeighbourTop_({nCells[0] + 2, nCells[1] + 2},   {meshWidth[0]/2.,   meshWidth[1]/2.}, meshWidth),
  hasFluidNeighbourBottom_({nCells[0] + 2, nCells[1] + 2},   {meshWidth[0]/2.,   meshWidth[1]/2.}, meshWidth),
  meshWidth_(meshWidth),
  nCells_(nCells)
{
    // init with zeros/false
    isObstacleCell_.setToZero();
    hasFluidNeighbourLeft_.setToZero();
    hasFluidNeighbourRight_.setToZero();
    hasFluidNeighbourTop_.setToZero();
    hasFluidNeighbourBottom_.setToZero();

    // init oldBoundaryValue
    // oldBoundaryValueTop_.resize(nCells[0], 0.0); // we go from 0 to n-1 on the boundary
    // oldBoundaryValueBottom_.resize(nCells[0], 0.0);
    // oldBoundaryValueLeft_.resize(nCells[1], 0.0);
    // oldBoundaryValueRight_.resize(nCells[1], 0.0);
}

void StaggeredGrid::setObstacleFlags(std::string pathToGeometry) 
{
    std::cout << "The path to the complex geometry: " << pathToGeometry << std::endl;

    // open file
    std::ifstream myFile(pathToGeometry);
    if(!myFile.is_open()) 
    {
        throw std::runtime_error("Could not open file");
    }

    // one row as a string with comma-separated values
    std::string rowString;

    // a single entry as a string
    std::string entryString;

    // read matrixDataFile row by row
    int i = 0; // column number
    int j = 0; // row number
    while (std::getline(myFile, rowString)) 
    {
        // convert rowString to a stream
        std::stringstream rowStringStream(rowString);

        // read rowStringStream entry by entry
        while (std::getline(rowStringStream, entryString, ',')) {
            int value = int(atof(entryString.c_str()));
            if (i < nCells_[0] && j < nCells_[1])
            {
                if(value == 1)
                {
                    // write top to bottom as the origin is in the bottom left corner
                    isObstacleCell(i, nCells_[1] - 1 - j) = 1.;
                    u(i, nCells_[1] - 1 - j) = NAN;
                    v(i, nCells_[1] - 1 - j) = NAN;
                    p(i, nCells_[1] - 1 - j) = NAN;
                    f(i, nCells_[1] - 1 - j) = NAN;
                    g(i, nCells_[1] - 1 - j) = NAN;
                    rhs(i, nCells_[1] - 1 - j) = NAN;
                }
            }
            ++i;
        }
        if (i != nCells_[0]) 
        {
            throw std::runtime_error("The image has too few columns!");
        }
    i = 0;
    ++j;
  }

  if (j != nCells_[1]) 
  {
    throw std::runtime_error("The image has too few rows!");
  }
}

void StaggeredGrid::setObstacleNeighbourFlags() 
{
    for ( int i = 0; i < nCells_[0]; i++)
    { 
        for (int j = 0; j < nCells_[1]; j++)
        {
            if (isObstacleCell(i,j)==1.)
            {
                // left neighbour
                if (isObstacleCell(i-1,j)!=1.)
                {
                    hasFluidNeighbourLeft(i,j)=1.;
                }
                // right neighbour
                if (isObstacleCell(i+1,j)!=1.)
                {
                    hasFluidNeighbourRight(i,j)=1.;
                }
                // top neighbour
                if (isObstacleCell(i,j+1)!=1.)
                {
                    hasFluidNeighbourTop(i,j)=1.;
                }
                // bottom neighbour
                if (isObstacleCell(i,j-1)!=1.)
                {
                    hasFluidNeighbourBottom(i,j)=1.;
                }
            }
            
        }
    }
}


// getter for parameters
const std::array<double,2> StaggeredGrid::meshWidth() const
{
    return meshWidth_;
}

const std::array<int,2> StaggeredGrid::nCells() const
{
    return nCells_;
}

// getter for data
const FieldVariable& StaggeredGrid::u() const
{
    return u_;
}

const FieldVariable& StaggeredGrid::v() const
{
    return v_;
}

const FieldVariable& StaggeredGrid::p() const
{
    return p_;
}

// obstacle flags
const FieldVariable& StaggeredGrid::isObstacleCell() const
{
    return isObstacleCell_;
}

const FieldVariable& StaggeredGrid::hasFluidNeighbourLeft() const
{
    return hasFluidNeighbourLeft_;
}

const FieldVariable& StaggeredGrid::hasFluidNeighbourRight() const
{
    return hasFluidNeighbourRight_;
}

const FieldVariable& StaggeredGrid::hasFluidNeighbourTop() const
{
    return hasFluidNeighbourTop_;
}

const FieldVariable& StaggeredGrid::hasFluidNeighbourBottom() const
{
    return hasFluidNeighbourBottom_;
}


// mesh width
double StaggeredGrid::dx() const
{
    return meshWidth_[0];
}

double StaggeredGrid::dy() const
{
    return meshWidth_[1];
}

// all the beginnings and endings
// u
int StaggeredGrid::uIBegin() const
{
    return -1;
}

int StaggeredGrid::uIEnd() const
{
    return nCells_[0];
}

int StaggeredGrid::uJBegin() const
{
    return -1;
}

int StaggeredGrid::uJEnd() const
{
    return nCells_[1]+1;
}

// v
int StaggeredGrid::vIBegin() const
{
    return -1;
}

int StaggeredGrid::vIEnd() const
{
    return nCells_[0]+1;
}

int StaggeredGrid::vJBegin() const
{
    return -1;
}

int StaggeredGrid::vJEnd() const
{
    return nCells_[1];
}

// p
int StaggeredGrid::pIBegin() const
{
    return -1;
}

int StaggeredGrid::pIEnd() const
{
    return nCells_[0]+1;
}

int StaggeredGrid::pJBegin() const
{
    return -1;
}

int StaggeredGrid::pJEnd() const
{
    return nCells_[1]+1;
}

// access
// u
double StaggeredGrid::u(int i, int j) const
{
    // check the validity of the indicies
    if (!(i >= uIBegin() && i < uIEnd()))
    {
        std::cout << "i-Index of u out of bounds error for index " << i << " in range " << uIBegin() << " to " << uIEnd() << std::endl;
        throw;
    }

    if (!(j >= uJBegin() && j < uJEnd()))
    {
        std::cout << "j-Index of u out of bounds error for index " << j << " in range " << uJBegin() << " to " << uJEnd() << std::endl;
        throw;
    }

    // transform to x, y
    int x = i +1;
    int y = j +1;
    return u_(x,y);
}

double& StaggeredGrid::u(int i, int j)
{
    // check the validity of the indicies
    if (!(i >= uIBegin() && i < uIEnd()))
    {
        std::cout << "i-Index of u out of bounds error for index " << i << " in range " << uIBegin() << " to " << uIEnd() << std::endl;
        throw;
    }

    if (!(j >= uJBegin() && j < uJEnd()))
    {
        std::cout << "i-Index of u out of bounds error for index " << j << " in range " << uJBegin() << " to " << uJEnd() << std::endl;
        throw;
    }

    // transform to x, y
    int x = i +1;
    int y = j +1;
    return u_(x,y);
}

// v
double StaggeredGrid::v(int i, int j) const
{
    // check the validity of the indicies
    if (!(i >= vIBegin() && i < vIEnd()))
    {
        std::cout << "i-Index of v out of bounds error for index " << i << " in range " << vIBegin() << " to " << vIEnd() << std::endl;
        throw;
    }

    if (!(j >= vJBegin() && j < vJEnd()))
    {
        std::cout << "j-Index of v out of bounds error for index " << j << " in range " << vJBegin() << " to " << vJEnd() << std::endl;
        throw;
    }

    // transform to x, y
    int x = i +1;
    int y = j +1;
    return v_(x,y);
}

double& StaggeredGrid::v(int i, int j)
{
    // check the validity of the indicies
    if (!(i >= vIBegin() && i < vIEnd()))
    {
        std::cout << "i-Index of v out of bounds error for index " << i << " in range " << vIBegin() << " to " << vIEnd() << std::endl;
        throw;
    }

    if (!(j >= vJBegin() && j < vJEnd()))
    {
        std::cout << "j-Index of v out of bounds error for index " << j << " in range " << vJBegin() << " to " << vJEnd() << std::endl;
        throw;
    }

    // transform to x, y
    int x = i +1;
    int y = j +1;
    return v_(x,y);
}

// p
double StaggeredGrid::p(int i, int j) const
{
    // check the validity of the indicies
    if (!(i >= pIBegin() && i < pIEnd()))
    {
        std::cout << "i-Index of p out of bounds error for index " << i << " in range " << pIBegin() << " to " << pIEnd() << std::endl;
        throw;
    }

    if (!(j >= pJBegin() && j < pJEnd()))
    {
        std::cout << "j-Index of p out of bounds error for index " << j << " in range " << pJBegin() << " to " << pJEnd() << std::endl;
        throw;
    }

    // transform to x, y
    int x = i +1;
    int y = j +1;
    return p_(x,y);
}

double& StaggeredGrid::p(int i, int j)
{
    // check the validity of the indicies
    if (!(i >= pIBegin() && i < pIEnd()))
    {
        std::cout << "i-Index of p out of bounds error for index " << i << " in range " << pIBegin() << " to " << pIEnd() << std::endl;
        throw;
    }

    if (!(j >= pJBegin() && j < pJEnd()))
    {
        std::cout << "j-Index of p out of bounds error for index " << j << " in range " << pJBegin() << " to " << pJEnd() << std::endl;
        throw;
    }

    // transform to x, y
    int x = i +1;
    int y = j +1;
    return p_(x,y);
}

// access intermediate results
double& StaggeredGrid::rhs(int i, int j)
{
    // check the validity of the indicies
    if (!(i >= pIBegin() && i < pIEnd()))
    {
        std::cout << "i-Index of rhs out of bounds error for index " << i << " in range " << pIBegin() << " to " << pIEnd() << std::endl;
        throw;  
    }

    if (!(j >= pJBegin() && j < pJEnd()))
    {
        std::cout << "j-Index of rhs out of bounds error for index " << j << " in range " << pJBegin() << " to " << pJEnd() << std::endl;
        throw;  
    }
  
    // transform to x, y
    int x = i +1;
    int y = j +1;
    return rhs_(x,y);
}

double& StaggeredGrid::f(int i, int j)
{
    // check the validity of the indicies
    if (!(i >= uIBegin() && i < uIEnd()))
    {
        std::cout << "i-Index of f out of bounds error for index " << i << " in range " << uIBegin() << " to " << uIEnd() << std::endl;
        throw;  
    }

    if (!(j >= uJBegin() && j < uJEnd()))
    {
        std::cout << "j-Index of f out of bounds error for index " << j << " in range " << uJBegin() << " to " << uJEnd() << std::endl;
        throw;  
    }

    // transform to x, y
    int x = i +1;
    int y = j +1;
    return f_(x,y);
}

double& StaggeredGrid::g(int i, int j)
{
    // check the validity of the indicies
    if (!(i >= vIBegin() && i < vIEnd()))
    {
        std::cout << "i-Index of g out of bounds error for index " << i << " in range " << vIBegin() << " to " << vIEnd() << std::endl;
        throw;  
    }

    if (!(j >= vJBegin() && j < vJEnd()))
    {
        std::cout << "j-Index of g out of bounds error for index " << j << " in range " << vJBegin() << " to " << vJEnd() << std::endl;
        throw;  
    }

    // transform to x, y
    int x = i +1;
    int y = j +1;
    return g_(x,y);
}

double& StaggeredGrid::isObstacleCell(int i, int j)
{
    // transform to x, y
    int x = i +1;
    int y = j +1;
    return isObstacleCell_(x,y);
}

double& StaggeredGrid::hasFluidNeighbourLeft(int i, int j)
{
    // transform to x, y
    int x = i +1;
    int y = j +1;
    return hasFluidNeighbourLeft_(x,y);
}
double& StaggeredGrid::hasFluidNeighbourRight(int i, int j)
{
    // transform to x, y
    int x = i +1;
    int y = j +1;
    return hasFluidNeighbourRight_(x,y);
}
double& StaggeredGrid::hasFluidNeighbourTop(int i, int j)
{
    // transform to x, y
    int x = i +1;
    int y = j +1;
    return hasFluidNeighbourTop_(x,y);
}
double& StaggeredGrid::hasFluidNeighbourBottom(int i, int j)
{
    // transform to x, y
    int x = i +1;
    int y = j +1;
    return hasFluidNeighbourBottom_(x,y);
}