#include "discretization/2_central_differences.h"

//! compute the 1st derivative ∂ u^2 / ∂x
double CentralDifferences::computeDu2Dx(int i, int j) const 
{
    const double uCenterRight = (u(i+1,j) + u(i,j))   / 2.; // u at the center of the right cell (i+1,j)
    const double uCenterCell  = (u(i,j)   + u(i-1,j)) / 2.; // u at the center of the current cell (i,j)

    // return central difference of squared values
    return  (uCenterRight*uCenterRight - uCenterCell*uCenterCell)/dx(); 
}

//! compute the 1st derivative ∂ v^2 / ∂x 
double CentralDifferences::computeDv2Dy(int i, int j) const
{
    const double vCenterAbove = (v(i,j+1) + v(i,j))   / 2.0; // v at the center of the cell above (i,j+1)
    const double vCenterCell  = (v(i,j)   + v(i,j-1)) / 2.0; // v at the center of the current cell (i,j) 

    // return central difference of squared values
    return (vCenterAbove*vCenterAbove - vCenterCell*vCenterCell)/dy() ;
}

//! compute the 1st derivative ∂ (uv) / ∂x at top edge of cell (for v)
double CentralDifferences::computeDuvDx(int i, int j) const 
{

    const double uTopRight    = (u(i,j+1)    + u(i,j))    / 2.0; // u at top right corner of cell i,j
    const double uTopLeft     = (u(i-1, j+1) + u(i-1,j))  / 2.0; // u at top left corner
    
    const double vTopRight    = (v(i+1,j)    + v(i,j))    / 2.0; // v at the top right corner of cell i,j
    const double vTopLeft     = (v(i,j)      + v(i-1, j)) / 2.0; // v at top left corner

    return (uTopRight * vTopRight - uTopLeft * vTopLeft)/dx(); // central difference at top edge of cell i,j
}

//! compute the 1st derivative ∂ (uv) / ∂y at right edge of cell (for u)
double CentralDifferences::computeDuvDy(int i, int j) const 
{
    const double uTopRight    = (u(i,j+1) + u(i,j))     / 2.0; // avg. u at top right corner of cell i,j
    const double uBottomRight = (u(i,j)   + u(i,j-1))   / 2.0; // avg. u at bottom  right corner of cell i,j

    const double vTopRight    = (v(i+1,j)   + v(i,j))   / 2.0; // avg. v at the top right corner of cell i,j
    const double vBottomRight = (v(i+1,j-1) + v(i,j-1)) / 2.0; // avg. v at the top right corner of cell i,j

    return (uTopRight * vTopRight - uBottomRight * vBottomRight)/dy(); // central difference at right edge
}