#include "discretization/2_donor_cell.h"

//! constructor 
DonorCell::DonorCell(std::array<int,2> nCells, std::array<double,2> meshWidth, double alpha) :
    Discretization(nCells, meshWidth),
    alpha_(alpha)
{
    
}

//! compute the 1st derivative ∂ u^2 / ∂x
double DonorCell::computeDu2Dx(int i, int j) const 
{   
    // (copied from central differences)
    const double uCenterRight = (u(i+1,j) + u(i,j))   / 2.; // u at the center of the right cell (i+1,j)
    const double uCenterCell  = (u(i,j)   + u(i-1,j)) / 2.; // u at the center of the current cell (i,j)

    const double du2dx_cdiff = (uCenterRight*uCenterRight - uCenterCell*uCenterCell)/ dx(); 

    // donor cell part
    const double uCenterRight_donor = std::abs(uCenterRight)*((u(i,   j) - u(i+1,j)) / 2.);
    const double uCenterCell_donor  = std::abs(uCenterCell) *((u(i-1, j) - u(i  ,j)) / 2.); 

    const double du2dx_donor = (uCenterRight_donor - uCenterCell_donor) / dx(); 

    // return combined and weighted sum
    return  du2dx_cdiff + alpha_*du2dx_donor; 
}

//! compute the 1st derivative ∂ v^2 / ∂x
double DonorCell::computeDv2Dy(int i, int j) const {

    // central diff part
    const double vCenterAbove = (v(i,j+1) + v(i,j))   / 2.0; // v at the center of the cell above (i,j+1)
    const double vCenterCell  = (v(i,j)   + v(i,j-1)) / 2.0; // v at the center of the current cell (i,j) 

    const double dv2dy_cdiff = (vCenterAbove*vCenterAbove - vCenterCell*vCenterCell)/dy() ;

    // donor cell part
    const double vCenterAbove_donor = std::abs(vCenterAbove) * ((v(i,j  ) - v(i,j+1)) / 2.0);
    const double vCenterCell_donor  = std::abs(vCenterCell)  * ((v(i,j-1) - v(i, j )) / 2.0);

    const double dv2dy_donor = (vCenterAbove_donor - vCenterCell_donor) / dy();

    // return combined and weighted sum
    return dv2dy_cdiff + alpha_*dv2dy_donor; 
}

//! compute the 1st derivative ∂ (uv) / ∂x
double DonorCell::computeDuvDx(int i, int j) const {
    
    // central diff part
    const double uTopRight    = (u(i,j+1)    + u(i,j))    / 2.0; // u at top right corner of cell i,j
    const double uTopLeft     = (u(i-1, j+1) + u(i-1,j))  / 2.0; // u at top left corner
    
    const double vTopRight    = (v(i+1,j)    + v(i,j))    / 2.0; // v at the top right corner of cell i,j
    const double vTopLeft     = (v(i,j)      + v(i-1, j)) / 2.0; // v at top left corner

    const double duvdx_cdiff = (uTopRight * vTopRight - uTopLeft * vTopLeft)/dx(); // central difference at top edge of cell i,j    

    // donor cell part
    const double uvTopRight_donor   = std::abs(uTopRight) * ((v(i,j)   - v(i+1, j)) / 2.0);
    const double uvTopLeft_donor    = std::abs(uTopLeft)  * ((v(i-1,j) - v(i,j))    / 2.0); 

    const double duvdx_donor = (uvTopRight_donor - uvTopLeft_donor) / dx();

    // weight and sum
    return duvdx_cdiff + alpha_ * duvdx_donor;
}

//! compute the 1st derivative ∂ (uv) / ∂y
double DonorCell::computeDuvDy(int i, int j) const {

    // central diff part
    const double uTopRight    = (u(i,j+1) + u(i,j))     / 2.0; // avg. u at top right corner of cell i,j
    const double uBottomRight = (u(i,j)   + u(i,j-1))   / 2.0; // avg. u at bottom  right corner of cell i,j

    const double vTopRight    = (v(i+1,j)   + v(i,j))   / 2.0; // avg. v at the top right corner of cell i,j
    const double vBottomRight = (v(i+1,j-1) + v(i,j-1)) / 2.0; // avg. v at the top right corner of cell i,j
    
    const double duvdy_cdiff  = (uTopRight * vTopRight - uBottomRight * vBottomRight) / dy(); // central difference at right edge 

    // donor cell part
    const double uvTopRight_donor    = std::abs(vTopRight)    * ((u(i,j)   - u(i,j+1)) / 2.0);
    const double uvBottomRight_donor = std::abs(vBottomRight) * ((u(i,j-1) - u(i,j))   / 2.0);

    const double duvdy_donor = (uvTopRight_donor - uvBottomRight_donor) / dy();

    return duvdy_cdiff + alpha_ * duvdy_donor; 

}
