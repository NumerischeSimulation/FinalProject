
#include "1_discretization.h"

// Constructor
Discretization::Discretization(std::array<int,2> nCells, std::array<double,2> meshWidth) :
    StaggeredGrid(nCells, meshWidth)
{
}


//! compute the 2nd derivative ∂^2 u / ∂x^2 at right edge
double Discretization::computeD2uDx2(int i, int j) const
{
    // return 2nd derivative ddu/dxx at right edge of cell i,j
    return 1./(dx()*dx()) * (u(i+1,j) - 2.*u(i,j) + u(i-1,j));

}

//! compute the 2nd derivative ∂^2 u / ∂y^2
double Discretization::computeD2uDy2(int i, int j) const 
{
    // return 2nd derivative ddu/dyy at right edge of cell i,j
    return 1./(dy()*dy()) * (u(i,j+1)-2.*u(i,j)+u(i,j-1));
}

//! compute the 2nd derivative ∂^2 v / ∂x^2
double Discretization::computeD2vDx2(int i, int j) const 
{

    // return 2nd derivative at top edge of cell
    return 1./(dx()*dx()) * (v(i+1,j)-2.*v(i,j)+v(i-1,j));
}

//! compute the 2nd derivative ∂^2 v / ∂y^2
double Discretization::computeD2vDy2(int i, int j) const 
{

    // return 2nd derivative at top edge of cell
    return 1./(dy()*dy()) * (v(i,j+1)-2.*v(i,j)+v(i,j-1));
}

//! compute the 1st derivative ∂p / ∂x
double Discretization::computeDpDx(int i, int j) const 
{

    // return 1st derivative at right edge of cell 
    return 1./dx() * (p(i+1,j) - p(i,j));
}

//! compute the 1st derivative ∂p / ∂y
double Discretization::computeDpDy(int i, int j) const 
{

    // return 1st derivative at top edge of cell i,j
    return 1./dy() * (p(i,j+1) - p(i,j));
}

