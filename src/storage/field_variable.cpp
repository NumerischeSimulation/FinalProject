#include "storage/field_variable.h"

FieldVariable::FieldVariable(std::array<int,2> size, std::array<double,2> origin, std::array<double,2> meshWidth) :
Array2D(size), 
origin_(origin), 
meshWidth_(meshWidth)
{};


double FieldVariable::interpolateAt(double x, double y) const // with x, y only in the computational domain 
{    
    const double dx = meshWidth_[0]; // mesh width in x dir. 
    const double dy = meshWidth_[1]; // mesh width in y dir. 

    // indicies of (cell (i,j) in which the point (x,y) lies
    int iLeftEdge  = (int) std::floor((x + origin_[0]) / dx); // -1 for  neg idx??
    int jLowerEdge = (int) std::floor((y + origin_[1]) / dy); 


    
    //  shift right and upper boundaries so that they don't use cells outside of the grid
    if (iLeftEdge >= (*this).size_[0] -1)
    {
        iLeftEdge = iLeftEdge -1; // shift it one column to the right
    }
    if (jLowerEdge >= (*this).size_[1] -1) 
    {
        jLowerEdge = jLowerEdge -1; // shift it down
    } 
    


    // relative position of x and y in the cell
    // one cell: |<-xr1-> x <-xr2->|
    //           |<--    dx      ->|
    const double xr1 = x  - ((meshWidth_[0]*iLeftEdge) - origin_[0]);   // relative position of x from left edge
    const double yr1 = y  - ((meshWidth_[1]*jLowerEdge) - origin_[1]);   // relative poistion of y from lower edge
    const double xr2 = dx - xr1;   // distance right_edge - x
    const double yr2 = dy - yr1;   // distance upper edge - y

    // transform to x, y coordinates when directly accessing the array2D
    int xLeftEdge = iLeftEdge;
    int yLowerEdge = jLowerEdge;

    // get values at corner points
    const double f_lowerLeft  = (*this)(xLeftEdge,     yLowerEdge);
    const double f_upperLeft  = (*this)(xLeftEdge,     yLowerEdge + 1); 
    const double f_lowerRight = (*this)(xLeftEdge + 1, yLowerEdge); 
    const double f_upperRight = (*this)(xLeftEdge + 1, yLowerEdge + 1);

    // bilinear interpolation
    const double f_intp = (f_lowerLeft * xr2 * yr2 
                        + f_lowerRight * xr1 * yr2
                        + f_upperLeft  * xr2 * yr1
                        + f_upperRight * xr1 * yr1) / (dx * dy);

return f_intp;

}
