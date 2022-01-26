#include "1_sor.h"
#include <math.h> 

SOR::SOR(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, double omega) :
    PressureSolver(discretization, epsilon, maximumNumberOfIterations),
    omega_(omega)
{
}

void SOR::solve()
{ 
    //int nCellsx = discretization_->nCells()[0]; // inner cells
    //int nCellsy = discretization_->nCells()[1]; // inner cells

    //sell size
    double dy = discretization_->dy();
    double dx = discretization_->dx();
    double dx2 = pow(dx,2);
    double dy2 = pow(dy,2);
    double factor = (dx2 * dy2) / ( 2.0 * (dx2 + dy2));

    int iteration = 0;
    
    //initial residual
    applyObstacleBoundaryValues();
    double res = calculateResidual();

   
    
    
     // setBoundaryValues();


    // iterate through grid 
     while( iteration < maximumNumberOfIterations_ && res > pow(epsilon_,2))
    {
        
        // one solver iteration
        for ( int i = discretization_->pIBegin() +1; i < discretization_->pIEnd() -1; i++)
        { 
            for ( int j = discretization_->pJBegin() +1; j < discretization_->pJEnd() -1; j++)
            {
                if (discretization_->isObstacleCell(i,j) != 1.)
                {
                    double sum_x = (discretization_->p(i+1, j) + discretization_->p(i-1, j)) / (dx2);
                    double sum_y = (discretization_->p(i, j+1) + discretization_->p(i, j-1)) / (dy2);
                    discretization_->p(i, j) = (1 - omega_) * discretization_->p(i, j) + omega_ * factor *(sum_x + sum_y - discretization_->rhs(i, j));
                }
            }
        }
        
       

        iteration +=1;
        
        //set new boundary values
        setBoundaryValues();
        

        res = calculateResidual();
    }
    
    std::cout << "SOR: " << iteration << " with a residuum of " << res << " from target " << std::pow(epsilon_,2) << std::endl;
       
}
