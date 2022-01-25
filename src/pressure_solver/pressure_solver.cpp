#include "pressure_solver.h"
#include <math.h> 

#include <memory>

PressureSolver::PressureSolver(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations):
  discretization_(discretization),
  epsilon_(epsilon),
  maximumNumberOfIterations_(maximumNumberOfIterations)
{
}

  void PressureSolver::setBoundaryValues()
  {
     
      //p_i,-1 = p_i,0. p_i,n+1 = p_i,n.
      for ( int i = discretization_->pIBegin() +1; i < discretization_->pIEnd() -1; i++)
      {
        // bottom
        discretization_->p(i, discretization_->pJBegin())  = discretization_->p(i, discretization_->pJBegin() +1);
        // top
        discretization_->p(i, discretization_->pJEnd() -1) = discretization_->p(i, discretization_->pJEnd() -2);
      }

       // prioritise left and right boundaries
       // p_-1,j = p_0,j. p_n+1,j = p_n,j.
      for ( int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++)
      {
        // left
        discretization_->p(discretization_->pIBegin(),j) =  discretization_->p(discretization_->pIBegin()+1,j);
        // right
        discretization_->p(discretization_->pIEnd() -1,j) =  discretization_->p(discretization_->pIEnd() -2,j);
      }

       
  }

  double PressureSolver::calculateResidual()
  {
    int nCellsx = discretization_->nCells()[0] -2; // inner cells
    int nCellsy = discretization_->nCells()[1] -2; // inner cells

    //sell size
    double dy = discretization_->dy();
    double dx = discretization_->dx();
    double dx2 = pow(dx,2);
    double dy2 = pow(dy,2);

   
    double res = 0.0;

        for ( int i = discretization_->pIBegin() +1; i < discretization_->pIEnd() -1; i++)
        { 
            for ( int j= discretization_->pJBegin() +1; j < discretization_->pJEnd() -1; j++)
            {
                // calculate residual 
                double pxx = (discretization_->p(i-1, j) - 2.0 *discretization_->p(i,j) + discretization_->p(i+1, j)) / (dx2);
                double pyy = (discretization_->p(i, j-1) - 2.0 *discretization_->p(i,j) + discretization_->p(i, j+1)) / (dy2);

                double resij = discretization_->rhs(i, j) - pxx - pyy;   
                res = res + (pow(resij,2));
            }
        }
        
        //calculate residual
         res = res/(nCellsx * nCellsy);
         return res;

       

     

  }
