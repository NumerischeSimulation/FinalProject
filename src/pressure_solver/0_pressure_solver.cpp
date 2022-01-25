#include "0_pressure_solver.h"
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
    setBoundaryValuesTop();
    setBoundaryValuesBottom();

    // prioritise left and right boundaries
    setBoundaryValuesLeft();
    setBoundaryValuesRight();      
  }

  void PressureSolver::setBoundaryValuesLeft()
  {
    // p_-1,j = p_0
    for ( int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++)
    {
      discretization_->p(-1,j) =  discretization_->p(0,j);
    }
  }

  void PressureSolver::setBoundaryValuesRight()
  {
    int n = discretization_->nCells()[0]; // number of cells in the computational domain in x-direction
    // p_n+1,j = p_n,j
    for ( int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++)
    {
      discretization_->p(n,j) =  discretization_->p(n-1 ,j);
    }
  }

  void PressureSolver::setBoundaryValuesTop()
  {
    int m = discretization_->nCells()[1]; // number of cells in the computational domain in y-direction
    // p_-1,j = p_0
    for ( int i = discretization_->pIBegin() +1; i < discretization_->pIEnd() -1; i++)
    {
      discretization_->p(i, m) = discretization_->p(i, m-1);
    }
  }

  void PressureSolver::setBoundaryValuesBottom()
  {
    // p_-1,j = p_0
    for ( int i = discretization_->pIBegin() +1; i < discretization_->pIEnd() -1; i++)
    {
      discretization_->p(i, -1)  = discretization_->p(i, 0);
    }
  }

  double PressureSolver::calculateResidual()
  {
    int nCellsx = discretization_->nCells()[0]; // inner cells in x direction
    int nCellsy = discretization_->nCells()[1]; // inner cells in y direction

    //cell size
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
