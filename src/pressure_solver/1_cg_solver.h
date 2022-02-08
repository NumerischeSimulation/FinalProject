#pragma once

#include "pressure_solver/0_pressure_solver.h"

#include <vector>
#include <cmath>

class CGSolver: public PressureSolver
{
public:
    CGSolver(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations);

    void solve();
    void calculateInitialResidual();

private:
    FieldVariable res;
};