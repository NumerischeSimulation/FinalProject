#pragma once

#include "pressure_solver/0_pressure_solver.h"

#include <cmath>

class SOR: public PressureSolver
{
public:
    SOR(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, double omega, bool outflowBottom, bool outflowTop, bool outflowLeft, bool outflowRight);

    void solve();
private:
  double omega_;
};
