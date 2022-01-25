#pragma once

#include "pressure_solver/0_pressure_solver.h"

class GaussSeidel: public PressureSolver{
public:
    using PressureSolver::PressureSolver;

    void solve();
};