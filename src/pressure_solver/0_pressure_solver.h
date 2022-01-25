#pragma once

#include <memory>

#include "discretization/1_discretization.h"

class PressureSolver {
public:
    PressureSolver(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations );
    virtual void solve() = 0;

    



protected:
    // set boundary values
    void setBoundaryValues();
    void setBoundaryValuesLeft();
    void setBoundaryValuesRight();
    void setBoundaryValuesTop();
    void setBoundaryValuesBottom();

    void applyObstacleBoundaryValues();

    double calculateResidual();

protected:
    std::shared_ptr<Discretization> discretization_;
    double epsilon_;
    int maximumNumberOfIterations_;
};