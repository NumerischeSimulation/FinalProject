#pragma once

#include <memory>

#include "discretization/1_discretization.h"

class PressureSolver {
public:
    PressureSolver(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, bool outflowBottom, bool outflowTop, bool outflowLeft, bool outflowRight);
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
    bool outflowBottom_;  //> prescribed state of the outflow at bottom of domain
    bool outflowTop_;  //> prescribed state of the outflow at top of domain
    bool outflowLeft_;  //> prescribed state of the outflow at left of domain
    bool outflowRight_;  //> prescribed state of the outflow at right of domain
};