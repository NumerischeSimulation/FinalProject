#pragma once

#include "settings.h"
#include "discretization/2_central_differences.h"
#include "discretization/2_donor_cell.h"
#include "pressure_solver/sor.h"
#include "pressure_solver/gauss_seidel.h"
#include "output_writer/output_writer_paraview.h"
#include "output_writer/output_writer_text.h"
#include "storage/field_variable.h"

#include <memory>
#include <iostream>  // for cout
#include <cmath>
#include <algorithm>

class Computation
{
public:    
    //! initialize the computation object, parse the settings from file that is given as the only command line argument
    void initialize(int argc, char *argv[]);

    //! run the whole simulation until t_end 
    void runSimulation();

    //! test the simulation, especially the pressure solver
    void runTest();

private:
    //! compute the time step width dt from maximum velocities 
    void computeTimeStepWidth();

    //!  set boundary values of u and v to correct values
    void applyBoundaryValues();

    //! compute the preliminary velocities, F and G
    void computePreliminaryVelocities();

    //! compute the right hand side of the Poisson equation for the pressure
    void computeRightHandSide();

    //! solve the Poisson equation for the pressure 
    void computePressure();

    //! compute the new velocities, u,v, from the preliminary velocities, F,G and the pressure, p
    void computeVelocities();

    Settings settings_;

    std::shared_ptr<Discretization> discretization_;
    std::unique_ptr<PressureSolver> pressureSolver_;

    std::unique_ptr<OutputWriterParaview> outputWriterParaview_;
    std::unique_ptr<OutputWriterText> outputWriterText_;

    std::array<double, 2> meshWidth_; //< h_x, h_y
    double dt_; //< time step
};
