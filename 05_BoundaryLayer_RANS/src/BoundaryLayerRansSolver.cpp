#include "BoundaryLayerRansSolver.hpp"

#include <cmath>
#include <iostream>

namespace bl_rans
{
BoundaryLayerRansSolver::BoundaryLayerRansSolver(Domain calc_domain,
                                                 unsigned int nx,
                                                 unsigned int ny)
    : calc_domain_(calc_domain), imax_(nx), jmax_(ny)
{
    u.resize(imax_ + 1);
    v.resize(imax_ + 1);
    kinematic_visc_eff.resize(imax_ + 1);
    for (unsigned int i = 0; i < imax_ + 1; ++i)
    {
        u[i].resize(jmax_ + 1);
        v[i].resize(jmax_ + 1);
        kinematic_visc_eff[i].resize(jmax_ + 1);
    }
}

BoundaryLayerRansSolver::~BoundaryLayerRansSolver() {}

void BoundaryLayerRansSolver::SetConfigs(double rho0, double kinematic_visc,
                                         double u0)
{
    rho0_           = rho0;
    u0_             = u0;
    kinematic_visc_ = kinematic_visc;
}

void BoundaryLayerRansSolver::SetInitialCondition()
{
    for (unsigned int j = 0; j < jmax_ + 1; ++j)
    {
        u[0][j] = u0_;
        v[0][j] = 0.0;
    }
}

void BoundaryLayerRansSolver::SetBoundaryCondition()
{
    for (unsigned int i = 0; i < imax_ + 1; ++i)
    {
        u[i][0]     = 0.0;
        v[i][0]     = 0.0;
        u[i][jmax_] = u0_;
    }
}

void BoundaryLayerRansSolver::SolveFlow() {}

void BoundaryLayerRansSolver::WriteResults() {}

}  // namespace bl_rans

int main()
{
    /* Configurations */
    double rho0        = 1.225;
    double temperature = 293.0;
    double kinematic_visc =
        1.458e-6 * std::pow(temperature, 1.5) / (temperature + 110.4);
    double u0   = 20.0;
    double xmax = 1.0;
    double Rex  = rho0 * u0 * xmax / kinematic_visc;
    std::cout << "Reynolds number: " << Rex << std::endl;
    // // Laminar BL analytical solution for a flat plate
    // double laminar_bl_thickness = xmax * 5.3 / std::sqrt(Rex);
    // std::cout << "Laminar Boundary Layer Thickness for Flat Plate: "
    //           << laminar_bl_thickness * 1000 << " [mm]" << std::endl;
    // Turnulent BL analytical solution for a flat plate
    double turb_bl_thickness = xmax * 0.37 / std::pow(Rex, 0.2);
    std::cout << "Turbulent Boundary Layer Thickness for Flat Plate: "
              << turb_bl_thickness * 1000 << " [mm]" << std::endl;

    double ymax        = turb_bl_thickness * 2;
    Bounds x_bounds    = {0.0, xmax};
    Bounds y_bounds    = {0.0, ymax};
    Domain calc_domain = {x_bounds, y_bounds};

    unsigned int nx = 50000;
    // At least 10 points for the thickness direction is necessary.
    unsigned int ny = 200;

    bl_rans::BoundaryLayerRansSolver blSolver(calc_domain, nx, ny);
    blSolver.SetConfigs(rho0, kinematic_visc, u0);
    blSolver.SetInitialCondition();
    blSolver.SetBoundaryCondition();

    blSolver.SolveFlow();
    blSolver.WriteResults();

    return 0;
}