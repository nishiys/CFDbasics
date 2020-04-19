#pragma once

#include <array>
#include <vector>

using Bounds = std::array<double, 2>;
using Domain = std::array<Bounds, 2>;

namespace bl_rans
{
class BoundaryLayerRansSolver
{
public:
    BoundaryLayerRansSolver(Domain calc_domain, unsigned int nx,
                            unsigned int ny);
    ~BoundaryLayerRansSolver();

    void SetConfigs(double rho0, double kinematic_viisc, double u0);
    void SetInitialCondition();
    void SetBoundaryCondition();
    void SolveFlow();
    void WriteResults();

private:
    Domain calc_domain_;
    double rho0_;
    double kinematic_visc_;
    double u0_;
    unsigned int imax_;
    unsigned int jmax_;

    std::vector<std::vector<double>> u;
    std::vector<std::vector<double>> v;
    std::vector<std::vector<double>> kinematic_visc_eff;
};

}  // namespace bl_rans
