#pragma once

#include <Eigen/Core>
#include <array>
#include <iostream>
#include <string>
#include <vector>

using Bounds = std::array<double, 2>;

class ShockTubeSolver
{
public:
    ShockTubeSolver(Bounds calc_bounds, double endtime, double dt, int n_cell);
    ~ShockTubeSolver();

    void SetGrid();
    void SetFlowField(double wall_position, double rhoL, double pL, double rhoR,
                      double pR);
    void SetSchemes(std::string convectiveflux_scheme, bool MUSCL, int k);
    void Solve();
    void WriteFlowFile(std::string filename);


    inline int Get_n_cells() const { return n_cells_; };
    inline double Get_dx() const { return dx_; };

private:
    std::vector<double> cells;
    std::vector<Eigen::Vector3d> flowVars;
    std::vector<Eigen::Vector3d> flowConservatives;

    Eigen::Vector3d PrimToCons(Eigen::Vector3d prim);
    Eigen::Vector3d ConsToPrim(Eigen::Vector3d cons);

    // Euler flux
    Eigen::Vector3d Roe(Eigen::Vector3d primL, Eigen::Vector3d primR);
    Eigen::Vector3d vanLeer(Eigen::Vector3d primL, Eigen::Vector3d primR);
    Eigen::Vector3d AUSM(Eigen::Vector3d primL, Eigen::Vector3d primR);

    // Limiter
    double minmod(double a, double b);


    int n_cells_;
    int n_timestep_;
    double dx_;
    Bounds calc_bounds_;
    double endtime_;
    double dt_;
    double CFL_;

    const double GAMMA = 1.40;

    std::string convectiveflux_scheme_;
    bool MUSCL_;
    //! A coefficient of MUSCL scheme
    //! k = -1: 2nd order upwind differencing
    //! k = 0: 2nd order upwind-biased differencing
    //! k = 1/3: 3rd order upwind-biased differencing
    double k_;
};
