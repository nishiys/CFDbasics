#pragma once

#include <Eigen/Core>
#include <array>
#include <iostream>
#include <string>
#include <vector>

using Vector = std::vector<double>;
using Bounds = std::array<double, 2>;

class ShockTubeSolver
{
public:
    ShockTubeSolver(Bounds calc_bounds, double endtime, int n_timestep,
                    int n_cell);
    ~ShockTubeSolver();

    void SetGrid();
    void SetFlowField(double wall_position, double rhoL, double pL, double rhoR,
                      double pR);
    void Solve();

    std::vector<double> cells;
    std::vector<Eigen::VectorXd> flowVars;
    Eigen::VectorXd primitives;
    Eigen::VectorXd conservatives;
    Eigen::VectorXd PrimToCons(Eigen::VectorXd prim);
    Eigen::VectorXd ConsToPrim(Eigen::VectorXd cons);

    Eigen::VectorXd Roe(Eigen::VectorXd primL, Eigen::VectorXd primR);
    Eigen::VectorXd vanLeer(Eigen::VectorXd primL, Eigen::VectorXd primR);

    inline int Get_n_cells() const { return n_cells_; };
    inline double Get_dx() const { return dx_; };

private:
    void WriteFlowFile(std::string filename);
    int n_cells_;
    int n_timestep_;
    double dx_;
    Bounds calc_bounds_;
    double endtime_;
    double dt_;

    const int N_EQ     = 3;
    const double GAMMA = 1.40;
};
