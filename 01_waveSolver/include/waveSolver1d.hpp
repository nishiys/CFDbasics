#pragma once

#include <array>
#include <iostream>
#include <string>
#include <vector>

using Vector = std::vector<double>;
using Bounds = std::array<double, 2>;


class waveSolver1d
{
public:
    waveSolver1d(Bounds calc_bounds, double endtime, int n_timestep, int n_cell,
                double wave_velocity);
    ~waveSolver1d();

    void solve();
    
    
    Vector U;
    Vector x;

    inline int get_n_cells() const {return n_cells_;};
    inline double get_dx() const {return dx_;};;


private:
    void writeFlowFile(std::string filename);
    int n_cells_;
    int n_timestep_;
    double dx_;
    Bounds calc_bounds_;
    double endtime_;
    double dt_;
    double c_;
    double CFL_;
};
