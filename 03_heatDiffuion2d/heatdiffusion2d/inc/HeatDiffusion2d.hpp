#pragma once

// #include <Eigen/Core>
#include <Eigen/Sparse>

#include "SU2meshparser.hpp"

namespace heatdiff
{
class HeatDiffusion2d
{
public:
    HeatDiffusion2d();
    ~HeatDiffusion2d();

    void SetConfiguration(std::string meshfilename);
    void GenerateGrid();
    void SetBoundaryConditions(std::string tag, std::string bc_type,
                               double temperature);
    void SetInitialConditions(double init_temperature);
    void Solve();

    // For debug
    void PrintDebug();

private:
    std::string meshfilename_;
    std::vector<CellQuad4> cells_;

    // Matrix initialization
    void InitializeMatrices();

    // For Simultaneous Linear Equation
    //! Coefficient matrix
    Eigen::SparseMatrix<double> A;
    //! Temperature vector
    Eigen::VectorXd T;
    //! Source vector
    Eigen::VectorXd B;
};
}  // namespace heatdiff