#pragma once

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
};
}  // namespace heatdiff