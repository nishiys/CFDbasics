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
    void SetBoundaryConditions();
    void Solve();

private:
    std::string meshfilename_;
};
}  // namespace heatdiff