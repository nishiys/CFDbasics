#include "HeatDiffusion2d.hpp"

int main()
{
    heatdiff::HeatDiffusion2d heat2dsolver;
    heat2dsolver.SetConfiguration("quad.su2");
    heat2dsolver.GenerateGrid();
    heat2dsolver.SetBoundaryConditions();
    heat2dsolver.Solve();
    return 0;
}