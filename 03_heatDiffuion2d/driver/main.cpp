#include "HeatDiffusion2d.hpp"

int main()
{
    /*--- Configuration ---*/
    // Heat Source Per Unit Volume (W/m3)
    double volumetric_source = 1000;
    // Thermal Conducivity (W/mK)
    double thermal_cond = 100;
    // Thickness of the plate (m)
    double thickness = 0.1;

    // BCs
    double top_temp    = 250;
    double left_temp   = 100;
    double bottom_temp = 150;
    double right_temp  = 200;
    // ICs.
    double init_temp = 100;

    /*--- Solve ---*/
    heatdiff::HeatDiffusion2d heat2dsolver;
    heat2dsolver.SetMeshConfiguration("quad.su2");
    heat2dsolver.GenerateGrid();
    heat2dsolver.SetFlowConfig(volumetric_source, thermal_cond, thickness);
    heat2dsolver.SetBoundaryConditions("top", "Dirichlet", top_temp);
    heat2dsolver.SetBoundaryConditions("left", "Dirichlet", left_temp);
    heat2dsolver.SetBoundaryConditions("bottom", "Dirichlet", bottom_temp);
    heat2dsolver.SetBoundaryConditions("right", "Dirichlet", right_temp);
    heat2dsolver.SetInitialConditions(init_temp);

    heat2dsolver.Solve();
    heat2dsolver.WriteResultsToVtk("quad_result.vtk");
    heat2dsolver.PrintDebug();

    return 0;
}