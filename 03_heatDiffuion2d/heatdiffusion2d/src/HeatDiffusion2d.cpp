#include "HeatDiffusion2d.hpp"

#include <iostream>

namespace heatdiff
{
HeatDiffusion2d::HeatDiffusion2d() {}

HeatDiffusion2d::~HeatDiffusion2d() {}

void HeatDiffusion2d::SetConfiguration(std::string meshfilename)
{
    meshfilename_ = meshfilename;
}

void HeatDiffusion2d::GenerateGrid()
{
    su2mesh::SU2meshparser meshparser(meshfilename_);
    meshparser.LoadData();
    meshparser.WriteVtkFile("grid.vtk");

    cells_ = std::move(meshparser.cellarray);
}

void HeatDiffusion2d::SetBoundaryConditions(std::string tag,
                                            std::string bc_type,
                                            double temperature)
{
    for (size_t iCell = 0; iCell < cells_.size(); ++iCell)
    {
        for (size_t iFace = 0; iFace < 4; ++iFace)
        {
            if (cells_[iCell].Face(iFace)->GetTag() == tag)
            {
                if (bc_type == "Dirichlet")
                {
                    cells_[iCell].var.temperature = temperature;
                }
            }
        }
    }
}
void HeatDiffusion2d::SetInitialConditions(double init_temperature)
{
    for (size_t iCell = 0; iCell < cells_.size(); ++iCell)
    {
        for (size_t iFace = 0; iFace < 4; ++iFace)
        {
            if (cells_[iCell].Face(iFace)->GetTag() == "interior")
            {
                cells_[iCell].var.temperature = init_temperature;
            }
        }
    }
}

void HeatDiffusion2d::Solve() {}

void HeatDiffusion2d::PrintDebug()
{
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    for (size_t iCell = 0; iCell < cells_.size(); ++iCell)
    {
        std::cout << "Cell " << cells_[iCell].GetID() << ": "  //
                  << cells_[iCell].var.temperature << " [deg]" << std::endl;
    }
}

}  // namespace heatdiff