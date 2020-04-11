#include "HeatDiffusion2d.hpp"

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

    cells_ = std::move(meshparser.cellarray_);
}

void HeatDiffusion2d::SetBoundaryConditions(std::string tag,
                                            std::string bc_type,
                                            double temperature)
{
    for (size_t iCell = 0; iCell < cells_.size(); ++iCell)
    {
        if (cells_[iCell].Face1()->GetTag() == tag)
        {
            if (bc_type == "Dirichlet")
            {
                cells_[iCell].var.temperature = temperature;
            }
        }
        if (cells_[iCell].Face2()->GetTag() == tag)
        {
            if (bc_type == "Dirichlet")
            {
                cells_[iCell].var.temperature = temperature;
            }
        }
        if (cells_[iCell].Face3()->GetTag() == tag)
        {
            if (bc_type == "Dirichlet")
            {
                cells_[iCell].var.temperature = temperature;
            }
        }
        if (cells_[iCell].Face4()->GetTag() == tag)
        {
            if (bc_type == "Dirichlet")
            {
                cells_[iCell].var.temperature = temperature;
            }
        }
    }
}
void HeatDiffusion2d::SetInitialConditions(double init_temperature)
{
    for (size_t iCell = 0; iCell < cells_.size(); ++iCell)
    {
        if (cells_[iCell].Face1()->GetTag() == "interior")
        {
            cells_[iCell].var.temperature = init_temperature;
        }
        if (cells_[iCell].Face2()->GetTag() == "interior")
        {
            cells_[iCell].var.temperature = init_temperature;
        }
        if (cells_[iCell].Face3()->GetTag() == "interior")
        {
            cells_[iCell].var.temperature = init_temperature;
        }
        if (cells_[iCell].Face4()->GetTag() == "interior")
        {
            cells_[iCell].var.temperature = init_temperature;
        }
    }
}

void HeatDiffusion2d::Solve() {}

}  // namespace heatdiff