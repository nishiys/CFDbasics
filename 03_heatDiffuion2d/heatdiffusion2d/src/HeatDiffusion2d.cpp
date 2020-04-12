#include "HeatDiffusion2d.hpp"

#include <fstream>
#include <iostream>

namespace heatdiff
{
HeatDiffusion2d::HeatDiffusion2d() {}

HeatDiffusion2d::~HeatDiffusion2d() {}

void HeatDiffusion2d::SetMeshConfiguration(std::string meshfilename)
{
    meshfilename_ = meshfilename;
}

void HeatDiffusion2d::GenerateGrid()
{
    su2mesh::SU2meshparser meshparser(meshfilename_);
    meshparser.LoadData();
    meshparser.WriteVtkFile("grid.vtk");

    cells_ = std::move(meshparser.cellarray);
    nodes_ = std::move(meshparser.nodearray);
}
void HeatDiffusion2d::SetFlowConfig(double volumetric_source,
                                    double thermal_cond, double thickness)
{
    for (size_t iCell = 0; iCell < cells_.size(); ++iCell)
    {
        cells_[iCell].cellvar.volumetric_source = volumetric_source;
        for (size_t iFace = 0; iFace < 4; ++iFace)
        {
            cells_[iCell].Face(iFace)->facevar.thermal_cond = thermal_cond;
        }
    }
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
                    cells_[iCell].Face(iFace)->facevar.temperature =
                        temperature;
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
                cells_[iCell].Face(iFace)->facevar.temperature =
                    init_temperature;
            }
        }
    }
}

void HeatDiffusion2d::CalcInitialValues()
{
    for (size_t iCell = 0; iCell < cells_.size(); ++iCell)
    {
        cells_[iCell].cellvar.temperature = 0;
        for (size_t iFace = 0; iFace < 4; ++iFace)
        {
            cells_[iCell].cellvar.temperature +=
                1.0 / 4.0 * cells_[iCell].Face(iFace)->facevar.temperature;
        }
    }
}

void HeatDiffusion2d::InitializeMatrices()
{
    CalcInitialValues();

    int nSize = cells_.size();
    A.resize(nSize, nSize);
    T.resize(nSize);
    B.resize(nSize);

    for (size_t iCell = 0; iCell < cells_.size(); ++iCell)
    {
        T(iCell) = cells_[iCell].cellvar.temperature;
    }
}

void HeatDiffusion2d::Solve() { InitializeMatrices(); }

void HeatDiffusion2d::WriteResultsToVtk(std::string vtkfilename)
{
    std::ofstream outfile;
    outfile.open(vtkfilename);

    std::cout << "Writing to " << vtkfilename << " ..." << std::endl;

    outfile << "# vtk DataFile Version 4.2" << std::endl;
    outfile << "Results on mesh created from " << meshfilename_ << std::endl;
    outfile << "ASCII" << std::endl;
    outfile << "DATASET UNSTRUCTURED_GRID" << std::endl;

    outfile << "POINTS " << nodes_.size() << " "
            << "double " << std::endl;
    for (unsigned int i = 0; i < nodes_.size(); ++i)
    {
        outfile << nodes_[i].GetX() << "\t"  //
                << nodes_[i].GetY() << "\t" << 0 << std::endl;
    }

    outfile << "CELLS " << cells_.size() << " " << 5 * cells_.size()
            << std::endl;
    for (unsigned int i = 0; i < cells_.size(); ++i)
    {
        outfile << "4\t"  // number of point
                << cells_[i].GetNode(0)->GetID() << "\t"
                << cells_[i].GetNode(1)->GetID() << "\t"
                << cells_[i].GetNode(2)->GetID() << "\t"
                << cells_[i].GetNode(3)->GetID() << std::endl;
    }

    outfile << "CELL_TYPES " << cells_.size() << std::endl;
    for (unsigned int i = 0; i < cells_.size(); ++i)
    {
        // VTK_QUAD cell type is defined as 9
        outfile << "9" << std::endl;
    }

    /* Set variables. You can see variables defined here from Paraview */
    outfile << "CELL_DATA " << cells_.size()
            << std::endl;  // Declare the variable below is set in cells
    outfile << "SCALARS ID int" << std::endl;  // give its KIND name type
    outfile << "LOOKUP_TABLE default" << std::endl;
    for (unsigned int i = 0; i < cells_.size(); i++)
    {
        outfile << cells_[i].GetID() << std::endl;
    }
    /* Set variables. You can see variables defined here from Paraview */
    outfile << "SCALARS Temperature float"
            << std::endl;  // give its KIND name type
    outfile << "LOOKUP_TABLE default" << std::endl;
    for (unsigned int i = 0; i < cells_.size(); i++)
    {
        outfile << cells_[i].cellvar.temperature << std::endl;
    }
    std::cout << "File writing finished!" << std::endl;
}

void HeatDiffusion2d::PrintDebug()
{
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    for (size_t iCell = 0; iCell < cells_.size(); ++iCell)
    {
        std::cout << "Cell " << cells_[iCell].GetID() << ": "
                  << cells_[iCell].cellvar.temperature << " [deg]" << std::endl;
        for (size_t iFace = 0; iFace < 4; ++iFace)
        {
            std::cout << "  Face " << iFace << ": "
                      << cells_[iCell].Face(iFace)->facevar.temperature
                      << " [deg]" << std::endl;
        }
    }
}

}  // namespace heatdiff