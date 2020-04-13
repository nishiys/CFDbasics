#include "HeatDiffusion2d.hpp"

#include <fstream>
#include <iostream>
#include <vector>

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
    thickness_ = thickness;
    for (size_t iCell = 0; iCell < cells_.size(); ++iCell)
    {
        cells_[iCell].cellvar.volumetric_source = volumetric_source;
        for (size_t iFace = 0; iFace < 4; ++iFace)
        {
            cells_[iCell].Face(iFace)->facevar.thermal_cond = thermal_cond;
        }
    }
}

void HeatDiffusion2d::SetBoundaryConditions(std::string marker,
                                            std::string bc_type,
                                            double temperature)
{
    // Set face & nodal values
    for (size_t iCell = 0; iCell < cells_.size(); ++iCell)
    {
        for (size_t iFace = 0; iFace < 4; ++iFace)
        {
            if (cells_[iCell].Face(iFace)->GetMarker() == marker)
            {
                if (bc_type == "Dirichlet")
                {
                    // Set face values
                    cells_[iCell].Face(iFace)->facevar.temperature =
                        temperature;
                    cells_[iCell].Face(iFace)->SetBcType(bc_type);
                    // Set nodal values
                    cells_[iCell].Face(iFace)->GetNode1()->nodevar.temperature =
                        temperature;
                    cells_[iCell].Face(iFace)->GetNode2()->nodevar.temperature =
                        temperature;
                }
                else if (bc_type == "Neumann")
                {
                    //TODO
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
            if (cells_[iCell].Face(iFace)->GetMarker() == "interior")
            {
                cells_[iCell].Face(iFace)->facevar.temperature =
                    init_temperature;
            }
            else if (cells_[iCell].Face(iFace)->GetBcType() == "Neumann")
            {
                // Set nodal values
                cells_[iCell].Face(iFace)->GetNode1()->nodevar.temperature =
                    init_temperature;
                cells_[iCell].Face(iFace)->GetNode2()->nodevar.temperature =
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

void HeatDiffusion2d::ConstructMatrices()
{
    CalcInitialValues();

    PrintDebug();

    int nSize = cells_.size();
    A.resize(nSize, nSize);
    T.resize(nSize);
    B.resize(nSize);

    /*--- B vector (Source terms) ---*/
    for (size_t iCell = 0; iCell < cells_.size(); ++iCell)
    {
        // volumetric contribution
        B(iCell) = cells_[iCell].cellvar.volumetric_source *
                   cells_[iCell].GetVolume() * thickness_;

        // boundary contribution
        for (size_t iFace = 0; iFace < 4; ++iFace)
        {
            if (cells_[iCell].Face(iFace)->GetBcType() == "Dirichlet")
            {
                double Twall = cells_[iCell].Face(iFace)->facevar.temperature;
                double source_from_dirichlet =
                    Twall * cells_[iCell].Face(iFace)->facevar.thermal_cond *
                    cells_[iCell].Face(iFace)->GetArea() /
                    cells_[iCell].GetNormalVectorToBoundary(iFace).norm();
                B(iCell) += source_from_dirichlet * thickness_;
            }
            else if (cells_[iCell].Face(iFace)->GetBcType() == "Neumann")
            {
                //TODO
            }
            else  // interior faces
            {
            }
        }
    }
    std::cout << "B vector:\n" << B << std::endl;

    /*--- A matrix ---*/
    // For eigen sparse matrix
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletvec;
    for (size_t iCell = 0; iCell < cells_.size(); ++iCell)
    {
        for (size_t iFace = 0; iFace < 4; ++iFace)
        {
            if (cells_[iCell].Face(iFace)->GetMarker() == "interior")
            {
                unsigned int neighbor_id =
                    cells_[iCell].GetNeighborPtr(iFace)->GetID();
                double ap = cells_[iCell].Face(iFace)->facevar.thermal_cond *
                            cells_[iCell].Face(iFace)->GetArea() * thickness_ *
                            cells_[iCell].GetMag_N1Vectors(iFace) /
                            cells_[iCell].GetVectorToNeighbor(iFace).norm();
                tripletvec.push_back(T(iCell, iCell, ap));
                tripletvec.push_back(T(iCell, neighbor_id, -ap));
            }
            else if (cells_[iCell].Face(iFace)->GetBcType() == "Dirichlet")
            {
                double ap =
                    cells_[iCell].Face(iFace)->facevar.thermal_cond *
                    cells_[iCell].Face(iFace)->GetArea() * thickness_ /
                    cells_[iCell].GetNormalVectorToBoundary(iFace).norm();
                tripletvec.push_back(T(iCell, iCell, ap));
            }
        }
    }
    A.setFromTriplets(tripletvec.begin(), tripletvec.end());
    std::cout << "A matrix: \n" << A << std::endl;
}
void HeatDiffusion2d::Solve()
{
    /*--- Set LSE ---*/
    ConstructMatrices();

    /*--- Solve LSE ---*/
    Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >
        solver;
    solver.compute(A);
    if (solver.info() != Eigen::Success)
    {
        std::cerr << "decomposition failed" << std::endl;
    }
    T = solver.solve(B);
    if (solver.info() != Eigen::Success)
    {
        std::cerr << "solving failed" << std::endl;
    }
    std::cout << "Temperature Vector: \n" << T << std::endl;

    /*--- Update variables ---*/
    for (size_t iCell = 0; iCell < cells_.size(); ++iCell)
    {
        // Cell value
        cells_[iCell].cellvar.temperature = T(iCell);
    }

    /*--- Calculate Nodal values ---*/
    CalcNodalValues();

    // Face value
    UpdateFaceValues();
}

void HeatDiffusion2d::CalcNodalValues()
{
    /*--- Calculate nodal value w/ weighted averaging ---*/
    for (size_t iNode = 0; iNode < nodes_.size(); ++iNode)
    {
        if (nodes_[iNode].IsInterior() == true)
        {
            double numerator   = 0;
            double denominator = 0;

            for (size_t i = 0; i < nodes_[iNode].GetCellPtrs().size(); ++i)
            {
                // search nodeID in this cell that is the same as iNode
                std::array<unsigned int, 4> nodeIDs_in_current_cell;
                nodeIDs_in_current_cell[0] =
                    nodes_[iNode].GetCellPtrs()[i]->GetNode(0)->GetID();
                nodeIDs_in_current_cell[1] =
                    nodes_[iNode].GetCellPtrs()[i]->GetNode(1)->GetID();
                nodeIDs_in_current_cell[2] =
                    nodes_[iNode].GetCellPtrs()[i]->GetNode(2)->GetID();
                nodeIDs_in_current_cell[3] =
                    nodes_[iNode].GetCellPtrs()[i]->GetNode(3)->GetID();

                unsigned int flag;
                for (unsigned int index = 0; index < 4; ++index)
                {
                    if (nodeIDs_in_current_cell[index] == iNode) flag = index;
                }

                numerator +=
                    nodes_[iNode].GetCellPtrs()[i]->cellvar.temperature /
                    nodes_[iNode]
                        .GetCellPtrs()[i]
                        ->GetVectorToNode(flag)
                        .norm();
                denominator += 1.0 / nodes_[iNode]
                                         .GetCellPtrs()[i]
                                         ->GetVectorToNode(flag)
                                         .norm();
            }

            nodes_[iNode].nodevar.temperature = numerator / denominator;
        }
        else
        {
            //TODO: if neumann bc exists, you have to add here
        }
    }
}

void HeatDiffusion2d::UpdateFaceValues()
{ /* w/ Nodal Avaraging scheme */
    for (size_t iCell = 0; iCell < cells_.size(); ++iCell)
    {
        for (size_t iFace = 0; iFace < 4; ++iFace)
        {
            cells_[iCell].Face(iFace)->facevar.temperature =
                (cells_[iCell].Face(iFace)->GetNode1()->nodevar.temperature +
                 cells_[iCell].Face(iFace)->GetNode2()->nodevar.temperature) /
                2.0;
        }
    }
}

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
        outfile << nodes_[i].GetCoords()(0) << "\t"  //
                << nodes_[i].GetCoords()(1) << "\t" << 0 << std::endl;
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

    outfile << "POINT_DATA " << nodes_.size() << std::endl;
    outfile << "SCALARS Temperature float" << std::endl;
    outfile << "LOOKUP_TABLE default" << std::endl;
    for (unsigned int i = 0; i < nodes_.size(); i++)
    {
        outfile << nodes_[i].nodevar.temperature << std::endl;
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
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
}

}  // namespace heatdiff