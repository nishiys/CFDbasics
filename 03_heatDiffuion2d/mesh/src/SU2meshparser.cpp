#include "SU2meshparser.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>

namespace su2mesh
{
SU2meshparser::SU2meshparser(std::string meshfilename)
    : meshfilename_(meshfilename)
{
}
SU2meshparser::~SU2meshparser() {}

void SU2meshparser::LoadData()
{
    ReadFile();
    CreateQuadArray();
    SetMarkersToFaces();
    SetNeighborCells();
    SetVectorsToNeighbors();
    SetN1Vectors();
    PrintDebug();
}

void SU2meshparser::ReadFile()
{
    std::ifstream infile;
    infile.open(meshfilename_, std::ios::in);
    if (!infile)
    {
        std::cerr << "File was not opend." << std::endl;
        std::exit(1);
    }
    else
    {
        std::cout << "File " << meshfilename_ << " was successfully opened."
                  << std::endl;
    }

    std::string line;

    // Get DIM_
    std::getline(infile, line);
    if (line.find("NDIME") != std::string::npos)
    {
        std::stringstream ss(line);
        // skip until '='
        ss.ignore(line.size(), '=');
        ss >> nDim_;
        std::cout << "Dimension: " << nDim_ << std::endl;
    }

    // Get NElement_
    std::getline(infile, line);
    if (line.find("NELEM") != std::string::npos)
    {
        std::stringstream ss(line);
        // skip until '='
        ss.ignore(line.size(), '=');
        ss >> nElement_;
        std::cout << "Number of Elements: " << nElement_ << std::endl;
    }

    // Get Element table
    element_table_.resize(nElement_);
    for (size_t i = 0; i < nElement_; ++i)
    {
        element_table_[i].resize(5);
    }

    for (size_t i = 0; i < nElement_; ++i)
    {
        std::getline(infile, line);
        std::stringstream ss(line);
        unsigned int celltype;
        ss >> celltype;
        if (celltype == VTK_QUAD4)
        {
            ss >> element_table_[i][1]    // face1
                >> element_table_[i][2]   // face2
                >> element_table_[i][3]   // face3
                >> element_table_[i][4]   // face4
                >> element_table_[i][0];  // id
        }
    }

    // Get NPoint_
    std::getline(infile, line);
    if (line.find("NPOIN") != std::string::npos)
    {
        std::stringstream ss(line);
        // skip until '='
        ss.ignore(line.size(), '=');
        ss >> nPoint_;
        std::cout << "Number of Points: " << nPoint_ << std::endl;
    }

    // Get nodearray
    nodearray.resize(nPoint_);
    for (size_t i = 0; i < nPoint_; ++i)
    {
        std::getline(infile, line);
        std::stringstream ss(line);
        double x;
        double y;
        unsigned int id;
        ss >> x >> y >> id;
        Node2d node_obj(id, x, y);
        nodearray[i] = node_obj;
    }

    // Get Markers
    std::getline(infile, line);
    if (line.find("NMARK") != std::string::npos)
    {
        std::stringstream ss(line);
        // skip until '='
        ss.ignore(line.size(), '=');
        ss >> nMarker_;
        std::cout << "Number of Markers: " << nMarker_ << std::endl;
    }

    for (size_t i = 0; i < nMarker_; ++i)
    {
        // Get tag
        std::string tag;
        std::getline(infile, line);
        if (line.find("MARKER_TAG") != std::string::npos)
        {
            std::stringstream ss(line);
            // skip until '='
            ss.ignore(line.size(), '=');
            ss >> tag;
        }
        // Get No. of elements which have this tag
        unsigned int nElem;
        std::getline(infile, line);
        if (line.find("MARKER_ELEMS") != std::string::npos)
        {
            std::stringstream ss(line);
            // skip until '='
            ss.ignore(line.size(), '=');
            ss >> nElem;
        }
        // Push back to marker table
        for (size_t i = 0; i < nElem; ++i)
        {
            std::getline(infile, line);
            std::stringstream ss(line);
            unsigned int vtk_type;
            ss >> vtk_type;
            if (vtk_type == VTK_LINE)
            {
                marker_table_.marker_array.push_back(tag);
                std::array<unsigned int, 2> edge;
                ss >> edge[0] >> edge[1];
                marker_table_.edge_array.push_back(edge);
            }
        }
    }

    // Skip unnecessary lines
    while (!infile.eof())
    {
        std::getline(infile, line);
    }

    std::cout << "------------------------------" << std::endl;
}

void SU2meshparser::CreateQuadArray()
{
    cellarray.resize(nElement_);
    for (size_t i = 0; i < element_table_.size(); ++i)
    {
        unsigned int id = element_table_[i][0];
        Node2d& node1   = nodearray[element_table_[i][1]];
        Node2d& node2   = nodearray[element_table_[i][2]];
        Node2d& node3   = nodearray[element_table_[i][3]];
        Node2d& node4   = nodearray[element_table_[i][4]];

        CellQuad4 quad(id, &node1, &node2, &node3, &node4);
        cellarray[i] = quad;
    }
}

void SU2meshparser::SetMarkersToFaces()
{
    /*--- Set read tags to boundary faces ---*/
    for (size_t iTable = 0; iTable < marker_table_.marker_array.size();
         ++iTable)
    {
        std::string mark = marker_table_.marker_array[iTable];
        std::array<unsigned int, 2> markered_edge =
            marker_table_.edge_array[iTable];
        for (size_t iCell = 0; iCell < cellarray.size(); ++iCell)
        {
            // Get nodes IDs that consist each cell edge
            std::array<std::array<unsigned int, 2>, 4> edges;
            edges[0][0] = cellarray[iCell].GetNode(0)->GetID();
            edges[0][1] = cellarray[iCell].GetNode(1)->GetID();
            edges[1][0] = cellarray[iCell].GetNode(1)->GetID();
            edges[1][1] = cellarray[iCell].GetNode(2)->GetID();
            edges[2][0] = cellarray[iCell].GetNode(2)->GetID();
            edges[2][1] = cellarray[iCell].GetNode(3)->GetID();
            edges[3][0] = cellarray[iCell].GetNode(3)->GetID();
            edges[3][1] = cellarray[iCell].GetNode(0)->GetID();

            // check if cell edges are the same as markered edges
            std::array<bool, 4> edge_flags;
            for (size_t iEdge = 0; iEdge < edges.size(); ++iEdge)
            {
                edge_flags[iEdge] = (edges[iEdge] == markered_edge);
            }

            // Set tags
            for (size_t iEdge = 0; iEdge < edges.size(); ++iEdge)
            {
                if (edge_flags[iEdge])
                {
                    cellarray[iCell].Face(iEdge)->SetTag(mark);
                }
            }
        }
    }

    /*--- Set interior tags to interior faces ---*/
    for (size_t iCell = 0; iCell < cellarray.size(); ++iCell)
    {
        for (size_t iFace = 0; iFace < 4; ++iFace)
        {
            if (cellarray[iCell].Face(iFace)->GetTag().empty())
            {
                cellarray[iCell].Face(iFace)->SetTag("interior");
            }
        }
    }
}

void SU2meshparser::SetNeighborCells()
{
    for (size_t iCell = 0; iCell < cellarray.size(); ++iCell)
    {
        // Search over the element table & Get candidate neighbor cells
        std::vector<unsigned int> candidateNeighborCellsIDs;
        unsigned int node1ID = cellarray[iCell].GetNode(0)->GetID();
        unsigned int node2ID = cellarray[iCell].GetNode(1)->GetID();
        unsigned int node3ID = cellarray[iCell].GetNode(2)->GetID();
        unsigned int node4ID = cellarray[iCell].GetNode(3)->GetID();

        std::cout << node1ID << " " << node2ID << " " << node3ID << " "
                  << node4ID << std::endl;

        for (size_t jCell = 0; jCell < element_table_.size(); ++jCell)
        {
            // add candidate neighbors IDs
            if (iCell != jCell)
            {
                if ((element_table_[jCell][1] == node1ID) ||
                    (element_table_[jCell][1] == node2ID) ||
                    (element_table_[jCell][1] == node3ID) ||
                    (element_table_[jCell][1] == node4ID) ||
                    (element_table_[jCell][2] == node1ID) ||
                    (element_table_[jCell][2] == node2ID) ||
                    (element_table_[jCell][2] == node3ID) ||
                    (element_table_[jCell][2] == node4ID) ||
                    (element_table_[jCell][3] == node1ID) ||
                    (element_table_[jCell][3] == node2ID) ||
                    (element_table_[jCell][3] == node3ID) ||
                    (element_table_[jCell][3] == node4ID) ||
                    (element_table_[jCell][4] == node1ID) ||
                    (element_table_[jCell][4] == node2ID) ||
                    (element_table_[jCell][4] == node3ID) ||
                    (element_table_[jCell][4] == node4ID))
                {
                    candidateNeighborCellsIDs.push_back(jCell);
                }
            }
        }

        // Loop over candidate neighbors & Set pointers to neighbors
        std::array<bool, 4> isInterior;
        for (size_t iFace = 0; iFace < isInterior.size(); ++iFace)
        {
            isInterior[iFace] =
                (cellarray[iCell].Face(iFace)->GetTag() == "interior");
        }

        // Get nodes IDs that consist each cell edge
        std::array<std::array<unsigned int, 2>, 4> edges;
        edges[0][0] = cellarray[iCell].GetNode(0)->GetID();
        edges[0][1] = cellarray[iCell].GetNode(1)->GetID();
        edges[1][0] = cellarray[iCell].GetNode(1)->GetID();
        edges[1][1] = cellarray[iCell].GetNode(2)->GetID();
        edges[2][0] = cellarray[iCell].GetNode(2)->GetID();
        edges[2][1] = cellarray[iCell].GetNode(3)->GetID();
        edges[3][0] = cellarray[iCell].GetNode(3)->GetID();
        edges[3][1] = cellarray[iCell].GetNode(0)->GetID();

        unsigned int nNeighbors = candidateNeighborCellsIDs.size();
        for (size_t iFace = 0; iFace < isInterior.size(); ++iFace)
        {
            if (isInterior[iFace])
            {
                for (size_t i = 0; i < nNeighbors; ++i)
                {
                    std::array<unsigned int, 2> candidate_edge1;
                    candidate_edge1[1] = cellarray[candidateNeighborCellsIDs[i]]
                                             .GetNode(0)
                                             ->GetID();
                    candidate_edge1[0] = cellarray[candidateNeighborCellsIDs[i]]
                                             .GetNode(1)
                                             ->GetID();
                    std::array<unsigned int, 2> candidate_edge2;
                    candidate_edge2[1] = cellarray[candidateNeighborCellsIDs[i]]
                                             .GetNode(1)
                                             ->GetID();
                    candidate_edge2[0] = cellarray[candidateNeighborCellsIDs[i]]
                                             .GetNode(2)
                                             ->GetID();
                    std::array<unsigned int, 2> candidate_edge3;
                    candidate_edge3[1] = cellarray[candidateNeighborCellsIDs[i]]
                                             .GetNode(2)
                                             ->GetID();
                    candidate_edge3[0] = cellarray[candidateNeighborCellsIDs[i]]
                                             .GetNode(3)
                                             ->GetID();
                    std::array<unsigned int, 2> candidate_edge4;
                    candidate_edge4[1] = cellarray[candidateNeighborCellsIDs[i]]
                                             .GetNode(3)
                                             ->GetID();
                    candidate_edge4[0] = cellarray[candidateNeighborCellsIDs[i]]
                                             .GetNode(0)
                                             ->GetID();

                    if (edges[iFace] == candidate_edge1 ||
                        edges[iFace] == candidate_edge2 ||
                        edges[iFace] == candidate_edge3 ||
                        edges[iFace] == candidate_edge4)
                    {
                        cellarray[iCell].SetNeighborPtr(
                            iFace, &cellarray[candidateNeighborCellsIDs[i]]);
                    }
                }
            }
        }
    }
}

void SU2meshparser::SetVectorsToNeighbors()
{
    for (size_t iCell = 0; iCell < cellarray.size(); ++iCell)
    {
        for (size_t iFace = 0; iFace < 4; ++iFace)
        {
            cellarray[iCell].SetVectorToNeighbor(iFace);
        }
    }
}

void SU2meshparser::SetN1Vectors()
{
    for (size_t iCell = 0; iCell < cellarray.size(); ++iCell)
    {
        for (size_t iFace = 0; iFace < 4; ++iFace)
        {
            cellarray[iCell].SetMag_N1Vectors(iFace);
        }
    }
}
void SU2meshparser::PrintDebug()
{
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

    std::cout << "Element table: " << std::endl;
    for (size_t i = 0; i < element_table_.size(); ++i)
    {
        std::cout << element_table_[i][0] << "\t"  //
                  << element_table_[i][1] << "\t"  //
                  << element_table_[i][2] << "\t"  //
                  << element_table_[i][3] << "\t"  //
                  << element_table_[i][4] << std::endl;
    }

    std::cout << "\nNode array: " << std::endl;
    for (size_t i = 0; i < nodearray.size(); ++i)
    {
        std::cout << nodearray[i].GetID() << "\t"  //
                  << nodearray[i].GetX() << "\t"   //
                  << nodearray[i].GetY() << std::endl;
    }

    std::cout << "\nCell array: " << std::endl;
    for (size_t i = 0; i < cellarray.size(); ++i)
    {
        std::cout << cellarray[i].GetID() << "\t"
                  << cellarray[i].GetNode(0)->GetID() << "\t"
                  << cellarray[i].GetNode(1)->GetID() << "\t"
                  << cellarray[i].GetNode(2)->GetID() << "\t"
                  << cellarray[i].GetNode(3)->GetID() << std::endl;
    }
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

    for (size_t iCell = 0; iCell < cellarray.size(); ++iCell)
    {
        for (size_t iFace = 0; iFace < 4; ++iFace)
        {
            if (!cellarray[iCell].Face(iFace)->GetTag().empty())
            {
                std::cout << "Cell " << cellarray[iCell].GetID()
                          << " has a face w/ tag "
                          << cellarray[iCell].Face(iFace)->GetTag()
                          << std::endl;
            }
        }
    }
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

    for (size_t iCell = 0; iCell < cellarray.size(); ++iCell)
    {
        std::array<CellQuad4*, 4> neighborPtrs;
        for (size_t iFace = 0; iFace < 4; ++iFace)
        {
            neighborPtrs[iFace] = cellarray[iCell].GetNeighborPtr(iFace);
            if (neighborPtrs[iFace] != nullptr)
            {
                std::cout << "Cell " << cellarray[iCell].GetID()
                          << " is ajacent to "
                          << cellarray[iCell].GetNeighborPtr(iFace)->GetID()
                          << std::endl;
            }
        }
    }
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

    for (size_t iCell = 0; iCell < cellarray.size(); ++iCell)
    {
        std::cout << "Cell " << iCell << ":" << std::endl;
        for (size_t iFace = 0; iFace < 4; ++iFace)
        {
            Eigen::IOFormat CleanFmt(3, 0, ", ", "\n", "[", "]");
            if (cellarray[iCell].GetNeighborPtr(iFace) != nullptr)
            {
                std::cout << "  vector to neighbor " << iFace << ": "
                          << cellarray[iCell]
                                 .GetVectorToNeighbor(iFace)
                                 .transpose()
                                 .format(CleanFmt)
                          << std::endl;
            }
            else
            {
                std::cout << "  normal vector to boundary " << iFace << ": "
                          << cellarray[iCell]
                                 .GetNormalVectorToBoundary(iFace)
                                 .transpose()
                                 .format(CleanFmt)
                          << std::endl;
            }
        }
    }
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
}

void SU2meshparser::WriteVtkFile(std::string vtkfilename)
{
    std::ofstream outfile;
    outfile.open(vtkfilename);

    std::cout << "Writing to " << vtkfilename << " ..." << std::endl;

    outfile << "# vtk DataFile Version 4.2" << std::endl;
    outfile << "Created from " << meshfilename_ << std::endl;
    outfile << "ASCII" << std::endl;
    outfile << "DATASET UNSTRUCTURED_GRID" << std::endl;

    outfile << "POINTS " << nPoint_ << " "
            << "double " << std::endl;
    for (unsigned int i = 0; i < nPoint_; ++i)
    {
        outfile << nodearray[i].GetX() << "\t"  //
                << nodearray[i].GetY() << "\t" << 0 << std::endl;
    }

    outfile << "CELLS " << nElement_ << " " << 5 * nElement_ << std::endl;
    for (unsigned int i = 0; i < nElement_; ++i)
    {
        outfile << "4\t"  // number of point
                << cellarray[i].GetNode(0)->GetID() << "\t"
                << cellarray[i].GetNode(1)->GetID() << "\t"
                << cellarray[i].GetNode(2)->GetID() << "\t"
                << cellarray[i].GetNode(3)->GetID() << std::endl;
    }

    outfile << "CELL_TYPES " << nElement_ << std::endl;
    for (unsigned int i = 0; i < nElement_; ++i)
    {
        // VTK_QUAD cell type is defined as 9
        outfile << "9" << std::endl;
    }

    /* Set variables. You can see variables defined here from Paraview */
    outfile << "CELL_DATA " << nElement_
            << std::endl;  // Declare the variable below is set in cells
    outfile << "SCALARS ID int" << std::endl;  // give its KIND name type
    outfile << "LOOKUP_TABLE default" << std::endl;
    for (unsigned int i = 0; i < nElement_; i++)
    {
        outfile << cellarray[i].GetID() << std::endl;
    }

    std::cout << "File writing finished!" << std::endl;
}

}  // namespace su2mesh
