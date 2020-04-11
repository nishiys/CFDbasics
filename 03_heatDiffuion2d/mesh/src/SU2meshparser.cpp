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
        ss >> DIM_;
        std::cout << "Dimension: " << DIM_ << std::endl;
    }

    // Get NElement_
    std::getline(infile, line);
    if (line.find("NELEM") != std::string::npos)
    {
        std::stringstream ss(line);
        // skip until '='
        ss.ignore(line.size(), '=');
        ss >> NElement_;
        std::cout << "Number of Elements: " << NElement_ << std::endl;
    }

    // Get Element table
    element_table_.resize(NElement_);
    for (size_t i = 0; i < NElement_; ++i)
    {
        element_table_[i].resize(5);
    }

    for (size_t i = 0; i < NElement_; ++i)
    {
        std::getline(infile, line);
        std::stringstream ss(line);
        unsigned int celltype;
        ss >> celltype;
        if (celltype == QUAD4)
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
        ss >> NPoint_;
        std::cout << "Number of Points: " << NPoint_ << std::endl;
    }

    // Get nodearray_
    nodearray_.resize(NPoint_);
    for (size_t i = 0; i < NPoint_; ++i)
    {
        std::getline(infile, line);
        std::stringstream ss(line);
        double x;
        double y;
        unsigned int id;
        ss >> x >> y >> id;
        Node2d node_obj(id, x, y);
        nodearray_[i] = node_obj;
    }

    // Get Markers
    std::getline(infile, line);
    if (line.find("NMARK") != std::string::npos)
    {
        std::stringstream ss(line);
        // skip until '='
        ss.ignore(line.size(), '=');
        ss >> NMarker_;
        std::cout << "Number of Markers: " << NMarker_ << std::endl;
    }

    for (size_t i = 0; i < NMarker_; ++i)
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
            unsigned int type;
            ss >> type;
            if (type == LINE)
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
    cellarray_.resize(NElement_);
    for (size_t i = 0; i < element_table_.size(); ++i)
    {
        unsigned int id = element_table_[i][0];
        Node2d& node1   = nodearray_[element_table_[i][1]];
        Node2d& node2   = nodearray_[element_table_[i][2]];
        Node2d& node3   = nodearray_[element_table_[i][3]];
        Node2d& node4   = nodearray_[element_table_[i][4]];

        CellQuad4 quad(id, &node1, &node2, &node3, &node4);
        cellarray_[i] = quad;
    }
}

void SU2meshparser::SetMarkersToFaces()
{
    for (size_t iTable = 0; iTable < marker_table_.marker_array.size();
         ++iTable)
    {
        std::string mark = marker_table_.marker_array[iTable];
        std::array<unsigned int, 2> markered_edge =
            marker_table_.edge_array[iTable];
        for (size_t iCell = 0; iCell < cellarray_.size(); ++iCell)
        {
            std::array<unsigned int, 2> edge1;
            edge1[0] = cellarray_[iCell].GetNode1()->GetID();
            edge1[1] = cellarray_[iCell].GetNode2()->GetID();
            std::array<unsigned int, 2> edge2;
            edge2[0] = cellarray_[iCell].GetNode2()->GetID();
            edge2[1] = cellarray_[iCell].GetNode3()->GetID();
            std::array<unsigned int, 2> edge3;
            edge3[0] = cellarray_[iCell].GetNode3()->GetID();
            edge3[1] = cellarray_[iCell].GetNode4()->GetID();
            std::array<unsigned int, 2> edge4;
            edge4[0] = cellarray_[iCell].GetNode4()->GetID();
            edge4[1] = cellarray_[iCell].GetNode1()->GetID();

            bool edge1_flag =
                (edge1[0] == markered_edge[0] && edge1[1] == markered_edge[1]);
            bool edge2_flag =
                (edge2[0] == markered_edge[0] && edge2[1] == markered_edge[1]);
            bool edge3_flag =
                (edge3[0] == markered_edge[0] && edge3[1] == markered_edge[1]);
            bool edge4_flag =
                (edge4[0] == markered_edge[0] && edge4[1] == markered_edge[1]);

            if (edge1_flag)
            {
                cellarray_[iCell].Face1()->SetTag(mark);
            }
            else if (edge2_flag)
            {
                cellarray_[iCell].Face2()->SetTag(mark);
            }
            else if (edge3_flag)
            {
                cellarray_[iCell].Face3()->SetTag(mark);
            }
            else if (edge4_flag)
            {
                cellarray_[iCell].Face4()->SetTag(mark);
            }
            else
            {
                //
            }
        }
    }

    /*--- Set interior tags to interior faces ---*/
    for (size_t iCell = 0; iCell < cellarray_.size(); ++iCell)
    {
        if (cellarray_[iCell].Face1()->GetTag().empty())
        {
            cellarray_[iCell].Face1()->SetTag("interior");
        }
        if (cellarray_[iCell].Face2()->GetTag().empty())
        {
            cellarray_[iCell].Face2()->SetTag("interior");
        }
        if (cellarray_[iCell].Face3()->GetTag().empty())
        {
            cellarray_[iCell].Face3()->SetTag("interior");
        }
        if (cellarray_[iCell].Face4()->GetTag().empty())
        {
            cellarray_[iCell].Face4()->SetTag("interior");
        }
    }
}

void SU2meshparser::SetNeighborCells()
{
    for (size_t iCell = 0; iCell < cellarray_.size(); ++iCell)
    {
        // Search over the element table & Get candidate neighbor cells
        std::vector<unsigned int> candidateNeighborCellsIDs;
        unsigned int node1ID = cellarray_[iCell].GetNode1()->GetID();
        unsigned int node2ID = cellarray_[iCell].GetNode2()->GetID();
        unsigned int node3ID = cellarray_[iCell].GetNode3()->GetID();
        unsigned int node4ID = cellarray_[iCell].GetNode4()->GetID();

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
        bool isInterior_face1 =
            (cellarray_[iCell].Face1()->GetTag() == "interior");
        bool isInterior_face2 =
            (cellarray_[iCell].Face2()->GetTag() == "interior");
        bool isInterior_face3 =
            (cellarray_[iCell].Face3()->GetTag() == "interior");
        bool isInterior_face4 =
            (cellarray_[iCell].Face4()->GetTag() == "interior");

        std::array<unsigned int, 2> edge1;
        edge1[0] = cellarray_[iCell].GetNode1()->GetID();
        edge1[1] = cellarray_[iCell].GetNode2()->GetID();
        std::array<unsigned int, 2> edge2;
        edge2[0] = cellarray_[iCell].GetNode2()->GetID();
        edge2[1] = cellarray_[iCell].GetNode3()->GetID();
        std::array<unsigned int, 2> edge3;
        edge3[0] = cellarray_[iCell].GetNode3()->GetID();
        edge3[1] = cellarray_[iCell].GetNode4()->GetID();
        std::array<unsigned int, 2> edge4;
        edge4[0] = cellarray_[iCell].GetNode4()->GetID();
        edge4[1] = cellarray_[iCell].GetNode1()->GetID();

        unsigned int nNeighbors = candidateNeighborCellsIDs.size();
        if (isInterior_face1)
        {
            for (size_t i = 0; i < nNeighbors; ++i)
            {
                std::array<unsigned int, 2> candidate_edge1;
                candidate_edge1[1] = cellarray_[candidateNeighborCellsIDs[i]]
                                         .GetNode1()
                                         ->GetID();
                candidate_edge1[0] = cellarray_[candidateNeighborCellsIDs[i]]
                                         .GetNode2()
                                         ->GetID();
                std::array<unsigned int, 2> candidate_edge2;
                candidate_edge2[1] = cellarray_[candidateNeighborCellsIDs[i]]
                                         .GetNode2()
                                         ->GetID();
                candidate_edge2[0] = cellarray_[candidateNeighborCellsIDs[i]]
                                         .GetNode3()
                                         ->GetID();
                std::array<unsigned int, 2> candidate_edge3;
                candidate_edge3[1] = cellarray_[candidateNeighborCellsIDs[i]]
                                         .GetNode3()
                                         ->GetID();
                candidate_edge3[0] = cellarray_[candidateNeighborCellsIDs[i]]
                                         .GetNode4()
                                         ->GetID();
                std::array<unsigned int, 2> candidate_edge4;
                candidate_edge4[1] = cellarray_[candidateNeighborCellsIDs[i]]
                                         .GetNode4()
                                         ->GetID();
                candidate_edge4[0] = cellarray_[candidateNeighborCellsIDs[i]]
                                         .GetNode1()
                                         ->GetID();

                if (edge1 == candidate_edge1 || edge1 == candidate_edge2 ||
                    edge1 == candidate_edge3 || edge1 == candidate_edge4)
                {
                    cellarray_[iCell].SetNeighbor1Ptr(
                        &cellarray_[candidateNeighborCellsIDs[i]]);
                }
            }
        }
        if (isInterior_face2)
        {
            for (size_t i = 0; i < nNeighbors; ++i)
            {
                std::array<unsigned int, 2> candidate_edge1;
                candidate_edge1[1] = cellarray_[candidateNeighborCellsIDs[i]]
                                         .GetNode1()
                                         ->GetID();
                candidate_edge1[0] = cellarray_[candidateNeighborCellsIDs[i]]
                                         .GetNode2()
                                         ->GetID();
                std::array<unsigned int, 2> candidate_edge2;
                candidate_edge2[1] = cellarray_[candidateNeighborCellsIDs[i]]
                                         .GetNode2()
                                         ->GetID();
                candidate_edge2[0] = cellarray_[candidateNeighborCellsIDs[i]]
                                         .GetNode3()
                                         ->GetID();
                std::array<unsigned int, 2> candidate_edge3;
                candidate_edge3[1] = cellarray_[candidateNeighborCellsIDs[i]]
                                         .GetNode3()
                                         ->GetID();
                candidate_edge3[0] = cellarray_[candidateNeighborCellsIDs[i]]
                                         .GetNode4()
                                         ->GetID();
                std::array<unsigned int, 2> candidate_edge4;
                candidate_edge4[1] = cellarray_[candidateNeighborCellsIDs[i]]
                                         .GetNode4()
                                         ->GetID();
                candidate_edge4[0] = cellarray_[candidateNeighborCellsIDs[i]]
                                         .GetNode1()
                                         ->GetID();

                if (edge2 == candidate_edge1 || edge2 == candidate_edge2 ||
                    edge2 == candidate_edge3 || edge2 == candidate_edge4)
                {
                    cellarray_[iCell].SetNeighbor2Ptr(
                        &cellarray_[candidateNeighborCellsIDs[i]]);
                }
            }
        }
        if (isInterior_face3)
        {
            for (size_t i = 0; i < nNeighbors; ++i)
            {
                std::array<unsigned int, 2> candidate_edge1;
                candidate_edge1[1] = cellarray_[candidateNeighborCellsIDs[i]]
                                         .GetNode1()
                                         ->GetID();
                candidate_edge1[0] = cellarray_[candidateNeighborCellsIDs[i]]
                                         .GetNode2()
                                         ->GetID();
                std::array<unsigned int, 2> candidate_edge2;
                candidate_edge2[1] = cellarray_[candidateNeighborCellsIDs[i]]
                                         .GetNode2()
                                         ->GetID();
                candidate_edge2[0] = cellarray_[candidateNeighborCellsIDs[i]]
                                         .GetNode3()
                                         ->GetID();
                std::array<unsigned int, 2> candidate_edge3;
                candidate_edge3[1] = cellarray_[candidateNeighborCellsIDs[i]]
                                         .GetNode3()
                                         ->GetID();
                candidate_edge3[0] = cellarray_[candidateNeighborCellsIDs[i]]
                                         .GetNode4()
                                         ->GetID();
                std::array<unsigned int, 2> candidate_edge4;
                candidate_edge4[1] = cellarray_[candidateNeighborCellsIDs[i]]
                                         .GetNode4()
                                         ->GetID();
                candidate_edge4[0] = cellarray_[candidateNeighborCellsIDs[i]]
                                         .GetNode1()
                                         ->GetID();

                if (edge3 == candidate_edge1 || edge3 == candidate_edge2 ||
                    edge3 == candidate_edge3 || edge3 == candidate_edge4)
                {
                    cellarray_[iCell].SetNeighbor3Ptr(
                        &cellarray_[candidateNeighborCellsIDs[i]]);
                }
            }
        }
        if (isInterior_face4)
        {
            for (size_t i = 0; i < nNeighbors; ++i)
            {
                std::array<unsigned int, 2> candidate_edge1;
                candidate_edge1[1] = cellarray_[candidateNeighborCellsIDs[i]]
                                         .GetNode1()
                                         ->GetID();
                candidate_edge1[0] = cellarray_[candidateNeighborCellsIDs[i]]
                                         .GetNode2()
                                         ->GetID();
                std::array<unsigned int, 2> candidate_edge2;
                candidate_edge2[1] = cellarray_[candidateNeighborCellsIDs[i]]
                                         .GetNode2()
                                         ->GetID();
                candidate_edge2[0] = cellarray_[candidateNeighborCellsIDs[i]]
                                         .GetNode3()
                                         ->GetID();
                std::array<unsigned int, 2> candidate_edge3;
                candidate_edge3[1] = cellarray_[candidateNeighborCellsIDs[i]]
                                         .GetNode3()
                                         ->GetID();
                candidate_edge3[0] = cellarray_[candidateNeighborCellsIDs[i]]
                                         .GetNode4()
                                         ->GetID();
                std::array<unsigned int, 2> candidate_edge4;
                candidate_edge4[1] = cellarray_[candidateNeighborCellsIDs[i]]
                                         .GetNode4()
                                         ->GetID();
                candidate_edge4[0] = cellarray_[candidateNeighborCellsIDs[i]]
                                         .GetNode1()
                                         ->GetID();

                if (edge4 == candidate_edge1 || edge4 == candidate_edge2 ||
                    edge4 == candidate_edge3 || edge4 == candidate_edge4)
                {
                    cellarray_[iCell].SetNeighbor4Ptr(
                        &cellarray_[candidateNeighborCellsIDs[i]]);
                }
            }
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
    for (size_t i = 0; i < nodearray_.size(); ++i)
    {
        std::cout << nodearray_[i].GetID() << "\t"  //
                  << nodearray_[i].GetX() << "\t"   //
                  << nodearray_[i].GetY() << std::endl;
    }

    std::cout << "\nCell array: " << std::endl;
    for (size_t i = 0; i < cellarray_.size(); ++i)
    {
        std::cout << cellarray_[i].GetID() << "\t"
                  << cellarray_[i].GetNode1()->GetID() << "\t"
                  << cellarray_[i].GetNode2()->GetID() << "\t"
                  << cellarray_[i].GetNode3()->GetID() << "\t"
                  << cellarray_[i].GetNode4()->GetID() << std::endl;
    }
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

    for (size_t iCell = 0; iCell < cellarray_.size(); ++iCell)
    {
        if (!cellarray_[iCell].Face1()->GetTag().empty())
        {
            std::cout << "Cell " << cellarray_[iCell].GetID()
                      << " has a face w/ tag "
                      << cellarray_[iCell].Face1()->GetTag() << std::endl;
        }
        if (!cellarray_[iCell].Face2()->GetTag().empty())
        {
            std::cout << "Cell " << cellarray_[iCell].GetID()
                      << " has a face w/ tag "
                      << cellarray_[iCell].Face2()->GetTag() << std::endl;
        }
        if (!cellarray_[iCell].Face3()->GetTag().empty())
        {
            std::cout << "Cell " << cellarray_[iCell].GetID()
                      << " has a face w/ tag "
                      << cellarray_[iCell].Face3()->GetTag() << std::endl;
        }
        if (!cellarray_[iCell].Face4()->GetTag().empty())
        {
            std::cout << "Cell " << cellarray_[iCell].GetID()
                      << " has a face w/ tag "
                      << cellarray_[iCell].Face4()->GetTag() << std::endl;
        }
    }
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

    for (size_t iCell = 0; iCell < cellarray_.size(); ++iCell)
    {
        auto neighbor1Ptr = cellarray_[iCell].GetNeighbors1Ptr();
        auto neighbor2Ptr = cellarray_[iCell].GetNeighbors2Ptr();
        auto neighbor3Ptr = cellarray_[iCell].GetNeighbors3Ptr();
        auto neighbor4Ptr = cellarray_[iCell].GetNeighbors4Ptr();

        if (neighbor1Ptr != nullptr)
        {
            std::cout << "Cell " << cellarray_[iCell].GetID()
                      << " is ajacent to "
                      << cellarray_[iCell].GetNeighbors1Ptr()->GetID()
                      << std::endl;
            ;
        }
        if (neighbor2Ptr != nullptr)
        {
            std::cout << "Cell " << cellarray_[iCell].GetID()
                      << " is ajacent to "
                      << cellarray_[iCell].GetNeighbors2Ptr()->GetID()
                      << std::endl;
        }
        if (neighbor3Ptr != nullptr)
        {
            std::cout << "Cell " << cellarray_[iCell].GetID()
                      << " is ajacent to "
                      << cellarray_[iCell].GetNeighbors3Ptr()->GetID()
                      << std::endl;
            ;
        }
        if (neighbor4Ptr != nullptr)
        {
            std::cout << "Cell " << cellarray_[iCell].GetID()
                      << " is ajacent to "
                      << cellarray_[iCell].GetNeighbors4Ptr()->GetID()
                      << std::endl;
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

    outfile << "POINTS " << NPoint_ << " "
            << "double " << std::endl;
    for (unsigned int i = 0; i < NPoint_; ++i)
    {
        outfile << nodearray_[i].GetX() << "\t"  //
                << nodearray_[i].GetY() << "\t" << 0 << std::endl;
    }

    outfile << "CELLS " << NElement_ << " " << 5 * NElement_ << std::endl;
    for (unsigned int i = 0; i < NElement_; ++i)
    {
        outfile << "4\t"  // number of point
                << cellarray_[i].GetNode1()->GetID() << "\t"
                << cellarray_[i].GetNode2()->GetID() << "\t"
                << cellarray_[i].GetNode3()->GetID() << "\t"
                << cellarray_[i].GetNode4()->GetID() << std::endl;
    }

    outfile << "CELL_TYPES " << NElement_ << std::endl;
    for (unsigned int i = 0; i < NElement_; ++i)
    {
        // VTK_QUAD cell type is defined as 9
        outfile << "9" << std::endl;
    }

    /* Set variables. You can see variables defined here from Paraview */
    outfile << "CELL_DATA " << NElement_
            << std::endl;  // Declare the variable below is set in cells
    outfile << "SCALARS ID int" << std::endl;  // give its KIND name type
    outfile << "LOOKUP_TABLE default" << std::endl;
    for (unsigned int i = 0; i < NElement_; i++)
    {
        outfile << cellarray_[i].GetID() << std::endl;
    }

    std::cout << "File writing finished!" << std::endl;
}

}  // namespace su2mesh
