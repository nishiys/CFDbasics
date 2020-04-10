#include "SU2meshparser.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

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
                std::vector<std::string> vec(3);
                vec[0] = tag;
                ss >> vec[1] >> vec[2];
                marker_table_.push_back(vec);
            }
        }
    }

    // Skip unnecessary lines
    while (!infile.eof())
    {
        std::getline(infile, line);
    }
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

void SU2meshparser::PrintDebug()
{
    std::cout << "\nElement table: " << std::endl;
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
