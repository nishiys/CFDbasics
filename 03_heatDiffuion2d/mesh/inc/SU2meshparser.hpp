#pragma once

#include <string>
#include <vector>

#include "CellQuad4.hpp"
#include "Face2d.hpp"
#include "Node2d.hpp"

namespace su2mesh
{
class SU2meshparser
{
public:
    SU2meshparser(std::string meshfilename);
    ~SU2meshparser();

    void LoadData();
    void WriteVtkFile(std::string vtkfilename);

    std::vector<CellQuad4> cellarray;

private:
    std::string meshfilename_;

    void ReadFile();
    void CreateQuadArray();
    void SetMarkersToFaces();
    void SetNeighborCells();

    std::vector<Node2d> nodearray_;

    unsigned int nDim_;
    unsigned int nElement_;
    unsigned int nPoint_;
    unsigned int nMarker_;

    const unsigned int VTK_LINE  = 3;
    const unsigned int VTK_QUAD4 = 9;
    std::vector<std::vector<unsigned int>> element_table_;
    struct MarkerTable
    {
    public:
        std::vector<std::string> marker_array;
        std::vector<std::array<unsigned int, 2>> edge_array;
    };
    MarkerTable marker_table_;

    void PrintDebug();
};
}  // namespace su2mesh
