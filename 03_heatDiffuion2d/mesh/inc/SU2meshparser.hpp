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

private:
    std::string meshfilename_;

    void ReadFile();
    void CreateQuadArray();

    std::vector<CellQuad4> cellarray_;
    // std::vector<Face2d> facearray_;
    std::vector<Node2d> nodearray_;

    unsigned int DIM_;
    unsigned int NElement_;
    unsigned int NPoint_;
    unsigned int NMarker_;

    const unsigned int LINE  = 3;
    const unsigned int QUAD4 = 9;
    std::vector<std::vector<unsigned int>> element_table_;
    std::vector<std::vector<std::string>> marker_table_;

    void PrintDebug();
};
}  // namespace su2mesh
