#include "SU2meshparser.hpp"

int main()
{
    su2mesh::SU2meshparser meshparser("quad.su2");
    meshparser.LoadData();
    meshparser.WriteVtkFile("quad.vtk");

    return 0;
}