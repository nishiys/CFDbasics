#include "Node2d.hpp"

#include "CellQuad4.hpp"

Node2d::Node2d() {}

Node2d::Node2d(unsigned int id, double x, double y) : id_(id)
{
    coords_(0) = x;
    coords_(1) = y;
}

Node2d::~Node2d() {}

/*----------------------------------------------------------------------------------*/
/* Global Methods */
/*----------------------------------------------------------------------------------*/

unsigned int Node2d::GetID() const { return id_; };

Eigen::Vector2d Node2d::GetCoords() const { return coords_; };

void Node2d::AddMarker(std::string marker) { markers_.push_back(marker); }

void Node2d::AddCellPointer(CellQuad4* p_cellquad4)
{
    pCellQuad4s_.push_back(p_cellquad4);
}

std::vector<CellQuad4*> Node2d::GetCellPtrs() { return pCellQuad4s_; }

bool Node2d::IsInterior() const
{
    if (markers_.size() == 0)
        return true;
    else
        return false;
}
bool Node2d::IsCorner() const
{
    if (markers_.size() > 1)
        return true;
    else
        return false;
}