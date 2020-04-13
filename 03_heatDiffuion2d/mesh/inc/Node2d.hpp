#pragma once

#include <Eigen/Core>
#include <array>
#include <string>
#include <vector>

// #include "CellQuad4.hpp"
#include "Variable.hpp"

class CellQuad4;
class Node2d
{
public:
    Node2d();
    Node2d(unsigned int id, double x, double y);
    ~Node2d();

    void AddCellPointer(CellQuad4* p_cellquad4);
    unsigned int GetID() const;
    Eigen::Vector2d GetCoords() const;
    std::vector<CellQuad4*> GetCellPtrs();
    void AddMarker(std::string marker);
    bool IsInterior() const;
    bool IsCorner() const;

    Variable nodevar;

private:
    unsigned int id_;
    Eigen::Vector2d coords_;
    std::vector<CellQuad4*> pCellQuad4s_;

    //! markers to detect the node status
    std::vector<std::string> markers_;
};