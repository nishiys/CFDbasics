#pragma once

#include <Eigen/Core>
#include <array>
#include <vector>

#include "Face2d.hpp"
// #include "Node2d.hpp"
#include "Variable.hpp"

class Node2d;
class CellQuad4
{
public:
    CellQuad4();
    CellQuad4(unsigned int id, Node2d* pNode1, Node2d* pNode2, Node2d* pNode3,
              Node2d* pNode4);
    ~CellQuad4();

    unsigned int GetID() const;
    Node2d* GetNode(unsigned int index) const;
    Face2d* Face(unsigned int index);
    void SetNeighborPtr(unsigned int index, CellQuad4* pCellQuad);
    CellQuad4* GetNeighborPtr(unsigned int index) const;
    Eigen::Vector2d GetCentroid() const;
    void SetVectorToNeighbor(unsigned int index);
    Eigen::Vector2d GetVectorToNeighbor(unsigned int index) const;
    Eigen::Vector2d GetNormalVectorToBoundary(unsigned int index) const;
    void SetMag_N1Vectors(unsigned int index);
    double GetMag_N1Vectors(unsigned int index) const;
    Eigen::Vector2d GetvectorToFace(unsigned int index) const;
    Eigen::Vector2d GetVectorToNode(unsigned int index) const;
    double GetVolume() const;

    Variable cellvar;

private:
    unsigned int id_;
    std::array<Node2d*, 4> pNodes_;
    std::array<Face2d, 4> faces_;
    std::array<CellQuad4*, 4> pNeighbors_;
    std::array<Eigen::Vector2d, 4> vectorToNeighbors_;
    std::array<Eigen::Vector2d, 4> normalvectorToBoundaries_;
    std::array<double, 4> mag_n1Vectors_;
    std::array<Eigen::Vector2d, 4> vectorToFaces_;
    std::array<Eigen::Vector2d, 4> vectorToNodes_;

    void CheckFaceNormals();
    void CalcCentroid();
    void CalcVolume();
    void CalcVectorToFaces();
    void CalcVectorToNodes();

    Eigen::Vector2d centroid_;
    double volume_;
};