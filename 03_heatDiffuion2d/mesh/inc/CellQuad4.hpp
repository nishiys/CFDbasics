#pragma once

#include <Eigen/Core>
#include <array>
#include <vector>

#include "Face2d.hpp"
#include "Node2d.hpp"
#include "Variable.hpp"

class CellQuad4
{
public:
    CellQuad4();
    CellQuad4(unsigned int id, Node2d* pNode1, Node2d* pNode2, Node2d* pNode3,
              Node2d* pNode4);
    ~CellQuad4();

    inline unsigned int GetID() const { return id_; }
    inline Node2d* GetNode(unsigned int index) const { return pNodes_[index]; }
    inline Face2d* Face(unsigned int index) { return &faces_[index]; }
    inline void SetNeighborPtr(unsigned int index, CellQuad4* pCellQuad)
    {
        pNeighbors_[index] = pCellQuad;
    }
    inline CellQuad4* GetNeighborPtr(unsigned int index) const
    {
        return pNeighbors_[index];
    }
    inline Eigen::Vector2d GetCentroid() const { return centroid_; }
    inline void SetVectorToNeighbor(unsigned int index)
    {
        if (pNeighbors_[index] != nullptr)
        {
            vectorToNeighbors_[index] =
                (pNeighbors_[index]->GetCentroid() - centroid_);
        }
        else
        {
            Eigen::Vector2d vectorToBoundary =
                faces_[index].GetFaceCenter() - centroid_;
            Eigen::Vector2d facenormal_vector = faces_[index].GetNormalVec();
            double y_magnitude = vectorToBoundary.dot(facenormal_vector);
            normalvectorToBoundaries_[index] = y_magnitude * facenormal_vector;
        }
    }
    inline Eigen::Vector2d GetVectorToNeighbor(unsigned int index) const
    {
        return vectorToNeighbors_[index];
    }
    inline Eigen::Vector2d GetNormalVectorToBoundary(unsigned int index) const
    {
        return normalvectorToBoundaries_[index];
    }
    inline void SetMag_N1Vectors(unsigned int index)
    {
        if (pNeighbors_[index] != nullptr)
        {
            mag_n1Vectors_[index] =
                1.0 / vectorToNeighbors_[index].norm() *
                faces_[index].GetNormalVec().dot(vectorToNeighbors_[index]);
        }
    }
    inline double GetMag_N1Vectors(unsigned int index) const
    {
        return mag_n1Vectors_[index];
    }

    inline double GetVolume() const { return volume_; }

    Variable cellvar;

private:
    unsigned int id_;
    std::array<Node2d*, 4> pNodes_;
    std::array<Face2d, 4> faces_;
    std::array<CellQuad4*, 4> pNeighbors_;
    std::array<Eigen::Vector2d, 4> vectorToNeighbors_;
    std::array<Eigen::Vector2d, 4> normalvectorToBoundaries_;
    std::array<double, 4> mag_n1Vectors_;

    void CheckFaceNormals();
    void CalcCentroid();
    void CalcVolume();

    Eigen::Vector2d centroid_;
    double volume_;
};