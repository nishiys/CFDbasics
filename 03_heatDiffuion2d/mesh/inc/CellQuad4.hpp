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

    inline unsigned int GetID() const { return id_; };

    inline Node2d* GetNode(unsigned int index) const { return pNodes_[index]; };
    // inline Node2d* GetNode1() const { return pNode1_; };
    // inline Node2d* GetNode2() const { return pNode2_; };
    // inline Node2d* GetNode3() const { return pNode3_; };
    // inline Node2d* GetNode4() const { return pNode4_; };

    inline Face2d* Face(unsigned int index) { return &faces_[index]; };
    // inline Face2d* Face1() { return &face1_; };
    // inline Face2d* Face2() { return &face2_; };
    // inline Face2d* Face3() { return &face3_; };
    // inline Face2d* Face4() { return &face4_; };

    inline void SetNeighborPtr(unsigned int index, CellQuad4* pCellQuad)
    {
        pNeighbors_[index] = pCellQuad;
    };

    // inline void SetNeighbor1Ptr(CellQuad4* pCellQuad)
    // {
    //     pNeighbor1_ = pCellQuad;
    // };
    // inline void SetNeighbor2Ptr(CellQuad4* pCellQuad)
    // {
    //     pNeighbor2_ = pCellQuad;
    // };
    // inline void SetNeighbor3Ptr(CellQuad4* pCellQuad)
    // {
    //     pNeighbor3_ = pCellQuad;
    // };
    // inline void SetNeighbor4Ptr(CellQuad4* pCellQuad)
    // {
    //     pNeighbor4_ = pCellQuad;
    // };

    inline CellQuad4* GetNeighborPtr(unsigned int index) const
    {
        return pNeighbors_[index];
    };
    // inline CellQuad4* GetNeighbors2Ptr() const { return pNeighbor2_; };t
    // inline CellQuad4* GetNeighbors1Ptr() const { return pNeighbor1_; };
    // inline CellQuad4* GetNeighbors2Ptr() const { return pNeighbor2_; };
    // inline CellQuad4* GetNeighbors3Ptr() const { return pNeighbor3_; };
    // inline CellQuad4* GetNeighbors4Ptr() const { return pNeighbor4_; };

    Variable var;

private:
    unsigned int id_;
    std::array<Node2d*, 4> pNodes_;
    // Node2d* pNode1_;
    // Node2d* pNode2_;
    // Node2d* pNode3_;
    // Node2d* pNode4_;

    std::array<Face2d, 4> faces_;
    // Face2d face1_;
    // Face2d face2_;
    // Face2d face3_;
    // Face2d face4_;

    std::array<CellQuad4*, 4> pNeighbors_;
    // CellQuad4* pNeighbor1_;
    // CellQuad4* pNeighbor2_;
    // CellQuad4* pNeighbor3_;
    // CellQuad4* pNeighbor4_;

    void CheckFaceNormals();
    void CalcCentroid();
    void CalcVolume();

    Eigen::VectorXd centroid_;
    double volume_;
};