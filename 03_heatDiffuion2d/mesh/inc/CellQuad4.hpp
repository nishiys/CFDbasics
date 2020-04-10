#pragma once

#include <Eigen/Core>
#include <array>
#include <vector>

#include "Face2d.hpp"
#include "Node2d.hpp"

class CellQuad4
{
public:
    CellQuad4();
    CellQuad4(unsigned int id, Node2d* pNode1, Node2d* pNode2, Node2d* pNode3,
              Node2d* pNode4);
    ~CellQuad4();

    inline unsigned int GetID() const { return id_; };
    inline Node2d* GetNode1() const { return pNode1_; };
    inline Node2d* GetNode2() const { return pNode2_; };
    inline Node2d* GetNode3() const { return pNode3_; };
    inline Node2d* GetNode4() const { return pNode4_; };

    inline Face2d* Face1() { return &face1_; };
    inline Face2d* Face2() { return &face2_; };
    inline Face2d* Face3() { return &face3_; };
    inline Face2d* Face4() { return &face4_; };

private:
    unsigned int id_;
    Node2d* pNode1_;
    Node2d* pNode2_;
    Node2d* pNode3_;
    Node2d* pNode4_;

    Face2d face1_;
    Face2d face2_;
    Face2d face3_;
    Face2d face4_;

    CellQuad4* neighbor1_;
    CellQuad4* neighbor2_;
    CellQuad4* neighbor3_;
    CellQuad4* neighbor4_;

    void CheckFaceNormals();
    void CalcCentroid();
    void CalcVolume();

    Eigen::VectorXd centroid_;
    double volume_;
};