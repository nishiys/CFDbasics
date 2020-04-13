#include "Face2d.hpp"

#include <iostream>

// #include "Node2d.hpp"

Face2d::Face2d() {}
Face2d::Face2d(Node2d* pNode1, Node2d* pNode2)
    : pNode1_(pNode1), pNode2_(pNode2), normalvec_(2), facecenter_(2)
{
    CalcFaceCenter();
    CalcNormalVec();
    CalcArea();
}

Face2d::~Face2d() {}

void Face2d::CalcFaceCenter()
{
    facecenter_ = 0.5 * (pNode1_->GetCoords() + pNode2_->GetCoords());
}

void Face2d::CalcNormalVec()
{
    Eigen::Vector2d edgevec = pNode2_->GetCoords() - pNode1_->GetCoords();

    Eigen::MatrixXd Rotm90(2, 2);
    Rotm90 << 0.0, 1.0, -1.0, 0.0;
    normalvec_ = Rotm90 * edgevec;

    normalvec_.normalize();
}

void Face2d::CalcArea()
{
    Eigen::Vector2d edgevec = pNode2_->GetCoords() - pNode1_->GetCoords();

    area_ = edgevec.norm() * 1.0;
}

void Face2d::FlipNormalVec()
{
    Eigen::MatrixXd Rot180(2, 2);
    Rot180 << -1.0, 0.0, 0.0, -1.0;
    normalvec_ = Rot180 * normalvec_;
}