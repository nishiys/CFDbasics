#include "Face2d.hpp"

#include <iostream>

#include "Node2d.hpp"

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
    facecenter_(0) = 0.5 * (pNode1_->GetX() + pNode2_->GetX());
    facecenter_(1) = 0.5 * (pNode1_->GetY() + pNode2_->GetY());
}

void Face2d::CalcNormalVec()
{
    Eigen::Vector2d edgevec(2);
    edgevec(0) = pNode2_->GetX() - pNode1_->GetX();
    edgevec(1) = pNode2_->GetY() - pNode1_->GetY();

    Eigen::MatrixXd Rotm90(2, 2);
    Rotm90 << 0.0, 1.0, -1.0, 0.0;
    normalvec_ = Rotm90 * edgevec;

    normalvec_.normalize();

    // std::cout << "Face Normal Vector: (" << normalvec_(0) << ", "
    //           << normalvec_(1) << ")" << std::endl;
}

void Face2d::CalcArea()
{
    Eigen::Vector2d edgevec(2);
    edgevec(0) = pNode2_->GetX() - pNode1_->GetX();
    edgevec(1) = pNode2_->GetY() - pNode1_->GetY();

    area_ = edgevec.norm() * 1.0;
}

void Face2d::FlipNormalVec()
{
    Eigen::MatrixXd Rot180(2, 2);
    Rot180 << -1.0, 0.0, 0.0, -1.0;
    normalvec_ = Rot180 * normalvec_;
}