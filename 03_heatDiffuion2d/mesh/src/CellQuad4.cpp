#include "CellQuad4.hpp"

#include <iostream>

CellQuad4::CellQuad4() {}
CellQuad4::CellQuad4(unsigned int id, Node2d* pNode1, Node2d* pNode2,
                     Node2d* pNode3, Node2d* pNode4)
    : pNode1_(pNode1),
      pNode2_(pNode2),
      pNode3_(pNode3),
      pNode4_(pNode4),
      id_(id),
      face1_(pNode1_, pNode2_),
      face2_(pNode2_, pNode3_),
      face3_(pNode3_, pNode4_),
      face4_(pNode4_, pNode1_),
      centroid_(2)
{
    CalcCentroid();
    CheckFaceNormals();
    CalcVolume();
}

CellQuad4::~CellQuad4() {}

void CellQuad4::CheckFaceNormals()
{
    Eigen::VectorXd cellcent_to_face1cent(2);
    Eigen::VectorXd cellcent_to_face2cent(2);
    Eigen::VectorXd cellcent_to_face3cent(2);
    Eigen::VectorXd cellcent_to_face4cent(2);

    cellcent_to_face1cent = face1_.GetFaceCenter() - centroid_;
    cellcent_to_face2cent = face2_.GetFaceCenter() - centroid_;
    cellcent_to_face3cent = face3_.GetFaceCenter() - centroid_;
    cellcent_to_face4cent = face4_.GetFaceCenter() - centroid_;

    if (cellcent_to_face1cent.dot(face1_.GetNormalVec()) < 0)
    {
        face1_.FlipNormalVec();
    }
    if (cellcent_to_face2cent.dot(face2_.GetNormalVec()) < 0)
    {
        face2_.FlipNormalVec();
    }
    if (cellcent_to_face3cent.dot(face3_.GetNormalVec()) < 0)
    {
        face3_.FlipNormalVec();
    }
    if (cellcent_to_face4cent.dot(face4_.GetNormalVec()) < 0)
    {
        face4_.FlipNormalVec();
    }
}

void CellQuad4::CalcCentroid()
{
    centroid_(0) = 0.25 * (pNode1_->GetX() + pNode2_->GetX() + pNode3_->GetX() +
                           pNode4_->GetX());
    centroid_(1) = 0.25 * (pNode1_->GetY() + pNode2_->GetY() + pNode3_->GetY() +
                           pNode4_->GetY());

    std::cout << "Cell " << id_ << " centroid: (" << centroid_(0) << ", "
              << centroid_(1) << ")" << std::endl;
}

void CellQuad4::CalcVolume()
{
    // using Gauss divergence theorem
    volume_ = 0.5 * face1_.GetNormalVec().dot(face1_.GetFaceCenter());
    volume_ += 0.5 * face2_.GetNormalVec().dot(face2_.GetFaceCenter());
    volume_ += 0.5 * face3_.GetNormalVec().dot(face3_.GetFaceCenter());
    volume_ += 0.5 * face4_.GetNormalVec().dot(face4_.GetFaceCenter());

    std::cout << "Cell " << id_ << " volume: " << volume_ << std::endl;
}