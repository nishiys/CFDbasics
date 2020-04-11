#include "CellQuad4.hpp"

#include <iostream>

CellQuad4::CellQuad4() {}
CellQuad4::CellQuad4(unsigned int id, Node2d* pNode1, Node2d* pNode2,
                     Node2d* pNode3, Node2d* pNode4)
    : pNodes_{pNode1, pNode2, pNode3, pNode4},
      id_(id),
      centroid_(2),
      pNeighbors_{nullptr, nullptr, nullptr, nullptr}
{
    Face2d face1(pNode1, pNode2);
    Face2d face2(pNode2, pNode3);
    Face2d face3(pNode3, pNode4);
    Face2d face4(pNode4, pNode1);
    faces_ = {face1, face2, face3, face4};

    /*--- Calculate fundamental values ---*/
    CalcCentroid();
    CheckFaceNormals();
    CalcVolume();
}

CellQuad4::~CellQuad4() {}

void CellQuad4::CheckFaceNormals()
{
    std::array<Eigen::Vector2d, 4> cellcent_to_facecent;

    for (size_t iFace = 0; iFace < faces_.size(); ++iFace)
    {
        cellcent_to_facecent[iFace] = faces_[iFace].GetFaceCenter() - centroid_;
        if (cellcent_to_facecent[iFace].dot(faces_[iFace].GetNormalVec()) < 0)
        {
            faces_[iFace].FlipNormalVec();
        }
    }
}

void CellQuad4::CalcCentroid()
{
    unsigned int nNode = pNodes_.size();
    centroid_(0)       = 0;
    centroid_(1)       = 0;
    for (unsigned int iNode = 0; iNode < nNode; ++iNode)
    {
        centroid_(0) +=
            1.0 / static_cast<double>(nNode) * pNodes_[iNode]->GetX();
        centroid_(1) +=
            1.0 / static_cast<double>(nNode) * pNodes_[iNode]->GetY();
    }

    std::cout << "Cell " << id_ << " centroid: (" << centroid_(0) << ", "
              << centroid_(1) << ")" << std::endl;
}

void CellQuad4::CalcVolume()
{
    // using Gauss divergence theorem
    unsigned int nDim = 2;
    volume_           = 0;
    for (unsigned int iFace = 0; iFace < faces_.size(); ++iFace)
    {
        volume_ +=
            1.0 / static_cast<double>(nDim) *
            faces_[iFace].GetNormalVec().dot(faces_[iFace].GetFaceCenter());
    }
    std::cout << "Cell " << id_ << " volume: " << volume_ << std::endl;
}