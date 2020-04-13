#include "CellQuad4.hpp"

#include <iostream>

#include "Node2d.hpp"

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
    CalcVectorToFaces();
    CalcVectorToNodes();
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
            1.0 / static_cast<double>(nNode) * pNodes_[iNode]->GetCoords()(0);
        centroid_(1) +=
            1.0 / static_cast<double>(nNode) * pNodes_[iNode]->GetCoords()(1);
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

void CellQuad4::CalcVectorToFaces()
{
    for (size_t iFace = 0; iFace < faces_.size(); ++iFace)
    {
        vectorToNeighbors_[iFace] = faces_[iFace].GetFaceCenter() - centroid_;
    }
}

void CellQuad4::CalcVectorToNodes()
{
    for (size_t iNode = 0; iNode < 4; ++iNode)
    {
        vectorToNodes_[iNode] = pNodes_[iNode]->GetCoords() - centroid_;
    }
}

/*---------------------------------------------------------*/
/* Global Methods */
/*---------------------------------------------------------*/
unsigned int CellQuad4::GetID() const { return id_; }
Node2d* CellQuad4::GetNode(unsigned int index) const { return pNodes_[index]; }
Face2d* CellQuad4::Face(unsigned int index) { return &faces_[index]; }
void CellQuad4::SetNeighborPtr(unsigned int index, CellQuad4* pCellQuad)
{
    pNeighbors_[index] = pCellQuad;
}
CellQuad4* CellQuad4::GetNeighborPtr(unsigned int index) const
{
    return pNeighbors_[index];
}
Eigen::Vector2d CellQuad4::GetCentroid() const { return centroid_; }
void CellQuad4::SetVectorToNeighbor(unsigned int index)
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
Eigen::Vector2d CellQuad4::GetVectorToNeighbor(unsigned int index) const
{
    return vectorToNeighbors_[index];
}
Eigen::Vector2d CellQuad4::GetNormalVectorToBoundary(unsigned int index) const
{
    return normalvectorToBoundaries_[index];
}
void CellQuad4::SetMag_N1Vectors(unsigned int index)
{
    if (pNeighbors_[index] != nullptr)
    {
        mag_n1Vectors_[index] =
            1.0 / vectorToNeighbors_[index].norm() *
            faces_[index].GetNormalVec().dot(vectorToNeighbors_[index]);
    }
}
double CellQuad4::GetMag_N1Vectors(unsigned int index) const
{
    return mag_n1Vectors_[index];
}
Eigen::Vector2d CellQuad4::GetvectorToFace(unsigned int index) const
{
    return vectorToFaces_[index];
}
Eigen::Vector2d CellQuad4::GetVectorToNode(unsigned int index) const
{
    return vectorToNodes_[index];
}

double CellQuad4::GetVolume() const { return volume_; }