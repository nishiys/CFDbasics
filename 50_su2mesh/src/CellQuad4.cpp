
#include "CellQuad4.hpp"

CellQuad4::CellQuad4() {}
CellQuad4::CellQuad4(unsigned int id, Node2d* pNode1, Node2d* pNode2,
                     Node2d* pNode3, Node2d* pNode4)
    : pNode1_(pNode1),
      pNode2_(pNode2),
      pNode3_(pNode3),
      pNode4_(pNode4),
      id_(id)
{
}

CellQuad4::~CellQuad4() {}
