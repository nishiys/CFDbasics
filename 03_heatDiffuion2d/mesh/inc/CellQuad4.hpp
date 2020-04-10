#pragma once

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

private:
    unsigned int id_;
    Node2d* pNode1_;
    Node2d* pNode2_;
    Node2d* pNode3_;
    Node2d* pNode4_;
};