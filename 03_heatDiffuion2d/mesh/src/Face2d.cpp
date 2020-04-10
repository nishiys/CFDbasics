#pragma once

#include "Face2d.hpp"

#include "Node2d.hpp"

Face2d::Face2d() {}
Face2d::Face2d(Node2d* pNode1, Node2d* pNode2)
    : pNode1_(pNode1), pNode2_(pNode2)
{
}
Face2d::~Face2d() {}