#pragma once

#include "Node2d.hpp"

class Face2d
{
public:
    Face2d();
    Face2d(Node2d* pNode1, Node2d* pNode2);
    ~Face2d();

private:
    unsigned int id_;
    Node2d* pNode1_;
    Node2d* pNode2_;
    std::string status_;
};