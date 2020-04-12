#pragma once

#include <array>
#include <string>

#include "Variable.hpp"

class Node2d
{
public:
    Node2d();
    Node2d(unsigned int id, double x, double y);
    ~Node2d();

    inline unsigned int GetID() const { return id_; };
    inline double GetX() const { return coords_[0]; };
    inline double GetY() const { return coords_[1]; };

    Variable nodevar;

private:
    unsigned int id_;
    std::array<double, 2> coords_;
};