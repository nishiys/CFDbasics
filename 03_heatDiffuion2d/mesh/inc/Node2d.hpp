#pragma once

#include <array>
#include <string>
class Node2d
{
public:
    Node2d();
    Node2d(unsigned int id, double x, double y);
    ~Node2d();

    void SetStatus(std::string status);

    inline unsigned int GetID() const { return id_; };
    inline double GetX() const { return coords_[0]; };
    inline double GetY() const { return coords_[1]; };

private:
    unsigned int id_;
    std::array<double, 2> coords_;
    //! boundary / interior
    std::string status_;
};