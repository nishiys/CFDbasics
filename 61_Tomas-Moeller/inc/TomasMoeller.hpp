#pragma once

#include <Eigen/Dense>

struct Triangle
{
public:
    Triangle(Eigen::Vector3d a_vertex0, Eigen::Vector3d a_vertex1,
             Eigen::Vector3d a_vertex2)
        : vertex0(a_vertex0), vertex1(a_vertex1), vertex2(a_vertex2)
    {
    }

    Eigen::Vector3d vertex0;
    Eigen::Vector3d vertex1;
    Eigen::Vector3d vertex2;
};

class RayCaster
{
public:
    RayCaster(Eigen::Vector3d ray_origin, Eigen::Vector3d direction_vector);
    ~RayCaster();

    /* Tomas Moeller algorithm */
    bool CheckTriangleIntersection(const Triangle& triangle);

private:
    Eigen::Vector3d ray_origin_;
    Eigen::Vector3d direction_vector_;
};
