#include "TomasMoeller.hpp"

#include <iostream>

RayCaster::RayCaster(Eigen::Vector3d ray_origin,
                     Eigen::Vector3d direction_vector)
    : ray_origin_(ray_origin), direction_vector_(direction_vector)
{
}
RayCaster::~RayCaster() {}

bool RayCaster::CheckTriangleIntersection(const Triangle& triangle)
{
    Eigen::Vector3d E1 = triangle.vertex1 - triangle.vertex0;
    Eigen::Vector3d E2 = triangle.vertex2 - triangle.vertex0;

    Eigen::Vector3d T;
    T = ray_origin_ - triangle.vertex0;

    /*--- Solve system of equation by using Cramer's rule ---*/
    // tmp vectors
    Eigen::Vector3d P = direction_vector_.cross(E2);
    Eigen::Vector3d Q = T.cross(E1);
    Eigen::Vector3d b;
    b(0) = Q.dot(E2);
    b(1) = P.dot(T);
    b(2) = Q.dot(direction_vector_);
    // solution vector: x = [t, u, v]
    Eigen::Vector3d x;
    x = 1.0 / (P.dot(E1)) * b;

    /*--- Detect if it is intersected ---*/
    bool is_intersect;
    if (x(1) >= 0 && x(1) <= 1)
    {
        if (x(2) >= 0 && x(2) <= 1)
        {
            if (x(0) > 0)
            {
                is_intersect = true;
            }
            else
            {
                is_intersect = false;
            }
        }
        else
        {
            is_intersect = false;
        }
    }
    else
    {
        is_intersect = false;
    }

    return is_intersect;
}

int main()
{
    Eigen::Vector3d v0;
    v0 << 1, 0, 0;
    Eigen::Vector3d v1;
    v1 << 0, 1, 0;
    Eigen::Vector3d v2;
    v2 << 0, 0, 1;
    Triangle triangle(v0, v1, v2);

    Eigen::Vector3d origin;
    origin << 0, 0, 0;
    Eigen::Vector3d direction_vector;
    direction_vector << 1, 1, 1;
    RayCaster raycaster(origin, direction_vector);

    bool isIntersect = raycaster.CheckTriangleIntersection(triangle);

    std::cout << isIntersect << std::endl;
}
