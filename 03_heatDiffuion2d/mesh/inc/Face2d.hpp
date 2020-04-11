#pragma once

#include <Eigen/Core>
#include <string>

#include "Node2d.hpp"
#include "Variable.hpp"

class Face2d
{
public:
    Face2d();
    Face2d(Node2d* pNode1, Node2d* pNode2);
    ~Face2d();

    void FlipNormalVec();

    inline Node2d* GetNode1() const { return pNode1_; };
    inline Node2d* GetNode2() const { return pNode2_; };
    // inline unsigned int GetID() const { return id_; };
    inline Eigen::VectorXd GetNormalVec() const { return normalvec_; };
    inline Eigen::VectorXd GetFaceCenter() const { return facecenter_; };
    inline double GetArea() const { return area_; };

    inline void SetTag(std::string tagname) { tag_ = tagname; };
    inline std::string GetTag() const { return tag_; };

    Variable facevar;

private:
    unsigned int id_;
    Node2d* pNode1_;
    Node2d* pNode2_;

    //! Tags set in .su2 file
    std::string tag_;
    //! Dirichlet, Neumann, Interior
    std::string status_;

    Eigen::VectorXd normalvec_;
    Eigen::VectorXd facecenter_;
    double area_;

    void CalcNormalVec();
    void CalcFaceCenter();
    void CalcArea();
};