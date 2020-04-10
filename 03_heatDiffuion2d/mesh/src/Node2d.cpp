#include "Node2d.hpp"

Node2d::Node2d() {}

Node2d::Node2d(unsigned int id, double x, double y) : id_(id), coords_{x, y} {}
Node2d::~Node2d() {}
