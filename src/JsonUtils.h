#pragma once

#include "Vector3.h"
#include "Triangle.h"

namespace utils {
namespace json {

std::vector<Vector3> parseVertices(const std::vector<double>& vertices_flat);
std::vector<Triangle> parseTriangles(const std::vector<int>& triangles_indices);

}//json
}//utils
