#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "Vector3.h"
#include "Triangle.h"
#include <vector>
#include <array>

class Geometry
{
public:
    static std::vector<Vector3> parseVertices(const std::vector<double>& vertices_flat);
    static std::vector<Triangle> parseTriangles(const std::vector<int>& triangles_indices);
    static std::array<double, 3> calculateFaceNormal(const Vector3& v0, const Vector3& v1, const Vector3& v2);
    static bool isPointInsideMesh(const Vector3& point, const std::vector<Vector3>& vertices, const std::vector<Triangle>& triangles);
};

#endif // GEOMETRY_H
