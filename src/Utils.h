#pragma once

#include "Vector3.h"
#include "Triangle.h"

#include <numeric>
#include <unordered_map>

namespace utils
{

bool nearlyEqual(double f1, double f2);
double calculateAngle(const Vector3& v0, const Vector3& v1);

std::vector<Vector3> calculateSmoothNormals(const std::vector<Vector3>& vertices,
                                            const std::vector<Triangle>& triangles);


double calculateTriangleArea(const Vector3& v0, const Vector3& v1, const Vector3& v2);

std::tuple<double, double, double> calculateTriangleAreas(const std::vector<Vector3>& vertices,
                                                       const std::vector<Triangle>& triangles);

inline Vector3 crossProduct(const Vector3& l, const Vector3& r)
{
    return Vector3(
        (l.y * r.z) - (l.z * r.y),
        (l.z * r.x) - (l.x * r.z),
        (l.x * r.y) - (l.y * r.x)
    );
}

inline bool isNearlyZero(double value) {
    return std::abs(value) < std::numeric_limits<double>::epsilon();
}

inline double length(const Vector3& v)
{
    double lengthSquared = v.x * v.x + v.y * v.y + v.z * v.z;
    if (isNearlyZero(lengthSquared)) {
        return 0.0;
    }
    double length = std::sqrt(lengthSquared);
    if (std::isnan(length) || std::isinf(length)) {
        return 0.0;
    }
    return length;
}

inline void normalize_r(Vector3& v)
{
    double l = length(v);
    if (l > 0) {
        v.x /= l;
        v.y /= l;
        v.z /= l;
    }
}

inline Vector3 normalize(Vector3 v)
{
    double l = length(v);
    if (l > 0) {
        v.x /= l;
        v.y /= l;
        v.z /= l;
    }
    return v;
}

struct Mesh
{
    std::vector<Vector3> vertices;
    std::vector<Triangle> triangles;
    std::vector<std::vector<int>> adjacencyList;
};

void buildAdjacencyList(Mesh& mesh);
std::vector<Vector3> calculateFaceNormals(const Mesh& mesh);
std::vector<Vector3> calculateSmoothNormals2(const Mesh& mesh, const std::vector<Vector3>& faceNormals);
std::vector<Vector3> calculateWeightedSmoothNormals(const Mesh& mesh, const std::vector<Vector3>& faceNormals);

struct MeshStatistics
{
    double minArea{std::numeric_limits<double>::max()};
    double maxArea{std::numeric_limits<double>::min()};
    double avgArea{0.0};
};

MeshStatistics calculateMeshStatistics(const std::vector<Vector3>& vertices,
                                       const std::vector<Triangle>& triangles);

} // utils
