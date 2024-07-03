#include "Geometry.h"
#include <cmath>

std::vector<Vector3> Geometry::parseVertices(const std::vector<double>& vertices_flat) {
    std::vector<Vector3> vertices;
    for (size_t i = 0; i < vertices_flat.size(); i += 3) {
        vertices.push_back({ vertices_flat[i], vertices_flat[i + 1], vertices_flat[i + 2] });
    }
    return vertices;
}

std::vector<Triangle> Geometry::parseTriangles(const std::vector<int>& triangles_indices) {
    std::vector<Triangle> triangles;
    for (size_t i = 0; i < triangles_indices.size(); i += 3) {
        triangles.push_back({static_cast<uint32_t>(triangles_indices[i]),
                             static_cast<uint32_t>(triangles_indices[i + 1]),
                             static_cast<uint32_t>(triangles_indices[i + 2]) });
    }
    return triangles;
}

std::array<double, 3> Geometry::calculateFaceNormal(const Vector3& v0, const Vector3& v1, const Vector3& v2) {
    std::array<double, 3> edge1 = { v1.x - v0.x, v1.y - v0.y, v1.z - v0.z };
    std::array<double, 3> edge2 = { v2.x - v0.x, v2.y - v0.y, v2.z - v0.z };
    std::array<double, 3> normal = {
        edge1[1] * edge2[2] - edge1[2] * edge2[1],
        edge1[2] * edge2[0] - edge1[0] * edge2[2],
        edge1[0] * edge2[1] - edge1[1] * edge2[0]
    };
    double length = std::sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
    return { normal[0] / length, normal[1] / length, normal[2] / length };
}

bool Geometry::isPointInsideMesh(const Vector3& point, const std::vector<Vector3>& vertices, const std::vector<Triangle>& triangles) {
    std::array<double, 3> ray_direction = { 1.0f, 0.0f, 0.0f };
    int intersection_count = 0;

    for (const auto& tri : triangles) {
        const Vector3& v0 = vertices[tri.vertices[0]];
        const Vector3& v1 = vertices[tri.vertices[1]];
        const Vector3& v2 = vertices[tri.vertices[2]];

        std::array<double, 3> edge1 = { v1.x - v0.x, v1.y - v0.y, v1.z - v0.z };
        std::array<double, 3> edge2 = { v2.x - v0.x, v2.y - v0.y, v2.z - v0.z };

        std::array<double, 3> h = {
            ray_direction[1] * edge2[2] - ray_direction[2] * edge2[1],
            ray_direction[2] * edge2[0] - ray_direction[0] * edge2[2],
            ray_direction[0] * edge2[1] - ray_direction[1] * edge2[0]
        };

        double a = edge1[0] * h[0] + edge1[1] * h[1] + edge1[2] * h[2];
        if (std::fabs(a) < 1e-8) continue;

        double f = 1.0f / a;
        std::array<double, 3> s = { point.x - v0.x, point.y - v0.y, point.z - v0.z };
        double u = f * (s[0] * h[0] + s[1] * h[1] + s[2] * h[2]);
        if (u < 0.0f || u > 1.0f) continue;

        std::array<double, 3> q = {
            s[1] * edge1[2] - s[2] * edge1[1],
            s[2] * edge1[0] - s[0] * edge1[2],
            s[0] * edge1[1] - s[1] * edge1[0]
        };

        double v = f * (ray_direction[0] * q[0] + ray_direction[1] * q[1] + ray_direction[2] * q[2]);
        if (v < 0.0f || u + v > 1.0f) continue;

        double t = f * (edge2[0] * q[0] + edge2[1] * q[1] + edge2[2] * q[2]);
        if (t > 1e-8) {
            intersection_count++;
        }
    }

    return (intersection_count % 2) == 1;
}
