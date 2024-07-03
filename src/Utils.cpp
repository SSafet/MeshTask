#include "Utils.h"

#include <mutex>
#include <future>
#include <thread>

namespace utils {

// double calculateAngle(const Vector3& v0, const Vector3& v1) {
//     double dot = (v0.x * v1.x) + (v0.y * v1.y) + (v0.z * v1.z);
//     double len0 = length(v0);
//     double len1 = length(v1);
//     return std::acos(dot / (length(v0) * length(v1)));
// }

bool nearlyEqual(double f1, double f2)
{
    static constexpr double epsilon = 1.0e-05f;
    if (std::abs(f1 - f2) <= epsilon) {
        return true;
    }
    return std::abs(f1 - f2) <= epsilon * std::max(std::abs(f1), std::abs(f2));
}

double calculateAngle(const Vector3& v0, const Vector3& v1)
{
    double dot = (v0.x * v1.x) + (v0.y * v1.y) + (v0.z * v1.z);
    double len0 = length(v0);
    double len1 = length(v1);

    if (nearlyEqual(len0, 0.0) || nearlyEqual(len1, 0.0)) {
        return 0.0;
    }

    double denominator = len0 * len1;
    double cosAngle = dot / denominator;
    
    cosAngle = std::max(-1.0, std::min(1.0, cosAngle));
    return std::acos(cosAngle);
}

std::vector<Vector3> calculateSmoothNormals(const std::vector<Vector3>& vertices,
                                            const std::vector<Triangle>& triangles)
{
    std::vector<Vector3> faceNormals(triangles.size());
    std::vector<Vector3> vertexNormals(vertices.size(), Vector3(0, 0, 0));

    std::vector<std::vector<int>> vertexToTriangles(vertices.size());
    //Vertex to Triangles = [vertex index] [triangle index] -> show which triangles share a vertex

    // Calculate face normals
    for (size_t i = 0; i < triangles.size(); ++i) {
        const Triangle& tri = triangles[i];
        const Vector3& v0 = vertices[tri.vertices[0]];
        const Vector3& v1 = vertices[tri.vertices[1]];
        const Vector3& v2 = vertices[tri.vertices[2]];

        Vector3 e0 = v1 - v0;//edges
        Vector3 e1 = v2 - v0;//edges
        
        faceNormals[i] = utils::normalize(utils::crossProduct(e0, e1));

        vertexToTriangles[tri.vertices[0]].push_back(i);
        vertexToTriangles[tri.vertices[1]].push_back(i);
        vertexToTriangles[tri.vertices[2]].push_back(i);
    }

    // Calculate vertex normals
    for (size_t currentVertexIndex = 0; currentVertexIndex < vertices.size(); ++currentVertexIndex) {
        Vector3 normal;
        for (int t : vertexToTriangles[currentVertexIndex]) {
            const auto& tri = triangles[t];
            int index = 0;
            for (int j = 0; j < 3; ++j) {
                if (tri.vertices[j] == currentVertexIndex) {
                    index = j;
                    break;
                }
            }
            Vector3 v0 = vertices[tri.vertices[(index + 1) % 3]] - vertices[tri.vertices[index]];
            Vector3 v1 = vertices[tri.vertices[(index + 2) % 3]] - vertices[tri.vertices[index]];
            double angle = calculateAngle(v0, v1);
            
            normal += faceNormals[t] * angle;
        }
        utils::normalize_r(normal);
        vertexNormals[currentVertexIndex] = std::move(normal);
    }

    return vertexNormals;
}

double calculateTriangleArea(const Vector3& v0, const Vector3& v1, const Vector3& v2)
{
    auto c = crossProduct(v1 - v0, v2 - v0);
    double area = 0.5 * std::sqrt(c.x * c.x + c.y * c.y + c.z * c.z);//TODO check error
    return area;
}

std::tuple<double, double, double> calculateTriangleAreas(const std::vector<Vector3>& vertices,
                                                          const std::vector<Triangle>& triangles)
{
    double maxArea{std::numeric_limits<double>::min()};
    double smallestArea{std::numeric_limits<double>::max()};
    double totalArea{0.f};
    for (const auto& tri : triangles) {
        double area = calculateTriangleArea(vertices[tri.vertices[0]], vertices[tri.vertices[1]], vertices[tri.vertices[2]]);
        totalArea += area;
        if (area > maxArea) {
            maxArea = area;
        }
        if (area > 0.f && area < smallestArea) {
            smallestArea = area;
        }
    }
    return {maxArea, smallestArea, totalArea / triangles.size()};
}

MeshStatistics calculateMeshStatistics(const std::vector<Vector3>& vertices,
                                       const std::vector<Triangle>& triangles)
{
    MeshStatistics result;

    auto computeArea = [&vertices, &triangles] (size_t start, size_t end) {
        MeshStatistics res;
        for (size_t i = start; i < end; ++i) {
            auto area = calculateTriangleArea(vertices[triangles[i].vertices[0]],
                                              vertices[triangles[i].vertices[1]],
                                              vertices[triangles[i].vertices[2]]);

            res.avgArea += area;
            if (area > res.maxArea) {
                res.maxArea = area;
            }
            if (area > 0. && area < res.minArea) {
                res.minArea = area;
            }
        }
        res.avgArea /= (end - start);
        return res;
    };

    size_t numThreads = std::min(std::thread::hardware_concurrency(), static_cast<uint32_t>(triangles.size()));
    std::vector<std::future<MeshStatistics>> futures;
    size_t chunkSize = triangles.size() / numThreads;

    for (size_t i = 0; i < numThreads; ++i) {
        size_t start = i * chunkSize;
        size_t end = (i == (numThreads - 1)) ? triangles.size() : start + chunkSize;
        futures.push_back(std::async(std::launch::async, computeArea, start, end));
    }

    for (auto& future : futures) {
        MeshStatistics res = future.get();
        result.avgArea += res.avgArea;
        result.minArea = std::min(result.minArea, res.minArea);
        result.maxArea = std::max(result.maxArea, res.maxArea);
    }

    result.avgArea /= numThreads;
    return result;
}

void buildAdjacencyList(Mesh& mesh) {
    mesh.adjacencyList.resize(mesh.vertices.size());
    for (size_t i = 0; i < mesh.triangles.size(); ++i) {
        const auto& tri = mesh.triangles[i];
        for (int v : tri.vertices) {
            mesh.adjacencyList[v].push_back(i);
        }
    }
}


std::vector<Vector3> calculateFaceNormals(const Mesh& mesh) {
    std::vector<Vector3> faceNormals(mesh.triangles.size());
    for (size_t i = 0; i < mesh.triangles.size(); ++i) {
        const auto& tri = mesh.triangles[i];
        Vector3 v0 = mesh.vertices[tri.vertices[0]];
        Vector3 v1 = mesh.vertices[tri.vertices[1]];
        Vector3 v2 = mesh.vertices[tri.vertices[2]];
        faceNormals[i] = normalize(crossProduct({v1.x - v0.x, v1.y - v0.y, v1.z - v0.z}, {v2.x - v0.x, v2.y - v0.y, v2.z - v0.z}));
    }
    return faceNormals;
}

std::vector<Vector3> calculateSmoothNormals2(const Mesh& mesh, const std::vector<Vector3>& faceNormals) {
    std::vector<Vector3> vertexNormals(mesh.vertices.size(), {0, 0, 0});

    for (size_t vertexIndex = 0; vertexIndex < mesh.adjacencyList.size(); ++vertexIndex) {
        Vector3 normalSum = std::accumulate(mesh.adjacencyList[vertexIndex].begin(), mesh.adjacencyList[vertexIndex].end(), Vector3{0, 0, 0},
                                            [&faceNormals](const Vector3& sum, int triIndex) {
                                                return sum + faceNormals[triIndex];
                                            });
        vertexNormals[vertexIndex] = normalize(normalSum);
    }
    
    return vertexNormals;
}

std::vector<Vector3> calculateWeightedSmoothNormals(const Mesh& mesh, const std::vector<Vector3>& faceNormals) {
    std::vector<Vector3> vertexNormals(mesh.vertices.size(), {0, 0, 0});

    for (size_t vertexIndex = 0; vertexIndex < mesh.adjacencyList.size(); ++vertexIndex) {
        Vector3 normalSum = {0, 0, 0};
        double totalWeight = 0.0;

        for (int triIndex : mesh.adjacencyList[vertexIndex]) {
            const Triangle& tri = mesh.triangles[triIndex];

            int v0Index = tri.vertices[0];
            int v1Index = tri.vertices[1];
            int v2Index = tri.vertices[2];

            Vector3 v0 = mesh.vertices[v0Index];
            Vector3 v1 = mesh.vertices[v1Index];
            Vector3 v2 = mesh.vertices[v2Index];

            Vector3 edge1, edge2;
            if (v0Index == vertexIndex) {
                edge1 = v1 - v0;
                edge2 = v2 - v0;
            } else if (v1Index == vertexIndex) {
                edge1 = v0 - v1;
                edge2 = v2 - v1;
            } else {
                edge1 = v0 - v2;
                edge2 = v1 - v2;
            }

            double angle = calculateAngle(edge1, edge2);
            normalSum += faceNormals[triIndex] * angle;
            totalWeight += angle;
        }

        if (totalWeight > 0) {
            vertexNormals[vertexIndex] = normalize(normalSum / totalWeight);
        }
    }

    return vertexNormals;
}


} //utils
