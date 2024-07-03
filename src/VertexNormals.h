#pragma once

#include "Utils.h"
#include <unordered_map>


class VertexNormals {
public:

    std::vector<Vector3> calculateSmoothNormals(const std::vector<Vector3> &vertices,
                                                const std::vector<Triangle> &triangles);

private:
    //Calculated edges by vertex index
    //[vertex1][vertex2] = calculated edge
    std::unordered_map<uint32_t, std::unordered_map<uint32_t, Vector3>> _edges;

    //[vertex index] [triangle index] ->
    //show which triangles share a vertex
    std::vector<std::vector<uint32_t>> _vertexToTriangles;


};
