#include "VertexNormals.h"


std::vector<Vector3> VertexNormals::calculateSmoothNormals(const std::vector<Vector3> &vertices,
                                                           const std::vector<Triangle> &triangles)
{
    //Calculate edges
    for (size_t i = 0; i < triangles.size(); ++i) {
        const Triangle& tri = triangles[i];
        if(auto it = _edges[tri.vertices[1]].find(tri.vertices[0]); it == _edges.begin()->second.end()) {
            it->second = vertices[tri.vertices[1]] - vertices[tri.vertices[0]];
        }
    
        if(auto it = _edges[tri.vertices[2]].find(tri.vertices[0]); it == _edges.begin()->second.end()) {
            it->second = vertices[tri.vertices[2]] - vertices[tri.vertices[0]];
        }
    }
    /////////////////



    return std::vector<Vector3>();
}


