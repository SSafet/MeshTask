#include "JsonUtils.h"

#include <vector>

namespace utils
{
    namespace json
    {

        std::vector<Vector3> parseVertices(const std::vector<double> &vertices_flat)
        {
            std::vector<Vector3> vertices;
            for (size_t i = 0; i < vertices_flat.size(); i += 3)
            {
                vertices.push_back({vertices_flat[i],
                                    vertices_flat[i + 1],
                                    vertices_flat[i + 2]});
            }
            return vertices;
        }

        std::vector<Triangle> parseTriangles(const std::vector<int> &triangles_indices)
        {
            std::vector<Triangle> triangles;
            for (size_t i = 0; i < triangles_indices.size(); i += 3)
            {
                triangles.push_back({static_cast<uint32_t>(triangles_indices[i]),
                                     static_cast<uint32_t>(triangles_indices[i + 1]),
                                     static_cast<uint32_t>(triangles_indices[i + 2])});
            }
            return triangles;
        }
    } // json
} // utils
