#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>
#include "Geometry.h"
#include "Vector3.h"
#include "Utils.h"

using json = nlohmann::json;


struct Timer {
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    std::chrono::high_resolution_clock::duration duration;
    std::string asd;

    Timer(std::string q) : asd(q) {
        start = std::chrono::high_resolution_clock::now();
    }

    ~Timer() {
        duration = std::chrono::high_resolution_clock::now() - start;
        std::cout << asd + " Elapsed time: " << std::chrono::duration_cast<std::chrono::microseconds>(duration).count() << " us." << std::endl;
    }
};

std::vector<Vector3> calculateVertexNormals(
                    const std::vector<Vector3>& vertices,
                    const std::vector<Triangle>& triangles)
{
    std::vector<Vector3> vertexNormals(vertices.size(), Vector3());

    for (const auto& tri : triangles) {
        Vector3 p1 = vertices[tri.vertices[0]];
        Vector3 p2 = vertices[tri.vertices[1]];
        Vector3 p3 = vertices[tri.vertices[2]];

        Vector3 normal = utils::normalize(utils::crossProduct(p2 - p1, p3 - p1));
        vertexNormals[tri.vertices[0]] += normal;
        vertexNormals[tri.vertices[1]] += normal;
        vertexNormals[tri.vertices[2]] += normal;
    }

    for(auto& v : vertexNormals) {
        utils::normalize_r(v);
    }
    return vertexNormals;
}

void exportToFile(const std::string& filename,
                  const std::vector<Vector3>& vertices,
                  const std::vector<Vector3>& normals)
{
    std::ofstream outFile(filename);
    for (size_t i = 0; i < vertices.size(); ++i) {
        outFile << vertices[i].x << " " << vertices[i].y << " " << vertices[i].z << " ";
        outFile << normals[i].x << " " << normals[i].y << " " << normals[i].z << "\n";
    }
    outFile.close();
}

int main() {
    //Pass the file from terminal or UI
    std::ifstream file("C:/Users/safet/Desktop/task_input/teapot.json");
    json mesh_json = json::parse(file);
    /////////////////////////////////////////////////////////////////////

    std::vector<double> vertices_flat = mesh_json["geometry_object"]["vertices"].get<std::vector<double>>();
    std::vector<int> triangles_indices = mesh_json["geometry_object"]["triangles"].get<std::vector<int>>();

    std::vector<Vector3> vertices = Geometry::parseVertices(vertices_flat);
    std::vector<Triangle> triangles = Geometry::parseTriangles(triangles_indices);

    ///////////////////////////////////////////////////////////////////////////////
    {
        utils::Mesh mesh;
        mesh.vertices = vertices;
        mesh.triangles = triangles;
        // Load mesh vertices and triangles here
        std::vector<Vector3> smoothNormals;
        {
            Timer t0("Adj list");
            utils::buildAdjacencyList(mesh);
            std::vector<Vector3> faceNormals = utils::calculateFaceNormals(mesh);
            smoothNormals = utils::calculateSmoothNormals2(mesh, faceNormals);
        }
        exportToFile("adj_smooth_normals.txt", vertices, smoothNormals);
    }
    ///////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////
    {
        utils::Mesh mesh;
        mesh.vertices = vertices;
        mesh.triangles = triangles;
        // Load mesh vertices and triangles here
        std::vector<Vector3> smoothNormals;
        {
            Timer t0("Adj WS list");
            utils::buildAdjacencyList(mesh);
            std::vector<Vector3> faceNormals = utils::calculateFaceNormals(mesh);
            smoothNormals = utils::calculateWeightedSmoothNormals(mesh, faceNormals);
        }
        exportToFile("adj_WS_smooth_normals.txt", vertices, smoothNormals);
    }
    ///////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////
    {
        std::vector<Vector3> vertexNormals;
        {
            Timer t1("Vertex Normals");
            vertexNormals = calculateVertexNormals(vertices, triangles);
        }
        exportToFile("normals.txt", vertices, vertexNormals);
    }
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    {
        std::vector<Vector3> smoothNormals;
        {
            Timer t2("Smooth Vertex Normals");
            smoothNormals = utils::calculateSmoothNormals(vertices, triangles);
        }
        exportToFile("smooth_normals.txt", vertices, smoothNormals);
    }
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    {
        std::tuple<double, double, double> m;
        {
            Timer t3("Areas");
            m = utils::calculateTriangleAreas(vertices, triangles);
        }
        std::cout << "Largest: " << std::get<0>(m) << " Smallest: " << std::get<1>(m) << " Average: " << std::get<2>(m) << std::endl;
    }
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    {
        utils::MeshStatistics ms;
        {
            Timer t4("Areas j12");
            ms = utils::calculateMeshStatistics(vertices, triangles);
        }
        std::cout << "Largest: " << ms.maxArea << " Smallest: " << ms.minArea << " Average: " << ms.avgArea << std::endl;
    }
    ///////////////////////////////////////////////////////////////////////////////

    // Vector3 random_point{ 0.1f, 0.1f, 1.f };
    // bool inside = Geometry::isPointInsideMesh(random_point, vertices, triangles);
    // std::cout << "Is the point (" << random_point.x << ", " << random_point.y << ", " << random_point.z << ") inside the mesh? "
    //     << (inside ? "Yes" : "No") << std::endl;

    return 0;
}
