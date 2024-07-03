cmake_minimum_required(VERSION 3.24)
project(MeshTask)

set(CMAKE_CXX_STANDARD 17)

# Add source files
file(GLOB_RECURSE SOURCES
     CONFIGURE_DEPENDS
     "${CMAKE_SOURCE_DIR}/src/*.cpp"
     "${CMAKE_SOURCE_DIR}/src/*.h"
     "${CMAKE_SOURCE_DIR}/src/*.hpp")

find_package(nlohmann_json REQUIRED)

# Add executable
add_executable(${PROJECT_NAME} ${SOURCES})

# Link nlohmann_json library
target_link_libraries(${PROJECT_NAME} PRIVATE nlohmann_json::nlohmann_json)
