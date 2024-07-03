

1. Open terminal
2. Go to project folder
3. mkdir build
4. cd build
5. conan install .. --output-folder=. --build=missing
6. cmake .. -DCMAKE_TOOLCHAIN_FILE="conan_toolchain.cmake" -DCMAKE_BUILD_TYPE=Debug
7. cmake --build .