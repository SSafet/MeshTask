#pragma once

#include <cmath>
#include <vector>

struct Vector3 {
    Vector3() = default;
    Vector3(double x, double y, double z);

    Vector3& operator+=(const Vector3& other);

    Vector3 operator+(const Vector3& other) const;
    Vector3 operator-(const Vector3& other) const;
    Vector3 operator*(double scalar) const;
    Vector3 operator/(double scalar) const;

    double x{0.0};
    double y{0.0};
    double z{0.0};
};
