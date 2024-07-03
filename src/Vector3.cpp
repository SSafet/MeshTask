#include "Vector3.h"

Vector3::Vector3(double x, double y, double z)
: x(x)
, y(y)
, z(z)
{}

Vector3 Vector3::operator+(const Vector3& other) const
{
    return Vector3(x + other.x, y + other.y, z + other.z);
}

Vector3& Vector3::operator+=(const Vector3& other)
{
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
}

Vector3 Vector3::operator-(const Vector3& other) const
{
    return Vector3(x - other.x, y - other.y, z - other.z);
}

Vector3 Vector3::operator*(double scalar) const
{
    return Vector3(x * scalar, y * scalar, z * scalar);
}

Vector3 Vector3::operator/(double scalar) const
{
    return {x / scalar, y / scalar, z / scalar};
}