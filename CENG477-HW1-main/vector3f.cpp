#include "vector3f.h"

vector3f::vector3f(float a, float b, float c) {x= a;   y = b; z = c;}

float vector3f::dotProduct(vector3f a)
{
    return (a.x*x)+(a.y*y)+(a.z*z);
}
vector3f vector3f::normalize()
{
    return *this*(1.0/sqrt(dotProduct(*this)));
}
vector3f vector3f::crossProduct(vector3f second){
    vector3f res;
    res.x = (y*second.z-z*second.y);
    res.y = (z*second.x-x*second.z);
    res.z = x*second.y-y*second.x;
    return res;
}