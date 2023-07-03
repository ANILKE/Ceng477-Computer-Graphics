#include "Scaling.h"
#include <iostream>
#include <iomanip>
// Our includes
#include "Matrix4.h"
using namespace std;

Scaling::Scaling() {}

Scaling::Scaling(int scalingId, double sx, double sy, double sz)
{
    this->scalingId = scalingId;
    this->sx = sx;
    this->sy = sy;
    this->sz = sz;
}
Matrix4 Scaling::getScaleMatrix()
{
    double tmp[4][4] = {
                            {this->sx,0,0,0},
                            {0,this->sy,0,0},
                            {0,0,this->sz,0},
                            {0,0,0,1}
                       };
    Matrix4 matrix = Matrix4(tmp);
    return matrix;
}


ostream &operator<<(ostream &os, const Scaling &s)
{
    os << fixed << setprecision(3) << "Scaling " << s.scalingId << " => [" << s.sx << ", " << s.sy << ", " << s.sz << "]";

    return os;
}
