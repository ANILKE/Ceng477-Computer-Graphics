#define _USE_MATH_DEFINES
#include "Rotation.h"
#include <iostream>
#include <iomanip>

#include "Matrix4.h"
#include "Vec3.h"
#include "Helpers.h"
#include "cmath"

using namespace std;

Rotation::Rotation() {}

Rotation::Rotation(int rotationId, double angle, double x, double y, double z)
{
    this->rotationId = rotationId;
    this->angle = angle;
    this->ux = x;
    this->uy = y;
    this->uz = z;
}


int findMinIdx(Vec3 u)
{
    double min = 1.7976931348623158e+308;
    int idx = -1;

    if(abs(u.getElementAt(0)) < min)
    {
        idx = 0;
        min = u.getElementAt(0);
    }
    if(abs(u.getElementAt(1)) < min)
    {
        idx = 1;
        min = u.getElementAt(1);
    }

    if(abs(u.getElementAt(2)) < min)
    {
        idx = 2;
        min = u.getElementAt(2);
    }

    return idx;
}

Vec3 calculateV(Vec3 u)
{
    int idx = findMinIdx(u);
    if(idx == 0)
    {
        double vy = -u.z;
        double vz = u.y;
        double vx = 0;

        Vec3 v(vx,vy,vz,-1);
        return v;
    }
    else if(idx == 1)
    {
        double vy = 0;
        double vx = -u.z;
        double vz = u.x;

        Vec3 v(vx,vy,vz,-1);
        return v;
    }
    else if(idx == 2)
    {
        double vz = 0;
        double vx = -u.y;
        double vy = u.x;

        Vec3 v(vx,vy,vz,-1);
        
        return v;
    }

}



Matrix4 Rotation::getRotationMatrix()
{
    // Orthonormal basis method.

    // Calculate u vector.
    Vec3 u(this->ux,this->uy,this->uz,-1);
    u = normalizeVec3(u);

    // Calculate v vector.
    Vec3 v = calculateV(u);


    // Calculate w vector.
    Vec3 w = crossProductVec3(u,v);
    v = normalizeVec3(v);
    w = normalizeVec3(w);
    

    
    //WARNING ROTATION ABOUT MAJOR AXIS MATRIS CAN BE WRONG, NOT SURE.
    double radian = this->angle * (M_PI / 180);
    double tmp1[4][4] = {
                            {1,          0,           0,0},
                            {0,cos(radian),-sin(radian),0},
                            {0,sin(radian), cos(radian),0},
                            {0,          0,           0,1}
                       };

    Matrix4 rotationMatrix(tmp1);


    double tmp2[4][4] = {
                            {u.x,u.y,u.z,0},
                            {v.x,v.y,v.z,0},
                            {w.x,w.y,w.z,0},
                            {  0,  0,  0,1}
                        };
    Matrix4 matrixM(tmp2);


    double tmp3[4][4] = {
                            {u.x,v.x,w.x,0},
                            {u.y,v.y,w.y,0},
                            {u.z,v.z,w.z,0},
                            {  0,  0,  0,1}
                        };
    Matrix4 inverseMatrixM(tmp3);


    Matrix4 returnValue = multiplyMatrixWithMatrix(inverseMatrixM,multiplyMatrixWithMatrix(rotationMatrix,matrixM));
    return returnValue;

}






ostream &operator<<(ostream &os, const Rotation &r)
{
    os << fixed << setprecision(3) << "Rotation " << r.rotationId << " => [angle=" << r.angle << ", " << r.ux << ", " << r.uy << ", " << r.uz << "]";

    return os;
}