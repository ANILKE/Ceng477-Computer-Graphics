#ifndef HWNEW_VECTOR3F_H
#define HWNEW_VECTOR3F_H


#include <math.h>
#include "parser.h"
class vector3f {
    public:


    float x;
    float y;
    float z;
    vector3f(float =0, float =0, float=0);
    vector3f  operator + (vector3f const &obj) {
        vector3f res;
        res.x= x+obj.x;
        res.y= y+obj.y;
        res.z= z+obj.z;
        return res;
    }
    vector3f operator - (vector3f const &obj) {
        vector3f res;
        res.x= x-obj.x;
        res.y= y-obj.y;
        res.z= z-obj.z;
        return res;
    }

    vector3f operator * (vector3f const &obj) {
        vector3f res;
        if(x!=0){
            res.x= x*obj.x;
        }
        if(y!=0){
            res.y= y*obj.y;
        }
        if(z!=0){
            res.z= z*obj.z;
        }

        return res;
    }
    vector3f operator * (float a) {
        vector3f res;
        res.x= x*a;
        res.y= y*a;
        res.z= z*a;
        return res;
    }
    vector3f operator / (float a) {
        vector3f res;
        res.x= x/a;
        res.y= y/a;
        res.z= z/a;
        return res;
    }

    float dotProduct(vector3f);
    vector3f normalize();
    vector3f crossProduct(vector3f);

};


#endif //HWNEW_VECTOR3F_H
