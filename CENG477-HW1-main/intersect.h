#ifndef HWNEW_INTERSECT_H
#define HWNEW_INTERSECT_H
#include "vector3f.h"

class intersect {
    public:
        float distance;
        vector3f normal;
        int material_id;
        vector3f intersection_point;

        intersect( float =-1, vector3f= {0,0,0}, int =-1,vector3f= {0,0,0} );
};


#endif //HWNEW_INTERSECT_H
