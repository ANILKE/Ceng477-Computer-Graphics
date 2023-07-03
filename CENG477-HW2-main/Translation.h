#ifndef __TRANSLATION_H__
#define __TRANSLATION_H__

#include <iostream>

using namespace std;
#include "Matrix4.h"
class Translation
{
public:
    int translationId;
    double tx, ty, tz;

    Translation();
    Translation(int translationId, double tx, double ty, double tz);
    friend ostream &operator<<(ostream &os, const Translation &t);

    // Our Helper Methods.
    Matrix4 getTranslationMatrix();




};

#endif