#include "Translation.h"
#include <iostream>
#include <iomanip>

// Our includes.
using namespace std;

Translation::Translation()
{
    this->translationId = -1;
    this->tx = 0.0;
    this->ty = 0.0;
    this->tz = 0.0;
}

Translation::Translation(int translationId, double tx, double ty, double tz)
{
    this->translationId = translationId;
    this->tx = tx;
    this->ty = ty;
    this->tz = tz;
}

Matrix4 Translation::getTranslationMatrix()
{
    double tmp[4][4] = {
                            {1,0,0,this->tx},
                            {0,1,0,this->ty},
                            {0,0,1,this->tz},
                            {0,0,0,1} 
                        };
    Matrix4 matrix = Matrix4(tmp);
    return matrix;
}   


ostream &operator<<(ostream &os, const Translation &t)
{
    os << fixed << setprecision(3) << "Translation " << t.translationId << " => [" << t.tx << ", " << t.ty << ", " << t.tz << "]";

    return os;
}