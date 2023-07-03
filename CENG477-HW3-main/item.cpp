#include "item.h"




item::item(int i,int j,glm::vec3 color,glm::mat4 T,glm::mat4 R,glm::mat4 S,float initScaleFact)
{
    this->color = color;
    this-> gridLocI = i;
    this-> gridLocJ = j;
    this->T = T;
    this->R = R;
    this->S = S;
    this->initialScaleFactor = initScaleFact;
    this->angle = 0.f;
    this->currScale = initScaleFact;

}

glm::mat4 item::getModelMat()
{
    return T*R*S;
}

glm::mat4 item::getModelMatInv()
{
    return glm::transpose(glm::inverse(getModelMat()));
}

void item::changeScaleMat()
{
    this->currScale += 0.01;
    this->S = glm::scale(glm::mat4(1.f), glm::vec3(this->currScale,this->currScale,this->currScale));
}

bool item::compareIdx(int i, int j)
{
    if(i == gridLocI && j == gridLocJ)
        return true;
    return false;
}
glm::vec3 item::getColor()
{
    return color;
}
void item::incrementAngle()
{
    this->angle += 0.5;
    this->R = glm::rotate(glm::mat4(1.f), glm::radians(this->angle), glm::vec3(0, 1, 0));
}

bool item::shouldDeleted()
{
    return  this->currScale / initialScaleFactor > 1.5;
}