#ifndef OBJECT_H
#define OBJECT_H
#include <cstdio>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <string>
#include <map>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <GL/glew.h>   // The GL Header File
#include <GL/gl.h>   // The GL Header File
#include <GLFW/glfw3.h> // The GLFW header
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <ft2build.h>


class item
{
    
    private:

        int gridLocI;
        int gridLocJ;
        glm::vec3 color;
        glm::mat4 T;
        glm::mat4 R;
        glm::mat4 S;
        float angle;
        float initialScaleFactor;
        float currScale;
    


    public:
        item(int i, int j, glm::vec3 color,glm::mat4 T,glm::mat4 R,glm::mat4 S,float initialScaleFactor);
        glm::mat4 getModelMat();
        void draw();
        glm::mat4 getModelMatInv();
        void changeScaleMat();
        bool compareIdx(int i, int j);
        glm::vec3 getColor();
        void incrementAngle();
        bool shouldDeleted();
        

};



#endif