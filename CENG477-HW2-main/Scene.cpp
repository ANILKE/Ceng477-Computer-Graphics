#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <cmath>

#include "Scene.h"
#include "Camera.h"
#include "Color.h"
#include "Mesh.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Triangle.h"
#include "Vec3.h"
#include "tinyxml2.h"
#include "Helpers.h"

using namespace tinyxml2;
using namespace std;

/*
	Transformations, clipping, culling, rasterization are done here.
	You may define helper functions.
*/

static bool is_visible(double den, double num, double& t_E, double& t_L) {
    if (den > 0) {
        double t = num / den;
        if (t > t_L){
            return false;
        }
        else if (t > t_E){
            t_E = t;
        }
    }
    else if (den < 0) {
        double t = num / den;
        if (t < t_E){
            return false;
        }
        else if (t < t_L){
            t_L = t;
        }

    }
    else if (num > 0) {
        return false;
    }
    return true;
}
bool Liang_Barsky(Vec4& p_0 , Vec4& p_1, Color& c_1, Color& c_2){
    bool isVisible = false;

    double t_E = 0, t_L = 1;
    double dx = p_1.x - p_0.x, dy = p_1.y - p_0.y, dz = p_1.z - p_0.z;
    double difR = c_2.r - c_1.r;
    double difG = c_2.g - c_1.g;
    double difB = c_2.b - c_1.b;
    double x_min = -1, y_min = -1, z_min = -1;

    double x_max = 1, y_max = 1, z_max = 1;
    if (is_visible(dx, x_min - p_0.x, t_E, t_L)){
        if( is_visible(-dx, p_0.x - x_max, t_E, t_L)) {
            if (is_visible(dy, y_min - p_0.y, t_E, t_L)) {
                if (is_visible(-dy, p_0.y - y_max, t_E, t_L)) {
                    if (is_visible(dz, z_min - p_0.z, t_E, t_L)) {
                        if (is_visible(-dz, p_0.z - z_max, t_E, t_L)) {
                            isVisible = true;
                            if (t_L < 1) {
                                p_1.x = p_0.x + (dx * t_L);
                                p_1.y = p_0.y + (dy * t_L);
                                p_1.z = p_0.z + (dz * t_L);
                                c_2.r = c_1.r + (difR * t_L);
                                c_2.g = c_1.g + (difG * t_L);
                                c_2.b = c_1.b + (difB * t_L);
                            }
                            if (t_E > 0) {
                                p_0.x = p_0.x + (dx * t_E);
                                p_0.y = p_0.y + (dy * t_E);
                                p_0.z = p_0.z + (dz * t_E);
                                c_1.r = c_1.r + (difR * t_E);
                                c_1.g = c_1.g + (difG * t_E);
                                c_1.b = c_1.b + (difB * t_E);
                            }
                        }
                    }
                }
            }
        }
    }
    return isVisible;
}

Matrix4 computeModelingMatrix(const Scene * scene,Mesh * mesh)
{
	Matrix4 matrix = getIdentityMatrix();
	int transformSize = mesh -> numberOfTransformations;
	for(int i = 0;i < transformSize;i++)
	{
		int transformIdx = mesh->transformationIds[i] - 1;
		char transformType = mesh->transformationTypes[i];

		if(transformType == 'r')
		{
			// WARNING: Rotation Matrix Not Implemented!
			Matrix4 rotationMatrix = scene->rotations[transformIdx]->getRotationMatrix();
			matrix = multiplyMatrixWithMatrix(rotationMatrix,matrix);
		}
		else if(transformType == 't')
		{
			Matrix4 translationMatrix = scene->translations[transformIdx] ->getTranslationMatrix();
			matrix = multiplyMatrixWithMatrix(translationMatrix,matrix);
		}
		else
		{
			Matrix4 scalingMatrix = scene->scalings[transformIdx] ->getScaleMatrix();
			matrix = multiplyMatrixWithMatrix(scalingMatrix,matrix);
		}
	}
	return matrix;
}

Matrix4 cameraTransMat(Camera * camera)
{
	double rot[4][4] = {
							{camera->u.x,camera->u.y,camera->u.z,0},
							{camera->v.x,camera->v.y,camera->v.z,0},
							{camera->w.x,camera->w.y,camera->w.z,0},
							{          0,          0,          0,1}
					   };

	double trans[4][4] = {
							{1,0,0,-camera->pos.x},
							{0,1,0,-camera->pos.y},
							{0,0,1,-camera->pos.z},
							{0,0,0,             1}
						 };

	Matrix4 rotationMat(rot);
	Matrix4 translationMat(trans);
	return multiplyMatrixWithMatrix(rotationMat,translationMat);
}
Matrix4 orthographicMat(Camera * camera)
{
	double orth[4][4] = {
							{2.0 / (camera->right - camera->left), 0,   0, - ((camera->right + camera->left) / (camera->right - camera->left))},
							{0, 2.0 / (camera->top - camera->bottom),   0, - ((camera->top + camera->bottom) / (camera->top - camera->bottom))},
							{0,0, - (2.0 / (camera ->far - camera->near)), -((camera->far + camera->near) / (camera->far - camera->near))},
							{0,0,0,1}
						};

	return Matrix4(orth);
}
Matrix4 perspectiveMat(Camera * camera)
{
	// Sikeyim bunlar ne .
	double pers[4][4] = {
							{(2 * camera->near) / (camera->right - camera->left),0,((camera->right + camera->left) / (camera->right - camera->left)),0},
							{0,(2 * camera->near) / (camera->top - camera->bottom),((camera->top + camera->bottom) / (camera->top - camera->bottom)),0},
							{0,0,-((camera->far + camera->near) / (camera->far - camera->near)),-((2*camera->far*camera->near) / (camera->far - camera->near))},
							{0,0,-1,0}
						};
	return Matrix4(pers);
}
Matrix4 calculateViewportMatrix(Camera * camera)
{
	double tmp[4][4] = {
							{double(camera->horRes) / 2, 0,0, double(camera->horRes - 1) / 2},
							{0 , double(camera->verRes) / 2,0,double(camera->verRes - 1) / 2},
							{0,0,0.5,0.5},
                            {0,0,0,1.0}
					   };
	return Matrix4(tmp);
}
Vec4 perspectiveDivide(Vec4 * vec)
{
	vec->x /=vec->t;
	vec->y /=vec->t;
	vec->z /=vec->t;
	vec->t /=vec->t;
	return *vec;
}
Matrix4 viewingMatrixCalculator(Camera * camera)
{

	switch (camera->projectionType)
	{
	case 0:
		return multiplyMatrixWithMatrix(orthographicMat(camera),cameraTransMat(camera));
		break;
	case 1:
		return multiplyMatrixWithMatrix(perspectiveMat(camera),cameraTransMat(camera));
		break;
	default:
		break;
	}
}
Vec3 calculateCenter(Vec4& vec1, Vec4& vec2, Vec4& vec3)
{
	double xLoc = (vec1.x + vec2.x + vec3.x) / 3;
    double yLoc = (vec1.y + vec2.y + vec3.y) / 3;
    double zLoc = (vec1.z + vec2.z + vec3.z) / 3;

	return Vec3(xLoc,yLoc,zLoc,-1);

}
Vec3 calculateNormal(Vec4& vec1, Vec4& vec2, Vec4& vec3)
{
	Vec3 BtoA = Vec3(vec2.x - vec1.x,vec2.y - vec1.y,vec2.z - vec1.z,-1);
	Vec3 CtoA = Vec3(vec3.x - vec1.x,vec3.y - vec1.y,vec3.z - vec1.z,-1);
	Vec3 norm = crossProductVec3(CtoA,BtoA);
    norm = normalizeVec3(norm);
	return norm;
}
double lineEquetionForEdge(double x, double y, double x_0, double x_1, double y_0, double y_1){
    double result;
    result = x*(y_0-y_1) + y*(x_1-x_0) + x_0*y_1 - y_0*x_1;
    return result;

}

void triangleRasterizer(vector<vector<Color>>& image,Vec3 vector1, Vec3 vector2, Vec3 vector3, Color& c_1, Color& c_2, Color& c_3,int limit_x,int limit_y){
    int x_0,x_1,x_2,y_0,y_1,y_2,x_min = 100000,x_max,y_min = 100000,y_max;
    x_0 = vector1.x;
    x_1 = vector2.x;
    x_2 = vector3.x;
    y_0 = vector1.y;
    y_1 = vector2.y;
    y_2 = vector3.y;
    x_min = min(min(x_0, x_1), x_2) >= 0 ? min(min(x_0, x_1), x_2) : 0;
    y_min = min(min(y_0, y_1), y_2) >= 0 ? min(min(y_0, y_1), y_2) : 0;
    x_max = max(max(x_0, x_1), x_2) < 0 ? 0 : max(max(x_0, x_1), x_2) ;
    y_max = max(max(y_0, y_1), y_2) < 0 ? 0 : max(max(y_0, y_1), y_2) ;

    if(x_min>limit_x){
        x_min = limit_x;
    }
    if(y_min>limit_y){
        y_min = limit_y;
    }
    if(x_max>limit_x){
        x_max = limit_x;
    }
    if(y_max>limit_y){
        y_max = limit_y;
    }

    double alpha_div = lineEquetionForEdge(x_0,y_0,x_1,x_2,y_1,y_2);
    double beta_div = lineEquetionForEdge(x_1,y_1,x_2,x_0,y_2,y_0);
    double gamma_div = lineEquetionForEdge(x_2,y_2,x_0,x_1,y_0,y_1);
    for(int y=y_min; y<=y_max; y++){
        for(int x=x_min; x<=x_max; x++){
            double alpha = lineEquetionForEdge(x,y,x_1,x_2,y_1,y_2) /alpha_div;
            double beta = lineEquetionForEdge(x,y,x_2,x_0,y_2,y_0) /beta_div;
            double gamma = lineEquetionForEdge(x,y,x_0,x_1,y_0,y_1) /gamma_div;
            if(alpha>=0 && beta>=0 && gamma>=0){
                Color col;
                col.r = int(((alpha*c_1.r)+(beta*c_2.r)+ (gamma*c_3.r))+0.5);
                col.g = int(((alpha*c_1.g)+(beta*c_2.g)+ (gamma*c_3.g))+0.5);
                col.b = int(((alpha*c_1.b)+(beta*c_2.b)+ (gamma*c_3.b))+0.5);
                image[x][y] = col;
            }
        }
        }

}

void lineRasterizer(vector<vector<Color>>& image,Vec3 vector1, Vec3 vector2, Color& c_1, Color& c_2)
{
    double slope_x ,slope_y, d ;
    Color dc;
    slope_x = vector2.x-vector1.x;
    slope_y = vector2.y-vector1.y;
    bool flag_for_minus = false;
    if (abs(slope_y) <= abs(slope_x)) {
        if(vector2.x< vector1.x){
            lineRasterizer(image,vector2,vector1,c_2,c_1);
            return;
        }
        if(vector2.y < vector1.y){
            flag_for_minus = true;
        }
        int y = vector1.y;
        if(flag_for_minus){
            d = (vector1.y - vector2.y) + (-1*0.5 * (vector2.x - vector1.x));
        }
        else {
            d = (vector1.y - vector2.y) + (0.5 * (vector2.x - vector1.x));
        }
        Color c = c_1;
        Color diff_c, div_c;
        diff_c.r = c_2.r-c_1.r;
        diff_c.b = c_2.b-c_1.b;
        diff_c.g = c_2.g-c_1.g;
        div_c.r = diff_c.r /(vector2.x - vector1.x);
        div_c.b = diff_c.b /(vector2.x - vector1.x);
        div_c.g = diff_c.g / (vector2.x - vector1.x);
        dc = div_c;
        for (int x = vector1.x; x <= vector2.x && y>=0; x++) {
            Color temp_c;
            temp_c.r = int(c.r+0.5);
            temp_c.g = int(c.g+0.5);
            temp_c.b = int(c.b+0.5);
            image[x][y] = temp_c;
            if(flag_for_minus){
                if (d *-1 < 0) { // choose NE
                    y -= 1;
                    d += (vector1.y - vector2.y) + -1*(vector2.x - vector1.x);
                }
                else // choose E
                    d += (vector1.y - vector2.y);
                c.r += dc.r;
                c.g += dc.g;
                c.b += dc.b;
            }
            else{
                if (d  < 0) { // choose NE
                    y += 1;
                    d += (vector1.y - vector2.y) + (vector2.x - vector1.x);
                }
                else // choose E
                    d += (vector1.y - vector2.y);
                c.r += dc.r;
                c.g += dc.g;
                c.b += dc.b;
            }

        }
    }
    else {
        if(vector2.y< vector1.y){
            lineRasterizer(image,vector2,vector1,c_2,c_1);
            return;
        }
        if(vector2.x< vector1.x){
            flag_for_minus = true;
        }
        int x = vector1.x;
        if(flag_for_minus){
            d = (vector2.x - vector1.x) + (-1*0.5 * (vector1.y - vector2.y));
        }
        else{
            d = (vector2.x - vector1.x) + (0.5 * (vector1.y - vector2.y));
        }
        Color c = c_1;

        Color diff_c, div_c;
        diff_c.r = c_2.r - c_1.r;
        diff_c.b = c_2.b - c_1.b;
        diff_c.g = c_2.g - c_1.g;
        div_c.r = diff_c.r / (vector2.y - vector1.y);
        div_c.b = diff_c.b / (vector2.y - vector1.y);
        div_c.g = diff_c.g / (vector2.y - vector1.y);
        dc = div_c;
        for (int y = vector1.y; y <= vector2.y && x>=0; y++) {
            Color temp_c;
            temp_c.r = int(c.r + 0.5);
            temp_c.g = int(c.g + 0.5);
            temp_c.b = int(c.b + 0.5);
            image[x][y] = temp_c;
            if(flag_for_minus){
                if (d*-1 > 0) {
                    x -= 1;
                    d += (vector2.x - vector1.x) + -1*(vector1.y - vector2.y);
                } else
                    d += (vector2.x - vector1.x);
            }
            else{
                if (d > 0) {
                    x += 1;
                    d += (vector2.x - vector1.x) + (vector1.y - vector2.y);
                } else
                    d += (vector2.x - vector1.x);
            }

            c.r += dc.r;
            c.g += dc.g;
            c.b += dc.b;
        }
    }
}

void Scene::forwardRenderingPipeline(Camera *camera)
{
	initializeImage(camera);
	// TODO: Implement this function.

	//Projection type handle.
	Matrix4 viewingMatrix = viewingMatrixCalculator(camera);

	//Viewport Matrix
	Matrix4 viewPortMatrix = calculateViewportMatrix(camera);

	// Modeling Transformation phase.



    for(Mesh * currMesh : this->meshes)
	{
		Matrix4 modelingMatrix = computeModelingMatrix(this,currMesh);

        vector<pair<Vec4,Color>> transformedPoints;
		for(Triangle triangle : currMesh->triangles)
		{
			// Get points indexes.
			int idx1 = triangle.getFirstVertexId() - 1;
			int idx2 = triangle.getSecondVertexId() - 1;
			int idx3 = triangle.getThirdVertexId() - 1;

			// Get points.
			Vec3 * point1 = this->vertices[idx1];
			Vec3 * point2 = this->vertices[idx2];
			Vec3 * point3 = this->vertices[idx3];

            // Get points color.
             const Color * c_1 = this->colorsOfVertices[point1->colorId-1];
             const Color * c_2 = this->colorsOfVertices[point2->colorId-1];
             const Color * c_3 = this->colorsOfVertices[point3->colorId-1];

			//Create Homogeneous Coordinates.
			Vec4 vec1(*point1);
			Vec4 vec2(*point2);
			Vec4 vec3(*point3);

			//Modeling Transforms.
			Vec4 modelingVec1 = multiplyMatrixWithVec4(modelingMatrix,vec1);
			Vec4 modelingVec2 = multiplyMatrixWithVec4(modelingMatrix,vec2);
			Vec4 modelingVec3 = multiplyMatrixWithVec4(modelingMatrix,vec3);



			// Camera Transforms.
			Vec4 cameraModelingVec1 = multiplyMatrixWithVec4(viewingMatrix,modelingVec1);
			Vec4 cameraModelingVec2 = multiplyMatrixWithVec4(viewingMatrix,modelingVec2);
			Vec4 cameraModelingVec3 = multiplyMatrixWithVec4(viewingMatrix,modelingVec3);


			//Perspective Divide.
			cameraModelingVec1 = perspectiveDivide(&cameraModelingVec1);
		    cameraModelingVec2 = perspectiveDivide(&cameraModelingVec2);
			cameraModelingVec3 = perspectiveDivide(&cameraModelingVec3);


			if(this->cullingEnabled == true)
			{
				// WARNING CAN BE WRONG, CHECK THE PRODUCT ORDERING.
				Vec3 centerOfMass = calculateCenter(cameraModelingVec1,cameraModelingVec2,cameraModelingVec3);
				Vec3 normal = calculateNormal(cameraModelingVec1,cameraModelingVec2,cameraModelingVec3);
				double facingDirection = dotProductVec3(centerOfMass,normal);

				if(facingDirection < 0)
				{
					// Front Facing triangle.
					transformedPoints.push_back(make_pair(cameraModelingVec1,*c_1));
					transformedPoints.push_back(make_pair(cameraModelingVec2,*c_2));
					transformedPoints.push_back(make_pair(cameraModelingVec3,*c_3));
				}
			}
			else
			{
                transformedPoints.push_back(make_pair(cameraModelingVec1,*c_1));
                transformedPoints.push_back(make_pair(cameraModelingVec2,*c_2));
                transformedPoints.push_back(make_pair(cameraModelingVec3,*c_3));
			}

        }
        //Liang_Barsky Clipping.
        vector<pair<pair<Vec4,bool>,Color>> clipped_lines;
        if(currMesh->type == 0 && transformedPoints.size()){
            for (int i = 0; i <= transformedPoints.size() - 3; i += 3) {
                Vec4 ver_12 = transformedPoints[i].first;
                Vec4 ver_13 = transformedPoints[i].first;
                Vec4 ver_21 = transformedPoints[i+1].first;
                Vec4 ver_31 = transformedPoints[i+2].first;
                Vec4 ver_23 = transformedPoints[i+1].first;
                Vec4 ver_32 = transformedPoints[i+2].first;

                Color c_12 = transformedPoints[i].second;
                Color c_13 = transformedPoints[i].second;
                Color c_21 = transformedPoints[i+1].second;
                Color c_23 = transformedPoints[i+1].second;
                Color c_31 = transformedPoints[i+2].second;
                Color c_32 = transformedPoints[i+2].second;

                bool visiblity1 = Liang_Barsky(ver_12,ver_21, c_12, c_21);
                bool visiblity2 = Liang_Barsky(ver_23,ver_32, c_23, c_32);
                bool visiblity3 = Liang_Barsky(ver_31,ver_13, c_31, c_13);


                clipped_lines.push_back(make_pair(make_pair(ver_12,visiblity1),c_12));
                clipped_lines.push_back(make_pair(make_pair(ver_13,visiblity3),c_13));
                clipped_lines.push_back(make_pair(make_pair(ver_21,visiblity1),c_21));
                clipped_lines.push_back(make_pair(make_pair(ver_23,visiblity2),c_23));
                clipped_lines.push_back(make_pair(make_pair(ver_31,visiblity3),c_31));
                clipped_lines.push_back(make_pair(make_pair(ver_32,visiblity2),c_32));


            }
        }


        if(currMesh->type == 0 && clipped_lines.size())
        {
            // Wireframe Mode.
            for(int i = 0;i <= clipped_lines.size() - 6; i+=6)
            {
                pair<pair<Vec4,bool>,Color> pair12 = clipped_lines[i];
                pair<pair<Vec4,bool>,Color> pair13 = clipped_lines[i + 1];
                pair<pair<Vec4,bool>,Color> pair21 = clipped_lines[i + 2];
                pair<pair<Vec4,bool>,Color> pair23 = clipped_lines[i + 3];
                pair<pair<Vec4,bool>,Color> pair31 = clipped_lines[i + 4];
                pair<pair<Vec4,bool>,Color> pair32 = clipped_lines[i + 5];
                if(pair12.first.second)
                {
                    Vec4 p1 = pair12.first.first;
                    Vec4 p2 = pair21.first.first;

                    Color c_1 = pair12.second;
                    Color c_2 = pair21.second;

                    Vec4 vec1 = multiplyMatrixWithVec4(viewPortMatrix,p1);
                    Vec4 vec2 = multiplyMatrixWithVec4(viewPortMatrix,p2);

                    Vec3 arg1(vec1.x,vec1.y,vec1.z,vec1.colorId);
                    Vec3 arg2(vec2.x,vec2.y,vec2.z,vec2.colorId);

                    lineRasterizer(this->image,arg1,arg2,c_1, c_2);

                }
                if(pair23.first.second)
                {
                    Vec4 p2 = pair23.first.first;
                    Vec4 p3 = pair32.first.first;

                    Color c_2 = pair23.second;
                    Color c_3 = pair32.second;

                    Vec4 vec2 = multiplyMatrixWithVec4(viewPortMatrix,p2);
                    Vec4 vec3 = multiplyMatrixWithVec4(viewPortMatrix,p3);

                    Vec3 arg1(vec2.x,vec2.y,vec2.z,vec2.colorId);
                    Vec3 arg2(vec3.x,vec3.y,vec3.z,vec3.colorId);

                    lineRasterizer(this->image,arg1,arg2,c_2,c_3);

                }
                if(pair31.first.second)
                {
                    Vec4 p3 = pair31.first.first;
                    Vec4 p1 = pair13.first.first;

                    Color c_3 = pair31.second;
                    Color c_1 = pair13.second;

                    Vec4 vec3 = multiplyMatrixWithVec4(viewPortMatrix,p3);
                    Vec4 vec1 = multiplyMatrixWithVec4(viewPortMatrix,p1);

                    Vec3 arg1(vec3.x,vec3.y,vec3.z,vec3.colorId);
                    Vec3 arg2(vec1.x,vec1.y,vec1.z,vec1.colorId);

                    lineRasterizer(this->image,arg1,arg2, c_3, c_1);
                }
                if(pair13.first.second)
                {
                    Vec4 p1 = pair13.first.first;
                    Vec4 p2 = pair31.first.first;

                    Color c_1 = pair13.second;
                    Color c_2 = pair31.second;

                    Vec4 vec1 = multiplyMatrixWithVec4(viewPortMatrix,p1);
                    Vec4 vec2 = multiplyMatrixWithVec4(viewPortMatrix,p2);

                    Vec3 arg1(vec1.x,vec1.y,vec1.z,vec1.colorId);
                    Vec3 arg2(vec2.x,vec2.y,vec2.z,vec2.colorId);

                    lineRasterizer(this->image,arg1,arg2,c_1, c_2);

                }
                if(pair21.first.second)
                {
                    Vec4 p2 = pair21.first.first;
                    Vec4 p3 = pair12.first.first;

                    Color c_2 = pair21.second;
                    Color c_3 = pair12.second;

                    Vec4 vec2 = multiplyMatrixWithVec4(viewPortMatrix,p2);
                    Vec4 vec3 = multiplyMatrixWithVec4(viewPortMatrix,p3);

                    Vec3 arg1(vec2.x,vec2.y,vec2.z,vec2.colorId);
                    Vec3 arg2(vec3.x,vec3.y,vec3.z,vec3.colorId);

                    lineRasterizer(this->image,arg1,arg2,c_2,c_3);

                }
                if(pair32.first.second)
                {
                    Vec4 p3 = pair32.first.first;
                    Vec4 p1 = pair23.first.first;

                    Color c_3 = pair32.second;
                    Color c_1 = pair23.second;

                    Vec4 vec3 = multiplyMatrixWithVec4(viewPortMatrix,p3);
                    Vec4 vec1 = multiplyMatrixWithVec4(viewPortMatrix,p1);

                    Vec3 arg1(vec3.x,vec3.y,vec3.z,vec3.colorId);
                    Vec3 arg2(vec1.x,vec1.y,vec1.z,vec1.colorId);

                    lineRasterizer(this->image,arg1,arg2, c_3, c_1);
                }
            }

        }
        else
        {
            // Solid Mode.
            for(int i = 0; i<= transformedPoints.size()-3;i+=3){
                Vec4 vec1 = multiplyMatrixWithVec4(viewPortMatrix,transformedPoints[i].first);
                Vec4 vec2 = multiplyMatrixWithVec4(viewPortMatrix,transformedPoints[i+1].first);
                Vec4 vec3 = multiplyMatrixWithVec4(viewPortMatrix,transformedPoints[i+2].first);
                Color c_1 = transformedPoints[i].second;
                Color c_2 = transformedPoints[i+1].second;
                Color c_3 = transformedPoints[i+2].second;
                Vec3 arg1(vec1.x,vec1.y,vec1.z,vec1.colorId);
                Vec3 arg2(vec2.x,vec2.y,vec2.z,vec2.colorId);
                Vec3 arg3(vec3.x,vec3.y,vec3.z,vec3.colorId);
                triangleRasterizer(this->image,arg1,arg2,arg3,c_1,c_2,c_3,camera->horRes-1,camera->verRes-1);
            }


        }

	}




	writeImageToPPMFile(camera);

}

/*
	Parses XML file
*/
Scene::Scene(const char *xmlPath)
{
	const char *str;
	XMLDocument xmlDoc;
	XMLElement *pElement;

	xmlDoc.LoadFile(xmlPath);

	XMLNode *pRoot = xmlDoc.FirstChild();

	// read background color
	pElement = pRoot->FirstChildElement("BackgroundColor");
	str = pElement->GetText();
	sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

	// read culling
	pElement = pRoot->FirstChildElement("Culling");
	if (pElement != NULL) {
		str = pElement->GetText();

		if (strcmp(str, "enabled") == 0) {
			cullingEnabled = true;
		}
		else {
			cullingEnabled = false;
		}
	}

	// read cameras
	pElement = pRoot->FirstChildElement("Cameras");
	XMLElement *pCamera = pElement->FirstChildElement("Camera");
	XMLElement *camElement;
	while (pCamera != NULL)
	{
		Camera *cam = new Camera();

		pCamera->QueryIntAttribute("id", &cam->cameraId);

		// read projection type
		str = pCamera->Attribute("type");

		if (strcmp(str, "orthographic") == 0) {
			cam->projectionType = 0;
		}
		else {
			cam->projectionType = 1;
		}

		camElement = pCamera->FirstChildElement("Position");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf", &cam->pos.x, &cam->pos.y, &cam->pos.z);

		camElement = pCamera->FirstChildElement("Gaze");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf", &cam->gaze.x, &cam->gaze.y, &cam->gaze.z);

		camElement = pCamera->FirstChildElement("Up");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf", &cam->v.x, &cam->v.y, &cam->v.z);

		cam->gaze = normalizeVec3(cam->gaze);
		cam->u = crossProductVec3(cam->gaze, cam->v);
		cam->u = normalizeVec3(cam->u);

		cam->w = inverseVec3(cam->gaze);
		cam->v = crossProductVec3(cam->u, cam->gaze);
		cam->v = normalizeVec3(cam->v);

		camElement = pCamera->FirstChildElement("ImagePlane");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
			   &cam->left, &cam->right, &cam->bottom, &cam->top,
			   &cam->near, &cam->far, &cam->horRes, &cam->verRes);

		camElement = pCamera->FirstChildElement("OutputName");
		str = camElement->GetText();
		cam->outputFileName = string(str);

		cameras.push_back(cam);

		pCamera = pCamera->NextSiblingElement("Camera");
	}

	// read vertices
	pElement = pRoot->FirstChildElement("Vertices");
	XMLElement *pVertex = pElement->FirstChildElement("Vertex");
	int vertexId = 1;

	while (pVertex != NULL)
	{
		Vec3 *vertex = new Vec3();
		Color *color = new Color();

		vertex->colorId = vertexId;

		str = pVertex->Attribute("position");
		sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

		str = pVertex->Attribute("color");
		sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

		vertices.push_back(vertex);
		colorsOfVertices.push_back(color);

		pVertex = pVertex->NextSiblingElement("Vertex");

		vertexId++;
	}

	// read translations
	pElement = pRoot->FirstChildElement("Translations");
	XMLElement *pTranslation = pElement->FirstChildElement("Translation");
	while (pTranslation != NULL)
	{
		Translation *translation = new Translation();

		pTranslation->QueryIntAttribute("id", &translation->translationId);

		str = pTranslation->Attribute("value");
		sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

		translations.push_back(translation);

		pTranslation = pTranslation->NextSiblingElement("Translation");
	}

	// read scalings
	pElement = pRoot->FirstChildElement("Scalings");
	XMLElement *pScaling = pElement->FirstChildElement("Scaling");
	while (pScaling != NULL)
	{
		Scaling *scaling = new Scaling();

		pScaling->QueryIntAttribute("id", &scaling->scalingId);
		str = pScaling->Attribute("value");
		sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

		scalings.push_back(scaling);

		pScaling = pScaling->NextSiblingElement("Scaling");
	}

	// read rotations
	pElement = pRoot->FirstChildElement("Rotations");
	XMLElement *pRotation = pElement->FirstChildElement("Rotation");
	while (pRotation != NULL)
	{
		Rotation *rotation = new Rotation();

		pRotation->QueryIntAttribute("id", &rotation->rotationId);
		str = pRotation->Attribute("value");
		sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

		rotations.push_back(rotation);

		pRotation = pRotation->NextSiblingElement("Rotation");
	}

	// read meshes
	pElement = pRoot->FirstChildElement("Meshes");

	XMLElement *pMesh = pElement->FirstChildElement("Mesh");
	XMLElement *meshElement;
	while (pMesh != NULL)
	{
		Mesh *mesh = new Mesh();

		pMesh->QueryIntAttribute("id", &mesh->meshId);

		// read projection type
		str = pMesh->Attribute("type");

		if (strcmp(str, "wireframe") == 0) {
			mesh->type = 0;
		}
		else {
			mesh->type = 1;
		}

		// read mesh transformations
		XMLElement *pTransformations = pMesh->FirstChildElement("Transformations");
		XMLElement *pTransformation = pTransformations->FirstChildElement("Transformation");

		while (pTransformation != NULL)
		{
			char transformationType;
			int transformationId;

			str = pTransformation->GetText();
			sscanf(str, "%c %d", &transformationType, &transformationId);

			mesh->transformationTypes.push_back(transformationType);
			mesh->transformationIds.push_back(transformationId);

			pTransformation = pTransformation->NextSiblingElement("Transformation");
		}

		mesh->numberOfTransformations = mesh->transformationIds.size();

		// read mesh faces
		char *row;
		char *clone_str;
		int v1, v2, v3;
		XMLElement *pFaces = pMesh->FirstChildElement("Faces");
        str = pFaces->GetText();
		clone_str = strdup(str);

		row = strtok(clone_str, "\n");
		while (row != NULL)
		{
			int result = sscanf(row, "%d %d %d", &v1, &v2, &v3);

			if (result != EOF) {
				mesh->triangles.push_back(Triangle(v1, v2, v3));
			}
			row = strtok(NULL, "\n");
		}
		mesh->numberOfTriangles = mesh->triangles.size();
		meshes.push_back(mesh);

		pMesh = pMesh->NextSiblingElement("Mesh");
	}
}

/*
	Initializes image with background color
*/
void Scene::initializeImage(Camera *camera)
{
	if (this->image.empty())
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			vector<Color> rowOfColors;

			for (int j = 0; j < camera->verRes; j++)
			{
				rowOfColors.push_back(this->backgroundColor);
			}

			this->image.push_back(rowOfColors);
		}
	}
	else
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			for (int j = 0; j < camera->verRes; j++)
			{
				this->image[i][j].r = this->backgroundColor.r;
				this->image[i][j].g = this->backgroundColor.g;
				this->image[i][j].b = this->backgroundColor.b;
			}
		}
	}
}

/*
	If given value is less than 0, converts value to 0.
	If given value is more than 255, converts value to 255.
	Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value)
{
	if (value >= 255.0)
		return 255;
	if (value <= 0.0)
		return 0;
	return (int)(value);
}

/*
	Writes contents of image (Color**) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera *camera)
{
	ofstream fout;

	fout.open(camera->outputFileName.c_str());

	fout << "P3" << endl;
	fout << "# " << camera->outputFileName << endl;
	fout << camera->horRes << " " << camera->verRes << endl;
	fout << "255" << endl;

	for (int j = camera->verRes - 1; j >= 0; j--)
	{
		for (int i = 0; i < camera->horRes; i++)
		{

			fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].g) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
		}
		fout << endl;
	}
	fout.close();
}

/*
	Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
	os_type == 1 		-> Ubuntu
	os_type == 2 		-> Windows
	os_type == other	-> No conversion
*/
void Scene::convertPPMToPNG(string ppmFileName, int osType)
{
	string command;

	// call command on Ubuntu
	if (osType == 1)
	{
		command = "convert " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// call command on Windows
	else if (osType == 2)
	{
		command = "magick convert " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// default action - don't do conversion
	else
	{
	}
}
