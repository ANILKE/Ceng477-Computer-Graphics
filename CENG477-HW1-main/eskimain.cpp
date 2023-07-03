
#include <iostream>
#include "parser.h"
#include "ppm.h"
#include <math.h>

parser::Scene scene;
float left, right, bottom, top,height_su,width_sv;;
parser::Vec3f m, q,e, u,v,w;
parser::Vec3i backgroundColor;
int height, width,num_of_sphere;
typedef unsigned char RGB[3];

typedef struct
{
    parser::Vec3f o,d;
} ray;


parser::Vec3f multS(parser::Vec3f a,float s)
{
    parser::Vec3f result;
    result.x = a.x*s;
    result.y = a.y*s;
    result.z = a.z*s;
    return result;
}

parser::Vec3f add(parser::Vec3f a, parser::Vec3f b)
{
    parser::Vec3f result;
    result.x = a.x+b.x;
    result.y = a.y+b.y;
    result.z = a.z+b.z;
    return result;
}


float dot(parser::Vec3f a,parser::Vec3f b)
{
    return a.x*b.x+a.y*b.y+a.z*b.z;
}

parser::Vec3f normalize(parser::Vec3f a)
{
    return multS(a,1.0/sqrt(dot(a,a)));
}


parser::Vec3f setCameraParams(parser::Vec3f vec){
    parser::Vec3f res;
    res.x = vec.x;
    res.y = vec.y;
    res.z = vec.z;

    return res;
}

parser::Vec3f setPosition(const parser::Vec3f & pos){
    parser::Vec3f res;
    res.x = pos.x;
    res.y = pos.y;
    res.z = pos.z;
    return res;
}
parser::Vec3f crossProduct(parser::Vec3f first,parser::Vec3f second){
    parser::Vec3f res;
    res.x = (first.y*second.z-first.z*second.y);
    res.y = (first.z*second.x-first.x*second.z);
    res.z = first.x*second.y-first.y*second.x;
    return res;
}

ray generateRay(int i, int j)
{
    ray result;
    float su,sv;
    parser::Vec3f s;
    su=(j+0.5)*height_su;
    sv=(i+0.5)*width_sv;
    s = add(q,add(multS(u,su),multS(v,-sv)));
    result.o = e;
    result.d = normalize(add(s,multS(e,-1)));
    return result;
}
parser::Vec3f multVects(parser::Vec3f a,parser::Vec3f b){
    parser::Vec3f result;
    result.x=a.x*b.x;
    result.y=a.y*b.y;
    result.z=a.z*b.z;
    return result;
}
float getDeterminant(const parser::Vec3f &first_column, const parser::Vec3f &second_column, const parser::Vec3f &third_column)
{
    float result;
    result=first_column.x * (second_column.y*third_column.z - third_column.y*second_column.z)
           + first_column.y * (third_column.x*second_column.z - second_column.x*third_column.z)
           + first_column.z * (second_column.x*third_column.y - second_column.y*third_column.x);
    return result;
}
float intersectTriangle(ray r,int index){
    parser::Vec3f point_A,point_B,point_C,A_to_B,A_to_C,A_to_O;
    int point_A_id=scene.triangles[index].indices.v0_id-1;
    int point_B_id=scene.triangles[index].indices.v1_id-1;
    int point_C_id=scene.triangles[index].indices.v2_id-1;
    float beta,gamma,t,A_determinant;
    point_A=scene.vertex_data[point_A_id];
    point_B=scene.vertex_data[point_B_id];
    point_C=scene.vertex_data[point_C_id];
    A_to_B= add(point_A, multS(point_B,-1));
    A_to_C= add(point_A, multS(point_C,-1));
    A_to_O= add(point_A, multS(r.o,-1));
    A_determinant= getDeterminant(A_to_B,A_to_C,r.d);
    if(A_determinant==0){
        return -1;
    }
    t= getDeterminant(A_to_B,A_to_C,A_to_O)/A_determinant;
    if(!(t>0)){
        return -1;
    }
    beta= getDeterminant(A_to_O,A_to_C,r.d)/A_determinant;
    if(beta<0){
        return -1;
    }
    gamma= getDeterminant(A_to_B,A_to_O,r.d)/A_determinant;
    if (gamma< 0)
        return -1;
    if (gamma +  beta >1)
        return -1;
    return t;
}
std::vector<float> intersectMesh(ray r, int meshIdx){
    std::vector<float> result = {99999999.0,-1};
    float t_new;
    parser::Vec3f A_to_B,A_to_C,A_to_O;
    int meshSize = scene.meshes.size();
    float beta,gamma,t=99999999.0,A_determinant;
    int last_face_index=-1;
    int faceSize = scene.meshes[meshIdx].faces.size();
    for(int faceNumber = 0; faceNumber < faceSize;faceNumber++)
    {
        parser::Vec3f point_A = scene.vertex_data[scene.meshes[meshIdx].faces[faceNumber].v0_id - 1];
        parser::Vec3f point_B = scene.vertex_data[scene.meshes[meshIdx].faces[faceNumber].v1_id - 1];
        parser::Vec3f point_C = scene.vertex_data[scene.meshes[meshIdx].faces[faceNumber].v2_id - 1];


        A_to_B= add(point_A, multS(point_B,-1));
        A_to_C= add(point_A, multS(point_C,-1));
        A_to_O= add(point_A, multS(r.o,-1));

        A_determinant= getDeterminant(A_to_B,A_to_C,r.d);
        if(A_determinant==0){
            t_new= -1;
        }
        t_new= getDeterminant(A_to_B,A_to_C,A_to_O)/A_determinant;
        if(!(t_new>0)){
            t_new= -1;
        }
        beta= getDeterminant(A_to_O,A_to_C,r.d)/A_determinant;
        if(beta<0){
            t_new= -1;
        }
        gamma= getDeterminant(A_to_B,A_to_O,r.d)/A_determinant;
        if (gamma< 0)
            t_new= -1;

        if (gamma +  beta >1)
            t_new = -1;

        if(t_new != -1 && t_new < t){
            result[0]=t_new;
            result[1]=faceNumber;

        }

    }

    if(result[0]==99999999.0)
        return {-1,-1};
    return result;
}
float intersectSphere(ray r,int index)
{
    parser::Sphere sphere=scene.spheres[index];

    std::vector<parser::Vec3f> vertex_data = scene.vertex_data;
    float A,B,C; //constants for the quadratic equation

    float delta;

    int c_ID;

    c_ID = sphere.center_vertex_id;
    parser::Vec3f center=vertex_data[c_ID-1];
    float t,t1,t2;

    C = (r.o.x-center.x)*(r.o.x-center.x)+(r.o.y-center.y)*(r.o.y-center.y)+(r.o.z-center.z)*(r.o.z-center.z)-sphere.radius*sphere.radius;

    B = 2*r.d.x*(r.o.x-center.x)+2*r.d.y*(r.o.y-center.y)+2*r.d.z*(r.o.z-center.z);

    A = r.d.x*r.d.x+r.d.y*r.d.y+r.d.z*r.d.z;

    delta = B*B-4*A*C;

    if (delta<0) return -1;
    else if (delta==0)
    {
        t = -B / (2*A);
    }
    else
    {
        delta = sqrt(delta);
        A = 2*A;
        t1 = (-B + delta) / A;
        t2 = (-B - delta) / A;

        if (t1<t2) t=t1; else t=t2;
    }

    return t;
}
float getLengthSquare(parser::Vec3f vec){
    float result=0;
    result+=pow(vec.y,2);
    result+=pow(vec.x,2);
    result+=pow(vec.z,2);
    return result;
}

parser::Vec3f computeSpecular_sphere(ray r,parser::Vec3f center,int minI,float minT)
{

    parser::Vec3f L,P,N,H;
    parser::Vec3f result = {0, 0, 0};
    P = multS(add(r.o, multS(r.d, minT)),1);
    int pointLightSize = scene.point_lights.size();
    for(int i = 0; i < pointLightSize;i++)
    {
        L = add(P, multS(scene.point_lights[i].position, -1));
        N = normalize(add(center, multS(P, -1)));
        float r_square = dot(L,L);
        L= normalize(L);
        H = normalize(add(L, normalize(add(P, multS(r.o,-1)))));
        float cos_alpha= fmax(0, dot(N, H));
        float cos_alpha_phong=pow(cos_alpha,scene.materials[scene.spheres[minI].material_id - 1].phong_exponent);
        parser::Vec3f irradiance = multS(scene.point_lights[i].intensity, 1 / r_square);
        float cos_theta= dot(L,N);
        float theta= acos(cos_theta);

        if(theta<90)
        {
            parser::Vec3f specular_const = multS(irradiance,cos_alpha_phong);
            parser::Vec3f specular = scene.materials[scene.spheres[minI].material_id - 1].specular;
            specular = multVects(specular,specular_const);
            result = add(result,specular);
        }
    }
    return result;
}

parser::Vec3f computeSpecular_triangle(ray r,int minI,float minT)
{
    parser::Vec3f L,P,N,H;
    parser::Vec3f result = {0, 0, 0};
    P = add(r.o, multS(r.d, minT));

    int pointLightSize = scene.point_lights.size();
    for(int i = 0; i < pointLightSize;i++)
    {
        int point_A_id=scene.triangles[minI].indices.v0_id-1;
        int point_B_id=scene.triangles[minI].indices.v1_id-1;
        int point_C_id=scene.triangles[minI].indices.v2_id-1;

        parser::Vec3f point_A=scene.vertex_data[point_A_id];
        parser::Vec3f point_B=scene.vertex_data[point_B_id];
        parser::Vec3f point_C=scene.vertex_data[point_C_id];
        N = normalize(crossProduct(add(point_B, multS(point_A,-1)),add(point_C, multS(point_A,-1))));
        H = normalize(add(normalize(L), multS(r.d,-1)));
        float cos_alpha= fmax(0, dot(N, H));
        float cos_alpha_phong=pow(cos_alpha,scene.materials[scene.triangles[minI].material_id - 1].phong_exponent);
        float r_square = getLengthSquare(L);
        parser::Vec3f irradiance = multS(scene.point_lights[i].intensity, 1 / r_square);
        float cos_theta= dot(normalize(L),N);
        float theta= acos(cos_theta);

        if(theta<90){
            parser::Vec3f specular_const = multS(irradiance,cos_alpha_phong);
            parser::Vec3f specular = scene.materials[scene.triangles[minI].material_id - 1].specular;
            specular = multVects(specular,specular_const);
            result = add(result,specular);
        }
    }
    return result;
}



parser::Vec3f computeSpecular_mesh(ray r,int minI,float minT,int mesh_index){

    parser::Vec3f L,P,N,result,H;
    result.x=result.y=result.z=0;
    P = add(r.o, multS(r.d, minT));

    for(int k=0;k<scene.point_lights.size();k++){
        L = add(scene.point_lights[k].position, multS(P, -1));
        int point_A_id=scene.meshes[minI].faces[mesh_index].v0_id - 1;
        int point_B_id=scene.meshes[minI].faces[mesh_index].v1_id - 1;
        int point_C_id=scene.meshes[minI].faces[mesh_index].v2_id - 1;
        parser::Vec3f point_A=scene.vertex_data[point_A_id];
        parser::Vec3f point_B=scene.vertex_data[point_B_id];
        parser::Vec3f point_C=scene.vertex_data[point_C_id];
        N= normalize(crossProduct(add(point_B, multS(point_A,-1)),add(point_C, multS(point_A,-1))));
        H= normalize(add(normalize(L), multS(r.d,-1)));
        float cos_alpha= fmax(0, dot(N, H));
        float cos_alpha_phong=pow(cos_alpha,scene.materials[scene.meshes[minI].material_id - 1].phong_exponent);
        float r_square = getLengthSquare(L);
        parser::Vec3f irradiance = multS(scene.point_lights[k].intensity, 1 / r_square);
        float cos_theta= dot(normalize(L),N);
        float theta= acos(cos_theta);

        if(theta<90){
            parser::Vec3f specular_const = multS(irradiance,cos_alpha_phong);
            parser::Vec3f specular = scene.materials[scene.meshes[minI].material_id - 1].specular;
            specular = multVects(specular,specular_const);
            result = add(result,specular);
        }

    }
    return result;
}

parser::Vec3f computeDiffuse_sphere(ray r,parser::Vec3f center,int minI,float minT){
    parser::Vec3f L,P,N,result;
    result.x=result.y=result.z=0;

    for(int k = 0;k<scene.point_lights.size();k++)
    {
        P = add(r.o, multS(r.d, minT));
        L = add(P, multS(scene.point_lights[k].position, -1));
        N = normalize(add(center, multS(P, -1)));
        float r_square = dot(L,L);
        parser::Vec3f irradiance = multS(scene.point_lights[k].intensity, 1 / r_square);
        L= normalize(L);
        float theta = fmax(0, dot(L, N));
        parser::Vec3f specular_const = multS(irradiance,theta);
        parser::Vec3f specular = scene.materials[scene.spheres[minI].material_id - 1].diffuse;
        specular = multVects(specular,specular_const);
        result = add(result,specular);
    }

    return result;
}
parser::Vec3f computeDiffuse_triangle(ray r,int minI,float minT){
    parser::Vec3f L,P,N,result;
    result.x=result.y=result.z=0;
    for(int k=0;k<scene.point_lights.size();k++){
        P = add(r.o, multS(r.d, minT));
        L = add(scene.point_lights[k].position, multS(P, -1));
        int point_A_id=scene.triangles[minI].indices.v0_id-1;
        int point_B_id=scene.triangles[minI].indices.v1_id-1;
        int point_C_id=scene.triangles[minI].indices.v2_id-1;
        parser::Vec3f point_A=scene.vertex_data[point_A_id];
        parser::Vec3f point_B=scene.vertex_data[point_B_id];
        parser::Vec3f point_C=scene.vertex_data[point_C_id];
        N= normalize(crossProduct(add(point_B, multS(point_A,-1)),add(point_C, multS(point_A,-1))));
        float r_square = getLengthSquare(L);
        parser::Vec3f irradiance = multS(scene.point_lights[k].intensity, 1 / r_square);
        float theta = fmax(0, dot(normalize(L), N));
        parser::Vec3f specular_const = multS(irradiance,theta);
        parser::Vec3f specular = scene.materials[scene.triangles[minI].material_id - 1].diffuse;
        specular = multVects(specular,specular_const);
        result = add(result,specular);
    }

    return result;
}
parser::Vec3f computeDiffuse_mesh(ray r,int minI,float minT,int mesh_index){
    parser::Vec3f L,P,N,result;
    result.x=result.y=result.z=0;
    for(int k = 0;k < scene.point_lights.size();k++){
        P = add(r.o, multS(r.d, minT));
        L = add(scene.point_lights[k].position, multS(P, -1));
        int point_A_id=scene.meshes[minI].faces[mesh_index].v0_id - 1;
        int point_B_id=scene.meshes[minI].faces[mesh_index].v1_id - 1;
        int point_C_id=scene.meshes[minI].faces[mesh_index].v2_id - 1;
        parser::Vec3f point_A=scene.vertex_data[point_A_id];
        parser::Vec3f point_B=scene.vertex_data[point_B_id];
        parser::Vec3f point_C=scene.vertex_data[point_C_id];
        N= normalize(crossProduct(add(point_B, multS(point_A,-1)),add(point_C, multS(point_A,-1))));
        float r_square = getLengthSquare(L);
        parser::Vec3f irradiance = multS(scene.point_lights[k].intensity, 1 / r_square);
        float theta = fmax(0, dot(normalize(L), N));
        parser::Vec3f specular_const = multS(irradiance,theta);
        parser::Vec3f specular = scene.materials[scene.meshes[minI].material_id - 1].diffuse;
        specular = multVects(specular,specular_const);
        result = add(result,specular);
    }
    return result;
}
parser::Vec3f computeColor(ray r,int max_rec_depth)
{
    std::vector<parser::Vec3f> vertex_data = scene.vertex_data;
    int i;
    int num_of_triangle=scene.triangles.size();
    int num_of_meshes=scene.meshes.size();
    parser::Vec3f c,center;
    float minT_sphere = 90000,minT_tri=90000,minT_mesh = 90000; // some large number
    float t_sphere,t_triangle;
    std::vector<float> t_mesh;
    int minI_sphere,minI_tri,minI_mesh;
    int mesh_face_index;
    c.x=backgroundColor.x;
    c.y=backgroundColor.y;
    c.z=backgroundColor.z;
    minI_sphere = -1;
    minI_tri=-1;
    minI_mesh = -1;
    for (i=0;i<num_of_sphere;i++)
    {
        t_sphere = intersectSphere(r,i);
        center = vertex_data[scene.spheres[i].center_vertex_id-1];
        if (t_sphere<minT_sphere && t_sphere>0)
        {
            minI_sphere = i;
            minT_sphere = t_sphere;
        }
    }
    for (i=0;i<num_of_triangle;i++)
    {

        t_triangle = intersectTriangle(r,i);
        if (t_triangle<minT_tri && t_triangle>=0)
        {
            minI_tri = i;
            minT_tri = t_triangle;
        }
    }
    for (int i=0;i<num_of_meshes;i++)
    {

        t_mesh = intersectMesh(r,i);
        if (t_mesh[0]<minT_mesh && t_mesh[0]>=0)
        {
            minI_mesh = i;
            minT_mesh = t_mesh[0];
            mesh_face_index=t_mesh[1];
        }
    }

    if (minI_sphere != -1 && minT_sphere < minT_tri &&  minT_sphere < minT_mesh )
    {
        c.x =  scene.materials[scene.spheres[minI_sphere].material_id - 1].ambient.x * scene.ambient_light.x;
        c.y =  scene.materials[scene.spheres[minI_sphere].material_id - 1].ambient.y * scene.ambient_light.y;
        c.z =  scene.materials[scene.spheres[minI_sphere].material_id - 1].ambient.z * scene.ambient_light.z;
        c=add(c,computeDiffuse_sphere(r,center,minI_sphere,minT_sphere));
        c=add(c,computeSpecular_sphere(r,center,minI_sphere,minT_sphere));

    }
    if (minI_tri != -1 && minT_tri < minT_sphere &&  minT_tri < minT_mesh)
    {
        c.x =  scene.materials[scene.triangles[minI_tri].material_id - 1].ambient.x * scene.ambient_light.x;
        c.y =  scene.materials[scene.triangles[minI_tri].material_id - 1].ambient.y * scene.ambient_light.y;
        c.z =  scene.materials[scene.triangles[minI_tri].material_id - 1].ambient.z * scene.ambient_light.z;
        c=add(c,computeDiffuse_triangle(r,minI_tri,minT_tri));
        c=add(c,computeSpecular_triangle(r,minI_tri,minT_tri));

    }

    if (minI_mesh!=-1 && minT_mesh<minT_sphere &&  minT_mesh<minT_tri)
    {
        c.x =  scene.materials[scene.meshes[minI_mesh].material_id - 1].ambient.x * scene.ambient_light.x;
        c.y =  scene.materials[scene.meshes[minI_mesh].material_id - 1].ambient.y * scene.ambient_light.y;
        c.z =  scene.materials[scene.meshes[minI_mesh].material_id - 1].ambient.z * scene.ambient_light.z;
        c=add(c,computeDiffuse_mesh(r,minI_mesh,minT_mesh,mesh_face_index));
        c=add(c,computeSpecular_mesh(r,minI_mesh,minT_mesh,mesh_face_index));

    }
    return c;
}


