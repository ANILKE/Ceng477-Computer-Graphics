#include <iostream>
#include "parser.h"
#include "ppm.h"
#include <math.h>
#include "vector3f.h"
#include "intersect.h"

parser::Scene scene;
float left, right, bottom, top,width_su,height_sv;;
vector3f m, q,e, u,v,w;
parser::Vec3i backgroundColor;
int num_of_sphere;

typedef struct
{
    vector3f o,d;
} ray;

float getDeterminant(const vector3f &first_column, const vector3f &second_column, const vector3f &third_column)
{
    float result;
    result=first_column.x * (second_column.y*third_column.z - third_column.y*second_column.z)
           + first_column.y * (third_column.x*second_column.z - second_column.x*third_column.z)
           + first_column.z * (second_column.x*third_column.y - second_column.y*third_column.x);
    return result;
}

ray  generateRay(int i, int j)
{
    ray result;
    float su,sv;
    vector3f s;
    su=(j+0.5)*height_sv;
    sv=(i+0.5)*width_su;
    s=q+(u*su)-(v*sv);
    result.o = e;
    result.d= s-e;
    result.d=result.d.normalize();
    return result;
}

intersect intersectSphere(ray r,int index)
{
    intersect result;
    result.distance=-1;
    parser::Sphere sphere=scene.spheres[index];
    vector3f P;
    std::vector<parser::Vec3f> vertex_data = scene.vertex_data;
    float A,B,C; //constants for the quadratic equation

    float delta;

    int c_ID;

    c_ID = sphere.center_vertex_id;
    vector3f center;
    center.x = (vertex_data[c_ID - 1]).x;
    center.y = (vertex_data[c_ID - 1]).y;
    center.z = (vertex_data[c_ID - 1]).z;
    float t,t1,t2;

    C = (r.o.x-center.x)*(r.o.x-center.x)+(r.o.y-center.y)*(r.o.y-center.y)+(r.o.z-center.z)*(r.o.z-center.z)-sphere.radius*sphere.radius;

    B = 2*r.d.x*(r.o.x-center.x)+2*r.d.y*(r.o.y-center.y)+2*r.d.z*(r.o.z-center.z);

    A = r.d.x*r.d.x+r.d.y*r.d.y+r.d.z*r.d.z;

    delta = B*B-4*A*C;

    if (delta<0) return result;
    else if (delta==0)
    {
        t = -B / (2*A);
        result.distance=t;
        result.material_id=scene.spheres[index].material_id - 1;
        P=r.o+(r.d*t);
        result.intersection_point= P;
        result.normal=P - center;
        result.normal=result.normal.normalize();
    }
    else
    {
        delta = sqrtf(delta);
        A = 2*A;
        t1 = (B*(-1) + delta) / A;
        t2 = (B*(-1) - delta) / A;

        if (t1<t2){
            result.distance=t1;
            result.material_id=scene.spheres[index].material_id - 1;
            P=r.o+(r.d*t1);
            result.intersection_point= P;
            result.normal= P-center;
            result.normal=result.normal.normalize();
        }
        else {
            result.distance=t2;
            result.material_id=scene.spheres[index].material_id - 1;
            P=r.o+(r.d*t2);
            result.intersection_point= P;
            result.normal= P - center;
            result.normal=result.normal.normalize();
        }
    }

    return result;
}

intersect intersectTriangle(ray r,int index){
    intersect result;
    result.distance=-1;
    vector3f point_A,point_B,point_C,A_to_B,A_to_C,A_to_O;
    int point_A_id=scene.triangles[index].indices.v0_id-1;
    int point_B_id=scene.triangles[index].indices.v1_id-1;
    int point_C_id=scene.triangles[index].indices.v2_id-1;
    float beta,gamma,t,A_determinant;
    point_A.x=(scene.vertex_data[point_A_id]).x;
    point_A.y=(scene.vertex_data[point_A_id]).y;
    point_A.z=(scene.vertex_data[point_A_id]).z;

    point_B.x=(scene.vertex_data[point_B_id]).x;
    point_B.y=(scene.vertex_data[point_B_id]).y;
    point_B.z=(scene.vertex_data[point_B_id]).z;

    point_C.x = (scene.vertex_data[point_C_id]).x;
    point_C.y = (scene.vertex_data[point_C_id]).y;
    point_C.z = (scene.vertex_data[point_C_id]).z;

    A_to_B=point_A-point_B;
    A_to_C=point_A-point_C;
    A_to_O = point_A - r.o;
    A_determinant= getDeterminant(A_to_B,A_to_C,r.d);
    if(A_determinant==0){
        return result;
    }
    t= getDeterminant(A_to_B,A_to_C,A_to_O)/A_determinant;
    if(!(t>0)){
        return result;
    }
    beta= getDeterminant(A_to_O,A_to_C,r.d)/A_determinant;
    if(beta<0){
        return result;
    }
    gamma= getDeterminant(A_to_B,A_to_O,r.d)/A_determinant;
    if (gamma< 0)
        return result;
    if (gamma +  beta >1)
        return result;
    result.distance=t;
    result.material_id=scene.triangles[index].material_id-1;
    vector3f b_to_a,c_to_a,norm;
    b_to_a=point_B-point_A;
    c_to_a=point_C-point_A;
    norm=b_to_a;
    norm = norm.crossProduct(c_to_a);
    result.normal= norm.normalize();
    vector3f P=r.o+(r.d*t);
    result.intersection_point = P;
    return result;
}


intersect intersectMesh(ray r, int index){
    intersect result;
    result.distance=100000000.0;
    float t_new;
    vector3f A_to_B,A_to_C,A_to_O;
    float beta,gamma,A_determinant;
    int faceSize = scene.meshes[index].faces.size();
    for(int faceNumber = 0; faceNumber < faceSize;faceNumber++)
    {
        vector3f point_A;
        point_A.x= (scene.vertex_data[scene.meshes[index].faces[faceNumber].v0_id - 1]).x;
        point_A.y= (scene.vertex_data[scene.meshes[index].faces[faceNumber].v0_id - 1]).y;
        point_A.z= (scene.vertex_data[scene.meshes[index].faces[faceNumber].v0_id - 1]).z;
        vector3f point_B;
        point_B.x= (scene.vertex_data[scene.meshes[index].faces[faceNumber].v1_id - 1]).x;
        point_B.y= (scene.vertex_data[scene.meshes[index].faces[faceNumber].v1_id - 1]).y;
        point_B.z= (scene.vertex_data[scene.meshes[index].faces[faceNumber].v1_id - 1]).z;
        vector3f point_C;
        point_C.x = (scene.vertex_data[scene.meshes[index].faces[faceNumber].v2_id - 1]).x;
        point_C.y = (scene.vertex_data[scene.meshes[index].faces[faceNumber].v2_id - 1]).y;
        point_C.z = (scene.vertex_data[scene.meshes[index].faces[faceNumber].v2_id - 1]).z;

        A_to_B = point_A-point_B;
        A_to_C = point_A- point_C;
        A_to_O = point_A-r.o;

        A_determinant= getDeterminant(A_to_B,A_to_C,r.d);
        if(A_determinant==0){
            t_new= -1;
            continue;
        }
        t_new= getDeterminant(A_to_B,A_to_C,A_to_O)/A_determinant;
        if(!(t_new>0)){
            t_new= -1;
            continue;
        }
        beta= getDeterminant(A_to_O,A_to_C,r.d)/A_determinant;
        if(beta<0){
            t_new= -1;
            continue;
        }
        gamma= getDeterminant(A_to_B,A_to_O,r.d)/A_determinant;
        if (gamma< 0){
            t_new= -1;
            continue;
        }

        if (gamma +  beta >1){
            t_new = -1;
            continue;
        }

        if(t_new != -1 && t_new <= result.distance){
            result.distance=t_new;
            result.material_id=scene.meshes[index].material_id-1;
            vector3f B_to_A,C_to_A,norm;
            B_to_A = point_B-point_A;
            C_to_A = point_C-point_A;
            norm = B_to_A;
            norm = norm.crossProduct(C_to_A);
            result.normal = norm.normalize();
            vector3f P=r.o+(r.d*t_new);
            result.intersection_point = P;
        }

    }

    if(result.distance==100000000.0){
        result.distance=-1;
        return result;
    }

    return result;
}
vector3f computeDiffuse(ray r,intersect intersection, int i){
    vector3f L,P,N,result;
    P=r.o+(r.d*intersection.distance);
    int pointLightSize = scene.point_lights.size();

    L.x = (scene.point_lights[i].position).x - P.x;
    L.y = (scene.point_lights[i].position).y - P.y;
    L.z = (scene.point_lights[i].position).z - P.z;

    N = intersection.normal;
    float r_square = L.dotProduct(L);
    L= L.normalize();
    vector3f intensity;
    intensity.x=(scene.point_lights[i].intensity).x;
    intensity.y=(scene.point_lights[i].intensity).y;
    intensity.z=(scene.point_lights[i].intensity).z;
    vector3f irradiance = intensity / (r_square );
    float cos_theta= fmax(0,L.dotProduct(N));
    vector3f specular_const = irradiance * cos_theta;
    vector3f specular;
    specular.x =  scene.materials[intersection.material_id].diffuse.x;
    specular.y =  scene.materials[intersection.material_id].diffuse.y;
    specular.z =  scene.materials[intersection.material_id].diffuse.z;
    specular = specular * specular_const;
    result = result+specular;

    return result;
}


vector3f computeSpecular(ray r,intersect intersection, int i)
{

    vector3f L,P,N,H,result;
    P=r.o+(r.d*intersection.distance);
    int pointLightSize = scene.point_lights.size();


    L.x = (scene.point_lights[i].position).x - P.x;
    L.y = (scene.point_lights[i].position).y - P.y;
    L.z = (scene.point_lights[i].position).z - P.z;

    N = intersection.normal;
    vector3f temp=L;
    float r_square = temp.dotProduct(L);
    L= L.normalize();
    H = L - r.d;
    H = H.normalize();
    temp=N;
    float cos_alpha = fmax(0, temp.dotProduct(H));
    float cos_alpha_phong = pow(cos_alpha,scene.materials[intersection.material_id].phong_exponent);
    vector3f intensity;
    intensity.x=(scene.point_lights[i].intensity).x;
    intensity.y=(scene.point_lights[i].intensity).y;
    intensity.z=(scene.point_lights[i].intensity).z;
    vector3f irradiance = intensity / (r_square);
    temp=L;
    float cos_theta= fmax(0,L.dotProduct(N));
    float theta= acos(cos_theta);
    if(theta<90)
    {
        vector3f specular_const = irradiance*cos_alpha_phong;
        vector3f specular;
        specular.x =  scene.materials[intersection.material_id].specular.x;
        specular.y =  scene.materials[intersection.material_id].specular.y;
        specular.z =  scene.materials[intersection.material_id].specular.z;
        specular = specular * specular_const;
        result = result+specular;
    }

    return result;
}

float getDistance(parser::Vec3f v_1, vector3f v_2){
    float result;
    result = pow(v_1.x-v_2.x,2) + pow(v_1.y-v_2.y,2) + pow(v_1.z-v_2.z,2);
    result = sqrt(result);
    return result;
}
vector3f computeColor(ray r,int max_rec_depth)
{
    std::vector<parser::Vec3f> vertex_data = scene.vertex_data;
    int i;
    int num_of_triangle=scene.triangles.size();
    int num_of_meshes=scene.meshes.size();
    vector3f c,center;
    intersect min_inter,shadow_min_inter;
    intersect inter_sphere,inter_triangle,inter_mesh;
    min_inter.distance = 90000; // some large number
    c.x=backgroundColor.x;
    c.y=backgroundColor.y;
    c.z=backgroundColor.z;
    for (i=0;i<num_of_sphere;i++)
    {
        inter_sphere = intersectSphere(r,i);
        if (inter_sphere.distance<min_inter.distance && inter_sphere.distance>0)
        {
            min_inter = inter_sphere;
        }
    }
    for (i=0;i<num_of_triangle;i++)
    {
        inter_triangle = intersectTriangle(r,i);
        if (inter_triangle.distance<min_inter.distance && inter_triangle.distance>0)
        {
            min_inter=inter_triangle;
        }
    }
    for (i=0;i<num_of_meshes;i++)
    {

        inter_mesh = intersectMesh(r,i);
        if (inter_mesh.distance<min_inter.distance && inter_mesh.distance>0)
        {
            min_inter = inter_mesh;
        }
    }
    if(min_inter.distance==90000){
        min_inter.distance=-1;
    }
    if (min_inter.distance != -1 )
    {
        c.x =  scene.materials[min_inter.material_id].ambient.x * scene.ambient_light.x;
        c.y =  scene.materials[min_inter.material_id].ambient.y * scene.ambient_light.y;
        c.z =  scene.materials[min_inter.material_id].ambient.z * scene.ambient_light.z;
        int light_count = scene.point_lights.size();
        bool flag= false;
        for(int a=0;a < light_count; a++){
            float min_distance_for_light = getDistance(scene.point_lights[a].position, min_inter.intersection_point);
            shadow_min_inter.distance = min_distance_for_light;
            ray shadow_ray;


            // NOT EPSİLON YUANLIŞ YERE EKLENİYOR.

            shadow_ray.d.x = scene.point_lights[a].position.x - min_inter.intersection_point.x;
            shadow_ray.d.y = scene.point_lights[a].position.y - min_inter.intersection_point.y;
            shadow_ray.d.z = scene.point_lights[a].position.z - min_inter.intersection_point.z;
            shadow_ray.d= shadow_ray.d.normalize();

            vector3f epsilon = min_inter.normal * scene.shadow_ray_epsilon;



            shadow_ray.o = min_inter.intersection_point + epsilon;
            for (i=0;i<num_of_sphere;i++)
            {
                inter_sphere = intersectSphere(shadow_ray,i);
                if (inter_sphere.distance<shadow_min_inter.distance && inter_sphere.distance>0 && inter_sphere.material_id!=min_inter.material_id)
                {
                    goto A;
                }
            }
            for (i=0;i<num_of_triangle;i++)
            {
                inter_triangle = intersectTriangle(shadow_ray,i);
                if (inter_triangle.distance< shadow_min_inter.distance && inter_triangle.distance>0 && inter_triangle.material_id!=min_inter.material_id)
                {
                    goto A;
                }
            }
            for (i=0;i<num_of_meshes;i++)
            {
                inter_mesh = intersectMesh(shadow_ray,i);
                if (inter_mesh.distance< shadow_min_inter.distance && inter_mesh.distance>0 && inter_mesh.material_id!=min_inter.material_id)
                {
                    goto A;
                }
            }
            c = c + computeDiffuse(r, min_inter,a);
            c = c + computeSpecular(r, min_inter,a);
            A:
            int b=5;
        }


        if(max_rec_depth && (scene.materials[min_inter.material_id].is_mirror)){

            // Shadowda kalanları almaması lazım.

            vector3f w0 = (r.d * -1).normalize();
            vector3f n = min_inter.normal;
            float n_dot_w0 = n.dotProduct(w0) * 2;

            vector3f epsilon = min_inter.normal*scene.shadow_ray_epsilon;
            vector3f _2nCosTheta = n * n_dot_w0;


            vector3f w_r = (w0 * -1);
            w_r = w_r + _2nCosTheta;
            w_r = w_r.normalize();


            ray mirror_ray;

            mirror_ray.o = min_inter.intersection_point + epsilon;
            mirror_ray.d  = w_r;



            vector3f mirror_color=computeColor(mirror_ray,max_rec_depth-1);
            c.x = c.x + mirror_color.x * scene.materials[min_inter.material_id].mirror.x;
            c.y = c.y + mirror_color.y * scene.materials[min_inter.material_id].mirror.y;
            c.z = c.z + mirror_color.z * scene.materials[min_inter.material_id].mirror.z;
        }
    }
    return c;
}


int main(int argc, char* argv[])
{
    // Sample usage for reading an XML scene file
    int height, width;
    scene.loadFromXml(argv[1]);
    num_of_sphere=scene.spheres.size();
    int camera_count= scene.cameras.size();
    backgroundColor=scene.background_color;
    for(int c = 0; c <camera_count;c++){
        parser::Camera currCamera = scene.cameras[c];
        width = currCamera.image_width;
        height = currCamera.image_height;

        unsigned char* image = new unsigned char [width * height * 3];

        left=currCamera.near_plane.x;
        right=currCamera.near_plane.y;
        bottom=currCamera.near_plane.z;
        top=currCamera.near_plane.w;

        width_su=(right-left)/width;
        height_sv =(top-bottom)/height;

        v.x = currCamera.up.x;
        v.y = currCamera.up.y;
        v.z = currCamera.up.z;
        v=v.normalize();

        e.x = currCamera.position.x;
        e.y = currCamera.position.y;
        e.z = currCamera.position.z;
        vector3f gaze;
        gaze.x = currCamera.gaze.x;
        gaze.y = currCamera.gaze.y;
        gaze.z = currCamera.gaze.z;
        gaze=gaze.normalize();
        w = gaze*(-1);
        u= gaze.crossProduct(v);
        u= u.normalize();
        m=(gaze*currCamera.near_distance)+e;
        q=m+(u*left)+(v*top);
        float ratio = (float) width/height;
        int idx = 0;

            for (int j = 0; j < height; j++) {
                for (int i=0;i<width;i++) {
                ray currRay = generateRay(j, i);
                vector3f currRayColor;
                currRayColor = computeColor(currRay, scene.max_recursion_depth);
                if (currRayColor.x > 255) {
                    image[idx++] = 255;
                } else {
                    image[idx++] = (int) (currRayColor.x);
                }
                if (currRayColor.y > 255) {
                    image[idx++] = 255;
                } else {
                    image[idx++] = (int) (currRayColor.y);
                }
                if (currRayColor.z > 255) {
                    image[idx++] = 255;
                } else {
                    image[idx++] = (int) (currRayColor.z);
                }
            }
        }
        const char* file_name = currCamera.image_name.c_str();
        write_ppm(file_name, image, width, height);
    }

    return 0;

}
