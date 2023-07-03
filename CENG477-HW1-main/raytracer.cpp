#include <iostream>
#include "parser.h"
#include "ppm.h"
#include <math.h>



typedef unsigned char RGB[3];



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



int main(int argc, char* argv[])
{
    // Sample usage for reading an XML scene file
    parser::Scene scene;

    scene.loadFromXml(argv[1]);

    // The code below creates a test pattern and writes
    // it to a PPM file to demonstrate the usage of the
    // ppm_write function.
    //
    // Normally, you would be running your ray tracing
    // code here to produce the desired image.
    const RGB BAR_COLOR[8] =
    {
        { 255, 255, 255 },  // 100% White
        { 255, 255,   0 },  // Yellow
        {   0, 255, 255 },  // Cyan
        {   0, 255,   0 },  // Green
        { 255,   0, 255 },  // Magenta
        { 255,   0,   0 },  // Red
        {   0,   0, 255 },  // Blue
        {   0,   0,   0 },  // Black
    };

    int width = 640, height = 480;
    int columnWidth = width / 8;

    unsigned char* image = new unsigned char [width * height * 3];

    int i = 0;
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            int colIdx = x / columnWidth;
            image[i++] = BAR_COLOR[colIdx][0];
            image[i++] = BAR_COLOR[colIdx][1];
            image[i++] = BAR_COLOR[colIdx][2];
        }
    }

    write_ppm("test.ppm", image, width, height);
    return 0;
}
