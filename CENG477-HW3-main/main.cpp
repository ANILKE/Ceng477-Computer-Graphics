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
#include "gridController.h"
#include "item.h"

#include FT_FREETYPE_H

#define BUFFER_OFFSET(i) ((char*)NULL + (i))

using namespace std;

GLuint gProgram[3];
GLint gIntensityLoc;
float gIntensity = 250;
int gWidth = 640, gHeight = 600;
bool mousePressed = false;
int moves = 0;
int score = 0;
bool flagger = true;
gridController * anil =NULL;
// Create 5 Colors

bool restart = false;
bool keyboardControl = false;


struct Vertex
{
    Vertex(GLfloat inX, GLfloat inY, GLfloat inZ) : x(inX), y(inY), z(inZ) { }
    GLfloat x, y, z;
};

struct Texture
{
    Texture(GLfloat inU, GLfloat inV) : u(inU), v(inV) { }
    GLfloat u, v;
};

struct Normal
{
    Normal(GLfloat inX, GLfloat inY, GLfloat inZ) : x(inX), y(inY), z(inZ) { }
    GLfloat x, y, z;
};

struct Face
{
	Face(int v[], int t[], int n[]) {
		vIndex[0] = v[0];
		vIndex[1] = v[1];
		vIndex[2] = v[2];
		tIndex[0] = t[0];
		tIndex[1] = t[1];
		tIndex[2] = t[2];
		nIndex[0] = n[0];
		nIndex[1] = n[1];
		nIndex[2] = n[2];
	}
    GLuint vIndex[3], tIndex[3], nIndex[3];
};

vector<Vertex> gVertices;
vector<Texture> gTextures;
vector<Normal> gNormals;
vector<Face> gFaces;

GLuint gVertexAttribBuffer, gTextVBO, gIndexBuffer;
GLint gInVertexLoc, gInNormalLoc;
int gVertexDataSizeInBytes, gNormalDataSizeInBytes;
double mouseX = -1 ,mouseY = - 1;
int mouseIndexI = -1, mouseIndexJ = -1;
double colInc,rowInc;


/// Holds all state information relevant to a character as loaded using FreeType
struct Character {
    GLuint TextureID;   // ID handle of the glyph texture
    glm::ivec2 Size;    // Size of glyph
    glm::ivec2 Bearing;  // Offset from baseline to left/top of glyph
    GLuint Advance;    // Horizontal offset to advance to next glyph
};

std::map<GLchar, Character> Characters;

bool ParseObj(const string& fileName)
{
    fstream myfile;

    // Open the input 
    myfile.open(fileName.c_str(), std::ios::in);

    if (myfile.is_open())
    {
        string curLine;

        while (getline(myfile, curLine))
        {
            stringstream str(curLine);
            GLfloat c1, c2, c3;
            GLuint index[9];
            string tmp;

            if (curLine.length() >= 2)
            {
                if (curLine[0] == '#') // comment
                {
                    continue;
                }
                else if (curLine[0] == 'v')
                {
                    if (curLine[1] == 't') // texture
                    {
                        str >> tmp; // consume "vt"
                        str >> c1 >> c2;
                        gTextures.push_back(Texture(c1, c2));
                    }
                    else if (curLine[1] == 'n') // normal
                    {
                        str >> tmp; // consume "vn"
                        str >> c1 >> c2 >> c3;
                        gNormals.push_back(Normal(c1, c2, c3));
                    }
                    else // vertex
                    {
                        str >> tmp; // consume "v"
                        str >> c1 >> c2 >> c3;
                        gVertices.push_back(Vertex(c1, c2, c3));
                    }
                }
                else if (curLine[0] == 'f') // face
                {
                    str >> tmp; // consume "f"
					char c;
					int vIndex[3],  nIndex[3], tIndex[3];
					str >> vIndex[0]; str >> c >> c; // consume "//"
					str >> nIndex[0]; 
					str >> vIndex[1]; str >> c >> c; // consume "//"
					str >> nIndex[1]; 
					str >> vIndex[2]; str >> c >> c; // consume "//"
					str >> nIndex[2]; 

					assert(vIndex[0] == nIndex[0] &&
						   vIndex[1] == nIndex[1] &&
						   vIndex[2] == nIndex[2]); // a limitation for now

					// make indices start from 0
					for (int c = 0; c < 3; ++c)
					{
						vIndex[c] -= 1;
						nIndex[c] -= 1;
						tIndex[c] -= 1;
					}

                    gFaces.push_back(Face(vIndex, tIndex, nIndex));
                }
                else
                {
                    cout << "Ignoring unidentified line in obj file: " << curLine << endl;
                }
            }

            //data += curLine;
            if (!myfile.eof())
            {
                //data += "\n";
            }
        }

        myfile.close();
    }
    else
    {
        return false;
    }

	assert(gVertices.size() == gNormals.size());

    return true;
}

bool ReadDataFromFile(
    const string& fileName, ///< [in]  Name of the shader file
    string&       data)     ///< [out] The contents of the file
{
    fstream myfile;

    // Open the input 
    myfile.open(fileName.c_str(), std::ios::in);

    if (myfile.is_open())
    {
        string curLine;

        while (getline(myfile, curLine))
        {
            data += curLine;
            if (!myfile.eof())
            {
                data += "\n";
            }
        }

        myfile.close();
    }
    else
    {
        return false;
    }

    return true;
}

void createVS(GLuint& program, const string& filename)
{
    string shaderSource;

    if (!ReadDataFromFile(filename, shaderSource))
    {
        cout << "Cannot find file name: " + filename << endl;
        exit(-1);
    }

    GLint length = shaderSource.length();
    const GLchar* shader = (const GLchar*) shaderSource.c_str();

    GLuint vs = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vs, 1, &shader, &length);
    glCompileShader(vs);

    char output[1024] = {0};
    glGetShaderInfoLog(vs, 1024, &length, output);
    printf("VS compile log: %s\n", output);

    glAttachShader(program, vs);
}

void createFS(GLuint& program, const string& filename)
{
    string shaderSource;

    if (!ReadDataFromFile(filename, shaderSource))
    {
        cout << "Cannot find file name: " + filename << endl;
        exit(-1);
    }

    GLint length = shaderSource.length();
    const GLchar* shader = (const GLchar*) shaderSource.c_str();

    GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fs, 1, &shader, &length);
    glCompileShader(fs);

    char output[1024] = {0};
    glGetShaderInfoLog(fs, 1024, &length, output);
    printf("FS compile log: %s\n", output);

    glAttachShader(program, fs);
}

void initShaders()
{
    gProgram[0] = glCreateProgram();
    gProgram[1] = glCreateProgram();
    gProgram[2] = glCreateProgram();

    createVS(gProgram[0], "vert0.glsl");
    createFS(gProgram[0], "frag0.glsl");

    createVS(gProgram[1], "vert1.glsl");
    createFS(gProgram[1], "frag1.glsl");

    createVS(gProgram[2], "vert_text.glsl");
    createFS(gProgram[2], "frag_text.glsl");

    glBindAttribLocation(gProgram[0], 0, "inVertex");
    glBindAttribLocation(gProgram[0], 1, "inNormal");
    glBindAttribLocation(gProgram[1], 0, "inVertex");
    glBindAttribLocation(gProgram[1], 1, "inNormal");
    glBindAttribLocation(gProgram[2], 2, "vertex");

    glLinkProgram(gProgram[0]);
    glLinkProgram(gProgram[1]);
    glLinkProgram(gProgram[2]);
    glUseProgram(gProgram[0]);

    gIntensityLoc = glGetUniformLocation(gProgram[0], "intensity");
    cout << "gIntensityLoc = " << gIntensityLoc << endl;
    glUniform1f(gIntensityLoc, gIntensity);
}

void initVBO()
{
    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);
    assert(glGetError() == GL_NONE);

    glGenBuffers(1, &gVertexAttribBuffer);
    glGenBuffers(1, &gIndexBuffer);

    assert(gVertexAttribBuffer > 0 && gIndexBuffer > 0);

    glBindBuffer(GL_ARRAY_BUFFER, gVertexAttribBuffer);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gIndexBuffer);

    gVertexDataSizeInBytes = gVertices.size() * 3 * sizeof(GLfloat);
    gNormalDataSizeInBytes = gNormals.size() * 3 * sizeof(GLfloat);
    int indexDataSizeInBytes = gFaces.size() * 3 * sizeof(GLuint);
    GLfloat* vertexData = new GLfloat [gVertices.size() * 3];
    GLfloat* normalData = new GLfloat [gNormals.size() * 3];
    GLuint* indexData = new GLuint [gFaces.size() * 3];

    float minX = 1e6, maxX = -1e6;
    float minY = 1e6, maxY = -1e6;
    float minZ = 1e6, maxZ = -1e6;

    for (int i = 0; i < gVertices.size(); ++i)
    {
        vertexData[3*i] = gVertices[i].x;
        vertexData[3*i+1] = gVertices[i].y;
        vertexData[3*i+2] = gVertices[i].z;

        minX = std::min(minX, gVertices[i].x);
        maxX = std::max(maxX, gVertices[i].x);
        minY = std::min(minY, gVertices[i].y);
        maxY = std::max(maxY, gVertices[i].y);
        minZ = std::min(minZ, gVertices[i].z);
        maxZ = std::max(maxZ, gVertices[i].z);
    }

    std::cout << "minX = " << minX << std::endl;
    std::cout << "maxX = " << maxX << std::endl;
    std::cout << "minY = " << minY << std::endl;
    std::cout << "maxY = " << maxY << std::endl;
    std::cout << "minZ = " << minZ << std::endl;
    std::cout << "maxZ = " << maxZ << std::endl;

    for (int i = 0; i < gNormals.size(); ++i)
    {
        normalData[3*i] = gNormals[i].x;
        normalData[3*i+1] = gNormals[i].y;
        normalData[3*i+2] = gNormals[i].z;
    }

    for (int i = 0; i < gFaces.size(); ++i)
    {
        indexData[3*i] = gFaces[i].vIndex[0];
        indexData[3*i+1] = gFaces[i].vIndex[1];
        indexData[3*i+2] = gFaces[i].vIndex[2];
    }


    glBufferData(GL_ARRAY_BUFFER, gVertexDataSizeInBytes + gNormalDataSizeInBytes, 0, GL_STATIC_DRAW);
    glBufferSubData(GL_ARRAY_BUFFER, 0, gVertexDataSizeInBytes, vertexData);
    glBufferSubData(GL_ARRAY_BUFFER, gVertexDataSizeInBytes, gNormalDataSizeInBytes, normalData);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indexDataSizeInBytes, indexData, GL_STATIC_DRAW);

    // done copying; can free now
    delete[] vertexData;
    delete[] normalData;
    delete[] indexData;

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(gVertexDataSizeInBytes));

}

void initFonts(int windowWidth, int windowHeight)
{
    // Set OpenGL options
    //glEnable(GL_CULL_FACE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glm::mat4 projection = glm::ortho(0.0f, static_cast<GLfloat>(windowWidth), 0.0f, static_cast<GLfloat>(windowHeight));
    glUseProgram(gProgram[2]);
    glUniformMatrix4fv(glGetUniformLocation(gProgram[2], "projection"), 1, GL_FALSE, glm::value_ptr(projection));

    // FreeType
    FT_Library ft;
    // All functions return a value different than 0 whenever an error occurred
    if (FT_Init_FreeType(&ft))
    {
        std::cout << "ERROR::FREETYPE: Could not init FreeType Library" << std::endl;
    }

    // Load font as face
    FT_Face face;
    if (FT_New_Face(ft, "/usr/share/fonts/truetype/liberation/LiberationSerif-Italic.ttf", 0, &face))
    {
        std::cout << "ERROR::FREETYPE: Failed to load font" << std::endl;
    }

    // Set size to load glyphs as
    FT_Set_Pixel_Sizes(face, 0, 48);

    // Disable byte-alignment restriction
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1); 

    // Load first 128 characters of ASCII set
    for (GLubyte c = 0; c < 128; c++)
    {
        // Load character glyph 
        if (FT_Load_Char(face, c, FT_LOAD_RENDER))
        {
            std::cout << "ERROR::FREETYTPE: Failed to load Glyph" << std::endl;
            continue;
        }
        // Generate texture
        GLuint texture;
        glGenTextures(1, &texture);
        glBindTexture(GL_TEXTURE_2D, texture);
        glTexImage2D(
                GL_TEXTURE_2D,
                0,
                GL_RED,
                face->glyph->bitmap.width,
                face->glyph->bitmap.rows,
                0,
                GL_RED,
                GL_UNSIGNED_BYTE,
                face->glyph->bitmap.buffer
                );
        // Set texture options
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        // Now store character for later use
        Character character = {
            texture,
            glm::ivec2(face->glyph->bitmap.width, face->glyph->bitmap.rows),
            glm::ivec2(face->glyph->bitmap_left, face->glyph->bitmap_top),
            face->glyph->advance.x
        };
        Characters.insert(std::pair<GLchar, Character>(c, character));
    }

    glBindTexture(GL_TEXTURE_2D, 0);
    // Destroy FreeType once we're finished
    FT_Done_Face(face);
    FT_Done_FreeType(ft);

    //
    // Configure VBO for texture quads
    //
    glGenBuffers(1, &gTextVBO);
    glBindBuffer(GL_ARRAY_BUFFER, gTextVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * 6 * 4, NULL, GL_DYNAMIC_DRAW);

    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), 0);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void init(std::string fileName) 
{
	ParseObj(fileName);
    glEnable(GL_DEPTH_TEST);
    initShaders();
    initFonts(gWidth, gHeight);
    initVBO();
}

void drawModel()
{
	glBindBuffer(GL_ARRAY_BUFFER, gVertexAttribBuffer);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gIndexBuffer);

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(gVertexDataSizeInBytes));

	glDrawElements(GL_TRIANGLES, gFaces.size() * 3, GL_UNSIGNED_INT, 0);
}

void renderText(const std::string& text, GLfloat x, GLfloat y, GLfloat scale, glm::vec3 color)
{
    // Activate corresponding render state	
    glUseProgram(gProgram[2]);
    glUniform3f(glGetUniformLocation(gProgram[2], "textColor"), color.x, color.y, color.z);
    glActiveTexture(GL_TEXTURE0);

    // Iterate through all characters
    std::string::const_iterator c;
    for (c = text.begin(); c != text.end(); c++) 
    {
        Character ch = Characters[*c];

        GLfloat xpos = x + ch.Bearing.x * scale;
        GLfloat ypos = y - (ch.Size.y - ch.Bearing.y) * scale;

        GLfloat w = ch.Size.x * scale;
        GLfloat h = ch.Size.y * scale;

        // Update VBO for each character
        GLfloat vertices[6][4] = {
            { xpos,     ypos + h,   0.0, 0.0 },            
            { xpos,     ypos,       0.0, 1.0 },
            { xpos + w, ypos,       1.0, 1.0 },

            { xpos,     ypos + h,   0.0, 0.0 },
            { xpos + w, ypos,       1.0, 1.0 },
            { xpos + w, ypos + h,   1.0, 0.0 }           
        };

        // Render glyph texture over quad
        glBindTexture(GL_TEXTURE_2D, ch.TextureID);

        // Update content of VBO memory
        glBindBuffer(GL_ARRAY_BUFFER, gTextVBO);
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices), vertices); // Be sure to use glBufferSubData and not glBufferData

        //glBindBuffer(GL_ARRAY_BUFFER, 0);

        // Render quad
        glDrawArrays(GL_TRIANGLES, 0, 6);
        // Now advance cursors for next glyph (note that advance is number of 1/64 pixels)

        x += (ch.Advance >> 6) * scale; // Bitshift by 6 to get value in pixels (2^6 = 64 (divide amount of 1/64th pixels by 64 to get amount of pixels))
    }

    glBindTexture(GL_TEXTURE_2D, 0);
}



bool isObjectExist(int mouseIndexI, int mouseIndexJ, std::vector<item *> & items)
{
    if (mouseIndexI < 0 || mouseIndexJ < 0)
    {
        return false;
    }
    for(auto iter = items.begin();iter != items.end();iter++)
    {
        if((*iter)->compareIdx(mouseIndexI,mouseIndexJ))
        {
            return true;
        }
    }
    return false;
}
bool isInVector(int i, int j, std::vector< pair<int, int>> vect){
    for(int k = 0; k<vect.size();k++){
        if(vect[k].first == i && vect[k].second == j){
            return true;
        }
    }
    return false;
} 
gridController * getMatchingColors(std::vector<item *> & items, int gridWidth, int gridHeight){
    gridController * result =new gridController(gridHeight,gridWidth);
    if(flagger){  //silme daha önceden başlamadıysa
        for(int i = 0; i<gridHeight ;i++){
            for(int j = 0; j<gridWidth ;j++){
                glm::vec3 curr_color = items[i*gridHeight+j]->getColor(); //i j color al
                int same_count = 0;
                for(int k = j+1; k<gridWidth;k++){ //j+1 den başla sağa doğru arka arkaya olan aynı renkleri say
                    glm::vec3 next_color = items[i*gridHeight+k]->getColor();
                    if(curr_color[0] == next_color[0] &&curr_color[1] == next_color[1] &&curr_color[2] == next_color[2] ){
                        same_count += 1;
                    }
                    else{
                        break;
                    }
                }
                if(same_count>=2){ //aynıı renk 2 den büyükse i,j den sonra listeye at indexlerini
                    for(int k = 0 ; k<=same_count;k++){
                        
                        if(! isInVector(j+k,i, *(result->getLocationList()))){
                            result->addLocation(j+k,i);
                        }
                    }
                }
                same_count = 0;
                for(int k = i+1; k<gridHeight;k++){  //i+1 den aşşağı doğru same color ara gerisi aynı
                    glm::vec3 next_color = items[k*gridHeight+j]->getColor();
                    if(curr_color[0] == next_color[0] &&curr_color[1] == next_color[1] &&curr_color[2] == next_color[2] ){
                        same_count += 1;
                    }
                    else{
                        break;
                    }
                }
                if(same_count>=2){

                    for(int k = 0 ; k<=same_count;k++){

                        if(! isInVector(j,i+k, *(result->getLocationList()))){
                            
                            result->addLocation(j,i+k);
                        }
                    }
                }
            }
        }
    }
    return result;
}

void setMatrixConfig(const glm::mat4 orthoMat, std::vector<item *>::iterator iter)
{
    glUniformMatrix4fv(glGetUniformLocation(gProgram[0], "modelingMat"), 1, GL_FALSE, glm::value_ptr((*iter)->getModelMat()));
    glUniformMatrix4fv(glGetUniformLocation(gProgram[0], "modelingMatInvTr"), 1, GL_FALSE, glm::value_ptr((*iter)->getModelMatInv()));
    glUniformMatrix4fv(glGetUniformLocation(gProgram[0], "orthoMat"), 1, GL_FALSE, glm::value_ptr(orthoMat));
    glUniform3fv(glGetUniformLocation(gProgram[0], "kd"), 1,glm::value_ptr((*iter)->getColor()));
}

void objectBomber(std::vector< pair<int, int>> * delete_list )
{






}




gridController * display(std::vector<item *> & items, int gridWidth, int gridHeight,gridController*  delete_color)
{   

    glClearColor(0, 0, 0, 1);
    glClearDepth(1.0f);
    glClearStencil(0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

    glm::mat4 orthoMat = glm::ortho(-10.0f, 10.0f, -10.0f, 10.0f, -20.0f, 20.0f);
    
    glUseProgram(gProgram[0]);
    if(delete_color == NULL || flagger ){ //silme başlamadıysa ya da hiç grid kontrolo yapılmadıysa aynı renkler için
        delete_color = getMatchingColors(items, gridWidth, gridHeight);  //grid kontrolü yap aynı renkler için
        
    }
    std::vector< pair<int, int>> * delete_list = NULL;
    if(delete_color != NULL){
        delete_list = delete_color->getLocationList();   //silinmesi gereken objelerin i j değerleri pair
    }


    if(!isObjectExist(mouseIndexI, mouseIndexJ, items) && mousePressed)
        mousePressed = false;
    


    for(auto iter = items.begin();iter != items.end();iter++)
    {   

        if((*iter)->checkIsAlive()){ //silinen ve görünmez olan itemleri drawlama

            if(mouseIndexI >= 0 && mouseIndexJ >= 0 && (*iter)->compareIdx(mouseIndexI,mouseIndexJ))
            {
                if((*iter) ->shouldDeleted())
                {
                    mouseIndexI = -1;
                    mouseIndexJ = -1;

                    cout << "Mouse set to false." << endl;
                    mousePressed = false;


                    items.erase(iter);                    
                    if(iter != items.begin())
                        iter--;
                    
                    score++;
                }
                else
                {
                    (*iter)->changeScaleMat();

                }
            }
            
            int tmpScore = score;
            //silinecek objeler var ise gir 

            for( auto colorIter = delete_list->begin(); colorIter != delete_list->end();colorIter++){ //iterate over same color vector list
                if((*iter)->compareIdx(colorIter->first,colorIter->second)) //if current iter is where we are looking for(silinecek objenin kordinatı)
                {
                    if((*iter) ->shouldDeleted())
                    {
                        (*iter)->setAlive(false);
                        
                        delete_list->erase(colorIter);  

                        if(colorIter != delete_list->begin()) 
                            colorIter--; 

                        tmpScore += 1;
                        break;
                    }
                    else
                    {
                        (*iter)->changeScaleMat(); //büyüt
                        flagger = false; //when deletetion of same objects start do not do new grid check before slide occude
                    }
                }
            }

            // Delete the items that are marked as dead.
            for (auto iter1 = items.begin(); iter1 != items.end(); iter1++)
            {
                if (!(*iter1)->checkIsAlive())
                {
                    items.erase(iter1);
                    if (iter1 != items.begin())
                    {
                        iter1--;
                        iter--;
                    }
                }
            }

            
            score = tmpScore;
            setMatrixConfig(orthoMat, iter);

            drawModel();
            (*iter)->incrementAngle();
        }
    }

    glUseProgram(gProgram[1]);
    assert(glGetError() == GL_NO_ERROR);
    std::string var = "Moves: " + std::to_string(moves) + "Score: " + std::to_string(score);  
    renderText(var, 0, 0, 1, glm::vec3(0, 1, 1));



    assert(glGetError() == GL_NO_ERROR);
    return delete_color;
}

void reshape(GLFWwindow* window, int w, int h)
{
    w = w < 1 ? 1 : w;
    h = h < 1 ? 1 : h;

    gWidth = w;
    gHeight = h;

    glViewport(0, 0, w, h);
}

void keyboard(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (keyboardControl)
        return;
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS )
    {
        glfwSetWindowShouldClose(window, GL_TRUE);
    }
    else if (key == GLFW_KEY_R && action == GLFW_PRESS)
    {
        cout << "Game is restarting." << endl;
        restart = true;
        score = 0;
        moves = 0;
        glUseProgram(gProgram[1]);

    }
    
}

void mouse(GLFWwindow* window, int button, int action, int mods)
{


    /* MOUSE A BASINCA İNDEXLERDEN ÖTÜRÜ BİR YANKDAİNE GİDİP BİR ADET BÜYÜTÜP SONRA SIFIRLIYOR*/
    // DÜZELTİLECEK (DÜZELTİLDİ.)

    // ŞİMDİ İSE MOSEU BASILDIKTAN SONRA MOUSEPRESSED FALSE OLMADIĞI İÇİN 
    // bazı yerlerde sıkıntı çıkartıyor.
    // DÜZELTİLECEK.

    if(mousePressed)
    {
        printf("mousePressed Before = True\n");
    }
    else
    {
        printf("mousePressed Before= False\n");
    }

    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS && !mousePressed)
    {
        cout << "Left mouse button pressed" << endl;
        printf("mousePressed is set to  True\n");

        //getting cursor position
        glfwGetCursorPos(window, &mouseX, &mouseY);

        cout << "xpos = " << mouseX << " ypos = " << mouseX << endl;
        mouseIndexI = mouseX / colInc;
        mouseIndexJ = mouseY / rowInc;
        mousePressed = true;
        moves += 1;
        // Orda bir obje olması lazım fakat yok bulamadığı için patlıyor.

    }
}


std::vector<item *> createItems(int gridHeight,int gridWidth)
{
    vector<glm::vec3> colors = {glm::vec3(1.0f,0.0f,0.f),
                            glm::vec3(0.0f,1.0f,0.f),
                            glm::vec3(0.0f,0.0f,1.f),
                            glm::vec3(1.0f,1.0f,0.f),
                            glm::vec3(1.0f,0.0f,1.f)
                        };

    const glm::mat4 initialT = glm::translate(glm::mat4(1.f), glm::vec3(10.0f/gridWidth - 10.0f, 10.0f - 9.5f / gridHeight + 3.f, 0));
    const float scaleFact = min(1.0f* ( 5.0f / gridWidth), 1.0f *(5.0f / gridHeight));
    glm::mat4 initS = glm::scale(glm::mat4(1.f), glm::vec3(scaleFact, scaleFact ,scaleFact));
    std::vector<item *> items;
    for(int i = 0; i < gridHeight;i++)
    {
        for(int j = 0; j < gridWidth;j++)
        {   
            glm::mat4 T = initialT * glm::translate(glm::mat4(1.f), glm::vec3(j * 20.0f / gridWidth, -i * (19.f / gridHeight), 0));
            int randomColorIdx = ((int) rand()) % 5;
            glm::mat4 R = glm::rotate(glm::mat4(1.f), glm::radians(0.0f), glm::vec3(0, 1, 0));
            item * newItem = new item(j,i,colors[randomColorIdx],T,R,initS,scaleFact,true);
            items.push_back(newItem);
        }
    }
    return items;
}


void resetTheGame(std::vector<item *> & items,int gridHeight,int gridWidth)
{
    for(int i = 0; i < items.size();i++)
        delete items[i];

    items = createItems(gridHeight,gridWidth);
    mousePressed = false;
    restart = false;
    mouseIndexI = -1;
    mouseIndexJ = -1;
    keyboardControl = false;
    flagger = true;
}

void mainLoop(GLFWwindow* window,int gridHeight,int gridWidth)
{

    std::vector<item *> items = createItems(gridHeight,gridWidth);

    while (!glfwWindowShouldClose(window))
    {
        if(restart == true)
            resetTheGame(items,gridHeight,gridWidth);


        anil = display(items, gridWidth, gridHeight,anil); //since when slide is not occure yet gridController needs to be same as before we are calling with it.
        glfwSwapBuffers(window);
        glfwPollEvents();
        
    }
}

int main(int argc, char** argv)   // Create Main Function For Bringing It All Together
{

    srand(time(NULL));

    if(argc < 4)
    {
        cout << "Usage: " << argv[0] << " <gridHeight> <gridWidth> <inputFile>" << endl;
        return 1;
    }

    int gridHeight = atoi(argv[1]);
    int gridWidth = atoi(argv[2]);
    string fileName = argv[3];


    rowInc = gHeight / gridHeight;
    colInc = gWidth / gridWidth;

    GLFWwindow* window;
    if (!glfwInit())
    {
        exit(-1);
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    //glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_COMPAT_PROFILE);
    //glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    window = glfwCreateWindow(gWidth, gHeight, "Simple Example", NULL, NULL);

    if (!window)
    {
        glfwTerminate();
        exit(-1);
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    // Initialize GLEW to setup the OpenGL Function pointers
    if (GLEW_OK != glewInit())
    {
        std::cout << "Failed to initialize GLEW" << std::endl;
        return EXIT_FAILURE;
    }

    char rendererInfo[512] = {0};
    strcpy(rendererInfo, (const char*) glGetString(GL_RENDERER));
    strcat(rendererInfo, " - ");
    strcat(rendererInfo, (const char*) glGetString(GL_VERSION));
    glfwSetWindowTitle(window, rendererInfo);

    init(fileName);

    glfwSetKeyCallback(window, keyboard);
    glfwSetMouseButtonCallback(window, mouse);
    // Disable mouse callback
    //glfwSetCursorPosCallback(window, mouse);



    mainLoop(window,gridHeight,gridWidth); // this does not return unless the window is closed

    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}

