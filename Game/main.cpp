
#define _USE_MATH_DEFINES // for C++
#include <cmath>
#include <iostream>
#include <fstream>

#include "InputManager.h"
// #include "../DisplayGLFW/display.h"
#include "game.h"
#include "../res/includes/glm/glm.hpp"
#include "stb_image.h"



void edgeDetect(const unsigned char* image, unsigned char* result);
void floydSteinberg(unsigned char* image, int w, int h, int org_colors, int new_colors);
unsigned char* classical_halftone(unsigned char* image, int w, int h);

int main(int argc,char *argv[])
{
    const int DISPLAY_WIDTH = 512;
    const int DISPLAY_HEIGHT = 512;
    const float CAMERA_ANGLE = 0.0f;
    const float NEAR = 1.0f;
    const float FAR = 100.0f;

    Game *scn = new Game(CAMERA_ANGLE,(float)DISPLAY_WIDTH/DISPLAY_HEIGHT,NEAR,FAR);

    Display display(DISPLAY_WIDTH, DISPLAY_HEIGHT, "OpenGL");

    Init(display);

    scn->Init();

    display.SetScene(scn);

    std::string lenaFile = "../res/textures/lena256.jpg";
    int width, height, numComponents;
    unsigned char* data1 = stbi_load((&lenaFile)->c_str(), &width, &height, &numComponents, 4);
    unsigned char* data2 = stbi_load((&lenaFile)->c_str(), &width, &height, &numComponents, 4);
    unsigned char* data3 = stbi_load((&lenaFile)->c_str(), &width, &height, &numComponents, 4);
    //implement edge detection
    //1 filter with derivative of gaussian
    //2 find magnitude and orientation of gradient
    //3 apply non-maximum suppression
    //linking and thresholding
    //  define 2 thresholds: high and low
    //  use high to start edge curves and low to continue them
    //      t1 - lower threshold value (fraction b/w 0-1)
    //      t2 - upper threshold value (fraction b/w 0-1)

    //regular lena
    scn->AddTexture(lenaFile, false);
    scn->SetShapeTex(0,0);
    scn->CustomDraw(1, 0, Game::BACK, true, false, 0);


    //edge detection
    unsigned char edgeLena[256 * 256 * 4];
    edgeDetect(data1, edgeLena);
    scn->AddTexture(width, height, edgeLena);
    scn->SetShapeTex(0,1);
    scn->CustomDraw(1, 0, Game::BACK, false, false, 1);

    //halftone
    //classical_halftone(data2, width, height);
    scn->AddTexture(width*2, height*2, classical_halftone(data2, width, height));
    scn->SetShapeTex(0,2);
    scn->CustomDraw(1, 0, Game::BACK, false, false, 2);

    //Floyd-Steinberg Algorithm
    floydSteinberg(data3, width, height, 256, 16);

    scn->AddTexture(256, 256, data3);
    scn->SetShapeTex(0,3);
    scn->CustomDraw(1, 0, Game::BACK, false, false, 3);

    scn->Motion();
    display.SwapBuffers();

    while(!display.CloseWindow())
    {
        //scn->Draw(1,0,scn->BACK,true,false);
        //scn->Motion();
        //display.SwapBuffers();
        display.PollEvents();

    }
    delete scn;
    return 0;
}

void setPixel(unsigned char* image, int w, int h, int i, int j, unsigned char val) {
    int index = 4 * (j + i * h);
    if (index >= 4 * w * h)return;
    image[index] = val;
    image[index+1] = val;
    image[index+2] = val;
    image[index+3] = 255;
}

unsigned char getPixel(unsigned char* image, int w, int h, int i, int j) {
    return image[4 * (j + i * h)];
}


void floydSteinberg(unsigned char* image, int w, int h, int org_colors, int new_colors) {
    float** fimage = (float**)malloc(sizeof(float*) * h);
    if (fimage == 0)return;
    for (int i = 0; i < h; i++) {
        fimage[i] = (float*)malloc(sizeof(float) * w);
        if (fimage[i] == 0)return;
    }
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            fimage[i][j] = getPixel(image,w,h,i,j);
        }
    }
    //printf("hell");

    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            float old_val = fimage[i][j];
            // unsigned char new_val = int(0.5 + (new_colors - 1) * old_val / (org_colors-1)) * floor((org_colors-1) / (new_colors-1))
            unsigned char new_val = int(old_val + 0.5)>=org_colors ? new_colors-1 : int(old_val+0.5)/(org_colors/new_colors);
            new_val *= org_colors / new_colors;
            fimage[i][j] = new_val;
            //printf("%f     %f    %d     %d\n", old_val, fimage[i][j], i, j);
            float error = old_val - new_val;

            if(i<h&&j+1<w)fimage[i][j+1] +=  error*7.0/16;
            if (i + 1 < h) {
                if(j>0)fimage[i + 1][j - 1] += error * 3.0 / 16;
                if(j<w)fimage[i + 1][j] += error * 5.0 / 16;
                if(j+1<w)fimage[i + 1][j + 1] += error * 1.0 / 16;
            }
        }
    }
    //printf("hell");
    std::ofstream outfile;
    outfile.open("../assignment/img6.txt", std::ios::out | std::ios::trunc);
    outfile.clear();
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            //printf("%f    %d\n", fimage[i][j], getPixel(image,w,h,i,j));
            setPixel(image, w, h, i, j, (unsigned char)fimage[i][j]);
            outfile << getPixel(image, w, h, i, j)/(org_colors/new_colors);
            if(i<h-1||j<w-1)outfile << ',';
        }
        outfile << 'n';
    }




}


unsigned char* classical_halftone(unsigned char* image, int w, int h) {
    unsigned char* himage = (unsigned char*)malloc(4 * sizeof(unsigned char) * w * h * 4);
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            float avg = getPixel(image, w, h, i, j)/255.0;
            //(i,j)->(2*i,2*j)
            i *= 2;
            j *= 2;
            w *= 2;
            h *= 2;
            if (avg >= 0.2) {
                setPixel(himage, w, h, i + 1, j, 255);
            }
            else {
                setPixel(himage, w, h, i + 1, j, 0);
            }
            if (avg >= 0.4) {
                setPixel(himage, w, h, i, j + 1, 255);
            }
            else {
                setPixel(himage, w, h, i, j + 1, 0);
            }
            if (avg >= 0.6) {
                setPixel(himage, w, h, i + 1, j + 1, 255);
            }
            else {
                setPixel(himage, w, h, i + 1, j + 1, 0);
            }
            if (avg >= 0.8) {
                setPixel(himage, w, h, i, j, 255);
            }
            else {
                setPixel(himage, w, h, i, j, 0);
            }
            i /= 2;
            j /= 2;
            w /= 2;
            h /= 2;
        }
    }

    std::ofstream outfile;
    outfile.open("../assignment/img5.txt", std::ios::out | std::ios::trunc);
    outfile.clear();
    h *= 2;
    w *= 2;
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            outfile << (255 == getPixel(himage, w, h, i, j)) ? '1' : '0';
            if (i < h - 1 || j < w - 1)outfile << ',';
        }
        outfile << 'n';
    }

    return himage;
}


void MaximumSuppression(const unsigned char* grayscale, const float* gradDirection, unsigned char* result) {
    const int width = 256;
    const int height = 256;
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {

            if (i == 0 || j == 0 || i == height - 1 || j == width - 1) {
                result[i * width + j] = 0;
                continue;
            }
            // Determine gradient directions
            float angle = gradDirection[i * width + j];

            // Determine neighboring pixels along the gradient direction
            float grad1, grad2;

            if (angle < 0)
                angle += M_PI;

            if ((angle >= 0 && angle < M_PI / 8) || (angle >= 15 * M_PI / 8 && angle < 2 * M_PI) || (angle >= 7 * M_PI / 8 && angle < 9 * M_PI / 8)) {
                grad1 = grayscale[i * width + j + 1];
                grad2 = grayscale[i * width + j - 1];
            }
            else if ((angle >= M_PI / 8 && angle < 3 * M_PI / 8) || (angle >= 9 * M_PI / 8 && angle < 11 * M_PI / 8)) {
                grad1 = grayscale[(i + 1) * width + j + 1];
                grad2 = grayscale[(i - 1) * width + j - 1];
            }
            else if ((angle >= 3 * M_PI / 8 && angle < 5 * M_PI / 8) || (angle >= 11 * M_PI / 8 && angle < 13 * M_PI / 8)) {
                grad1 = grayscale[(i + 1) * width + j];
                grad2 = grayscale[(i - 1) * width + j];
            }
            else if ((angle >= 5 * M_PI / 8 && angle < 7 * M_PI / 8) || (angle >= 13 * M_PI / 8 && angle < 15 * M_PI / 8)) {
                grad1 = grayscale[(i - 1) * width + j + 1];
                grad2 = grayscale[(i + 1) * width + j - 1];
            }
            else {
                grad1 = 0;
                grad2 = 0;
            }

            // Suppress non-maxima
            if (grayscale[i * width + j] >= grad1 && grayscale[i * width + j] >= grad2) {
                if (grayscale[i * width + j] < 28) result[i * width + j] = 0;
                else result[i * width + j] = 255;
                //result[i * width + j] = grayscale[i * width + j];
            }
            else {
                result[i * width + j] = 0;
            }
        }
    }
}

void convolution(const unsigned char* grayscale, const float* kernel, unsigned char* result) {
    const int width = 256;
    const int height = 256;
    int count = 0;
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            int sum = 0;
            for (int k = 0; k < 3; ++k) { //row i
                for (int l = 0; l < 3; ++l) { //element ij
                    int imagei = (i + k - 1);
                    int imagej = (j + l - 1);
                    if (imagei < 0) imagei = 1;
                    if (imagej < 0) imagej = 1;
                    if (imagei >= height) imagei = height - 1;
                    if (imagej >= width) imagej = width - 1;
                    sum += grayscale[imagei * width + imagej] * kernel[k * 3 + l];
                }
            }
            count++;
            result[i * (width)+j] = sum;
        }
    }

}


void calculateDerivative(const unsigned char* grayscale, unsigned char* result, float* gradDirection) {
    const int width = 256;
    const int height = 256;
    unsigned char grayscaledx[width * height];
    unsigned char grayscaledy[width * height];
    //calculate derivative on both axis
    for (int i = 1; i < height; ++i) {
        for (int j = 1; j < width; ++j) {
            grayscaledx[i * width + j] = abs(grayscale[i * width + j] - grayscale[i * width + j - 1]);
            grayscaledy[i * width + j] = abs(grayscale[i * width + j] - grayscale[(i - 1) * width + j]);
            gradDirection[i * width + j] = atan2(grayscaledy[i * width + j], grayscaledx[i * width + j]);
        }
    }

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; ++j) {
            if (i == 0 || j == 0)  result[i * width + j] = 0;
            result[i * width + j] = abs((int)std::round(sqrt(grayscaledx[i * width + j] * grayscaledx[i * width + j] + grayscaledy[i * width + j] * grayscaledy[i * width + j]))) % 256;
        }
    }

}

void edgeDetect(const unsigned char* image, unsigned char* result) {
    const int width = 256;
    const int height = 256;
    unsigned char grayscale[height * width];
    for (int i = 0; i < height * width * 4; i += 4) {
        grayscale[i / 4] = image[i];
    }

    unsigned char rConvolution[width * height];
    unsigned char rDerivative[width * height];
    float gradDirection[width * height];
    unsigned char rMaximumSuppression[height * width];

    float gaussian_kernel[] = { 0.125,0.125,0.125,0.125,0.25, 0.125,0.0625,0.125,0.0625 };

    convolution(grayscale, gaussian_kernel, rConvolution);

    calculateDerivative(rConvolution, rDerivative, gradDirection);
    MaximumSuppression(rDerivative, gradDirection, rMaximumSuppression);

    std::ofstream outfile;
    outfile.open("../assignment/img4.txt", std::ios::out | std::ios::trunc);
    outfile.clear();

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            outfile << (255 == rMaximumSuppression[i * width + j]) ? '1' : '0';
            if (i < height - 1 || j < width - 1)outfile << ',';
        }
        outfile << "\n";
    }
    outfile.close();

    //todo: turn the image back into a RGBA
    for (int i = 0; i < height * width; i++) {
        result[4 * i] = rMaximumSuppression[i];
        result[4 * i + 1] = rMaximumSuppression[i];
        result[4 * i + 2] = rMaximumSuppression[i];
        result[4 * i + 3] = 255;
    }
}