#include <iostream>
#include "InputManager.h"
// #include "../DisplayGLFW/display.h"
#include "game.h"
#include "../res/includes/glm/glm.hpp"
#include "stb_image.h"
#include "cmath"

void calculateDerivative(const unsigned char* grayscale, unsigned char* result, float* gradDirection, int height, int width){

    unsigned char grayscaledx[width*height];
    unsigned char grayscaledy[width*height];
    //calculate derivative on both axis
    for (int i = 1; i < height; ++i) {
        for (int j = 1; j < width; ++j) {
            grayscaledx[i*width + j] = abs(grayscale[i*width + j] - grayscale[i*width + j-1]);
            grayscaledy[i*width + j] = abs(grayscale[i*width + j] - grayscale[(i-1)*width + j]);
            gradDirection[i * width + j] = atan2(grayscaledy[i * width + j], grayscaledx[i * width + j]);
        }
    }

    for(int i = 0; i < height; i++){
        for (int j = 0; j < width; ++j) {
            if(i ==0 || j==0)  result[i*width+j] = 0;
            result[i*width+j] = abs((int) std::round(sqrt(grayscaledx[i*width+j] * grayscaledx[i*width+j] + grayscaledy[i*width+j] * grayscaledy[i*width+j])))%256;
        }
    }

}

void MaximumSuppression(const unsigned char* grayscale,  const float* gradDirection, unsigned char* result, int height, int width){

    for (int i = 1;i < height - 1; ++i) {
        for (int j = 1; j < width - 1; ++j) {
            // Determine gradient directions
            float angle = gradDirection[i * width + j];

            // Determine neighboring pixels along the gradient direction
            float grad1, grad2;

            if (angle < 0)
                angle += M_PI;

            if ((angle >= 0 && angle < M_PI / 8) || (angle >= 15 * M_PI / 8 && angle < 2 * M_PI) || (angle >= 7*M_PI / 8 && angle < 9* M_PI / 8)) {
                grad1 = grayscale[i * width + j + 1];
                grad2 = grayscale[i * width + j - 1];
            } else if ((angle >= M_PI / 8 && angle < 3 * M_PI / 8) || (angle >= 9*M_PI / 8 && angle < 11* M_PI / 8)) {
                grad1 = grayscale[(i + 1) * width + j + 1];
                grad2 = grayscale[(i - 1) * width + j - 1];
            } else if ((angle >= 3 * M_PI / 8 && angle < 5 * M_PI / 8) || (angle >= 11*M_PI / 8 && angle < 13* M_PI / 8)) {
                grad1 = grayscale[(i + 1) * width + j];
                grad2 = grayscale[(i - 1) * width + j];
            } else if ((angle >= 5 * M_PI / 8 && angle < 7 * M_PI / 8) || (angle >= 13*M_PI / 8 && angle < 15* M_PI / 8)) {
                grad1 = grayscale[(i - 1) * width + j + 1];
                grad2 = grayscale[(i + 1) * width + j - 1];
            } else {
                grad1 = 0;
                grad2 = 0;
            }

            // Suppress non-maxima
            if (grayscale[i * width + j] >= grad1 && grayscale[i * width + j] >= grad2) {
                if(grayscale[i * width + j] < 28) result[i * width + j] = 0;
                else result[i * width + j] = 255;
                    //result[i * width + j] = grayscale[i * width + j];
            } else {
                result[i * width + j] = 0;
            }
        }
    }


    /*for(int i = 1; i < height-1; i++){
        for(int j = 1; j < width-1; j++){
            if(grayscale[i*width+j] > grayscale[(i-1)*width+j] && grayscale[i*width+j] > grayscale[(i+1)*width+j] && grayscale[i*width+j] > grayscale[i*width+j-1] &&  grayscale[i*width+j] > grayscale[i*width+j+1]){
                result[i*width+j] = grayscale[i*width+j];
            } else {
                result[i*width+j] = 0; //suppress
            }
        }
    }*/
}

void convolution(const unsigned char* grayscale, int height, int width, const float* kernel, unsigned char* result) {

    int count = 0;
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            int sum = 0;
            for (int k = 0; k < 3; ++k) { //row i
                for (int l = 0; l < 3; ++l) { //element ij
                    int imagei = (i + k-1);
                    int imagej = (j + l-1);
                    if(imagei<0) imagei = 1;
                    if(imagej<0) imagej = 1;
                    if(imagei>=height) imagei = height-1;
                    if(imagej>=width) imagej = width-1;
                    sum += grayscale[imagei * width + imagej] * kernel[k * 3 + l];
                }
            }
            count++;
            result[i * (width) + j] = sum;
        }
    }

}

void edgeDetect(const unsigned char* image, int height, int width, unsigned char* result) {

    unsigned char grayscale[height*width];
    for(int i = 0; i<height*width*4; i+=4){
        grayscale[i/4] = image[i];
    }

    unsigned char rConvolution[width*height];
    unsigned char rDerivative[width*height];
    float gradDirection[width*height];
    unsigned char rMaximumSuppression[height*width];

    float gaussian_kernel[]={0.125,0.125,0.125,0.125,0.25, 0.125,0.0625,0.125,0.0625};
    //float gaussian_kernel_d[9];
    //calculateDerivative(gaussian_kernel, gaussian_kernel_d, 3, 3);
    //for(float i : gaussian_kernel_d) printf("%f\n", i);
    convolution(grayscale, height, width, gaussian_kernel, rConvolution);

    calculateDerivative(rConvolution, rDerivative, gradDirection, height, width);
    MaximumSuppression(rDerivative, gradDirection, rMaximumSuppression, height, width);

    //todo: turn the image back into a RGBA
    for(int i = 0; i<height*width; i++){
        result[4*i] = rMaximumSuppression[i];
        result[4*i+1] = rMaximumSuppression[i];
        result[4*i+2] = rMaximumSuppression[i];
        result[4*i+3] = 255;
    }
}

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

    auto* edgeLena = new unsigned char[width*height*4];
    edgeDetect(data1, height, width, edgeLena);

    for(int i = 0 ; i < width*height*4; i++){
        printf("%d\n", edgeLena[i]);
    }

    //regular lena
    scn->AddTexture(width, height, data1);
    scn->SetShapeTex(0,0);
    scn->CustomDraw(1, 0, Game::BACK, true, false, 0);

    //edge detection
    scn->AddTexture(width, height, edgeLena);
    scn->SetShapeTex(0,1);
    scn->CustomDraw(1, 0, Game::BACK, false, false, 1);

    //halftone
    scn->AddTexture(width, height, data1);
    scn->SetShapeTex(0,2);
    scn->CustomDraw(1, 0, Game::BACK, false, false, 2);

    //Floyd-Steinberg Algorithm
    scn->AddTexture(width, height, data1);
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
    delete[] edgeLena;
    stbi_image_free(data1);
	return 0;
}

/*void convolution(const unsigned char* image, int height, int width, const float* kernel, unsigned char* result) {

    //todo: turn the image into a different gray-scale array
    unsigned char grayscale[height*width];
    for(int i = 0; i<height*width*4; i+=4){
        grayscale[i/4] = image[i];
    }

    unsigned char grayscaleResult[height*width];

    int count = 0;
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            int sum = 0;
            for (int k = 0; k < 3; ++k) { //row i
                for (int l = 0; l < 3; ++l) { //element ij
                    int imagei = (i + k-1);
                    int imagej = (j + l-1);
                    if(imagei<0) imagei = 1;
                    if(imagej<0) imagej = 1;
                    if(imagei>=height) imagei = height-1;
                    if(imagej>=width) imagej = width-1;
                    sum += grayscale[imagei * width + imagej] * kernel[k * 3 + l];
                }
            }
            count++;
            grayscaleResult[i * (width) + j] = sum;
        }
    }

    //calculateDerivative(grayscaleResult, grayscaleResultDer, height, width);
    //MaximumSuppression(grayscaleResultDer, grayscaleResultMax, height, width);

    //todo: turn the image back into a RGBA
    for(int i = 0; i<height*width; i++){
        result[4*i] = grayscaleResult[i];
        result[4*i+1] = grayscaleResult[i];
        result[4*i+2] = grayscaleResult[i];
        result[4*i+3] = 255;
    }

}*/
