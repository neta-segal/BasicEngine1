#include <iostream>
#include "InputManager.h"
// #include "../DisplayGLFW/display.h"
#include "game.h"
#include "../res/includes/glm/glm.hpp"
#include "stb_image.h"

void convolution(const unsigned char* image, int height, int width, float* kernel, unsigned char* result) {
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
                    sum += image[imagei * width + imagej] * kernel[k * 3 + l];
                }
            }
            count++;
            result[i * (width) + j] = sum;
        }
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

    std::string lenaFile = "../res/textures/lena256.png";
    int width, height, numComponents;
    unsigned char* data1 = stbi_load((&lenaFile)->c_str(), &width, &height, &numComponents, 4);
    //printf("%d\n%d\n", height, width);
    //implement edge detection
    //1 filter with derivative of gaussian
    //float gaussian_kernel[]={0.01,0.08,0.01,0.08,0.64, 0.08,0.01,0.08,0.01}; //small theta
    unsigned char* data2 = new unsigned char[65536];
    //convolution(data1, height, width, gaussian_kernel, data2);
    //2 find magnitude and orientation of gradient
    //3 apply non-maximum suppression
    //linking and thresholding
    //  define 2 thresholds: high and low
    //  use high to start edge curves and low to continue them
    //      t1 - lower threshold value (fraction b/w 0-1)
    //      t2 - upper threshold value (fraction b/w 0-1)

    //regular lena
    scn->AddTexture(256, 256, data1);
    scn->SetShapeTex(0,0);
    scn->CustomDraw(1, 0, Game::BACK, true, false, 0);

    //edge detection
    scn->AddTexture(256, 256, data2);
    scn->SetShapeTex(0,0);
    scn->CustomDraw(1, 0, Game::BACK, false, false, 1);

    //halftone
    scn->AddTexture(256, 256, data1);
    scn->SetShapeTex(0,1);
    scn->CustomDraw(1, 0, Game::BACK, false, false, 2);

    //Floyd-Steinberg Algorithm
    scn->AddTexture(256, 256, data1);
    scn->SetShapeTex(0,2);
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
    stbi_image_free(data1);
    delete[] data2;
	return 0;
}
