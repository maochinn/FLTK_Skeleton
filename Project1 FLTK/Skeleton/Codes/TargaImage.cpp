///////////////////////////////////////////////////////////////////////////////
//
//      TargaImage.cpp                          Author:     Stephen Chenney
//                                              Modified:   Eric McDaniel
//                                              Date:       Fall 2004
//
//      Implementation of TargaImage methods.  You must implement the image
//  modification functions.
//
///////////////////////////////////////////////////////////////////////////////

#include "Globals.h"
//#include "TargaImage.h"
#include "Header.h"
#include "libtarga.h"
#include <stdlib.h>
#include <assert.h>
#include <memory.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <list>


using namespace std;

// constants
const int           RED             = 0;                // red channel
const int           GREEN           = 1;                // green channel
const int           BLUE            = 2;                // blue channel
const int			ALPHA			= 3;				// alpha channel
const unsigned char BACKGROUND[3]   = { 0, 0, 0 };      // background color


// Computes n choose s, efficiently
double Binomial(int n, int s)
{
    double        res;

    res = 1;
    for (int i = 1 ; i <= s ; i++)
        res = (n - i + 1) * res / i ;

    return res;
}// Binomial
glm::vec4* toFloat(unsigned char* u_rgba, int w, int h)
{
	glm::vec4* f_rgba;

	int size = w * h;
	f_rgba = new glm::vec4[size];

	for (int i = 0; i < size; i++)
	{
		int offset_rgba = i * 4;

		f_rgba[i].r = (float)u_rgba[offset_rgba + RED] / (float)255;
		f_rgba[i].g = (float)u_rgba[offset_rgba + GREEN] / (float)255;
		f_rgba[i].b = (float)u_rgba[offset_rgba + BLUE] / (float)255;
		f_rgba[i].a = (float)u_rgba[offset_rgba + ALPHA] / (float)255;
	}
	return f_rgba;
}
unsigned char* toUnsignedChar(glm::vec4* f_rgba, int w, int h)
{
	unsigned char* u_rgba;
	int size = w * h;
	u_rgba = new unsigned char[size * 4];

	for (int i = 0; i < size; i++)
	{
		int offset_rgba = i * 4;

		f_rgba[i] *= 255;
		u_rgba[offset_rgba + RED] = std::max(0, std::min((int)f_rgba[i].r, 255));
		u_rgba[offset_rgba + GREEN] = std::max(0, std::min((int)f_rgba[i].g, 255));
		u_rgba[offset_rgba + BLUE] = std::max(0, std::min((int)f_rgba[i].b, 255));
		u_rgba[offset_rgba + ALPHA] = std::max(0, std::min((int)f_rgba[i].a, 255));
	}
	return u_rgba;
}
double* difference(const glm::vec4* A, const glm::vec4* B, int size)
{
	double* diff = new double[size];
	for (int i = 0; i < size; i++)
	{
		diff[i] = sqrt
		(
			(A[i].r - B[i].r) * (A[i].r - B[i].r) +
			(A[i].g - B[i].g) * (A[i].g - B[i].g) +
			(A[i].b - B[i].b) * (A[i].b - B[i].b)
		);
	}
	return diff;
}

glm::vec3 convertTo332(glm::vec3 color)
{
	unsigned char R[8] = { 0, 36, 73, 109, 146, 182, 219, 255 };
	unsigned char G[8] = { 0, 36, 73, 109, 146, 182, 219, 255 };
	unsigned char B[4] = { 0, 85, 170, 255 };

	color.r = std::max(std::min(color.r, 1.0f), 0.0f);
	color.g = std::max(std::min(color.g, 1.0f), 0.0f);
	color.b = std::max(std::min(color.b, 1.0f), 0.0f);

	//convert to unsigned char
	unsigned char r = color.r * 255.0f;
	unsigned char g = color.g * 255.0f;
	unsigned char b = color.b * 255.0f;

	unsigned char min_error = INT_MAX;
	unsigned char min_r, min_g, min_b;
	for (int i = 0; i < 8; i++)
	{
		unsigned char error = abs(R[i] - r);
		if (error < min_error)
		{
			min_error = error;
			min_r = R[i];
		}
	}
	min_error = INT_MAX;
	for (int i = 0; i < 8; i++)
	{
		unsigned char error = abs(G[i] - g);
		if (error < min_error)
		{
			min_error = error;
			min_g = G[i];
		}
	}
	min_error = INT_MAX;
	for (int i = 0; i < 4; i++)
	{
		unsigned char error = abs(B[i] - b);
		if (error < min_error)
		{
			min_error = error;
			min_b = B[i];
		}
	}

	glm::vec3 new_color;
	new_color.r = min_r / 255.0f;
	new_color.g = min_g / 255.0f;
	new_color.b = min_b / 255.0f;
	return new_color;
}

///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage() : width(0), height(0), data(NULL)
{}// TargaImage

///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(int w, int h) : width(w), height(h)
{
   data = new unsigned char[width * height * 4];
   ClearToBlack();
}// TargaImage



///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables to values given.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(int w, int h, unsigned char *d)
{
    int i;

    width = w;
    height = h;
    data = new unsigned char[width * height * 4];

    for (i = 0; i < width * height * 4; i++)
	    data[i] = d[i];
}// TargaImage

///////////////////////////////////////////////////////////////////////////////
//
//      Copy Constructor.  Initialize member to that of input
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(const TargaImage& image) 
{
   width = image.width;
   height = image.height;
   data = NULL; 
   if (image.data != NULL) {
      data = new unsigned char[width * height * 4];
      memcpy(data, image.data, sizeof(unsigned char) * width * height * 4);
   }
}


///////////////////////////////////////////////////////////////////////////////
//
//      Destructor.  Free image memory.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::~TargaImage()
{
    if (data)
        delete[] data;
}// ~TargaImage


///////////////////////////////////////////////////////////////////////////////
//
//      Converts an image to RGB form, and returns the rgb pixel data - 24 
//  bits per pixel. The returned space should be deleted when no longer 
//  required.
//
///////////////////////////////////////////////////////////////////////////////
unsigned char* TargaImage::To_RGB(void)
{
    unsigned char   *rgb = new unsigned char[width * height * 3];
    int		    i, j;

    if (! data)
	    return NULL;

    // Divide out the alpha
    for (i = 0 ; i < height ; i++)
    {
	    int in_offset = i * width * 4;
	    int out_offset = i * width * 3;

	    for (j = 0 ; j < width ; j++)
        {
	        RGBA_To_RGB(data + (in_offset + j*4), rgb + (out_offset + j*3));
	    }
    }

    return rgb;
}// TargaImage


///////////////////////////////////////////////////////////////////////////////
//
//      Save the image to a targa file. Returns 1 on success, 0 on failure.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Save_Image(const char *filename)
{
    TargaImage	*out_image = Reverse_Rows();

    if (! out_image)
	    return false;

    if (!tga_write_raw(filename, width, height, out_image->data, TGA_TRUECOLOR_32))
    {
	    cout << "TGA Save EREDor: %s\n", tga_error_string(tga_get_last_error());
	    return false;
    }

    delete out_image;

    return true;
}// Save_Image


///////////////////////////////////////////////////////////////////////////////
//
//      Load a targa image from a file.  Return a new TargaImage object which 
//  must be deleted by caller.  Return NULL on failure.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage* TargaImage::Load_Image(char *filename)
{
    unsigned char   *temp_data;
    TargaImage	    *temp_image;
    TargaImage	    *result;
    int		        width, height;

    if (!filename)
    {
        cout << "No filename given." << endl;
        return NULL;
    }// if

    temp_data = (unsigned char*)tga_load(filename, &width, &height, TGA_TRUECOLOR_32);
    if (!temp_data)
    {
        cout << "TGA EREDor: %s\n", tga_error_string(tga_get_last_error());
	    width = height = 0;
	    return NULL;
    }
    temp_image = new TargaImage(width, height, temp_data);
    free(temp_data);

    result = temp_image->Reverse_Rows();

    delete temp_image;

    return result;
}// Load_Image


///////////////////////////////////////////////////////////////////////////////
//
//      Convert image to grayscale.  Red, green, and blue channels should all 
//  contain grayscale value.  Alpha channel shoould be left unchanged.  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::To_Grayscale()
{
	glm::vec4* img = toFloat(this->data, this->width, this->height);
	
	int size = this->width * this->height;
	for (int i = 0; i < size; i++)
	{
		float value = 0.299f * img[i].r + 0.587f * img[i].g + 0.114f * img[i].b;
		img[i].r = value;
		img[i].g = value;
		img[i].b = value;
	}

	this->data = toUnsignedChar(img, this->width, this->height);

    return true;
}// To_Grayscale


///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using uniform quantization.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Uniform()
{
	unsigned char* rgb = this->To_RGB();

	for (int y = 0; y < this->height; y++)
		for (int x = 0; x < this->width; x++)
		{
			int offset_rgba = x * 4 + (y * 4) * this->width;
			int offset_rgb = x * 3 + (y * 3) * this->width;

			unsigned char r = rgb[offset_rgb + RED];
			unsigned char g = rgb[offset_rgb + GREEN];
			unsigned char b = rgb[offset_rgb + BLUE];

			//convert to 3 or 2 bit;
			r = (r >> 5) << 5;
			g = (g >> 5) << 5;
			b = (b >> 6) << 6;

			data[offset_rgba + RED] = r;
			data[offset_rgba + GREEN] = g;
			data[offset_rgba + BLUE] = b;
			data[offset_rgba + ALPHA] = 255;
		}

	delete[] rgb;
    return true;
}// Quant_Uniform


///////////////////////////////////////////////////////////////////////////////
//
//      Convert the image to an 8 bit image using populosity quantization.  
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Populosity()
{
	unsigned char* rgb = this->To_RGB();

	//Populosity
	//count color that use RGB = (5, 5, 5)bit
	unsigned int color_size = 32 * 32 * 32;
	unsigned int* histogram = new unsigned int[color_size];
	for (int i = 0; i < color_size; i++)
		histogram[i] = 0;

	for (int y = 0; y < this->height; y++)
		for (int x = 0; x < this->width; x++)
		{
			int offset_rgba = x * 4 + (y * 4) * this->width;
			int offset_rgb = x * 3 + (y * 3) * this->width;

			unsigned char r = rgb[offset_rgb + RED];
			unsigned char g = rgb[offset_rgb + GREEN];
			unsigned char b = rgb[offset_rgb + BLUE];

			//convert to 5 bit(32 shade)
			r = (r >> 3);
			g = (g >> 3);
			b = (b >> 3);

			unsigned int idx = 0;
			idx = (idx | r) << 5;
			idx = (idx | g) << 5;
			idx = (idx | b);

			histogram[idx]++;
		}

	//find 256 most popular colors
	unsigned int popular = 256;
	unsigned char* pupular_rgbs = new unsigned char[popular * 3];
	for (int i = 0; i < popular; i++)
	{
		unsigned int max_count = 0;
		unsigned int max_idx = 0;
		//search i-th popular color
		for (int j = 0; j < color_size; j++)
		{
			if (histogram[j] > max_count)
			{
				max_count = histogram[j];
				max_idx = j;
			}
		}
		//record color
		{
			//set to zero to label already choose
			histogram[max_idx] = 0;

			//decode rgb color
			unsigned char b = ((unsigned char)max_idx << 3);
			max_idx = max_idx >> 5;
			unsigned char g = ((unsigned char)max_idx << 3);
			max_idx = max_idx >> 5;
			unsigned char r = ((unsigned char)max_idx << 3);

			pupular_rgbs[i * 3 + RED] = r;
			pupular_rgbs[i * 3 + GREEN] = g;
			pupular_rgbs[i * 3 + BLUE] = b;
		}
	}
	//set color by closest popular colors
	for (int y = 0; y < this->height; y++)
		for (int x = 0; x < this->width; x++)
		{
			int offset_rgb = x * 3 + (y * 3) * this->width;
			int offset_rgba = x * 4 + (y * 4) * this->width;

			float min_distance = FLT_MAX;
			unsigned char min_rgb[3];
			unsigned char now_rgb[3] =
			{
				rgb[offset_rgb + RED],
				rgb[offset_rgb + GREEN],
				rgb[offset_rgb + BLUE],
			};
			for (unsigned int j = 0; j < popular; j++)
			{
				unsigned char pupular_rgb[3] =
				{
					pupular_rgbs[j * 3 + RED],
					pupular_rgbs[j * 3 + GREEN],
					pupular_rgbs[j * 3 + BLUE],
				};

				float distance =
					((float)now_rgb[0] - pupular_rgb[0]) * ((float)now_rgb[0] - pupular_rgb[0]) +
					((float)now_rgb[1] - pupular_rgb[1]) * ((float)now_rgb[1] - pupular_rgb[1]) +
					((float)now_rgb[2] - pupular_rgb[2]) * ((float)now_rgb[2] - pupular_rgb[2]);
				if (distance < min_distance)
				{
					min_distance = distance;
					min_rgb[RED] = pupular_rgb[RED];
					min_rgb[GREEN] = pupular_rgb[GREEN];
					min_rgb[BLUE] = pupular_rgb[BLUE];
				}
			}
			//set color
			{
				this->data[offset_rgba + RED] = min_rgb[RED];
				this->data[offset_rgba + GREEN] = min_rgb[GREEN];
				this->data[offset_rgba + BLUE] = min_rgb[BLUE];
				this->data[offset_rgba + ALPHA] = 255;
			}
		}
	delete[] histogram;
	delete[] pupular_rgbs;

	delete[] rgb;
    return false;
}// Quant_Populosity


///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image using a threshold of 1/2.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Threshold()
{
	this->To_Grayscale();

	unsigned char* rgb = this->To_RGB();

	unsigned int size = this->width * this->height;
	for (unsigned int i = 0; i < size; i++)
	{
		int offset_rgb = i * 3;
		int offset_rgba = i * 4;
		if (rgb[offset_rgb + RED] < 128)
			this->data[offset_rgba + RED] = this->data[offset_rgba + GREEN] = this->data[offset_rgba + BLUE] = 0;
		else
			this->data[offset_rgba + RED] = this->data[offset_rgba + GREEN] = this->data[offset_rgba + BLUE] = 255;
		this->data[offset_rgba + ALPHA] = 255;
	}

	delete[] rgb;
    return true;
}// Dither_Threshold


///////////////////////////////////////////////////////////////////////////////
//
//      Dither image using random dithering.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Random()
{
	this->To_Grayscale();

	unsigned char* rgb = this->To_RGB();

	unsigned int size = this->width * this->height;
	for (unsigned int i = 0; i < size; i++)
	{
		int offset_rgb = i * 3;
		int offset_rgba = i * 4;

		//generate random number [-0.2f, 0.2f]
		float random = (0.2f - -0.2f) * rand() / (RAND_MAX + 1.0) - 0.2f;
		char ran = random * 255;

		if (rgb[offset_rgb + RED] + ran < 128)
			this->data[offset_rgba + RED] = this->data[offset_rgba + GREEN] = this->data[offset_rgba + BLUE] = 0;
		else
			this->data[offset_rgba + RED] = this->data[offset_rgba + GREEN] = this->data[offset_rgba + BLUE] = 255;
		this->data[offset_rgba + ALPHA] = 255;
	}

	delete[] rgb;
    return true;
}// Dither_Random


///////////////////////////////////////////////////////////////////////////////
//
//      Perform Floyd-Steinberg dithering on the image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_FS()
{
	this->To_Grayscale();

	unsigned char* rgb = this->To_RGB();

	//convert to float
	int size = this->width * this->height;
	glm::vec3* img = new glm::vec3[size];

	for (int i = 0; i < size; i++)
	{
		int offset_rgb = i * 3;

		img[i].r = (float)rgb[offset_rgb + RED] / 255.0f;
		img[i].g = (float)rgb[offset_rgb + GREEN] / 255.0f;
		img[i].b = (float)rgb[offset_rgb + BLUE] / 255.0f;
	}



	float threshold = 0.5;

	for (int y = 0; y < this->height; y++)
	{
		//zigzag
		if (y % 2)
		{
			for (int x = 0; x < this->width; x++)
			{
				float old_value = img[x + y * this->width].r;
				float new_value;
				if (old_value < threshold)
					new_value = 0.0f;
				else
					new_value = 1.0f;
				img[x + y * this->width] = glm::vec4(glm::vec3(new_value), 1.0f);
				float error = old_value - new_value;

				if (x + 1 < this->width)
					img[(x + 1) + (y)*this->width] += glm::vec3(error * 0.4375f);
				if (x - 1 >= 0 && y + 1 < this->height)
					img[(x - 1) + (y + 1) * this->width] += glm::vec3(error * 0.1875f);
				if (y + 1 < this->height)
					img[(x)+(y + 1) * this->width] += glm::vec3(error * 0.3125f);
				if (x + 1 < this->width && y + 1 < this->height)
					img[(x + 1) + (y + 1) * this->width] += glm::vec3(error * 0.0625f);
			}
		}
		else
		{
			for (int x = this->width - 1; x >= 0; x--)
			{
				float old_value = img[x + y * this->width].r;
				float new_value;
				if (old_value < threshold)
					new_value = 0.0f;
				else
					new_value = 1.0f;
				img[x + y * this->width] = glm::vec4(glm::vec3(new_value), 1.0f);
				float error = old_value - new_value;

				if (x - 1 >= 0)
					img[(x - 1) + (y)*this->width] += glm::vec3(error * 0.4375f);
				if (x + 1 < this->width && y + 1 < this->height)
					img[(x + 1) + (y + 1) * this->width] += glm::vec3(error * 0.1875f);
				if (y + 1 < this->height)
					img[(x)+(y + 1) * this->width] += glm::vec3(error * 0.3125f);
				if (x - 1 >= 0 && y + 1 < this->height)
					img[(x - 1) + (y + 1) * this->width] += glm::vec3(error * 0.0625f);


			}
		}
	}

	for (unsigned int i = 0; i < size; i++)
	{
		int offset_rgb = i * 3;
		int offset_rgba = i * 4;

		glm::clamp(img[i], glm::vec3(0.0f), glm::vec3(1.0f));

		this->data[offset_rgba + RED] = img[i].r * 255;
		this->data[offset_rgba + GREEN] = img[i].g * 255;
		this->data[offset_rgba + BLUE] = img[i].b * 255;
		this->data[offset_rgba + ALPHA] = 255;
	}

	delete[] rgb;
    return true;
}// Dither_FS


///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image while conserving the average brightness.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Bright()
{
	this->To_Grayscale();

	unsigned char* rgb = this->To_RGB();

	std::vector<unsigned char> values;

	float average = 0.0f;
	unsigned int size = this->width * this->height;
	for (unsigned int i = 0; i < size; i++)
	{
		values.push_back(rgb[i * 3 + RED]);
		average += values[i];
	}
	average /= (float)size;
	std::sort(values.begin(), values.end());
	unsigned int idx = ((255 - average) / 255) * size;
	unsigned char threshold = values[idx];

	for (unsigned int i = 0; i < size; i++)
	{
		if (rgb[i * 3 + RED] < threshold)
			this->data[i * 4 + RED] = this->data[i * 4 + GREEN] = this->data[i * 4 + BLUE] = 0;
		else
			this->data[i * 4 + RED] = this->data[i * 4 + GREEN] = this->data[i * 4 + BLUE] = 255;

		this->data[i * 4 + ALPHA] = 255;
	}

	values.clear();

	delete[] rgb;
    return true;
}// Dither_Bright


///////////////////////////////////////////////////////////////////////////////
//
//      Perform clustered differing of the image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Cluster()
{
	this->To_Grayscale();

	unsigned char* rgb = this->To_RGB();

	//convert to float
	int size = this->width * this->height;
	glm::vec3* img = new glm::vec3[size];

	for (int i = 0; i < size; i++)
	{
		int offset_rgb = i * 3;

		img[i].r = (float)rgb[offset_rgb + RED] / 255.0f;
		img[i].g = (float)rgb[offset_rgb + GREEN] / 255.0f;
		img[i].b = (float)rgb[offset_rgb + BLUE] / 255.0f;
	}

	glm::mat4 threshold;
	threshold[0] = glm::vec4(0.7059f, 0.3529f, 0.5882f, 0.2353f);
	threshold[1] = glm::vec4(0.0588f, 0.9412f, 0.8235f, 0.4118f);
	threshold[2] = glm::vec4(0.4706f, 0.7647f, 0.8824f, 0.1176f);
	threshold[3] = glm::vec4(0.1765f, 0.5294f, 0.2941f, 0.6471f);

	for (unsigned int j = 0; j < this->height; j++)
		for (unsigned int i = 0; i < this->width; i++)
		{
			if (img[i + j * this->width].r >= threshold[i % 4][j % 4])
				img[i + j * this->width] = glm::vec4(1.0f);
			else
				img[i] = glm::vec4(0.0f);
		}

	//set back
	for (unsigned int i = 0; i < size; i++)
	{
		int offset_rgb = i * 3;
		int offset_rgba = i * 4;

		this->data[offset_rgba + RED] = img[i].r * 255;
		this->data[offset_rgba + GREEN] = img[i].g * 255;
		this->data[offset_rgba + BLUE] = img[i].b * 255;
		this->data[offset_rgba + ALPHA] = 255;
	}



	delete[] rgb;
    return true;
}// Dither_Cluster


///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using Floyd-Steinberg dithering over
//  a uniform quantization - the same quantization as in Quant_Uniform.
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Color()
{
	unsigned char* rgb = this->To_RGB();

	//convert to float
	int size = this->width * this->height;
	glm::vec3* img = new glm::vec3[size];

	for (int i = 0; i < size; i++)
	{
		int offset_rgb = i * 3;

		img[i].r = (float)rgb[offset_rgb + RED] / 255.0f;
		img[i].g = (float)rgb[offset_rgb + GREEN] / 255.0f;
		img[i].b = (float)rgb[offset_rgb + BLUE] / 255.0f;
	}


	for (int y = 0; y < this->height; y++)
	{
		//zigzag
		if (y % 2)
		{
			for (int x = 0; x < this->width; x++)
			{
				glm::vec3 old_rgb = img[x + y * this->width];
				glm::vec3 new_rgb = convertTo332(old_rgb);

				img[x + y * this->width] = new_rgb;
				glm::vec3 error = old_rgb - new_rgb;

				if (x + 1 < this->width)
					img[(x + 1) + (y)*this->width] += error * 0.4375f;
				if (x - 1 >= 0 && y + 1 < this->height)
					img[(x - 1) + (y + 1) * this->width] += error * 0.1875f;
				if (y + 1 < this->height)
					img[(x)+(y + 1) * this->width] += error * 0.3125f;
				if (x + 1 < this->width && y + 1 < this->height)
					img[(x + 1) + (y + 1) * this->width] += error * 0.0625f;

			}
		}
		else
		{
			for (int x = this->width - 1; x >= 0; x--)
			{
				glm::vec3 old_rgb = img[x + y * this->width];
				glm::vec3 new_rgb = convertTo332(old_rgb);

				img[x + y * this->width] = new_rgb;
				glm::vec3 error = old_rgb - new_rgb;

				if (x - 1 >= 0)
					img[(x - 1) + (y)*this->width] += error * 0.4375f;
				if (x + 1 < this->width && y + 1 < this->height)
					img[(x + 1) + (y + 1) * this->width] += error * 0.1875f;
				if (y + 1 < this->height)
					img[(x)+(y + 1) * this->width] += error * 0.3125f;
				if (x - 1 >= 0 && y + 1 < this->height)
					img[(x - 1) + (y + 1) * this->width] += error * 0.0625f;

			}
		}
	}

	//set back
	for (unsigned int i = 0; i < size; i++)
	{
		int offset_rgb = i * 3;
		int offset_rgba = i * 4;

		this->data[offset_rgba + RED] = img[i].r * 255;
		this->data[offset_rgba + GREEN] = img[i].g * 255;
		this->data[offset_rgba + BLUE] = img[i].b * 255;
		this->data[offset_rgba + ALPHA] = 255;
	}

	delete[] rgb;
    return true;
}// Dither_Color


///////////////////////////////////////////////////////////////////////////////
//
//      Composite the cuREDent image over the given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Over(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout <<  "Comp_Over: Images not the same size\n";
        return false;
    }

	this->compImage(pImage, _over);

    return true;
}// Comp_Over


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "in" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_In(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_In: Images not the same size\n";
        return false;
    }

	this->compImage(pImage, _in);

    return true;
}// Comp_In


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "out" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Out(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Out: Images not the same size\n";
        return false;
    }

	this->compImage(pImage, _out);

    return true;
}// Comp_Out


///////////////////////////////////////////////////////////////////////////////
//
//      Composite cuREDent image "atop" given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Atop(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Atop: Images not the same size\n";
        return false;
    }

	this->compImage(pImage, _atop);

    return true;
}// Comp_Atop


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image with given image using exclusive or (XOR).  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Xor(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Xor: Images not the same size\n";
        return false;
    }

	this->compImage(pImage, _xor);

    return true;
}// Comp_Xor


///////////////////////////////////////////////////////////////////////////////
//
//      Calculate the difference bewteen this imag and the given one.  Image 
//  dimensions must be equal.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Difference(TargaImage* pImage)
{
    if (!pImage)
        return false;

    if (width != pImage->width || height != pImage->height)
    {
        cout << "Difference: Images not the same size\n";
        return false;
    }// if

    for (int i = 0 ; i < width * height * 4 ; i += 4)
    {
        unsigned char        rgb1[3];
        unsigned char        rgb2[3];

        RGBA_To_RGB(data + i, rgb1);
        RGBA_To_RGB(pImage->data + i, rgb2);

        data[i] = abs(rgb1[0] - rgb2[0]);
        data[i+1] = abs(rgb1[1] - rgb2[1]);
        data[i+2] = abs(rgb1[2] - rgb2[2]);
        data[i+3] = 255;
    }

    return true;
}// Difference


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 box filter on this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Box()
{
	int size = 5;

	double** mask = new double* [size];
	for (int i = 0; i < size; i++)
		mask[i] = new double[size];

	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
			mask[i][j] = 1.0f / (size * size);

	this->filtering(mask, size);

	for (int i = 0; i < size; i++)
		delete[] mask[i];
	delete[] mask;

    return true;
}// Filter_Box


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Bartlett filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Bartlett()
{
	double bartlett[5][5] =
	{
		1.0f, 2.0f, 3.0f, 2.0f, 1.0f,
		2.0f, 4.0f, 6.0f, 4.0f, 2.0f,
		3.0f, 6.0f, 9.0f, 6.0f, 3.0f,
		2.0f, 4.0f, 6.0f, 4.0f, 2.0f,
		1.0f, 2.0f, 3.0f, 2.0f, 1.0f,
	};

	int size = 5;

	double** mask = new double* [size];
	for (int i = 0; i < size; i++)
		mask[i] = new double[size];

	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
			mask[i][j] = bartlett[i][j] / 81.0;

	this->filtering(mask, size);

	for (int i = 0; i < size; i++)
		delete[] mask[i];
	delete[] mask;

    return true;
}// Filter_Bartlett


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Gaussian()
{
	double gaussian[5][5] =
	{
		1.0f, 4.0f, 6.0f, 4.0f, 1.0f,
		4.0f, 16.0f, 24.0f, 16.0f, 4.0f,
		6.0f, 24.0f, 36.0f, 24.0f, 6.0f,
		4.0f, 16.0f, 24.0f, 16.0f, 4.0f,
		1.0f, 4.0f, 6.0f, 4.0f, 1.0f,
	};
	int size = 5;

	double** mask = new double* [size];
	for (int i = 0; i < size; i++)
		mask[i] = new double[size];

	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
			mask[i][j] = gaussian[i][j] / 256.0;

	this->filtering(mask, size);

	for (int i = 0; i < size; i++)
		delete[] mask[i];
	delete[] mask;

    return true;
}// Filter_Gaussian

///////////////////////////////////////////////////////////////////////////////
//
//      Perform NxN Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////

bool TargaImage::Filter_Gaussian_N(unsigned int N)
{
	if (N % 2 == 0)
		N--;

	double** mask = new double* [N];
	for (int i = 0; i < N; i++)
		mask[i] = new double[N];

	for (int i = 0; i <= N / 2; i++)
	{
		int n = (N - 1) + 2 * i;
		for (int j = 0; j < N; j++)
		{
			int s = i + j;
			mask[i][j] = (float)Binomial(n, s);
		}
	}
	for (int i = N / 2 + 1; i < N; i++)
		for (int j = 0; j < N; j++)
			mask[i][j] = mask[N - i - 1][j];

	float sum = 0.0f;
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			sum += mask[i][j];

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			mask[i][j] /= sum;

	this->filtering(mask, N);
	return false;
}// Filter_Gaussian_N


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 edge detect (high pass) filter on this image.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Edge()
{
	unsigned char* old_rgb = this->To_RGB();
	this->Filter_Gaussian();
	unsigned char* new_rgb = this->To_RGB();

	int size = this->width * this->height;
	for (int i = 0; i < size; i++)
	{
		int offset_rgb = i * 3;
		int offset_rgba = i * 4;

		glm::ivec3 diff_rgb
		{
			old_rgb[offset_rgb + RED] - new_rgb[offset_rgb + RED],
			old_rgb[offset_rgb + GREEN] - new_rgb[offset_rgb + GREEN],
			old_rgb[offset_rgb + BLUE] - new_rgb[offset_rgb + BLUE]
		};

		diff_rgb.r = std::max(diff_rgb.r, 0);
		diff_rgb.g = std::max(diff_rgb.g, 0);
		diff_rgb.b = std::max(diff_rgb.b, 0);

		this->data[offset_rgba + RED] = diff_rgb.r;
		this->data[offset_rgba + GREEN] = diff_rgb.g;
		this->data[offset_rgba + BLUE] = diff_rgb.b;
		this->data[offset_rgba + ALPHA] = 255;
	}

	delete[] old_rgb;
	delete[] new_rgb;
    return true;
}// Filter_Edge


///////////////////////////////////////////////////////////////////////////////
//
//      Perform a 5x5 enhancement filter to this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Enhance()
{
	unsigned char* old_rgb = this->To_RGB();
	this->Filter_Edge();
	unsigned char* new_rgb = this->To_RGB();

	int size = this->width * this->height;
	for (int i = 0; i < size; i++)
	{
		int offset_rgb = i * 3;
		int offset_rgba = i * 4;

		glm::ivec3 add_rgb
		{
			old_rgb[offset_rgb + RED] + new_rgb[offset_rgb + RED],
			old_rgb[offset_rgb + GREEN] + new_rgb[offset_rgb + GREEN],
			old_rgb[offset_rgb + BLUE] + new_rgb[offset_rgb + BLUE]
		};

		add_rgb.r = std::min(add_rgb.r, 255);
		add_rgb.g = std::min(add_rgb.g, 255);
		add_rgb.b = std::min(add_rgb.b, 255);

		this->data[offset_rgba + RED] = add_rgb.r;
		this->data[offset_rgba + GREEN] = add_rgb.g;
		this->data[offset_rgba + BLUE] = add_rgb.b;
		this->data[offset_rgba + ALPHA] = 255;
	}

	delete[] old_rgb;
	delete[] new_rgb;
    return true;
}// Filter_Enhance


///////////////////////////////////////////////////////////////////////////////
//
//      Run simplified version of Hertzmann's painterly image filter.
//      You probably will want to use the Draw_Stroke funciton and the
//      Stroke class to help.
// Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::NPR_Paint()
{
	//canvas is rgba image
	int size = this->width * this->height;
	glm::vec4* canvas = new glm::vec4[size];
	//canvas is constant image
	for (int i = 0; i < size; i++)
	{
		//起始畫布顏色要與任何顏色都差很多(>threshold)
		canvas[i] = glm::vec4(10.0f, 10.0f, 10.0f, 1.0f);
	}
	unsigned char* source_img = new unsigned char[size * 4];
	for (int i = 0; i < size; i++)
	{
		int offset = i * 4;
		source_img[offset + RED]	= this->data[offset + RED];
		source_img[offset + GREEN]	= this->data[offset + GREEN];
		source_img[offset + BLUE]	= this->data[offset + BLUE];
		source_img[offset + ALPHA]	= this->data[offset + ALPHA];
	}

	//paint the canvas
	//brush radius
	float R[3] = { 7.0f, 3.0f, 1.0f };
	for (int i = 0; i < 3; i++)
	{
		//apply Gaussian blur
		{
			//put source_img to img_data
			for (int j = 0; j < size; j++)
			{
				int offset = j * 4;
				this->data[offset + RED]	= source_img[offset + RED];
				this->data[offset + GREEN]	= source_img[offset + GREEN];
				this->data[offset + BLUE]	= source_img[offset + BLUE];
				this->data[offset + ALPHA]	= source_img[offset + ALPHA];
			}
			int coeffient = 2 * R[i] + 1;
			//int coeffient = R[i];
			if (coeffient > 1)
				this->Filter_Gaussian_N(coeffient);
		}
		glm::vec4* reference_img = toFloat(this->data, this->width, this->height);
		//paint a layer

		NPR_Paint_Layer(canvas, reference_img, R[i]);

		delete[] reference_img;
	}

	//set canvas to img_data

	for (int i = 0; i < this->width * this->height; i++)
	{
		int offset = i * 4;

		canvas[i].r = std::max(std::min(canvas[i].r, 1.0f), 0.0f);
		canvas[i].g = std::max(std::min(canvas[i].g, 1.0f), 0.0f);
		canvas[i].b = std::max(std::min(canvas[i].b, 1.0f), 0.0f);
		canvas[i].a = std::max(std::min(canvas[i].a, 1.0f), 0.0f);

		this->data[offset + RED] = canvas[i].r * 255;
		this->data[offset + GREEN] = canvas[i].g * 255;
		this->data[offset + BLUE] = canvas[i].b * 255;
		this->data[offset + ALPHA] = canvas[i].a * 255;
	}


	delete[] source_img;
	delete[] canvas;
    return false;
}
void TargaImage::NPR_Paint_Layer(glm::vec4* tCanvas, const glm::vec4* tReferenceImage, int tBrushSize)
{
	float f_g = 0.8f;
	float T = 25.0f / 255.0f;

	std::list<Stroke> S;

	//create a pointwise difference image
	double* D = difference(tCanvas, tReferenceImage, this->width * this->height);

	int grid = std::max(1, (int)(f_g * tBrushSize));
	int half_grid = grid / 2;

	for (int x = 0; x < this->width; x += grid)
		for (int y = 0; y < this->height; y += grid)
		{

			//sum the error near(x, y)
			int M_size = (2 * half_grid + 1) * (2 * half_grid + 1);
			glm::uvec2* M = new glm::uvec2[M_size];
			int idx = 0;
			for (int i = x - half_grid; i <= x + half_grid; i++)
				for (int j = y - half_grid; j <= y + half_grid; j++)
				{
					int u = std::max(std::min(i, this->width - 1), 0);
					int v = std::max(std::min(j, this->height - 1), 0);

					M[idx].x = u;
					M[idx].y = v;
					idx++;
				}

			float areaError = 0.0f;
			for (int i = 0; i < M_size; i++)
			{
				areaError += D[M[i].x + M[i].y * this->width];
			}
			areaError /= (float)(grid * grid);

			if (areaError > T)
			{

				//find the largest error point
				float max_error = FLT_MIN;
				glm::uvec2 max_point;
				for (int i = 0; i < M_size; i++)
				{
					int offset = M[i].x + M[i].y * this->width;
					if (D[offset] > max_error)
					{
						max_error = D[offset];
						max_point = M[i];
					}
				}
				glm::vec4 refer_color =
					tReferenceImage[max_point.x + max_point.y * this->width];
				unsigned char r = refer_color.r * 255;
				unsigned char g = refer_color.g * 255;
				unsigned char b = refer_color.b * 255;
				unsigned char a = 255;


				Stroke s(tBrushSize, max_point.x, max_point.y, r, g, b, a);
				S.push_back(s);
			}

			delete[] M;
		}

	//paint all strokes in S on the canvas, in random order
	//put canvas to img_data
	for (int i = 0; i < this->width * this->height; i++)
	{
		int offset = i * 4;
		this->data[offset + RED] = tCanvas[i].r * 255;
		this->data[offset + GREEN] = tCanvas[i].g * 255;
		this->data[offset + BLUE] = tCanvas[i].b * 255;
		this->data[offset + ALPHA] = tCanvas[i].a * 255;
	}

	while (S.empty() != true)
	{
		int random = rand() % (S.size());
		auto it = S.begin();
		std::advance(it, random);
		Paint_Stroke(*it);
		S.erase(it);
	}

	//put img_data to canvas
	for (int i = 0; i < this->width * this->height; i++)
	{
		int offset = i * 4;
		tCanvas[i].r = this->data[offset + RED] / 255.0f;
		tCanvas[i].g = this->data[offset + GREEN] / 255.0f;
		tCanvas[i].b = this->data[offset + BLUE] / 255.0f;
		tCanvas[i].a = this->data[offset + ALPHA] / 255.0f;
	}

	delete[] D;
}


///////////////////////////////////////////////////////////////////////////////
//
//      Halve the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Half_Size()
{
	int new_w = this->width / 2;
	int new_h = this->height / 2;
	int new_size = new_w * new_h;

	glm::vec4* new_img = new glm::vec4[new_size];

	for (int y = 0; y < new_h; y++)
		for (int x = 0; x < new_w; x++)
		{
			int u = 2 * x;
			int v = 2 * y;
			new_img[x + y * new_w] = resample3X3(u, v);
		}

	delete[] this->data;
	this->data = new unsigned char[new_size * 4];

	//set back
	for (int i = 0; i < new_size; i++)
	{
		int offset_rgba = i * 4;
		data[offset_rgba + RED] = new_img[i].r * 255;
		data[offset_rgba + GREEN] = new_img[i].g * 255;
		data[offset_rgba + BLUE] = new_img[i].b * 255;
		data[offset_rgba + ALPHA] = new_img[i].a * 255;
	}
	this->width = new_w;
	this->height = new_h;

	delete[] new_img;
    return true;
}// Half_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Double the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Double_Size()
{
	int new_w = this->width * 2;
	int new_h = this->height * 2;
	int new_size = new_w * new_h;

	glm::vec4* new_img = new glm::vec4[new_size];


	for (int y = 0; y < new_h; y++)
		for (int x = 0; x < new_w; x++)
		{
			int u = x / 2;
			int v = y / 2;

			if (x % 2 == 0 && y % 2 == 0)
				new_img[x + y * new_w] = resample3X3(u, v);
			else if (x % 2 == 1 && y % 2 == 1)
				new_img[x + y * new_w] = resample4X4(u, v);
			else if (x % 2 == 0 && y % 2 == 1)
				new_img[x + y * new_w] = resample3X4(u, v);
			else if (x % 2 == 1 && y % 2 == 0)
				new_img[x + y * new_w] = resample4X3(u, v);
		}

	delete[] this->data;
	this->data = new unsigned char[new_size * 4];

	//set back
	for (int i = 0; i < new_size; i++)
	{
		int offset_rgba = i * 4;
		this->data[offset_rgba + RED] = new_img[i].r * 255;
		this->data[offset_rgba + GREEN] = new_img[i].g * 255;
		this->data[offset_rgba + BLUE] = new_img[i].b * 255;
		this->data[offset_rgba + ALPHA] = new_img[i].a * 255;
	}
	this->width = new_w;
	this->height = new_h;

	delete[] new_img;
    return true;
}// Double_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Scale the image dimensions by the given factor.  The given factor is 
//  assumed to be greater than one.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Resize(float scale)
{
	int new_w = this->width * scale;
	int new_h = this->height * scale;
	int new_size = new_w * new_h;

	glm::vec4* new_img = new glm::vec4[new_size];

	for (int y = 0; y < new_h; y++)
		for (int x = 0; x < new_w; x++)
		{
			int u = (int)(x / scale);
			int v = (int)(y / scale);

			if (x % 2 == 0 && y % 2 == 0)
				new_img[x + y * new_w] = resample3X3(u, v);
			else if (x % 2 == 1 && y % 2 == 1)
				new_img[x + y * new_w] = resample4X4(u, v);
			else if (x % 2 == 0 && y % 2 == 1)
				new_img[x + y * new_w] = resample3X4(u, v);
			else if (x % 2 == 1 && y % 2 == 0)
				new_img[x + y * new_w] = resample4X3(u, v);
		}
	delete[] this->data;
	this->data = new unsigned char[new_size * 4];

	//set back
	for (int i = 0; i < new_size; i++)
	{
		int offset_rgba = i * 4;
		data[offset_rgba + RED] = new_img[i].r * 255;
		data[offset_rgba + GREEN] = new_img[i].g * 255;
		data[offset_rgba + BLUE] = new_img[i].b * 255;
		data[offset_rgba + ALPHA] = new_img[i].a * 255;
	}
	this->width = new_w;
	this->height = new_h;

	delete[] new_img;
    return true;
}// Resize


//////////////////////////////////////////////////////////////////////////////
//
//      Rotate the image clockwise by the given angle.  Do not resize the 
//  image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Rotate(float angleDegrees)
{
	int new_w = this->width;
	int new_h = this->height;
	int new_size = new_w * new_h;

	float radians = glm::radians(angleDegrees);

	glm::vec4* new_img = new glm::vec4[new_size];


	for (int y = 0; y < new_h; y++)
		for (int x = 0; x < new_w; x++)
		{
			int u = (int)((x - new_w / 2) * cos(-radians) - (y - new_h / 2) * sin(-radians)) + (width / 2);
			int v = (int)((x - new_w / 2) * sin(-radians) + (y - new_h / 2) * cos(-radians)) + (height / 2);

			if (x % 2 == 0 && y % 2 == 0)
				new_img[x + y * new_w] = resample3X3(u, v);
			else if (x % 2 == 1 && y % 2 == 1)
				new_img[x + y * new_w] = resample4X4(u, v);
			else if (x % 2 == 0 && y % 2 == 1)
				new_img[x + y * new_w] = resample3X4(u, v);
			else if (x % 2 == 1 && y % 2 == 0)
				new_img[x + y * new_w] = resample4X3(u, v);
		}

	//set back
	for (int i = 0; i < new_size; i++)
	{
		int offset_rgba = i * 4;
		this->data[offset_rgba + RED] = new_img[i].r * 255;
		this->data[offset_rgba + GREEN] = new_img[i].g * 255;
		this->data[offset_rgba + BLUE] = new_img[i].b * 255;
		this->data[offset_rgba + ALPHA] = new_img[i].a * 255;
	}
	this->width = new_w;
	this->height = new_h;

	delete[] new_img;
    return false;
}// Rotate


//////////////////////////////////////////////////////////////////////////////
//
//      Given a single RGBA pixel return, via the second argument, the RGB
//      equivalent composited with a black background.
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::RGBA_To_RGB(unsigned char *rgba, unsigned char *rgb)
{
    const unsigned char	BACKGROUND[3] = { 0, 0, 0 };

    unsigned char  alpha = rgba[3];

    if (alpha == 0)
    {
        rgb[0] = BACKGROUND[0];
        rgb[1] = BACKGROUND[1];
        rgb[2] = BACKGROUND[2];
    }
    else
    {
	    float	alpha_scale = (float)255 / (float)alpha;
	    int	val;
	    int	i;

	    for (i = 0 ; i < 3 ; i++)
	    {
	        val = (int)floor(rgba[i] * alpha_scale);
	        if (val < 0)
		    rgb[i] = 0;
	        else if (val > 255)
		    rgb[i] = 255;
	        else
		    rgb[i] = val;
	    }
    }
}// RGA_To_RGB
glm::vec4 TargaImage::pixel(int x, int y)
{
	//if (x < 0 || x >= this->width ||
	//	y < 0 || y >= this->height)
	//	return glm::vec4(0.0f);
	
	x = std::max(0, std::min(x, this->width - 1));
	y = std::max(0, std::min(y, this->height - 1));

	int offset = (x * 4) + (y * 4) * this->width;

	glm::vec4 rgba
	{
		this->data[offset + RED] / 255.0f,
		this->data[offset + GREEN] / 255.0f,
		this->data[offset + BLUE] / 255.0f,
		this->data[offset + ALPHA] / 255.0f,
	};

	return rgba;
}
void TargaImage::filtering(double** filter, int n)
{
	for (unsigned int y = 0; y < this->height; y++)
		for (unsigned int x = 0; x < this->width; x++)
		{
			glm::vec4 sum(0.0f);
			for (int i = 0; i < n; i++)
				for (int j = 0; j < n; j++)
				{
					int idx_x = x - n / 2 + i;
					int idx_y = y - n / 2 + j;
					sum += this->pixel(idx_x, idx_y) * (float)filter[i][j];
				}

			int offset_rgba = x * 4 + (y * 4) * this->width;

			this->data[offset_rgba + RED] = sum.r * 255;
			this->data[offset_rgba + GREEN] = sum.g * 255;
			this->data[offset_rgba + BLUE] = sum.b * 255;
			this->data[offset_rgba + ALPHA] = sum.a * 255;
		}
}
glm::vec4 TargaImage::resample3X3(int u, int v)
{
	//use filter to resasmple
	float mask[3][3] =
	{
		 0.0625f, 0.125f, 0.0625f,
		 0.125f,  0.25f,  0.125f,
		 0.0625f, 0.125f, 0.0625f
	};

	glm::vec4 sum(0.0f);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
		{
			int idx_x = u - 1 + i;
			int idx_y = v - 1 + j;
			sum += pixel(idx_x, idx_y) * mask[i][j];
		}

	return sum;
}
glm::vec4 TargaImage::resample4X4(int u, int v)
{
	float mask[4][4] =
	{
		0.015625, 0.046875, 0.046875, 0.015625,
		0.046875, 0.140625, 0.140625, 0.046875,
		0.046875, 0.140625, 0.140625, 0.046875,
		0.015625, 0.046875, 0.046875, 0.015625
	};
	glm::vec4 sum(0.0f);
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
		{
			int idx_x = u - 1 + i;
			int idx_y = v - 1 + j;
			sum += pixel(idx_x, idx_y) * mask[i][j];
		}

	return sum;
}
glm::vec4 TargaImage::resample3X4(int u, int v)
{
	float mask[3][4] =
	{
		0.03125, 0.0625, 0.03125,
		0.09375, 0.1875, 0.09375,
		0.09375, 0.1875, 0.09375,
		0.03125, 0.0625, 0.03125,
	};
	glm::vec4 sum(0.0f);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 4; j++)
		{
			int idx_x = u - 1 + i;
			int idx_y = v - 1 + j;
			sum += pixel(idx_x, idx_y) * mask[i][j];
		}

	return sum;
}
glm::vec4 TargaImage::resample4X3(int u, int v)
{
	float mask[4][3] =
	{
		0.03125, 0.09375, 0.09375, 0.03125,
		0.0625, 0.1875, 0.1875, 0.0625,
		0.03125, 0.09375, 0.09375, 0.03125,
	};
	glm::vec4 sum(0.0f);
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 3; j++)
		{
			int idx_x = u - 1 + i;
			int idx_y = v - 1 + j;
			sum += pixel(idx_x, idx_y) * mask[i][j];
		}

	return sum;
}
void TargaImage::compImage(TargaImage* img, int method)
{
	int size = this->width * this->height;

	//unsigned char to float
	glm::vec4* img_A = toFloat(this->data, this->width, this->height);
	glm::vec4* img_B = toFloat(img->data, img->width, img->height);

	//pre-multiplied pixels
	for (int i = 0; i < size; i++)
	{
		img_A[i].r *= img_A[i].a;
		img_A[i].g *= img_A[i].a;
		img_A[i].b *= img_A[i].a;

		img_B[i].r *= img_B[i].a;
		img_B[i].g *= img_B[i].a;
		img_B[i].b *= img_B[i].a;
	}

	//Composite
	for (int i = 0; i < size; i++)
	{
		float F_A, F_B;

		if (method == _over)
		{
			F_A = 1.0f;
			F_B = 1.0f - img_A[i].a;
		}
		else if (method == _in)
		{
			F_A = img_B[i].a;
			F_B = 0.0f;
		}
		else if (method == _out)
		{
			F_A = 1.0f - img_B[i].a;
			F_B = 0.0f;
		}
		else if (method == _atop)
		{
			F_A = img_B[i].a;
			F_B = 1.0f - img_A[i].a;
		}
		else if (method == _xor)
		{
			F_A = 1.0f - img_B[i].a;
			F_B = 1.0f - img_A[i].a;
		}
		else
		{
			//clear
			F_A = 0.0f;
			F_B = 0.0f;
		}

		glm::vec3 c =
			F_A * glm::vec3(img_A[i]) +
			F_B * glm::vec3(img_B[i]);
		float a = img_A[i].a + img_B[i].a;

		c.r = std::min(std::max(c.r, 0.0f), 1.0f);
		c.g = std::min(std::max(c.g, 0.0f), 1.0f);
		c.b = std::min(std::max(c.b, 0.0f), 1.0f);
		a = std::min(std::max(a, 0.0f), 1.0f);

		//set back
		int offset_rgba = i * 4;
		this->data[offset_rgba + RED] = c.r * 255;
		this->data[offset_rgba + GREEN] = c.g * 255;
		this->data[offset_rgba + BLUE] = c.b * 255;
		this->data[offset_rgba + ALPHA] = a * 255;
	}
	delete[] img_A;
	delete[] img_B;
}
///////////////////////////////////////////////////////////////////////////////
//
//      Copy this into a new image, reversing the rows as it goes. A pointer
//  to the new image is returned.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage* TargaImage::Reverse_Rows(void)
{
    unsigned char   *dest = new unsigned char[width * height * 4];
    TargaImage	    *result;
    int 	        i, j;

    if (! data)
    	return NULL;

    for (i = 0 ; i < height ; i++)
    {
	    int in_offset = (height - i - 1) * width * 4;
	    int out_offset = i * width * 4;

	    for (j = 0 ; j < width ; j++)
        {
	        dest[out_offset + j * 4] = data[in_offset + j * 4];
	        dest[out_offset + j * 4 + 1] = data[in_offset + j * 4 + 1];
	        dest[out_offset + j * 4 + 2] = data[in_offset + j * 4 + 2];
	        dest[out_offset + j * 4 + 3] = data[in_offset + j * 4 + 3];
        }
    }

    result = new TargaImage(width, height, dest);
    delete[] dest;
    return result;
}// Reverse_Rows


///////////////////////////////////////////////////////////////////////////////
//
//      Clear the image to all black.
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::ClearToBlack()
{
    memset(data, 0, width * height * 4);
}// ClearToBlack


///////////////////////////////////////////////////////////////////////////////
//
//      Helper function for the painterly filter; paint a stroke at
// the given location
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::Paint_Stroke(const Stroke& s) {
   int radius_squared = (int)s.radius * (int)s.radius;
   for (int x_off = -((int)s.radius); x_off <= (int)s.radius; x_off++) {
      for (int y_off = -((int)s.radius); y_off <= (int)s.radius; y_off++) {
         int x_loc = (int)s.x + x_off;
         int y_loc = (int)s.y + y_off;
         // are we inside the circle, and inside the image?
         if ((x_loc >= 0 && x_loc < width && y_loc >= 0 && y_loc < height)) {
            int dist_squared = x_off * x_off + y_off * y_off;
            if (dist_squared <= radius_squared) {
               data[(y_loc * width + x_loc) * 4 + 0] = s.r;
               data[(y_loc * width + x_loc) * 4 + 1] = s.g;
               data[(y_loc * width + x_loc) * 4 + 2] = s.b;
               data[(y_loc * width + x_loc) * 4 + 3] = s.a;
            } else if (dist_squared == radius_squared + 1) {
               data[(y_loc * width + x_loc) * 4 + 0] = 
                  (data[(y_loc * width + x_loc) * 4 + 0] + s.r) / 2;
               data[(y_loc * width + x_loc) * 4 + 1] = 
                  (data[(y_loc * width + x_loc) * 4 + 1] + s.g) / 2;
               data[(y_loc * width + x_loc) * 4 + 2] = 
                  (data[(y_loc * width + x_loc) * 4 + 2] + s.b) / 2;
               data[(y_loc * width + x_loc) * 4 + 3] = 
                  (data[(y_loc * width + x_loc) * 4 + 3] + s.a) / 2;
            }
         }
      }
   }
}


///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke() {}

///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke(unsigned int iradius, unsigned int ix, unsigned int iy,
               unsigned char ir, unsigned char ig, unsigned char ib, unsigned char ia) :
   radius(iradius),x(ix),y(iy),r(ir),g(ig),b(ib),a(ia)
{
}

