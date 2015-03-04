/* 
 * rotate.cc
 * Joshua Enns
 * CPSC 4310 Assignment 1 test
 */

#include <pam.h>
#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <math.h>

using namespace std;

const double pi = 3.14159;

tuple **read_image(char *filename, pam &inpam);
void write_image(char *filename, const pam &inpam, tuple **array, const int degree);
tuple **rotate_image(tuple **array, const pam &inpam, const int degree);
int minx(const pam &temppam, double radian);
int miny(const pam &temppam, double radian);
int maxx(const pam &temppam, double radian);
int maxy(const pam &temppam, double radian);
double interpolate(tuple **array, float row, float col, int height, int width);

int main(int argc, char *argv[])
{
  int degree = atoi(argv[3]);
  /* structures for input image */
  pam inpam;

  /* a dynamic two-dimensional array to store the pixels...note that
     pnm uses a tuple (for color images with multiple planes) for
     each pixel.  For PGM files it will only be one plane. */
  tuple **array;
  tuple **newarray;
  /* initializes the library */
  pm_init(argv[0], 0);

  /* read the image */
  array = read_image(argv[1], inpam);
  
  /* rotate the image */
  newarray = rotate_image(array, inpam, degree);

  /* write the output */
  write_image(argv[2], inpam, newarray, degree);

  /* clean up */
  pnm_freepamarray(array, &inpam);

  return 0;
}

/* reads the image into the netpbm structure */
tuple **read_image(char *filename, pam &inpam)
{
  FILE *f;
  tuple **A;

  if ((f = pm_openr(filename)) == NULL) {
    cerr << "Cannot open file \"" << filename << "\" for reading." << endl;
    exit(1);
  }

  if ((A = pnm_readpam(f, &inpam, PAM_STRUCT_SIZE(tuple_type))) == NULL) {
    cerr << "Cannot read image \"" << filename << "\"." << endl;
    exit(1);
  }
  
  pm_close(f);
  return A;
}

/* writes the image to the given file */
void write_image(char *filename, const pam &inpam, tuple **array, const int degree)
{
  FILE *f;
  pam outpam = inpam;
  
  if ((f = pm_openw(filename)) == NULL) {
    cerr << "Cannot open file \"" << filename << "\" for writing." << endl;
    exit(1);
  }

  outpam.file = f;

  //Calculate and set output image dimensions
  double radian = degree * pi/180 * 1;
  int xmin, ymin, xmax, ymax;
  xmin = minx(inpam, radian);
  ymin = miny(inpam, radian);
  xmax = maxx(inpam, radian);
  ymax = maxy(inpam, radian);

  outpam.width = ymax - ymin;
  outpam.height = xmax - xmin;
  pnm_writepam(&outpam, array);

  pm_close(f);
}

/* rotates the image */
tuple **rotate_image(tuple **array, const pam &inpam, const int degree)
{
  //Create variables for function
  double radian;
  int xmin, ymin, xmax, ymax, row, col;
  float prevrow, prevcol;
  int pixelval = 0;
  
  //Convert degrees to radians and invert to implement inverse mapping
  radian = degree * pi/180 * -1;
  
  //set dimensions for after rotation
  xmin = minx(inpam, radian);
  ymin = miny(inpam, radian);
  xmax = maxx(inpam, radian);
  ymax = maxy(inpam, radian);

  //Create a temp container to place rotated image into
  pam temppam = inpam;
  temppam.height =  ymax - ymin;
  temppam.width = xmax - xmin;
  tuple **newarray;
  newarray = pnm_allocpamarray(&temppam);

  for (row = ymin; row < temppam.height + ymin; row++) {

    for (col = xmin; col < temppam.width + xmin; col++) {
      //Calculate coordinates of inversely rotated image
      prevrow = (row*cos(radian)) - (col*sin(radian));
      prevcol = (col*cos(radian)) + (row*sin(radian));

      if ((prevrow >= 0)&&(prevcol >= 0)&&(prevrow < inpam.height)&&(prevcol < inpam.width)){
        pixelval = interpolate (array, prevrow, prevcol, inpam.height, inpam.width);
      }
      else {
        pixelval = 0;
      }
      //Set pixel value of rotated image
      newarray[row+abs(ymin)][col+abs(xmin)][0] = pixelval;

    }
  }
  return newarray;
}

double interpolate(tuple **array, float row, float col, int height, int width){
  //set x and y as floor of row and col
  int x1 = (int)row;
  int x2 = x1+1;
  int y1 = (int)col;
  int y2 = y1+1;

  //Calculate intensities of neighbors
  double q11, q12, q21, q22 = 0;
  q11 = array[x1][y1][0];
  if (y2 < width) {
    q12 = array[x1][y2][0];
  }
  if(x2 < height) {
    q21 = array[x2][y1][0];
  }
  if((x2 < height) && (y2 < width)){
    q22 = array[x2][y2][0];
  }

  //Perform linear interpolation in the x-direction
  double r1 = ((x2-row)/(x2-x1) * q11) + ((row-x1)/(x2-x1) * q21);
  double r2 = ((x2-row)/(x2-x1) * q12) + ((row-x1)/(x2-x1) * q22);

  //Perform linear interpolation in the y-direction
  double i = (y2-col)/(y2-y1)*r1 + (col-y1)/(y2-y1)*r2;
  //std::cout << i << "\n";

  //return sum of all neighboring pixel intensities
  return i;
}

int minx(const pam &temppam, double radian)
{
  int x0y0, x0y1, x1y0, x1y1, currmin, ct;
  //calculate rotated coordinates of all four corners 
  x0y0 = 0;
  x0y1 = 0 - (temppam.height*sin(radian));
  x1y0 = (temppam.width*cos(radian));
  x1y1 = (temppam.width*cos(radian)) - (temppam.height*sin(radian));

  int minarray[4] = {x0y0, x0y1, x1y0, x1y1};

  //Find lowest coordinate and return
  currmin = x0y0;
  for (ct = 0; ct < 4; ct++){
    currmin = min (currmin, minarray[ct]);
  }
  return floor(currmin);

}

int maxx(const pam &temppam, double radian)
{
  int x0y0, x0y1, x1y0, x1y1, currmax, ct;
  //calculate rotated coordinates of all four corners 
  x0y0 = 0;
  x0y1 = 0 - (temppam.height*sin(radian));
  x1y0 = (temppam.width*cos(radian));
  x1y1 = (temppam.width*cos(radian)) - (temppam.height*sin(radian));
  
  int maxarray[4] = {x0y0, x0y1, x1y0, x1y1};
  //Find highest coordinate and return
  currmax = x0y0;
  for (ct = 0; ct < 4; ct++){
    currmax = max (currmax, maxarray[ct]);
  }
  return ceil(currmax);
}

int miny(const pam &temppam, double radian)
{

  int x0y0, x0y1, x1y0, x1y1, currmin, ct;
  //calculate rotated coordinates of all four corners
  x0y0 = 0;
  x0y1 = (temppam.height*sin(radian));
  x1y0 = (temppam.width*cos(radian));
  x1y1 = (temppam.width*cos(radian)) + (temppam.height*sin(radian));

  int minarray[4] = {x0y0, x0y1, x1y0, x1y1};
  //Find lowest coordinate and return
  currmin = x0y0;
  for (ct = 0; ct < 4; ct++){
    currmin = min (currmin, minarray[ct]);
  }
  return floor(currmin);
}


int maxy(const pam &temppam, double radian)
{
  int x0y0, x0y1, x1y0, x1y1, currmax, ct;
  //calculate rotated coordinates of all four corners
  x0y0 = 0;
  x0y1 = (temppam.height*sin(radian));
  x1y0 = (temppam.width*cos(radian));
  x1y1 = (temppam.width*cos(radian)) + (temppam.height*sin(radian));
  int maxarray[4] = {x0y0, x0y1, x1y0, x1y1};
  //Find lowest coordinate and return
  currmax = x0y0;
  for (ct = 0; ct < 4; ct++){
    currmax = max (currmax, maxarray[ct]);
  }
  return ceil(currmax);
}

