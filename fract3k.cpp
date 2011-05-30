#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <strstream>
#include <cv.h>
#include <highgui.h>

using namespace std;

#define LEVELS 8
#define SIZE (1 + (1 << LEVELS))
#define f3(delta,x0,x1,x2) (x0+x1+x2)/3+delta*Gauss()
#define f4(delta,x0,x1,x2,x3) (x0+x1+x2+x3)/4+delta*Gauss()
#define GetPixel(img, u, v) &(img->imageData[v*img->widthStep+u*img->depth/8*img->nChannels])

#if defined (_WIN32)
    #define SEED srand
    #define RANDNUM rand
    #define RANDMAX RAND_MAX
#else
    #define SEED srandom
    #define RANDNUM random
    #define RANDMAX (1u<<31)-1
#endif

typedef double grid[SIZE][SIZE];

// Globals
int Arand, Nrand;
double GaussAdd, GaussFac, delta, h_mid, h_diff, h_low, h_high;
grid X;

// Declarations
void InitGauss (int seed);
double Gauss ();
void MidPointFM2D (grid X, int maxlevel, double sigma, double H, int addition,
                   int seed);
void h_values ();
void MidPointFM2D_bad (IplImage* img, int maxlevel, double sigma, double H, int addition,
                   int seed); // Don't use this one.

int main(int argc, char **argv)
{
	int additions = 1, pixvalue, userhi, userlo;
	double H = .7, p = 1.0;
	unsigned int x, y;
	IplImage *mapimg = 0, *colorimg = 0;
	char *pixel;
	cvvInitSystem(argc, argv);

	cout << "Jaggedness of mountains (0.0-1.0, >=.7 recommended): ";
	cin >> H;
	cout << "Steepness of mountains (0<x<3): ";
	cin >> p;
	cout << "Maximum mountain height (32-255): ";
	cin >> userhi;
	cout << "Minimum ocean depth (0-31): ";
	cin >> userlo;
/*	if (userlo > userhi)
	{
		pixvalue = userlo;
		userlo = userhi;
		userhi = pixvalue;
	}
*/	
	cvvNamedWindow("Map", 0);
	cvvNamedWindow("Colored Map", 0);
	mapimg = cvCreateImage(cvSize(SIZE, SIZE), IPL_DEPTH_8U, 3);
	colorimg = cvCreateImage(cvSize(SIZE,SIZE), IPL_DEPTH_8U, 3);

	MidPointFM2D(X, LEVELS, 1.0, H, additions, (int)time(NULL));
//	MidPointFM2D_bad(mapimg, LEVELS, 1.0, H, additions, (int)time(NULL));
	h_values();

	cout << "h_mid:  " << h_mid << endl;
	cout << "h_diff: " << h_diff << endl;

	/* Converting grid format to image format */
	for (y=0; y<SIZE; y++)
		for (x=0; x<SIZE; x++)
		{
			// Maybe use powers to get more water vs. mtns: pow((X-l)/(h-l), p)
			pixvalue = (int)(userlo+((userhi-userlo)*pow((X[x][y]-h_low)/(h_high-h_low),p)+.5));

			// For gray image
			pixel = GetPixel(mapimg, x, y);
			pixel[0] = pixvalue;
			pixel[1] = pixvalue;
			pixel[2] = pixvalue;

			// For color image.. [0]=Blue, [1]=Green, [2]=Red 
			pixel = GetPixel(colorimg, x, y);
			pixel[2] = 0;
			if (pixvalue < 31)
			{
				pixel[0] = 200;
				pixel[1] = 0;
			}
			else if (pixvalue < 32)
			{	//tan for beach
				pixel[0] = 145;
				pixel[1] = 200;
				pixel[2] = 210;
			}
			else if (pixvalue < 86) // green to yellow
			{	// lo+(hi-lo)*(curr-logray)/steps
				pixel[0] = 0+(25-0)*(pixvalue-32)/53;
				pixel[1] = 33+(145-33)*(pixvalue-32)/53;
				pixel[2] = 0+(155-0)*(pixvalue-32)/53;
			}
			else if (pixvalue < 126) // yellow to lt.brown
			{	// lo+(hi-lo)*(curr-logray)/steps
				pixel[0] = 25+(25-25)*(pixvalue-86)/39;
				pixel[1] = 145+(60-145)*(pixvalue-86)/39;
				pixel[2] = 155+(85-155)*(pixvalue-86)/39;
			}
			else if (pixvalue < 216) // lt.brown to dk.brown
			{	// lo+(hi-lo)*(curr-logray)/steps
				pixel[0] = 25+(3-25)*(pixvalue-126)/89;
				pixel[1] = 60+(10-60)*(pixvalue-126)/89;
				pixel[2] = 85+(30-85)*(pixvalue-126)/89;
			}
			else
			{	// whites
				pixel[0] = pixvalue-25;
				pixel[1] = pixvalue-25;
				pixel[2] = pixvalue-25;
			}
			
		}

	cvvResizeWindow("Map", SIZE, SIZE);
	cvvShowImage("Map", mapimg);
	cvvResizeWindow("Colored Map", 2*SIZE, 2*SIZE);
	cvvShowImage("Colored Map", colorimg);
	cvvWaitKey(0);
	cvvSaveImage("map.bmp", mapimg);

	return(0);
}

void InitGauss (int seed)
/* Routine for initializing the Gaussian random number generator.  This is an
 * implementation of algorithm InitGauss on page 77 of "The Science of Fractal
 * Images".
 */
{
   Nrand = 4;
   Arand = RANDMAX;
   GaussAdd = sqrt (3 * Nrand);
   GaussFac = 2.0 * GaussAdd / ((double)Nrand * (double)Arand);
   SEED (seed);
}

double Gauss ()
/* Routine to generate a Gaussian random number.  This is an implementation of
 * algorithm Gauss on page 77 of "The Science of Fractal Images."
 */
{
   double sum;
   int i;
   
   sum = 0.0;
   for (i=1; i<=Nrand; i++) sum += (double) RANDNUM();
   return (GaussFac * sum - GaussAdd);
}

void MidPointFM2D (grid X, int maxlevel, double sigma, double H, 
   int addition, int seed)
/* Routine for computing the points that simulate fractional Brownian motion
 * via midpoint displacement.  This is an implementation of algorithm
 * MidPointFM2D on page 100 of "The Science of Fractal Images."
 */
{
   int N, stage;
   double delta, hpow;
   int x, y, D, d;
   
   InitGauss (seed);
   N = (int) pow (2.0, (double)maxlevel);
   
   /* Set the initial random corners. */
   
   delta = sigma;
   X[0][0] = delta * Gauss();
   X[0][N] = delta * Gauss();
   X[N][0] = delta * Gauss();
   X[N][N] = delta * Gauss();
   D = N;
   d = N >> 1;
   hpow = pow (0.5, 0.5*H);
   for (stage=1; stage<=maxlevel; stage++)
   {
      cout << "Computing level " << stage << " of " << maxlevel << " total...\r";
      cout.flush();
  
      /* Going from grid type I to type II. */
      
      delta *= hpow;
      
      /* Interpolate and offset points. */
      
      for (x=d; x<=N-d; x+=D)
         for (y=d; y<=N-d; y+=D)
            X[x][y] = f4(delta,X[x+d][y+d],X[x+d][y-d],X[x-d][y+d],X[x-d][y-d]);

      /* Displace other points also if needed. */
      
      if (addition)
      {
         for (x=0; x<=N; x+=D)
            for (y=0; y<=N; y+=D)
               X[x][y] += delta * Gauss();
      }
      
      /* Going from grid type II to type I. */
      
      delta *= hpow;
      
      /* Interpolate and offset boundary grid points. */
      
      for (x=d; x<=N-d; x+=D)
      {
         X[x][0] = f3(delta, X[x+d][0], X[x-d][0], X[x][d]);
         X[x][N] = f3(delta, X[x+d][N], X[x-d][N], X[x][N-d]);
         X[0][x] = f3(delta, X[0][x+d], X[0][x-d], X[d][x]);
         X[N][x] = f3(delta, X[N][x+d], X[N][x-d], X[N-d][x]);
      }
      
      /* Interpolate and offset interior grid points. */
      
      for (x=d; x<=N-d; x+=D)
         for (y=D; y<=N-d; y+=D)
            X[x][y] = f4(delta, X[x][y+d], X[x][y-d], X[x+d][y], X[x-d][y]);
      for (x=D; x<=N-d; x+=D)
         for (y=d; y<=N-d; y+=D)
            X[x][y] = f4(delta, X[x][y+d], X[x][y-d], X[x+d][y], X[x-d][y]);
      
      /* Displace other points also if needed. */
      
      if (addition)
      {
         for (x=0; x<=N; x+=D)
            for (y=0; y<=N; y+=D)
               X[x][y] += delta * Gauss();
         for (x=d; x<=N-d; x+=D)
            for (y=d; y<=N-d; y+=D)
               X[x][y] += delta * Gauss();
      }
      D >>= 1;
      d >>= 1;
   }
   cout << endl;
}

void h_values ()
/* This routine determines the h_mid and h_diff values needed for the color
   map. */
{
   double *x;
   int i;
   
   h_high = -(h_low = 10000.0);
   for (x=(double *)X, i=0; i<SIZE*SIZE; i++, x++)
      if (*x > h_high)
         h_high = *x;
      else if (*x < h_low)
         h_low = *x;
   h_mid = (h_low + h_high) / 2.0;
   h_diff = h_high - h_mid;
   cout << "h_high: " << h_high << endl;
   cout << "h_low:  " << h_low << endl;
}



// ** The following function is a failed attempt **
void MidPointFM2D_bad (IplImage* img, int maxlevel, double sigma, double H, 
   int addition, int seed)
/* Routine for computing the points that simulate fractional Brownian motion
 * via midpoint displacement.  This is an implementation of algorithm
 * MidPointFM2D on page 100 of "The Science of Fractal Images."
 */
{
   int N, stage;
   double delta, hpow;
   int x, y, D, d;
   char *X, *pix0, *pix1, *pix2, *pix3;
   
   InitGauss (seed);
   N = (int) pow (2.0, (double)maxlevel);
   
   /* Set the initial random corners. */
   
   delta = sigma;
   X = GetPixel(img, 0, 0);
   X[0] = delta * Gauss();
   X = GetPixel(img, 0, N);
   X[0] = delta * Gauss();
   X = GetPixel(img, 0, N);
   X[0] = delta * Gauss();
   X = GetPixel(img, 0, N);
   X[0] = delta * Gauss();
   D = N;
   d = N >> 1;
   hpow = pow (0.5, 0.5*H);
   for (stage=1; stage<=maxlevel; stage++)
   {
      cout << "Computing level " << stage << " of " << maxlevel << " total...\r";
      cout.flush();
  
      /* Going from grid type I to type II. */
      
      delta *= hpow;
      
      /* Interpolate and offset points. */
      
      for (x=d; x<=N-d; x+=D)
         for (y=d; y<=N-d; y+=D)
		 {
			 X = GetPixel(img, x, y);
			 pix0 = GetPixel(img, x+d, y+d);
			 pix1 = GetPixel(img, x+d, y-d);
			 pix2 = GetPixel(img, x-d, y+d);
			 pix3 = GetPixel(img, x-d, y-d);
			 X[0] = f4(delta,pix0[0],pix1[0],pix2[0],pix3[0]);
		 }

      /* Displace other points also if needed. */
      
      if (addition)
      {
         for (x=0; x<=N; x+=D)
            for (y=0; y<=N; y+=D)
			{
				X = GetPixel(img, x, y);
				X[0] += delta * Gauss();
			}
      }
      
      /* Going from grid type II to type I. */
      
      delta *= hpow;
      
      /* Interpolate and offset boundary grid points. */
      
      for (x=d; x<=N-d; x+=D)
      {
		  X = GetPixel(img, x, 0);
		  pix0 = GetPixel(img, x+d, 0);
		  pix1 = GetPixel(img, x-d, 0);
		  pix2 = GetPixel(img, x, d);
		  X[0] = f3(delta, pix0[0], pix1[0], pix2[0]);

		  X = GetPixel(img, x, N);
		  pix0 = GetPixel(img, x+d, N);
		  pix1 = GetPixel(img, x-d, N);
		  pix2 = GetPixel(img, x, N-d);
          X[0] = f3(delta, pix0[0], pix1[0], pix2[0]);

		  X = GetPixel(img, 0, x);
		  pix0 = GetPixel(img, 0, x+d);
		  pix1 = GetPixel(img, 0, x-d);
		  pix2 = GetPixel(img, d, x);
          X[0] = f3(delta, pix0[0], pix1[0], pix2[0]);

		  X = GetPixel(img, N, x);
		  pix0 = GetPixel(img, N, x+d);
		  pix1 = GetPixel(img, N, x-d);
		  pix2 = GetPixel(img, N-d, x);
          X[0] = f3(delta, pix0[0], pix1[0], pix2[0]);
      }
      
      /* Interpolate and offset interior grid points. */
      
      for (x=d; x<=N-d; x+=D)
         for (y=D; y<=N-d; y+=D)
		 {
			 X = GetPixel(img, x, y);
			 pix0 = GetPixel(img, x, y+d);
			 pix1 = GetPixel(img, x, y-d);
			 pix2 = GetPixel(img, x+d, y);
			 pix3 = GetPixel(img, x-d, y);
             X[0] = f4(delta, pix0[0], pix1[0], pix2[0], pix3[0]);
		 }
      for (x=D; x<=N-d; x+=D)
         for (y=d; y<=N-d; y+=D)
		 {
			 X = GetPixel(img, x, y);
			 pix0 = GetPixel(img, x, y+d);
			 pix1 = GetPixel(img, x, y-d);
			 pix2 = GetPixel(img, x+d, y);
			 pix3 = GetPixel(img, x-d, y);
             X[0] = f4(delta, pix0[0], pix1[0], pix2[0], pix3[0]);
		 }
      
      /* Displace other points also if needed. */
      
      if (addition)
      {
         for (x=0; x<=N; x+=D)
            for (y=0; y<=N; y+=D)
			{
				X = GetPixel(img, x, y);
                X[0] += delta * Gauss();
			}
         for (x=d; x<=N-d; x+=D)
            for (y=d; y<=N-d; y+=D)
			{
				X = GetPixel(img, x, y);
                X[0] += delta * Gauss();
			}
      }
      D >>= 1;
      d >>= 1;
   }
   cout << endl;
} 