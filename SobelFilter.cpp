#include "ReadWrite.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cv.h>
#include "highgui.h"
#include <conio.h>

void main()
{
	unsigned char *in = NULL;
	int ht = 0;
	int wd = 0;
	int a, b;
	int	GMASK[5][5];
	int Dxx[3][3];
	int Dyy[3][3];
	int Dxy[3][3];
	unsigned int X, Y;
	int I, J, t;
	double  sumX, sumY, sumZ, sumXT, sumYT, sumZT;
	double pi=3.141593, SUM;
	long double Xf,Yf,Zf;
	double theta;
	Xf=Yf=Zf=0;

	printf("Enter The Value For Theta\n");
	scanf("%d",&t);
	theta = (-t*(pi/180));
	
	/*gaussian GMASK*/
   GMASK[0][0] = 1; GMASK[0][1] = 4; GMASK[0][2] = 7; GMASK[0][3] = 4; GMASK[0][4] = 1;
   GMASK[1][0] = 4; GMASK[1][1] = 16; GMASK[1][2] = 26; GMASK[1][3] = 16; GMASK[1][4] = 4;
   GMASK[2][0] = 7; GMASK[2][1] = 26; GMASK[2][2] = 41; GMASK[2][3] = 26; GMASK[2][4] = 7;
   GMASK[3][0] = 4; GMASK[3][1] = 16; GMASK[3][2] = 26; GMASK[3][3] = 16; GMASK[3][4] = 4;
   GMASK[4][0] = -1; GMASK[4][1] = 4; GMASK[4][2] = 7; GMASK[4][3] = 4; GMASK[4][4] = 1;

				
	
	Dxx[0][0] = 1;	Dxx[0][1] = -2;	Dxx[0][2] = 1;
	Dxx[1][0] = 0;	Dxx[1][1] = 0;	Dxx[1][2] = 0;
	Dxx[2][0] = 0;	Dxx[2][1] = 0;	Dxx[2][2] = 0;

	Dyy[0][0] = 1;	Dyy[0][1] = 0;	Dyy[0][2] = 0;
	Dyy[1][0] = -2;	Dyy[1][1] = 0;	Dyy[1][2] = 0;
	Dyy[2][0] = 1;	Dyy[2][1] = 0;	Dyy[2][2] = 0;

	Dxy[0][0] = 1;	Dxy[0][1] = 0;	Dxy[0][2] = 1;
	Dxy[1][0] = 0;	Dxy[1][1] = -4;	Dxy[1][2] = 0;
	Dxy[2][0] = 1;	Dxy[2][1] = 0;	Dxy[2][2] = 1;	
			
			
	//reading the PGM file format
	read_pgm_image("grandeur29-NIR-000231.pgm", &in, &ht, &wd);

	//allocating the memory to get the output image
	unsigned char *out = (unsigned char*)calloc((ht)*(wd), sizeof(unsigned char));
	

	//Gauss Algo Starts Here.
	for(Y=0; Y<=(ht-1); Y++)  {
		for(X=0; X<=(wd-1); X++)  {
			SUM = 0;


			/* image boundaries */
			if(Y==0 || Y==1 || Y==ht-2 || Y==ht-1)
				SUM = 0;
			else if(X==0 || X==1 || X==wd-2 || X==wd-1)
				SUM = 0;

			/* Convolution starts here */
			else   {

				for(I=-2; I<=2; I++)  {
					for(J=-2; J<=2; J++)  {
						SUM = SUM + (int)( (*(in + X + I + (Y + J)*wd)) * GMASK[I+2][J+2]);
					}
				}
			}

			*(out + X + Y*wd) = (unsigned char)((SUM)/273);
		}
	}
	write_pgm_image("Gauss.pgm", out, ht, wd, "", 255);

	
	
	//reading the PGM file format
	read_pgm_image("Gauss.pgm", &in, &ht, &wd);
	
	//allocating the memory to get the output image
	unsigned char *out1 = (unsigned char*)calloc((ht)*(wd), sizeof(unsigned char));
		
	/* operation starts here*/

	for(Y=0; Y<=(ht-1); Y++)  {
		for(X=0; X<=(wd-1); X++)  {
			sumX = 0;
			sumY = 0;
			sumZ = 0;

			/* image boundaries */
			if(Y==0 || Y==ht-1)
				SUM = 0;
			else if(X==0 || X==wd-1)
				SUM = 0;

			/* Convolution starts here */
			else   {

				/*-------X GRADIENT APPROXIMATION------*/
				for(I=-1; I<=1; I++)  {
					for(J=-1; J<=1; J++)  {
						sumX = sumX + (int)( (*(in + X + I + (Y + J)*wd)) * Dxx[I+1][J+1]);
						sumXT = sumX * cos(theta) * cos(theta);
					}
				}

				/*-------Y GRADIENT APPROXIMATION-------*/
				for(I=-1; I<=1; I++)  {
					for(J=-1; J<=1; J++)  {
						sumY = sumY + (int)( (*(in + X + I + (Y + J)*wd)) * Dyy[I+1][J+1]);
						sumYT = sumY * sin(theta) * sin(theta);
					}
				}

				/*-------XY GRADIENT APPROXIMATION-------*/
				for(I=-1; I<=1; I++)  {
					for(J=-1; J<=1; J++)  {
						sumZ = sumZ + (int)( (*(in + X + I + (Y + J)*wd)) * Dxy[I+1][J+1]);
						sumZT = sumZ * cos(theta) * sin(theta);
					}
				}
				
				
				/*---GRADIENT MAGNITUDE APPROXIMATION (Myler p.218)----*/
				
				SUM = abs(sumXT) + abs(sumYT) + abs(sumZT);
				Xf= Xf + sumX; Yf= Yf + sumY; Zf= Zf + sumZ; 
			} 
			
			if(SUM>10) SUM=255;
            else SUM=0;
		
			
			*(out1 + X + Y*wd) = 255 - (unsigned char)((SUM));
		
		}
	
	}
	

	printf("values of x y z are %Lf %Lf %Lf \n",Xf,Yf,Zf);
	
	float A = sqrt((Xf * Xf) - (2 * Xf * Yf) + (Yf * Yf) + (4 * Zf));

	long double Tmin = atan((Xf-Yf-A)/(2*Zf));
	printf("\nTmin = %Lf\n",Tmin);
	long double Tmax = atan((Xf-Yf+A)/(2*Zf));
	printf("\nTmax = %Lf\n",Tmax);
	write_pgm_image("Out360.pgm", out1, ht, wd, "", 255);
	getch();
}
