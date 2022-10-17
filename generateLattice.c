//Program to generate 3D lattice

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<unistd.h>

double generateLattice(int np,double **r,double h, double w,int lp)
{
 
//npl is number of particles
 int i,j,npl=640;
 double a=6.0;
 
//pe is potential energy
 double pe=0.0;

 //Find the nearest square lattice
 i=0;

 //Lattice coordinates
 double **l;

 //Allocate memory
 l=(double**)calloc(np,sizeof(double*));
 for(i=0;i<np;i++)l[i]=(double*)calloc(3,sizeof(double));

//x coordinate
 l[0][0]=0.0;
 
//y coordinate
 l[0][1]=h;
 
//z coordinate
 l[0][2]=0.0;
 

 printf("0  ");
 printf("%lf ", l[0][0]);
 printf("%lf ", l[0][1]);
 printf("%lf\n", l[0][2]);
 
 int li=1;

 for (i = 1; i < np; i++) // number of layer
 {
   if (i < 4)
   {
     l[i][0] = l[i - 1][0] + r[i - 1][9] + 2 * r[i][9] + 0.00022; //[i][9];                           //increment in x
     l[i][1] = l[i - 1][1];
     l[i][2] = l[i - 1][2];
   }

   else if (i < 16)
   {
     l[i][0] = l[i - 1][0] + r[i - 1][9] + 2 * r[i][9] + 0.00022; //[i][9];                           //increment in x
     l[i][1] = l[i - 4][1] + r[i - 4][9] + 2 * r[i][9];           //+0.00001;
     l[i][2] = l[i - 1][2];
   }

   else
   {
     l[i][0] = l[i - 1][0] + r[i - 1][9] + 2 * r[i][9] + 0.00022; //[i][9];                           //increment in x
     l[i][1] = l[i - 4][1] + r[i - 4][9] + 2 * r[i][9];           //+0.00001;
     if (i % 16 == 0)
     {
      l[i][2] = l[i - 1][2] + r[i - 1][9] + 2 * r[i][9] + 0.00022;
     }
     else
     {
      l[i][2] = l[i - 1][2];
     }
   }

   // If one line is filled start, start a new line
   if (i % 4 == 0)
   {
     li++;
     if (li % 2 == 0)
     {
       l[i][0] = l[i - 4][0] + 1.55 * r[i][9] + 0.00022;  //[i][9];  +0.525*r[i][9]                                 //increment in x
       l[i][1] = l[i - 4][1] + r[i - 4][9] + 2 * r[i][9]; //+0.00001;//[i][9];                    //increment in y
     }
     else
     {
       l[i][0] = l[0][0];
       l[i][1] = l[i - 4][1] + r[i - 4][9] + 2 * r[i][9];
     }
   }

   printf("%d  ", i);
   printf("%lf ", l[i][0]);
   printf("%lf ", l[i][1]);
   printf("%lf\n", l[i][2]);
 }

 FILE *out=fopen("generate3d.dat","w");     //do we place fopen and fclose in such a way that any results generated between the 2 statements get written in the .dat file?
   for(i=0;i<np;i++)
      {
      r[i][0]=l[i][0];                      //why multiplying with scaling factor?
      r[i][1]=l[i][1];

      pe+=w*r[i][1];
      fprintf(out,"create_atoms   %d single  %lf %lf %lf %lf units box\n",1,r[i][0],0.0,r[i][1], l[i][2]);   //fprintf?
      }

   for(i=1;i<=np;i++)
      {
      fprintf(out,"set  atom %d diameter %lf density %lf vx %lf vy %lf vz %lf\n",i,2*r[i][9],920.0,0.0,0.0,0.0);   //fprintf?
      }

   fprintf(out,"%lf\n",pe);
   fclose(out);
   return pe;
 }
