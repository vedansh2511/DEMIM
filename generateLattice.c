//Program to generate a square lattice
//
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<unistd.h>

double plfvalidation(int np,double **r,double h, double w,int lp)
{
 //printf("hello");
 int i,j,npl=640;
 double a=6.0;
 double pe=0.0;

 //Find the nearest square lattice
 i=0;

 //Lattice coordinates
 double **l;

 //Allocate memory
 l=(double**)calloc(np,sizeof(double*));
 for(i=0;i<np;i++)l[i]=(double*)calloc(3,sizeof(double));

 l[0][0]=0.0;//x
 l[0][1]=h;//y
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


int main()
{

  int i, j, x, ko, ki, xo, xi, cnt = 0;
  double z1, z2;
  int lp;

  // max number of particles
  int n = 2100;

  int np;

  printf("Enter no. of particles: ");
  scanf("%d", &np);

  // Data structure of coordinates
  // Positions r[1----np][0--12]
  // 0--1 positions
  // 2--3 velocities
  // 4--5 forces
  // 6--7 old forces
  // 8    mas
  // 9    rad
  // 10 W_min
  // 11 xold
  // 12 yold
  // 13 broken=0 or not=1
  // 14 density
  // 15 new_S
  // 16 particle interacting with base
  // 17 stores the index of the particle from which the kinetic energy has to be removed in order to cause dispersion of progeny
  // 18 mass specific collision energy that is used to cause breakage (the maximum of all mass specific collision energies (because of multiple interactions))
  // 19 old torque
  // 20 new torque
  // 21 moment of inertia
  // 22 old angular velocity
  // 23 new angular velocity
  // 24 rke
  // 25 lke
  // 26 tangential overlap
  // 30 above particles interacting with each other or not (store the tangential overlap)//if one particle breaks into 2 (20-420)
  double **r;

  // Allocate memory for the variables
  r = (double **)calloc(n, sizeof(double *)); // syntax?
  for (i = 0; i < n; i++)
    r[i] = (double *)calloc(3000, sizeof(double)); // cant we directly write this line?

  double **cw;
  // contact with wall
  cw = (double **)calloc(n, sizeof(double *));
  for (i = 0; i < n; i++)
    cw[i] = (double *)calloc(2, sizeof(double)); // cw[0] stores 1 and 0 for contact ; cw[1] stores total tangential displacement

  double urad, lrad, h, w, mas, t, density;
  double g = 9.81;

  // printf("enter the initial minimum height of drop\n");
  // scanf("%lf",&h);

  //    printf("enter upper dia of particles\n");
  //    scanf("%lf",&urad);
  //
  //    printf("enter lower dia of particles\n");
  //    scanf("%lf",&lrad);

  printf("Enter density of particles: ");
  scanf("%lf", &density);

  // initialize radius of particles
  for (i = 0; i < np; i++)
  {
    // srand(time(NULL));
    // int upper=1000000000;
    // int lower=979795897;
    // double rnd=(rand() % (upper-lower+1))+lower;
    // r[i][9] = rnd*0.5*0.0025/upper;//250*22*(r[i][9]*r[i][9])/7;//for iron//mass for circular area
    r[i][9] = 0.125;
    //printf("%lf\n", r[i][9]);
  }

  for (i = 0; i < np; i++)
  {
    r[i][8] = 1.33333333333333333333 * r[i][9] * r[i][9] * r[i][9] * M_PI * density; // 0.001*((rand() % (umas - lmas + 1)) + lmas);
    //printf("%lf\n", r[i][8]);
  }

  // mas=r[i-1][8];

  printf("Enter simulation time: ");
  scanf("%lf", &t);

  printf("Enter height: ");
  scanf("%lf", &h);

  // weight
  // w=mas*g;

  // double wmin = 0.297/(2*rad);//for glass spheres

  // initialize W_min(minimum mass specific energy for breakage) for particles
  for (i = 0; i < np; i++)
  {
    r[i][10] = 0.000001 / 2 * r[i][9]; // 250*22*(r[i][9]*r[i][9])/7;//for iron//mass for circular area
  }

  // initializs broken or not
  for (i = 0; i < np; i++)
  {
    r[i][13] = 0; // 0 means not broken
  }

  // initialize probability of not breaking
  for (i = 0; i < np; i++)
  {
    r[i][15] = 1;
  }

  // initialize interactions
  for (i = 0; i < np; i++)
  {
    for (j = i + 1; j < np; j++)
    {
      r[i][30 + j] = 0;
      r[j][30 + i] = 0;
    }
  }

  // Initialize coordinates
  double po = plfvalidation(np, r, h, w, lp);


}
