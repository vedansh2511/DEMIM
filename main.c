#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

// Define value of Pi
#define M_PI 3.142857

// generate a lattice
#include "generateLattice.c"

// Calculate forces and Potential energy
#include "pforcef.c"

// Integrate the equation of motion
#include "pmovef.c"

// placement of particle
#include "replace.c"

int main()
{
  int i, j, x, ko, ki, xo, xi, cnt = 0;
  double z1, z2;
  int lp;

  // max number of particles
  int n = 2100;

  int np;

  printf("enter no. of particles");
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

  printf("enter density of particles\n");
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
    printf("%lf\n", r[i][9]);
  }

  for (i = 0; i < np; i++)
  {
    r[i][8] = 1.33333333333333333333 * r[i][9] * r[i][9] * r[i][9] * M_PI * density; // 0.001*((rand() % (umas - lmas + 1)) + lmas);
    printf("%lf\n", r[i][8]);
  }

  // mas=r[i-1][8];

  printf("enter simulation time\n");
  scanf("%lf", &t);

  printf("enter height\n");
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

  for (i = 0; i < n; i++)
  {
    if (r[i][0] - r[i][9] < -0.125)
    {
      cw[i][0] = 1; // in contact with wall
    }
    else
    {
      cw[i][0] = 0; // not in contact with wall
    }
  }

  for (i = 0; i < n; i++)
  {
    cw[i][1] = 0; // initial tangential displacement
  }

  double **cb;
  // contact with base
  cb = (double **)calloc(n, sizeof(double *));
  for (i = 0; i < n; i++)
    cb[i] = (double *)calloc(2, sizeof(double)); // cw[0] stores 1 and 0 for contact ; cw[1] stores total tangential displacement

  for (i = 0; i < n; i++)
  {
    cb[i][0] = 0; // not in contact with base
  }

  for (i = 0; i < n; i++)
  {
    cb[i][1] = 0; // initial tangential displacement
  }

  double *csw; // contact with sidewall
  // contact with side-wall
  csw = (double *)calloc(n, sizeof(double));
  // for(i=0;i<n;i++)csw[i]=(double*)calloc(1,sizeof(double));  //cw[0] stores 1 and 0 for contact ; cw[1] stores total tangential displacement

  for (i = 0; i < n; i++)
  {
    csw[i] = 0; // initial tangential displacement
  }

  // Initialize velocities and density and mass specific collision energy
  for (i = 0; i < np; i++)
  {
    r[i][2] = 0;
    r[i][3] = 0;
    //       r[0][0]=0;
    //       r[0][1]=0.1;
    ////
    //       r[0][2]=10;
    //       r[0][3]=0;
    ////
    ////       r[0][11]=0;
    ////       r[0][12]=0.1;
    ////
    //       r[1][0]=0.02;
    //       r[1][1]=0.1;
    ////
    //       r[1][2]=-10;
    //       r[1][3]=-5;
    //
    //       r[1][11]=0.1;
    //       r[1][12]=0.1;

    r[i][14] = density;
    r[i][18] = 0; // mass specific collision energy is initially zero since no collisions happening
    r[i][16] = 0; // no interactions with the base initially
    r[i][19] = 0; // old torque
    r[i][20] = 0; // new torque
    r[i][22] = 0; // old angular velocity
    r[i][23] = 0; // new angular velocity
    r[i][26] = 0; // tangential overlap
  }

  //      r[0][0]=0;
  //      r[0][1]=2.01;
  //      r[0][2]=5;
  //      r[0][3]=0;
  //
  //      r[1][0]=0.032;
  //      r[1][1]=2.0;
  //      r[1][2]=-5;
  //      r[1][3]=0;

  // double E = 70000000000.0;//young's modulus for Al-alloy

  // Ep = 109251528702;//E/(2*(1 - 0.291*0.291)) poisson's ratio = 0.291

  // double kb = 1.33333333*E*1.098901099*(sqrt(rad));//E'=E/(1-nu^2)             500000000000000000

  // double k = 1.33333333*E*0.549450549*(sqrt(rad/2.0));//E'=E/2*(1-nu^2)

  // Time step dt
  double dt = 0.000018;

  // Number of steps
  int nst = (int)(t / dt);

  // potential energy
  double pe;
  pe = 0.0;

  // Kinetic Energy
  double ke;
  ke = 0.0;

  double *b; // force from base
  b = (double *)calloc(n, sizeof(double));
  for (i = 0; i < n; i++)
  {
    b[i] = 0.0;
  }

  double *uold;
  uold = (double *)calloc(n, sizeof(double));

  double *vold;
  vold = (double *)calloc(n, sizeof(double));

  // initialize old velocities(needed in the calculation of friction)
  for (i = 0; i < n; i++)
  {
    uold[i] = 0.0;
    vold[i] = 0.0;
  }

  int flag = 0;

  FILE *files[60000];
  char filename[200];

  FILE *out = fopen("generallattice2d_10poblique_breakage probability_replacement_lessdense_nodamping_velocities_Al-alloy_3.4e-8_g9.81_r.01_t1.csv", "w"); // general2d_2p_nodamping_Al-alloy_5e-9_x1_0_x2_0.1_y1_0.2_y2_0.41_u1_0_u2_0_v1_0_v2_0_r.1_t.24.csv
  // fprintf(out,"ke+pe,ke,pe,time,x0,y0,x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7,x8,y8,x9,y9,x10,y10,x11,y11,x12,y12,x13,y13,x14,y14,x15,y15,x16,y16,x17,y17,x18,y18,x19,y19,x20,y20,x21,y21,x22,y22,x23,y23,x24,y24,x25,y25,x40,x49,x50,x55,x56,x60,x64,x69,x74,x79,x130,x169,z1,z2,fx0,fy0,fx1,fy1,fx2,fy2,fx3,fy3,fx4,fy4,fx5,fy5,fx6,fy6,fx7,fy7,fx8,fy8,fx9,fy9,fx10,fy10,fx11,fy11,fx12,fy12,fx13,fy13,fx14,fy14,fx15,fy15,fx16,fy16,fx17,fy17,fx18,fy18,fx19,fy19,fx20,fy20,fx21,fy21,fx22,fy22,fx23,fy23,fx24,fy24,fx25,fy25,u1,v1,u2,v2,u3,v3,u4,v4,u5,v5,u6,v6,u7,v7,u8,v8,u9,v9,u10,v10,u11,v11,u12,v12,u13,v13,u14,v14,u15,v15,u16,v16,u17,v17,u18,v18,u19,v19,u20,v20,u21,v21,u22,v22,u23,v23,u24,v24,u25,v25,u26,v26,k,s1,s2,ib1,ib2,i1,i2\n");
  // fprintf(out,"overlap,time,u1,v1,u2,v2,x1,y1,x2,y2\n");
  for (xo = 0; xo < np; xo++)
  {
    // for(xi=0;xi<np;xi++)
    //{
    // if(xo!=xi)
    //{
    fprintf(out, "r[%d][15],", xo);
    //}
    //}
  }

  fprintf(out, "\n");

  for (i = 0; i < nst; i++)
  {
    pe = pforcef(np, r, g, dt, uold, vold, b, i, cw, cb, csw);
    ke = pmovef(r, np, dt, uold, vold);

    if (i % 1500 == 0)
    {
      z1 = sqrt((r[1][0] - r[0][0]) * (r[1][0] - r[0][0]) + (r[1][1] - r[0][1]) * (r[1][1] - r[0][1]));
      z2 = sqrt((r[1][0] - r[2][0]) * (r[1][0] - r[2][0]) + (r[1][1] - r[2][1]) * (r[1][1] - r[2][1]));
      // fprintf(out,"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",pe,i*dt,r[0][2],r[0][3],r[1][2],r[1][3],r[0][0],r[0][1],r[1][0],r[1][1]);
      fprintf(out, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", pe + ke, ke, pe, 1 - r[0][15], 1 - r[1][15], 1 - r[2][15], 1 - r[3][15], 1 - r[4][15], 1 - r[5][15], 1 - r[6][15], 1 - r[7][15], 1 - r[8][15], i * dt);
      // for(xo=0;xo<np;xo++)
      //{
      //               //for(xi=0;xi<np;xi++)
      //                 // {
      //                   //if(xo!=xi)
      //                     //{
      // fprintf(out,"%lf,",r[xo][15]);
      //                     //}
      //                  //}
      //}
      fprintf(out, "\n");
    }

    np = replace(r, np, g);

    if (i % 2000 == 0)
    {
      sprintf(filename, "rsn%d.atom", cnt);
      files[cnt] = fopen(filename, "w");

      fprintf(files[cnt], "ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS ff ff ff\n0 2\n0 15\n0 2\nITEM: ATOMS id type x y z ix iy iz vx vy vz fx fy fz omegax omegay omegaz radius density\n", i, np);

      for (x = 0; x < np; x++)
      {
        fprintf(files[cnt], "%d %d %lf %lf %d %d %d %d %lf %lf %d %lf %lf %d %lf %d %d %lf %lf\n", x, 1, r[x][0], r[x][1], 0, 0, 0, 0, r[x][2], r[x][3], 0, r[x][6], r[x][7], 0, r[x][23], 0, 0, r[x][9], r[x][14]);
      }

      fclose(files[cnt]);
      cnt++;
    }

    /*else if(i%8333 == 0)
      {
       sprintf(filename, "rsn%d.atom", cnt);
       files[cnt] = fopen(filename, "w");

       fprintf(files[cnt],"ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS ff ff ff\n0 2\n0 15\n0 2\nITEM: ATOMS id type x y z ix iy iz vx vy vz fx fy fz omegax omegay omegaz radius density\n",i,np);

       for(x=0;x<np;x++)
         {
          fprintf(files[cnt],"%d %d %lf %lf %d %d %d %d %lf %lf %d %lf %lf %d %lf %d %d %lf %lf\n",x,1,r[x][0],r[x][1],0,0,0,0,r[x][2],r[x][3],0,r[x][6],r[x][7],0,r[x][23],0,0,r[x][9],r[x][14]);
         }

       fclose(files[cnt]);
       cnt++;
      }*/

    if (i % 1000 == 0)
    {
      printf("%d\n", i);
    }
  }
  fclose(out);

  printf("number of particles = %d", np);
}
