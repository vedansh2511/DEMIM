#include <stdlib.h>
#include <math.h>
#include <time.h>
int replace(double **r, int np, double g)
{
   double p_theta, alpha, dsq = 0, rsq = 0, s, ss, x, y, kebefore, keafter, rnd; // angle of line perpendicular to velocity of parent particle
   // s=sum of new u-velocities of progeny after dispersion  ss=sum of squares of new u-velocities of progeny after dispersion
   int flag = 0;
   int b = 0, k, c, iter; // b=number of particles that have broken in this iteration  k=index of particle with which wmax is achieved  c=variable to break loop of dispersion
   // srand(time(NULL));
   for (int i = 0; i < np; i++)
   {
      // srand(time(NULL));
      if ((1 - r[i][15]) > 0 && r[i][16] == 0 && r[i][13] == 0) // if probability of breakage is greater than zero and if not interacting with the base and if not already broken
      {
         for (int j = 0; j < np; j++) // check if not interacting with any other particle
         {
            if (i != j && r[i][30 + j] == 1)
            {
               flag = 1; // particle interacting with atleast one other particle
               break;
            }
         }
         // printf("%d",flag);
         if (flag == 0)
         {
            // srand(time(NULL));
            rnd = ((double)rand() / (RAND_MAX));
            printf("%lf\n", rnd);
            if ((1 - r[i][15]) > rnd)
            {
               b++;
               k = r[i][17];
               c = 1, alpha = 0.9;
               kebefore = 0, keafter = 0;

               r[np + b - 1][9] = 0.5 * r[i][9];                                                                               // radius
               r[np + b - 1][8] = 0.5 * r[i][8];                                                                               // mass
               r[np + b - 1][14] = 3 * r[np + b - 1][8] / (4 * M_PI * r[np + b - 1][9] * r[np + b - 1][9] * r[np + b - 1][9]); // density

               r[i][9] = 0.5 * r[i][9];                                           // radius
               r[i][8] = 0.5 * r[i][8];                                           // mass
               r[i][14] = 3 * r[i][8] / (4 * M_PI * r[i][9] * r[i][9] * r[i][9]); // density

               r[np + b - 1][0] = r[i][0] + r[np + b - 1][9];
               r[np + b - 1][1] = r[i][1];
               r[np + b - 1][2] = r[i][2];
               r[np + b - 1][3] = r[i][3];
               r[np + b - 1][4] = 0;
               r[np + b - 1][5] = -r[np + b - 1][8] * g;
               r[np + b - 1][6] = 0;
               r[np + b - 1][7] = r[np + b - 1][5];
               r[np + b - 1][10] = 0.000001 / (2 * r[np + b - 1][9]);
               r[np + b - 1][11] = r[np + b - 1][0];
               r[np + b - 1][12] = r[np + b - 1][1];
               r[np + b - 1][13] = 1; // broken
               r[np + b - 1][15] = 1; // probability of not breaking
               r[np + b - 1][16] = 0; // not interacting with the base

               r[i][0] = r[i][0] - r[np + b - 1][9];
               r[i][1] = r[i][1];
               r[i][2] = r[i][2];
               r[i][3] = r[i][3];
               r[i][4] = 0;
               r[i][5] = -r[i][8] * g;
               r[i][6] = 0;
               r[i][7] = r[i][5];
               r[i][10] = 0.000001 / (2 * r[i][9]);
               r[i][11] = r[i][0];
               r[i][12] = r[i][1];
               r[i][13] = 1;
               r[i][15] = 1;
               r[i][16] = 0;

               // initialize interactions
               for (iter = 30; iter < 3000; iter++)
               {
                  r[np + b - 1][iter] = 0;
                  r[i][iter] = 0;
               }

               kebefore = 0.5 * r[np + b - 1][8] * r[np + b - 1][2] * r[np + b - 1][2] + 0.5 * r[i][8] * r[i][2] * r[i][2] + 0.5 * r[k][8] * r[k][2] * r[k][2];

               // now we come to dispersion
               /*do
               {
                s = (r[k][8]*r[k][2]*(1-alpha))/r[i][8] + 2*r[i][2];
                dsq = 0.5*s*s;
                ss = (r[k][8]*r[k][2]*r[k][2]*(1 - alpha*alpha))/r[i][8] + 2*r[i][2]*r[i][2];

                printf("dsq-ss = %lf\n",dsq-ss);
                if(dsq<=ss)
                {
                 x = -2*s;
                 y = s*s - ss;

                 if(x*x>4*2*y)
                 {
                  c=0;
                  printf("c=%d",c);
                 }
                }
                alpha-=0.1;
                printf("inside dispersion alpha = %lf\n",alpha);
               }
               while(c==1 && alpha>=0);*/

               /*if(c==0)
               {
               r[np+b-1][2] = 0.25*((-x + sqrt(x*x - 4*2*y)));
               r[i][2] = 0.25*((-x - sqrt(x*x - 4*2*y)));
               r[k][2] = (alpha+0.1)*r[k][2];

               keafter = 0.5*r[np+b-1][8]*r[np+b-1][2]*r[np+b-1][2] + 0.5*r[i][8]*r[i][2]*r[i][2] + 0.5*r[k][8]*r[k][2]*r[k][2];

               printf("change in ke after dispersion = %lf",keafter-kebefore);
               }*/
            }
            // else r[i][15]=0;
         }
      }
   }

   np = np + b;
   return np;
}
