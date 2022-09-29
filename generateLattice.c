//Program to generate a square lattice
//
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<unistd.h>

double plfvalidation(int np,double **r,double h, double w,int lp)
{
 int i,j,npl=640;
 double a=6.0;
 double pe=0.0;

 //Find the nearest square lattice
 i=0;
// do{
//    j=a*a*i*i;
//    i++;
//    }while(j<np);
//
// //Total number of lattice points
// int npl;
// npl=j;
//
// //lattice parameter
// lp=i-1;

 //Lattice coordinates
 double **l;

 //Allocate memory
 l=(double**)calloc(np,sizeof(double*));
 for(i=0;i<np;i++)l[i]=(double*)calloc(2,sizeof(double));

 l[0][0]=0.0;//x
 l[0][1]=h;//y

 int li=1;

// for(i=1;i<npl;i++)
//       {
//        l[i][0]=l[i-1][0]+2.8*rad;//[i][9];                           //increment in x
//        l[i][1]=l[i-1][1];
//
//        //If one line is filled start, start a new line
//        if(i%lp==0)
//                {
//                li++;
//                if(li%2==0){
//                            l[i][0]=l[i-lp][0]+1.5*rad;//[i][9];                                        //increment in x
//                            l[i][1]=l[i-lp][1]+2.1*rad;//[i][9];                    //increment in y
//                           }
//                else{
//                     l[i][0]=l[0][0];
//                     l[i][1]=l[i-lp][1]+2.1*rad;
//                    }
//                }
//       }

 for(i=1;i<np;i++)//number of layer
       {
        if(i<8)
        {
        l[i][0]=l[i-1][0]+r[i-1][9]+2*r[i][9]+0.00022;//[i][9];                           //increment in x
        l[i][1]=l[i-1][1];
        }
        else
        {
        l[i][0]=l[i-1][0]+r[i-1][9]+2*r[i][9]+0.00022;//[i][9];                           //increment in x
        l[i][1]=l[i-8][1]+r[i-8][9]+2*r[i][9];//+0.00001;
        }

        //If one line is filled start, start a new line
        if(i%8==0)
                {
                li++;
                if(li%2==0){
                            l[i][0]=l[i-8][0]+1.55*r[i][9]+0.00022;//[i][9];       +0.525*r[i][9]                                 //increment in x
                            l[i][1]=l[i-8][1]+r[i-8][9]+2*r[i][9];//+0.00001;//[i][9];                    //increment in y
                           }
                else{
                     l[i][0]=l[0][0];
                     l[i][1]=l[i-8][1]+r[i-8][9]+2*r[i][9];//+0.00001;
                    }
                }
       }


 FILE *out=fopen("trylat.dat","w");     //do we place fopen and fclose in such a way that any results generated between the 2 statements get written in the .dat file?
   for(i=0;i<np;i++)
      {
      r[i][0]=l[i][0];                      //why multiplying with scaling factor?
      r[i][1]=l[i][1];

      pe+=w*r[i][1];
      fprintf(out,"create_atoms   %d single  %lf %lf %lf units box\n",1,r[i][0],0.0,r[i][1]);   //fprintf?
      }

   for(i=1;i<=np;i++)
      {
      fprintf(out,"set  atom %d diameter %lf density %lf vx %lf vy %lf vz %lf\n",i,2*r[i][9],920.0,0.0,0.0,0.0);   //fprintf?
      }

   fprintf(out,"%lf\n",pe);
   fclose(out);
   return pe;
 }
