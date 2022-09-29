//Function to integrate the equation of motion

double pmovef(double **r, int np, double dt, double *uold, double *vold)
{
 //Index
 int j;

 //Constants
 double dt2b2=dt*dt*0.5;

 //Kinetic energy
 double ke;
 ke=0.0;

 //Compute the kinetic energy of the system
 double rke,lke;
 rke=0.0,lke=0.0;

 //Compute the kinetic energy of the system
 for(j=0;j<np;j++)
    {
     lke+=(r[j][2]*r[j][2]+r[j][3]*r[j][3])*r[j][8]*0.5;//linear
     rke+=(r[j][23]*r[j][23])*(0.5*r[j][8]*r[j][9]*r[j][9])*0.5;//rke
    }
    ke=lke+rke;

 //update velocity and position
 for(j=0;j<np;j++)
    {
     uold[j]=r[j][2];//needed to calculate friction
     vold[j]=r[j][3];//needed to calculate friction
     r[j][22]=r[j][23];//store angular velocity to calculate friction

     r[j][11]=r[j][0];//xold - to calculate damping co-efficient
     r[j][12]=r[j][1];//yold - to calculate damping co-efficient

     //Update velocities
     r[j][2]+=0.5*(r[j][6]+r[j][4])*dt/r[j][8];
     r[j][3]+=0.5*(r[j][7]+r[j][5])*dt/r[j][8];

     //update angular velocities
     r[j][23]+=(r[j][19]+r[j][20])*dt/(r[j][8]*r[j][9]*r[j][9]);

     //Update position
     r[j][0]+=r[j][2]*dt+r[j][4]*dt2b2/r[j][8];
     r[j][1]+=r[j][3]*dt+r[j][5]*dt2b2/r[j][8];

     //Store new forces in old
     r[j][6]=r[j][4];
     r[j][7]=r[j][5];

     //store new torque in old
     r[j][19]=r[j][20];
    }

    return ke;

}
