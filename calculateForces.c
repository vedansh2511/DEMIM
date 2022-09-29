// Function to calculate the force
#include <stdlib.h>
#include <math.h>
#include <time.h>
double pforcef(int np, double **r, double g, double dt, double *uold, double *vold, double *b, double tis, double **cw, double **cb, double *csw)
{
  srand(time(NULL));
  double mu = 0.35;
  // double rnd=((double) rand() / (RAND_MAX));
  double eff_tang_mu = mu;
  double gstar, kt, c2;

  double nu = 0.42; // polypropylene

  // Potential energy and radii
  double pe, fx, fy, r1, r2, contactforce, stsw, ftsw, ftspsw, ftdsw, stw, ftw, ftspw, ftdw, st, ft, ftsp, ftd, ftspb, stb, ftdb, ftb, deltatb, deltat, deltatw, deltatsw, vn_i, vn_j, vn_ij, wm, wmbase, S, Sbase, tangentialunitvector;

  // co-efficient of restitution and damping ratio
  double e = 0.3;
  double alpha = -(log(e)) * sqrt(5 / ((log(e)) * (log(e)) + M_PI * M_PI));
  double beta = (log(e)) / (M_PI * M_PI + (log(e)) * (log(e)));

  // Index variables
  int i, j;
  double rcut2;

  // initialize forces
  for (i = 0; i < np; i++)
  {
    r[i][4] = r[i][5] = r[i][20] = 0.0;
  }

  // Initialize potential energy
  pe = 0.0;

  // Distances
  double rx, ry, ds2, lap, lapb, theta, theta1, lapold, lap_leftwall, vit, vjt, vijt;

  // fmat characterises the resistance of particulate material against fracture in impact comminution
  double fmat = 0.944; // kg/Jm

  // Forces and potential
  double ff, ff_wall, ffz, ffzx, ffzy, fnormalwall; // ground friction and wall friction

  // spring stiffness
  double k, kb, E, Ep; // k is for inter-particle and kb is for base

  // tangential unit vector tcap
  // tcap = (tx)icap + (ty)jcap
  double tx, ty;
  double pt = dt;

  E = 1500000; // young's modulus for polypropylene

  // ffz=2*0.33*42927.8*(9*pow(10,-7))*sqrt(9*pow(10,-7));//bcoz of tolerance in diameter of polypropylene spheres

  // shear modulus=E/2(1+nu)

  // Ep = 109251528702;//E/(2*(1 - 0.291*0.291)) poisson's ratio = 0.291

  // kb = 1.33333333*E*1.098901099*(sqrt(rad));//E'=E/(1-nu^2)

  // k = 1.33333333*E*0.549450549*(sqrt(rad/2.0));//E'=E/2*(1-nu^2)

  // ffz=2*0.35*42927.8*(10^(-5))*sqrt(10^(-5));

  // calculate force between particles and base
  for (i = 0; i < np; i++)
  {
    ffz = 2 * mu * (1.33333333333333 * E * 1.214181641 * 0.5 * sqrt(r[i][9])) * (0.0000011895) * sqrt(0.0000011895); // 2*mu*Fn, Fn=(4/3)*E'*R'
    ffzx = 0, ffzy = 0;
    gstar = E / (2 * (2 - nu) * (2 + 2 * nu));

    csw[i] += (r[i][2] + r[i][3]) * dt;
    stsw = 8 * gstar * (0.0000011895 * sqrt(r1));
    deltatsw = fabs(csw[i]); // in +x direction
    ftspsw = stsw * deltatsw;
    ftdsw = 2 * sqrt(5 / 6) * beta * sqrt(stsw * (r[i][8])) * fabs(sqrt(r[i][2] * r[i][2] + r[i][3] * r[i][3]));
    ftsw = 2 * (ftdsw + ftspsw);
    if (fabs(ftsw) > ffz)
    {
      ftsw = ffz;
    }

    //      if(((double) rand() / (RAND_MAX))>0.1)
    //        {
    if (r[i][2] == 0 && r[i][3] == 0)
    {
      ffzx = 0;
      ffzy = 0;
      // printf("ffz=%lf\n",ffz);
      // printf("ffzx1=%lf\n",ffzx);
      // printf("ffzy1=%lf\n",ffzy);
    }
    else if (r[i][2] == 0 && r[i][3] != 0)
    {
      // printf("%lf\n",ffz);
      ffzx = 0;
      ffzy = ftsw * (-r[i][3] / fabs(r[i][3]));
      // printf("ffzx2=%lf\n",ffzx);
      // printf("ffzy2=%lf\n",ffzy);
    }
    else if (r[i][2] != 0 && r[i][3] == 0)
    {
      ffzx = ftsw * (-r[i][2] / fabs(r[i][2]));
      ffzy = 0;
      // printf("ffzx3=%lf\n",ffzx);
      // printf("ffzy3=%lf\n",ffzy);
    }
    else
    {
      theta1 = atan(fabs(r[i][3] / r[i][2]));
      // printf("theta1=%lf",theta1);
      ffzx = ftsw * (-r[i][2] / fabs(r[i][2])) * cos(theta1);
      ffzy = ftsw * (-r[i][3] / fabs(r[i][3])) * sin(theta1);
      // printf("ffzx4=%lf\n",ffzx);
      // printf("ffzy4=%lf\n",ffzy);
    }
    r[i][4] += ffzx;
    r[i][5] += ffzy;

    //        }
    //        else
    //          {
    //          ffzx=0;
    //          ffzy=0;
    //          }

    r1 = r[i][9];
    if (r[i][1] <= r1) // i.e. particle is being compressed at the base
    {
      cb[i][1] += r[i][2] * dt;
      kb = 1.33333333333333333333 * E * 1.214181641 * 0.5 * (sqrt(r1));
      lapb = r1 - r[i][1];
      b[i] = kb * lapb * sqrt(lapb) + sqrt(sqrt(lapb)) * (r[i][12] - r[i][1]) * alpha * sqrt(r[i][8] * kb) / dt;
      stb = 8 * gstar * (lapb * sqrt(r1));

      deltatb = fabs(cb[i][1]); // in +x direction
      ftspb = stb * deltatb;
      ftdb = 2 * sqrt(5 / 6) * beta * sqrt(stb * (r[i][8])) * fabs(r[i][2]);
      ftb = ftdb + ftspb;
      if (fabs(ftb) > mu * b[i])
      {
        ftb = mu * b[i];
      }
      // r[i][16] = 1;//interacting with base

      if (r[i][12] > r1)
      {
        if (r[i][13] == 0)
        {
          wmbase = 0.5 * r[i][3] * r[i][3];
          if (wmbase > r[i][10]) // if mass specific collision energy is greater than mass specific breakage energy
          {
            Sbase = 1 - exp(-2 * fmat * r[i][9] * (wmbase - r[i][10]));
            // if(Sbase > r[i][15])//if current probabibility of breakage is greater than old probability of breakage
            //   {
            r[i][15] = r[i][15] * (1 - Sbase); // update new probability of breakage
            //}
            if (wmbase > r[i][18])
            {
              r[i][18] = wmbase;
            }
          }
        }
        r[i][16] = 1;
      }
      if (r[i][2] != 0)
      {
        if (fabs(uold[i] - r[i][2]) > (10) ^ (-4))
        {
          if (r[i][2] > 0)
          {
            r[i][4] += -ftb;
            r[i][20] += -ftb * r[i][9]; // calculate new torque
          }

          else
          {
            r[i][4] += ftb;
            r[i][20] += ftb * r[i][9]; // calculate new torque
          }
        }
      }

      else
      {
        r[i][4] += 0.0;
        r[i][20] += 0.0;
      }
      pe += 0.4 * kb * lapb * lapb * sqrt(lapb);
    }

    else
    {
      b[i] = 0.0;
      cb[i][0] = 0.0;
      cb[i][1] = 0.0;

      if (r[i][12] < r1)
      {
        r[i][16] = 0;
      }
    }

    // if particle interacting with the left wall
    if (r[i][0] - r[i][9] < -0.125) // radius max=0.00125
    {
      cw[i][1] = 1;
      cw[i][1] += r[i][3] * dt;
      lap_leftwall = fabs(r[i][9] - (r[i][0] + 0.125)); // radius max=0.00125
      fnormalwall = kb * lap_leftwall * sqrt(lap_leftwall) + sqrt(sqrt(lap_leftwall)) * (r[i][11] - r[i][0]) * alpha * sqrt(r[i][8] * kb) / dt;
      r[i][4] += fnormalwall;
      stw = 8 * gstar * (lap_leftwall * sqrt(r1));

      deltatw = fabs(cw[i][1]); // in +y direction
      ftspw = stw * deltatw;
      ftdw = 2 * sqrt(5 / 6) * beta * sqrt(stw * (r[i][8])) * fabs(r[i][3]);
      ftw = ftdw + ftspw;
      if (fabs(ftw) > mu * fnormalwall)
      {
        ftw = mu * fnormalwall;
      }
      if (r[i][3] != 0)
      {
        if (fabs(vold[i] - r[i][3]) > (10) ^ (-5))
        {
          if (r[i][3] > 0)
          {
            r[i][5] += -ftw;
            r[i][20] += r[i][9] * ftw;
          }

          else if (r[i][3] < 0)
          {
            r[i][5] += ftw;
            r[i][20] += -r[i][9] * ftw;
          }

          //                 else
          //                 {
          //                 r[i][5]+=0;
          //                 r[i][20]+=0;
          //                 }
        }
        else
        {
          r[i][5] += 0;
          r[i][20] += 0;
        }
      }

      else
      {
        r[i][5] += 0.0;
        r[i][20] += 0.0;
      }
    }
    else
    {
      cw[i][1] = 0;
      cw[i][0] = 0;
    }

    for (j = i + 1; j < np; j++)
    {
      r2 = r[j][9]; // radius of j'th particle
      // printf("r2=%lf\n",r2);
      rx = r[j][0] - r[i][0];
      ry = r[j][1] - r[i][1];
      // printf("ry=%lf\n",ry);
      // Square of the distance
      ds2 = rx * rx + ry * ry;

      // check the cut-off distance for interaction
      rcut2 = 2 * r1 * r2 + r1 * r1 + r2 * r2;

      // Check cut-off-distance
      if (ds2 < rcut2)
      {
        lap = r1 + r2 - sqrt(ds2);                                                                                              // overlap according to the current position
        lapold = r1 + r2 - sqrt((r[j][12] - r[i][12]) * (r[j][12] - r[i][12]) + (r[j][11] - r[i][11]) * (r[j][11] - r[i][11])); // overlap according to old position

        k = 1.33333333333333333333 * E * 1.214181641 * 0.5 * (sqrt(r1 * r2 / (r1 + r2)));

        st = 8 * gstar * (sqrt(r1 * r2 * lap / (r1 + r2)));

        // printf("lapold = %lf\n",lapold);

        if (lapold < 0)
        {
          lapold = 0;

          if (r[i][13] == 0) // if particle is not broken
          {
            vn_i = (r[i][2] * rx + r[i][3] * ry) / (sqrt(ds2));  // normal velocity of particle i
            vn_j = (-r[j][2] * rx - r[j][3] * ry) / (sqrt(ds2)); // normal velocity of particle j
            vn_ij = vn_i + vn_j;                                 // normal relative velocity of i and j

            // printf("vni = %lf\n",vn_i);
            // printf("vnj = %lf\n",vn_j);
            // printf("vnij = %lf\n",vn_ij);

            wm = 0.5 * vn_ij * vn_ij; // mass specific collision energy

            // printf("%lf\n",wm);

            if (wm > r[i][10]) // if mass specific collision energy is greater than mass specific breakage energy
            {
              S = 1 - exp(-2 * fmat * (r1 * r2 / (r1 + r2)) * (wm - r[i][10])); // probability of breakage

              //                    if(r[i][15] == 1)//if current probabibility of breakage is greater than old probability of breakage
              //                    {
              //                     r[i][15] = 1 - S;//update new probability of breakage
              //                     //printf("S = %lf\n",r[i][15]);
              //                    }
              //                    else
              //                    {
              r[i][15] = r[i][15] * (1 - S); // updated probability of not breaking after nth impact
              //                    }

              if (wm > r[i][18])
              {
                r[i][18] = wm;
                r[i][17] = j;
              }
            }
          }
          r[i][30 + j] = 1; // ith particle interacting with jth particle
          r[j][30 + i] = 1; // jth particle interacting with ith particle

          // printf("zero in one %lf\n",r[0][21]);
        }

        pe += 0.4 * k * lap * lap * sqrt(lap);

        // calculate angle of collision
        if (rx == 0)
        {
          theta = 1.57075;
        }

        else if (ry == 0)
        {
          theta = 0;
        }

        else
        {
          theta = atan(fabs(ry / rx));
        }

        contactforce = (k * lap * sqrt(lap)) + sqrt(sqrt(lap)) * (lap - lapold) * alpha * (sqrt((r[i][8] * r[j][8] / (r[i][8] + r[j][8])) * k)) / dt;

        // assign forces
        if (r[j][1] > r[i][1]) // check relative y-co-ordinate
        {
          r[j][5] += (contactforce)*sin(theta);
          r[i][5] += -(contactforce)*sin(theta);

          if (r[j][0] > r[i][0]) // check relative x-co-ordinate
          {
            tx = 1 / (sqrt(1 + pow(rx / ry, 2)));
            ty = (-rx / ry) / (sqrt(1 + pow(rx / ry, 2)));

            vit = (tx * r[i][2] + ty * r[i][3]) * sqrt(r[i][2] * r[i][2] + r[i][3] * r[i][3]);
            vjt = (tx * r[j][2] + ty * r[j][3]) * sqrt(r[j][2] * r[j][2] + r[j][3] * r[j][3]);

            vijt = (vit - r[i][9] * r[i][23]) - (vjt - r[j][9] * r[j][23]);
            // r[i][30+j]+=vijt*pt;
            if (r[i][13] > r[j][13] || r[j][13] > r[i][13])
            {
              deltat = 0.07 * lap;
            }
            else
              deltat = 0.15 * lap;
            ftsp = st * deltat;
            ftd = 2 * sqrt(5 / 6) * beta * sqrt(st * (r[i][8] * r[j][8] / (r[i][8] + r[j][8]))) * vijt;
            ft = ftd + ftsp;
            if (ft > mu * contactforce)
            {
              ft = mu * contactforce;
            }

            if (fabs(sqrt(r[i][2] * r[i][2] + r[i][3] * r[i][3]) - sqrt(uold[i] * uold[i] + vold[i] * vold[i])) > 10 ^ (-5)) // since the tangential direction can change frequently, hence we use total velocity
            {
              if (vijt > 0)
              {
                r[i][4] += -(ft)*tx;      // eff_tang_mu*contactforce
                r[i][5] += -(ft)*ty;      // eff_tang_mu*contactforce
                r[i][20] += (ft)*r[i][9]; // eff_tang_mu*contactforce

                r[j][4] += (ft)*tx;       // eff_tang_mu*contactforce
                r[j][5] += (ft)*ty;       // eff_tang_mu*contactforce
                r[j][20] += (ft)*r[i][9]; // eff_tang_mu*contactforce
              }
              else if (vijt < 0)
              {
                r[i][4] += (ft)*tx;        // eff_tang_mu*contactforce
                r[i][5] += (ft)*ty;        // eff_tang_mu*contactforce
                r[i][20] += -(ft)*r[i][9]; // eff_tang_mu*contactforce

                r[j][4] += -(ft)*tx;       // eff_tang_mu*contactforce
                r[j][5] += -(ft)*ty;       // eff_tang_mu*contactforce
                r[j][20] += -(ft)*r[i][9]; // eff_tang_mu*contactforce
              }
              else
              {
                if (vit - vjt > 0)
                {
                  r[i][4] += -(ft)*tx;      // eff_tang_mu*contactforce
                  r[i][5] += -(ft)*ty;      // eff_tang_mu*contactforce
                  r[i][20] += (ft)*r[i][9]; // eff_tang_mu*contactforce

                  r[j][4] += (ft)*tx;       // eff_tang_mu*contactforce
                  r[j][5] += (ft)*ty;       // eff_tang_mu*contactforce
                  r[j][20] += (ft)*r[i][9]; // eff_tang_mu*contactforce
                }
                else if (vit - vjt < 0)
                {
                  r[i][4] += (ft)*tx;        // eff_tang_mu*contactforce
                  r[i][5] += (ft)*ty;        // eff_tang_mu*contactforce
                  r[i][20] += -(ft)*r[i][9]; // eff_tang_mu*contactforce

                  r[j][4] += -(ft)*tx;       // eff_tang_mu*contactforce
                  r[j][5] += -(ft)*ty;       // eff_tang_mu*contactforce
                  r[j][20] += -(ft)*r[i][9]; // eff_tang_mu*contactforce
                }
                else
                {
                  r[i][4] += 0;
                  r[i][5] += 0;
                  r[i][20] += 0;

                  r[j][4] += 0;
                  r[j][5] += 0;
                  r[j][20] += 0;
                }
              }
            }
            r[j][4] += (contactforce)*cos(theta);
            r[i][4] += -(contactforce)*cos(theta);
          }

          else if (r[i][0] > r[j][0]) // check relative x-co-ordinate
          {
            tx = 1 / (sqrt(1 + pow(rx / ry, 2)));
            ty = (-rx / ry) / (sqrt(1 + pow(rx / ry, 2)));

            vit = (tx * r[i][2] + ty * r[i][3]) * sqrt(r[i][2] * r[i][2] + r[i][3] * r[i][3]);
            vjt = (tx * r[j][2] + ty * r[j][3]) * sqrt(r[j][2] * r[j][2] + r[j][3] * r[j][3]);

            vijt = (vit - r[i][9] * r[i][23]) - (vjt + r[j][9] * r[j][23]);
            // r[i][30+j]+=vijt*pt;
            if (r[i][13] > r[j][13] || r[j][13] > r[i][13])
            {
              deltat = 0.07 * lap;
            }
            else
              deltat = 0.15 * lap;
            ftsp = st * deltat;
            ftd = 2 * sqrt(5 / 6) * beta * sqrt(st * (r[i][8] * r[j][8] / (r[i][8] + r[j][8]))) * vijt;
            ft = ftd + ftsp;
            if (ft > mu * contactforce)
            {
              ft = mu * contactforce;
            }

            if (fabs(sqrt(r[i][2] * r[i][2] + r[i][3] * r[i][3]) - sqrt(uold[i] * uold[i] + vold[i] * vold[i])) > 10 ^ (-5)) // since the tangential direction can change frequently, hence we use total velocity
            {
              if (vijt > 0)
              {
                r[i][4] += -(ft)*tx;      // eff_tang_mu*contactforce
                r[i][5] += -(ft)*ty;      // eff_tang_mu*contactforce
                r[i][20] += (ft)*r[i][9]; // eff_tang_mu*contactforce

                r[j][4] += (ft)*tx;       // eff_tang_mu*contactforce
                r[j][5] += (ft)*ty;       // eff_tang_mu*contactforce
                r[j][20] += (ft)*r[i][9]; // eff_tang_mu*contactforce
              }
              else if (vijt < 0)
              {
                r[i][4] += (ft)*tx;        // eff_tang_mu*contactforce
                r[i][5] += (ft)*ty;        // eff_tang_mu*contactforce
                r[i][20] += -(ft)*r[i][9]; // eff_tang_mu*contactforce

                r[j][4] += -(ft)*tx;       // eff_tang_mu*contactforce
                r[j][5] += -(ft)*ty;       // eff_tang_mu*contactforce
                r[j][20] += -(ft)*r[i][9]; // eff_tang_mu*contactforce
              }
              else
              {
                if (vit - vjt > 0)
                {
                  r[i][4] += -(ft)*tx;      // eff_tang_mu*contactforce
                  r[i][5] += -(ft)*ty;      // eff_tang_mu*contactforce
                  r[i][20] += (ft)*r[i][9]; // eff_tang_mu*contactforce

                  r[j][4] += (ft)*tx;       // eff_tang_mu*contactforce
                  r[j][5] += (ft)*ty;       // eff_tang_mu*contactforce
                  r[j][20] += (ft)*r[i][9]; // eff_tang_mu*contactforce
                }
                else if (vit - vjt < 0)
                {
                  r[i][4] += (ft)*tx;        // eff_tang_mu*contactforce
                  r[i][5] += (ft)*ty;        // eff_tang_mu*contactforce
                  r[i][20] += -(ft)*r[i][9]; // eff_tang_mu*contactforce

                  r[j][4] += -(ft)*tx;       // eff_tang_mu*contactforce
                  r[j][5] += -(ft)*ty;       // eff_tang_mu*contactforce
                  r[j][20] += -(ft)*r[i][9]; // eff_tang_mu*contactforce
                }
                else
                {
                  r[i][4] += 0;
                  r[i][5] += 0;
                  r[i][20] += 0;

                  r[j][4] += 0;
                  r[j][5] += 0;
                  r[j][20] += 0;
                }
              }
            }
            r[i][4] += (contactforce)*cos(theta);
            r[j][4] += -(contactforce)*cos(theta);
          }

          else if (r[i][0] = r[j][0]) // check relative x-co-ordinate
          {
            // check relative tangential velocity in x-direction
            if ((r[i][2] - r[i][9] * r[i][23] - (r[j][2] + r[j][9] * r[j][23])) != 0)
            {
              tangentialunitvector = (r[i][2] - r[i][9] * r[i][23] - (r[j][2] + r[j][9] * r[j][23])) / fabs((r[i][2] - r[i][9] * r[i][23] - (r[j][2] + r[j][9] * r[j][23])));
            }
            else
            {
              if ((r[i][2] - r[j][2]) != 0)
              {
                tangentialunitvector = (r[i][2] - r[j][2]) / fabs(r[i][2] - r[j][2]);
              }
              else
                tangentialunitvector = 0.0;
            }
            // apply a tangential force accordingly(mu*Fn)
            // apply an appropriate torque
            r[i][4] += -tangentialunitvector * (mu * contactforce); // eff_tang_mu*contactforce
            r[j][4] += tangentialunitvector * (mu * contactforce);

            r[i][20] += -tangentialunitvector * (mu * contactforce) * r[i][9]; // eff_tang_mu*contactforce
            r[j][20] += tangentialunitvector * (mu * contactforce) * r[j][9];  // eff_tang_mu*contactforce
          }
        }

        else if (r[j][1] < r[i][1]) // check relative y-co-ordinate
        {
          r[i][5] += (contactforce)*sin(theta);
          r[j][5] += -(contactforce)*sin(theta);

          if (r[j][0] > r[i][0]) // check relative x-co-ordinate
          {
            tx = 1 / (sqrt(1 + pow(rx / ry, 2)));
            ty = (-rx / ry) / (sqrt(1 + pow(rx / ry, 2)));

            vit = (tx * r[i][2] + ty * r[i][3]) * sqrt(r[i][2] * r[i][2] + r[i][3] * r[i][3]);
            vjt = (tx * r[j][2] + ty * r[j][3]) * sqrt(r[j][2] * r[j][2] + r[j][3] * r[j][3]);

            vijt = (vit + r[i][9] * r[i][23]) - (vjt - r[j][9] * r[j][23]);
            // r[i][30+j]+=vijt*pt;
            if (r[i][13] > r[j][13] || r[j][13] > r[i][13])
            {
              deltat = 0.07 * lap;
            }
            else
              deltat = 0.15 * lap;
            ftsp = st * deltat;
            ftd = 2 * sqrt(5 / 6) * beta * sqrt(st * (r[i][8] * r[j][8] / (r[i][8] + r[j][8]))) * vijt;
            ft = ftd + ftsp;
            if (ft > mu * contactforce)
            {
              ft = mu * contactforce;
            }

            if (fabs(sqrt(r[i][2] * r[i][2] + r[i][3] * r[i][3]) - sqrt(uold[i] * uold[i] + vold[i] * vold[i])) > 10 ^ (-5)) // since the tangential direction can change frequently, hence we use total velocity
            {
              if (vijt > 0)
              {
                r[i][4] += -(ft)*tx;       // eff_tang_mu*contactforce
                r[i][5] += -(ft)*ty;       // eff_tang_mu*contactforce
                r[i][20] += -(ft)*r[i][9]; // eff_tang_mu*contactforce

                r[j][4] += (ft)*tx;        // eff_tang_mu*contactforce
                r[j][5] += (ft)*ty;        // eff_tang_mu*contactforce
                r[j][20] += -(ft)*r[i][9]; // eff_tang_mu*contactforce
              }
              else if (vijt < 0)
              {
                r[i][4] += (ft)*tx;       // eff_tang_mu*contactforce
                r[i][5] += (ft)*ty;       // eff_tang_mu*contactforce
                r[i][20] += (ft)*r[i][9]; // eff_tang_mu*contactforce

                r[j][4] += -(ft)*tx;      // eff_tang_mu*contactforce
                r[j][5] += -(ft)*ty;      // eff_tang_mu*contactforce
                r[j][20] += (ft)*r[i][9]; // eff_tang_mu*contactforce
              }
              else
              {
                if (vit - vjt > 0)
                {
                  r[i][4] += -(ft)*tx;       // eff_tang_mu*contactforce
                  r[i][5] += -(ft)*ty;       // eff_tang_mu*contactforce
                  r[i][20] += -(ft)*r[i][9]; // eff_tang_mu*contactforce

                  r[j][4] += (ft)*tx;        // eff_tang_mu*contactforce
                  r[j][5] += (ft)*ty;        // eff_tang_mu*contactforce
                  r[j][20] += -(ft)*r[i][9]; // eff_tang_mu*contactforce
                }
                else if (vit - vjt < 0)
                {
                  r[i][4] += (ft)*tx;       // eff_tang_mu*contactforce
                  r[i][5] += (ft)*ty;       // eff_tang_mu*contactforce
                  r[i][20] += (ft)*r[i][9]; // eff_tang_mu*contactforce

                  r[j][4] += -(ft)*tx;      // eff_tang_mu*contactforce
                  r[j][5] += -(ft)*ty;      // eff_tang_mu*contactforce
                  r[j][20] += (ft)*r[i][9]; // eff_tang_mu*contactforce
                }
                else
                {
                  r[i][4] += 0;
                  r[i][5] += 0;
                  r[i][20] += 0;

                  r[j][4] += 0;
                  r[j][5] += 0;
                  r[j][20] += 0;
                }
              }
            }
            r[j][4] += (contactforce)*cos(theta);
            r[i][4] += -(contactforce)*cos(theta);
          }

          else if (r[i][0] > r[j][0]) // check relative x-co-ordinate
          {
            tx = 1 / (sqrt(1 + pow(rx / ry, 2)));
            ty = (-rx / ry) / (sqrt(1 + pow(rx / ry, 2)));

            vit = (tx * r[i][2] + ty * r[i][3]) * sqrt(r[i][2] * r[i][2] + r[i][3] * r[i][3]);
            vjt = (tx * r[j][2] + ty * r[j][3]) * sqrt(r[j][2] * r[j][2] + r[j][3] * r[j][3]);

            vijt = (vit + r[i][9] * r[i][23]) - (vjt - r[j][9] * r[j][23]);
            // r[i][30+j]+=vijt*pt;
            if (r[i][13] > r[j][13] || r[j][13] > r[i][13])
            {
              deltat = 0.07 * lap;
            }
            else
              deltat = 0.15 * lap;
            ftsp = st * deltat;
            ftd = 2 * sqrt(5 / 6) * beta * sqrt(st * (r[i][8] * r[j][8] / (r[i][8] + r[j][8]))) * vijt;
            ft = ftd + ftsp;
            if (ft > mu * contactforce)
            {
              ft = mu * contactforce;
            }

            if (fabs(sqrt(r[i][2] * r[i][2] + r[i][3] * r[i][3]) - sqrt(uold[i] * uold[i] + vold[i] * vold[i])) > 10 ^ (-5)) // since the tangential direction can change frequently, hence we use total velocity
            {
              if (vijt > 0)
              {
                r[i][4] += -(ft)*tx;       // eff_tang_mu*contactforce
                r[i][5] += -(ft)*ty;       // eff_tang_mu*contactforce
                r[i][20] += -(ft)*r[i][9]; // eff_tang_mu*contactforce

                r[j][4] += (ft)*tx;        // eff_tang_mu*contactforce
                r[j][5] += (ft)*ty;        // eff_tang_mu*contactforce
                r[j][20] += -(ft)*r[i][9]; // eff_tang_mu*contactforce
              }
              else if (vijt < 0)
              {
                r[i][4] += (ft)*tx;       // eff_tang_mu*contactforce
                r[i][5] += (ft)*ty;       // eff_tang_mu*contactforce
                r[i][20] += (ft)*r[i][9]; // eff_tang_mu*contactforce

                r[j][4] += -(ft)*tx;      // eff_tang_mu*contactforce
                r[j][5] += -(ft)*ty;      // eff_tang_mu*contactforce
                r[j][20] += (ft)*r[i][9]; // eff_tang_mu*contactforce
              }
              else
              {
                if (vit - vjt > 0)
                {
                  r[i][4] += -(ft)*tx;       // eff_tang_mu*contactforce
                  r[i][5] += -(ft)*ty;       // eff_tang_mu*contactforce
                  r[i][20] += -(ft)*r[i][9]; // eff_tang_mu*contactforce

                  r[j][4] += (ft)*tx;        // eff_tang_mu*contactforce
                  r[j][5] += (ft)*ty;        // eff_tang_mu*contactforce
                  r[j][20] += -(ft)*r[i][9]; // eff_tang_mu*contactforce
                }
                else if (vit - vjt < 0)
                {
                  r[i][4] += (ft)*tx;       // eff_tang_mu*contactforce
                  r[i][5] += (ft)*ty;       // eff_tang_mu*contactforce
                  r[i][20] += (ft)*r[i][9]; // eff_tang_mu*contactforce

                  r[j][4] += -(ft)*tx;      // eff_tang_mu*contactforce
                  r[j][5] += -(ft)*ty;      // eff_tang_mu*contactforce
                  r[j][20] += (ft)*r[i][9]; // eff_tang_mu*contactforce
                }
                else
                {
                  r[i][4] += 0;
                  r[i][5] += 0;
                  r[i][20] += 0;

                  r[j][4] += 0;
                  r[j][5] += 0;
                  r[j][20] += 0;
                }
              }
            }
            r[i][4] += (contactforce)*cos(theta);
            r[j][4] += -(contactforce)*cos(theta);
          }

          else if (r[i][0] = r[j][0]) // check relative x-co-ordinate
          {
            // check relative tangential velocity in x-direction
            if ((r[i][2] + r[i][9] * r[i][23] - (r[j][2] - r[j][9] * r[j][23])) != 0)
            {
              tangentialunitvector = (r[i][2] + r[i][9] * r[i][23] - (r[j][2] - r[j][9] * r[j][23])) / fabs((r[i][2] + r[i][9] * r[i][23] - (r[j][2] - r[j][9] * r[j][23])));
            }
            else
            {
              if ((r[i][2] - r[j][2]) != 0)
              {
                tangentialunitvector = (r[i][2] - r[j][2]) / fabs(r[i][2] - r[j][2]);
              }
              else
                tangentialunitvector = 0.0;
            }
            // apply a tangential force accordingly(mu*Fn)
            // apply an appropriate torque
            r[i][4] += -tangentialunitvector * (mu * contactforce); // eff_tang_mu*contactforce
            r[j][4] += tangentialunitvector * (mu * contactforce);  // eff_tang_mu*contactforce

            r[i][20] += -tangentialunitvector * (mu * contactforce) * r[i][9]; // eff_tang_mu*contactforce
            r[j][20] += tangentialunitvector * (mu * contactforce) * r[j][9];  // eff_tang_mu*contactforce
          }
        }

        else // same y co-ordinates
        {

          /*vijt = (vit + r[i][9]*r[i][23]) - (vjt - r[j][9]*r[j][23]);
          deltat = vijt*6*dt;
          ftsp = st*vijt*deltat;
          ftd = 2*sqrt(5/6)*beta*sqrt(st*(r[i][8]*r[j][8]/(r[i][8]+r[j][8])))*vijt;
          if(ft>mu*contactforce)
          {
          ft=mu*contactforce;
          }*/

          if (r[j][0] > r[i][0]) // check relative x-co-ordinate
          {
            // check relative tangential velocity in x-direction
            if ((r[i][3] + r[i][9] * r[i][23] - (r[j][3] - r[j][9] * r[j][23])) != 0)
            {
              tangentialunitvector = ((r[i][3] + r[i][9] * r[i][23] - (r[j][3] - r[j][9] * r[j][23]))) / fabs((r[i][3] + r[i][9] * r[i][23] - (r[j][3] - r[j][9] * r[j][23])));
            }
            else
            {
              if ((r[i][3] - r[j][3]) != 0)
              {
                tangentialunitvector = (r[i][3] - r[j][3]) / fabs(r[i][3] - r[j][3]);
              }
              else
                tangentialunitvector = 0.0;
            }
            // apply a tangential force accordingly(mu*Fn)
            r[i][5] += -tangentialunitvector * (mu * contactforce); // eff_tang_mu*contactforce
            r[j][5] += tangentialunitvector * (mu * contactforce);  // eff_tang_mu*contactforce
            // apply an appropriate torque
            r[i][20] += -tangentialunitvector * (mu * contactforce) * r[i][9]; // eff_tang_mu*contactforce
            r[j][20] += tangentialunitvector * (mu * contactforce) * r[j][9];  // eff_tang_mu*contactforce
            // contactforce normal
            r[j][4] += (contactforce);
            r[i][4] += -(contactforce);
          }

          else if (r[i][0] > r[j][0]) // check relative x-co-ordinate
          {
            // check relative tangential velocity in x-direction
            if ((r[i][3] - r[i][9] * r[i][23] - (r[j][3] + r[j][9] * r[j][23])) != 0)
            {
              tangentialunitvector = ((r[i][3] - r[i][9] * r[i][23] - (r[j][3] + r[j][9] * r[j][23]))) / fabs((r[i][3] - r[i][9] * r[i][23] - (r[j][3] + r[j][9] * r[j][23])));
            }
            else
            {
              if ((r[i][3] - r[j][3]) != 0)
              {
                tangentialunitvector = (r[i][3] - r[j][3]) / fabs(r[i][3] - r[j][3]);
              }
              else
                tangentialunitvector = 0.0;
            }
            // apply a tangential force accordingly(mu*Fn)
            r[i][5] += -tangentialunitvector * (mu * contactforce); // eff_tang_mu*contactforce
            r[j][5] += tangentialunitvector * (mu * contactforce);  // eff_tang_mu*contactforce
            // apply an appropriate torque
            r[i][20] += -tangentialunitvector * (mu * contactforce) * r[i][9]; // eff_tang_mu*contactforce
            r[j][20] += tangentialunitvector * (mu * contactforce) * r[j][9];  // eff_tang_mu*contactforce
            // contactforce normal
            r[i][4] += (contactforce);
            r[j][4] += -(contactforce);
          }
        }
      }
      else // lap=0;
      {
        r[i][30 + j] = 0;                                                                                                       // interaction has stopped
        r[j][30 + i] = 0;                                                                                                       // interaction has stopped
        lapold = r1 + r2 - sqrt((r[j][12] - r[i][12]) * (r[j][12] - r[i][12]) + (r[j][11] - r[i][11]) * (r[j][11] - r[i][11])); // overlap according to old position
        /*if(lapold>0)//i.e. if the ith particle is not interacting with jth particle, but they were interacting in the previous time-step
        {
         r[i][30+j] = 0;//interaction has stopped
         r[j][30+i] = 0;//interaction has stopped
         //printf("zero in one%lf\n",r[0][21]);
        }*/
      }
    }

    r[i][5] += b[i] - r[i][8] * g;

    pe += r[i][8] * g * r[i][1];
  }

  return pe;
}
