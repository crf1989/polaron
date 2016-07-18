#ifndef MAGNETIZATION_HPP
#define MAGNETIZATION_HPP 1

#include <cmath>
#include <cstdlib>
#include <cstdio>
using namespace std;

double  brillouin_function(double j, double y)
{
  if(fabs(y) < 1.0e-5)
    return y*(j+1)/(3*j); 
  else
    return 1./tanh(y*(2*j+1)/(2*j)) * (2*j+1)/(2*j) - 1./tanh(y/(2*j))/(2*j);
}

double magnetization (double T, double B)
{
  const double J1 = -0.037;
  const double J2 = 0.069;
  const double g = 2;
  const double u = 0.6717139;
  const double j = 3.5;

  double sax, saz, sbx, sbz;
  double nsax, nsaz, nsbx, nsbz;
  double Hax, Haz, Hbx, Hbz;
  double diff = 0;

  sax = -0.2; saz = 2.5;
  sbx = +0.2; sbz = -2.54;

  int counter = 0;  
  do 
    {
      ++counter;
      Hax = 6*J1*sbx + 12*J2*sax;
      Haz = 6*J1*sbz + 12*J2*saz + g*u*B;
      Hbx = 6*J1*sax + 12*J2*sbx;
      Hbz = 6*J1*saz + 12*J2*sbz + g*u*B;
      
      nsax = j*Hax*brillouin_function(j, j*sqrt(Hax*Hax+Haz*Haz)/T) 
	/ sqrt(Hax*Hax+Haz*Haz);
      nsaz = j*Haz*brillouin_function(j, j*sqrt(Hax*Hax+Haz*Haz)/T) 
	/ sqrt(Hax*Hax+Haz*Haz);
      nsbx = j*Hbx*brillouin_function(j, j*sqrt(Hbx*Hbx+Hbz*Hbz)/T) 
	/ sqrt(Hbx*Hbx+Hbz*Hbz);
      nsbz = j*Hbz*brillouin_function(j, j*sqrt(Hbx*Hbx+Hbz*Hbz)/T) 
	/ sqrt(Hbx*Hbx+Hbz*Hbz);

      diff = fabs(nsax-sax) + fabs(nsaz-saz) + fabs(nsbx-sbx) + fabs(nsbz-sbz);

      sax = nsax; saz = nsaz;
      sbx = nsbx; sbz = nsbz;
    } while (diff > 1e-6);
  
  return (saz+sbz)/2.;
}

#endif //MAGNETIZATION_HPP
