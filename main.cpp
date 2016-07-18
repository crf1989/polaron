#include <iostream>
#include <cstdio>
#include <cstdlib>
#include "density_of_states.hpp"
#include "green_function.hpp"
#include "dmft.hpp"
#include "magnetization.hpp"
using namespace std;

double kubo_greenwood (green_function& G, double t, double u, double T)
{
  int L = 40;
  double h = 2*M_PI/L;
  double result = 0;

  for (int i = 0; i < G.size(); ++i)
    {
      if (fabs(G.sigma(i).imag()) < 1e-4)
	continue;
      if (G.omega(i) > 0)
	continue;
      for (int x = 0; x < L; ++x)
	for (int y = 0; y < L; ++y)
	  for (int z = 0; z < L; ++z)
	    {
	      double epsilon = 2*t*(cos(-M_PI+x*h) + 
				    cos(-M_PI+y*h) + cos(-M_PI+z*h));
	      double ImG = (1./(G.omega(i)-G.sigma(i)-epsilon)).imag();
	      double v = -sin(-M_PI+x*h);
	      result += exp((u-G.omega(i))/T)/T*v*v*ImG*ImG;
	    }
    }
  return result/(2*M_PI*2*M_PI*2*M_PI)*h*h*h;
}
  

int main (int argc, char* argv[])
{
  double t = 1;
  double g = 3;
  double w0 = 3;

  density_of_states rho (10000, t);
  green_function G (10000, -20, 20);
  
  dmft (G, rho, g, w0);
  FILE* fp = fopen ("dat", "w");
  for (int i = 0; i < G.size(); ++i)
    fprintf (fp, "%g\t%g\t%g\t%g\t%g\n", G.omega(i),
  	     G.sigma(i).real(), G.sigma(i).imag(),
  	     G.value(i).real(), G.value(i).imag());
  fclose (fp);
  fp = fopen ("conduct.dat", "w");

  //printf ("conductivity = %g\n", kubo_greenwood(G, t, -10, 0.1));
  for (double T = 0.1; T <= 10; T += 1e-2)
    {
      double conductivity = kubo_greenwood (G, t, -10, T);
      fprintf (fp, "%g\t%g\n", T, conductivity);
    }
}
