#include <iostream>
#include <cstdio>
#include <cstdlib>
#include "density_of_states.hpp"
#include "green_function.hpp"
#include "dmft.hpp"
using namespace std;

int main (int argc, char* argv[])
{
  double t = 1;
  double g = 2;
  double w0 = 2;

  density_of_states rho (10000, t);
  green_function G (10000, -20, 20);
  
  dmft (G, rho, g, w0);
  FILE* fp = fopen ("dat", "w");
  for (int i = 0; i < G.size(); ++i)
    {
        int L = 20;
	double h = 2*M_PI/L;
	double result = 0;

	if (fabs(G.sigma(i).imag()) < 0.1)
	  goto print;
	
	for (int x = 0; x < L; ++x)
	  for (int y = 0; y < L; ++y)
	    for (int z = 0; z < L; ++z)
	      {
		double epsilon = 2*t*(cos(-M_PI+x*h) + 
				      cos(-M_PI+y*h) + cos(-M_PI+z*h));
		double ImG = (1./(G.omega(i)-G.sigma(i)-epsilon)).imag();
		double v = -sin(-M_PI+x*h);
		result += v*v*ImG*ImG;
	      }
	result *= h*h*h;
	result /= (2*M_PI*2*M_PI*2*M_PI);
	
    print:
	fprintf (fp, "%g\t%g\t%g\t%g\t%g\t%g\n", G.omega(i),
		 G.sigma(i).real(), G.sigma(i).imag(),
		 G.value(i).real(), G.value(i).imag(),
		 result);
    }
  fclose (fp);
}
