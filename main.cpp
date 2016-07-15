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
  double g = 10;
  double w0 = 10;

  density_of_states rho (10000, t);
  green_function G (10000, -20, 20);
  
  dmft (G, rho, g, w0);
  FILE* fp = fopen ("dat", "w");
  for (int i = 0; i < G.size(); ++i)
    fprintf (fp, "%g\t%g\t%g\t%g\t%g\n", G.omega(i),
  	     G.sigma(i).real(), G.sigma(i).imag(),
  	     G.value(i).real(), G.value(i).imag());
}
