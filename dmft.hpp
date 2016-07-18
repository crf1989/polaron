#ifndef DMFT_HPP
#define DMFT_HPP 1

#include "density_of_states.hpp"
#include "green_function.hpp"

void dmft (green_function& G, density_of_states& rho,
	   double g, double w0)
{
  for (int i = 0; i < G.size(); ++i)
    G.sigma(i) = 0.;
  for (int i = 0; i < G.size(); ++i)
    G.value(i) = rho.hilbert_transform (G.omega(i), G.sigma(i));

  int counter = 0;
  double diff = 0;
  do
    {
      diff = 0;
      ++counter;

      //get new G0
      for (int i = 0; i < G.size(); ++i)
  	G[i] = 1./(1./G[i]+G.sigma(i));

      //get new sigma
      for (int i = 0; i < G.size(); ++i)
  	{
  	  complex<double> z = G.CFE (G.omega(i), g, w0, 1);
  	  diff += abs(z-G.sigma(i));
  	  G.sigma(i) = z;
  	}

      for (int i = 0; i < G.size(); ++i)
  	G.value(i) = rho.hilbert_transform (G.omega(i), G.sigma(i));

      printf ("counter = %d, diff = %g\n", counter, diff);
    } while (diff > 1e-2);
}

double kubo_greenwood (green_function& G, double t)
{
  double result = 0;
  double h = G.omega(1)-G.omega(0);	
}
  
#endif //DMFT_HPP
