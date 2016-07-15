#ifndef GREEN_FUNCTION_HPP
#define GREEN_FUNCTION_HPP 1

#include <vector>
#include <complex>
#include <cmath>
#include "density_of_states.hpp"
using namespace std;

class green_function
{
public:
  green_function (int, double, double);
  int size ();
  double& omega (int);
  complex<double>& sigma (int);
  complex<double>& value (int);
  complex<double> interpolate (double);
  complex<double>& operator[] (int);
  complex<double> operator() (double);
  complex<double> CFE (double, double, double, int);
private:
  vector<double> _omega;
  vector<complex<double> > _sigma;
  vector<complex<double> > _value;
};

green_function::green_function (int size, double a, double b)
{
  _omega.resize (size);
  _sigma.resize (size);
  _value.resize (size);
  for (int i = 0; i < size; ++i)
    _omega[i] = a + i*(b-a)/(size-1);
}
  
int green_function::size ()
{
  return _omega.size ();
}

double& green_function::omega (int i)
{
  return _omega[i];
}

complex<double>& green_function::sigma (int i)
{
  return _sigma[i];
}

complex<double>& green_function::value (int i)
{
  return _value[i];
}

complex<double> green_function::interpolate (double w)
{
  if (w < _omega.front() || w > _omega.back ())
    return 1/w;
  else
    {
      int i = 0;
      while (!(_omega[i] <= w && w <= _omega[i+1]))
	++i;
      double h = _omega[1]-_omega[0];
      return _value[i] + (_value[i+1]-_value[i])*(w-_omega[i])/h;
    }
}

complex<double>& green_function::operator[] (int i)
{
  return _value[i];
}

complex<double> green_function::operator() (double w)
{
  return interpolate (w);
}

complex<double> green_function::CFE (double w, double g, double w0, int n)
{
  if (n >= 10)
    return n*g*g*interpolate(w-n*w0);
  else
    return n*g*g/(1./interpolate(w-n*w0)-CFE(w,g,w0,n+1));
}

#endif //GREEN_FUNCTION_HPP
