#ifndef DENSITY_OF_STATES
#define DENSITY_OF_STATES 1

#include <vector>
#include <complex>
#include <cmath>
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_integration.h>
using namespace std;

class density_of_states
{
public:
  density_of_states (int, double);
  int size ();
  double begin ();
  double end ();
  double& omega (int);
  double& value (int);
  double interpolate (double);
  double& operator[] (int);
  double operator() (double);
  complex<double> hilbert_transform (double, complex<double>);

private:
  vector<double> _omega;
  vector<double> _value;
};

double K (double x, void* params)
{
  double w = fabs(((double*)params)[0]);
  double t = fabs(((double*)params)[1]);
  double k = (w-2*t*cos(x))/(4*t);  
  if (fabs (k) < 1e-6)
    return 0;
  else
    return gsl_sf_ellint_Kcomp (sqrt(1-k*k), GSL_PREC_DOUBLE);
}

density_of_states::density_of_states (int size, double t)
{
  _omega.resize (size);
  _value.resize (size);
  for (int i = 0; i < size; ++i)
    _omega[i] = -6*t + 12.*i*t/(size-1);

  double epsabs = 1e-6;
  double epsrel = 1e-6;
  double abserr;
  size_t nevals;
  double params[2];
  gsl_function F;
  F.function = K;
  F.params = params;
  gsl_integration_cquad_workspace* workspace = 
    gsl_integration_cquad_workspace_alloc (size);

  for (int i = 0; i < size; ++i)
    {
      params[0] = _omega[i];
      params[1] = t;

      if (2*t <= fabs (_omega[i]) && fabs(_omega[i]) <= 6*t)
        gsl_integration_cquad (&F, 0, acos((fabs(_omega[i])-4*t)/(2*t)),
			       epsabs, epsrel, workspace, 
			       &(_value[i]), &abserr, &nevals);
      else if (0 <= fabs(_omega[i]) && fabs(_omega[i]) < 2*t)
        gsl_integration_cquad (&F, 0, M_PI,
			       epsabs, epsrel, workspace, 
			       &(_value[i]), &abserr, &nevals);
      _value[i] /= (2*t*M_PI*M_PI*M_PI);      
    }
  gsl_integration_cquad_workspace_free (workspace);
}

int density_of_states::size ()
{
  return _omega.size ();
}

double density_of_states::begin ()
{
  return _omega.front ();
}

double density_of_states::end ()
{
  return _omega.back ();
}

double& density_of_states::omega (int i)
{
  return _omega[i];
}

double& density_of_states::value (int i)
{
  return _value[i];
}

double density_of_states::interpolate (double w)
{
  if (w < _omega.front() || w > _omega.back())
    return 0;
  else
    {
      int i = 0;
      while (!(_omega[i] <= w && w <= _omega[i+1]))
	++i;
      double h = _omega[1]-_omega[0];
      return _value[i] + 
	(_value[i+1]-_value[i])*(w-_omega[i])/h;
    }
}

double& density_of_states::operator[] (int i)
{
  return _value[i];
}

double density_of_states::operator() (double w)
{
  return interpolate (w);
}

complex<double>
density_of_states::hilbert_transform (double w, complex<double> sigma)
{
  double h = _omega[1]-_omega[0];
  sigma -= complex<double> (0,1e-3);
  complex<double> result = (_value.front()/(w-_omega.front()-sigma)+
			    _value.back()/(w-_omega.back()-sigma))/2.;
  for (int i = 1; i < size()-1; ++i)
    result += _value[i]/(w-_omega[i]-sigma);
  
  return result*h;
}
    


#endif //DENSITY_OF_STATES
