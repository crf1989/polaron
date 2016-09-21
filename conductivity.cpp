#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include "magnetization.hpp"
using namespace std;

int main ()
{
  double T=0.01;
  double u=-7;
  
  double a,b,c,d,e,f;
  
  vector<double> omega (0);
  vector<double> ksum (0);

  FILE* fp = fopen ("dat", "r");
  while (fscanf (fp, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", &a, &b, &c, &d, &e, &f) != EOF)
    {
      omega.push_back (a);
      ksum.push_back (f);
    }
  fclose (fp);

  double result = 0;
  double h = omega[1]-omega[0];
  for (int i = 0; i < omega.size(); ++i)
    result += (u-omega[i])/T*exp((u-omega[i])/T)
  
}
  
  
  
