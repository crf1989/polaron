#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include "magnetization.hpp"
using namespace std;

int main ()
{
  double T;
  double u;
  double a,b,c,d,e,f;
  
  vector<double> omega (0);
  vector<double> ksum (0);

  FILE* fp = fopen ("dat", "r");
  while (fscanf (fp, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", &a, &b, &c, &d, &e, &f) != EOF)
    {
      omega.push_back (a);
      ksum.push_back (f);
    }
}
  
  
  