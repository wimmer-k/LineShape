#ifndef __RANGE_H 
#define __RANGE_H

#include <iostream>
#include <fstream>
#include "TSpline.h"
using namespace std;

class Range
{
public:
  Range(double mass, char* filename);
  double GetRange(double beta);
  double GetBetaAfter(double beta, double thick);

  void CreateSplines(char* filename);


  double beta2gamma(double beta);
  double gamma2beta(double gamma);
  double beta2E(double beta);
  double E2beta(double E);
 


private:
  TSpline3* frange2E;
  TSpline3* fE2range;
  double fMass;
};
#endif
