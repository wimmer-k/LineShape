#include "Range.hh"

Range::Range(double mass, char* filename){
  fMass = mass;
  CreateSplines(filename);
}

double Range::GetRange(double beta){
  double e = beta2E(beta);
  //cout << "beta = " << beta << ", energy = "<< e << " MeV, "<< e/fMass << "MeV/u";
  double r = fE2range->Eval(e/fMass);
  //cout << ", range = " << r << endl;
  return r;
}

double Range::GetBetaAfter(double beta, double thick){
  double r = GetRange(beta);
  double rafter = r - thick;
  //cout << "thickness " << thick << ", range after " << rafter << endl;
  double e = frange2E->Eval(rafter);
  double betaafter = E2beta(e*fMass);
  //cout << "betaafter = " << betaafter << ", energy = "<< e*fMass << " MeV, "<< e << "MeV/u" << endl;
  return betaafter;
}

void Range::CreateSplines(char *filename){
  
  ifstream data;
  data.open(filename);
  double e[30];
  double r[2][30];
  
  for(int i=0;i<30;i++){
    data >> e[i] >> r[0][i] >> r[1][i];
  }
  fE2range = new TSpline3("E2range",e,r[1],30);
  frange2E = new TSpline3("range2E",r[1],e,30);
}

double Range::beta2gamma(double beta){
  return 1./sqrt(1-beta*beta);
}
double Range::gamma2beta(double gamma){
  //cout << "gamma = " << gamma << endl;
  return sqrt(1-1./gamma/gamma);
}
double Range::beta2E(double beta){
  double gamma = beta2gamma(beta);
  return (gamma-1)*fMass*931.4940954;
}
double Range::E2beta(double E){
  //cout << __PRETTY_FUNCTION__ << endl;
  //cout << "E = "<< E << ", gamma2beta("<<E<<"/("<<fMass<<"*931.4940954)+1) = gamma2beta("<<E/(fMass*931.4940954)+1<<")"<< endl;
  return gamma2beta(E/(fMass*931.4940954)+1);
}
