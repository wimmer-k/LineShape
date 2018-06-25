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
  double e[600];
  double r[2][600];
  int n =0;
  data >> n;
  data.ignore(1000,'\n');
  for(int i=0;i<n;i++){
    data >> e[i] >> r[0][i] >> r[1][i];
    data.ignore(1000,'\n');
  }
  fE2range = new TSpline3("E2range",e,r[0],n);
  frange2E = new TSpline3("range2E",r[0],e,n);
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
