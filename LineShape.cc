#include <sys/time.h>
#include <signal.h>
#include "TVector3.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TEnv.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TStopwatch.h"
#include "CommandLineInterface.hh"
#include "Range.hh"

TRandom3 *myrandom;
#ifndef PI
#define PI                       (TMath::Pi())
#endif

#ifndef SPEEDOFLIGHT
#define SPEEDOFLIGHT 299.792458 // mm/ns
#endif

TVector3* randdirection();
double doppler(double lab, double beta, TVector3* pg);
bool signal_received = false;
void signalhandler(int sig);
double get_time();
int main(int argc, char* argv[]){
  double time_start = get_time();  
  TStopwatch timer;
  timer.Start();
  signal(SIGINT,signalhandler);
  char *SetFile = NULL;
  char *OutFile = NULL;
  int NEvent =-1;
  double clifetime =-1;
  //Read in the command line arguments
  CommandLineInterface* interface = new CommandLineInterface();
  interface->Add("-o", "output file", &OutFile);    
  interface->Add("-s", "settings file", &SetFile);    
  interface->Add("-n", "number of events", &NEvent);  
  interface->Add("-l", "lifetime", &clifetime);  
  interface->CheckFlags(argc, argv);
  if(OutFile == NULL){
    cout << "No output ROOT file given " << endl;
    OutFile = (char*)"test.root";
  }
  if(SetFile == NULL){
    cout << "No settings file given " << endl;
    return 2;
  }
  if(NEvent <0){
    NEvent = 1000;
  }
  

  myrandom = new TRandom3();
  myrandom->ReadRandom("random.dat");

  double reactionZ;
  double reactionB;

  double emissionT;
  double emissionZ;
  double emissionB;
  TVector3 *emissionD_cm = new TVector3(0,0,1);
  double emissionE_cm;
  TVector3 *emissionD = new TVector3(0,0,1);
  double emissionE;

  double detectionE;
  TVector3 *detectionD = new TVector3(0,0,1);

  double correctionE;

  TFile *file = new TFile(OutFile, "RECREATE");

  TTree* tr = new TTree("tr","Simulated Gammas");
  tr->Branch("reactionZ",&reactionZ,32000,0);
  tr->Branch("reactionB",&reactionB,32000,0);

  tr->Branch("emissionT",&emissionT,32000,0);
  tr->Branch("emissionZ",&emissionZ,32000,0);
  tr->Branch("emissionB",&emissionB,32000,0);
  tr->Branch("emissionE_cm",&emissionE_cm,32000,0); 
  tr->Branch("emissionD_cm",&emissionD_cm,32000,0);
  tr->Branch("emissionE",&emissionE,32000,0); 
  tr->Branch("emissionD",&emissionD,32000,0);

  tr->Branch("detectionE",&detectionE,32000,0);
  tr->Branch("detectionD",&detectionD,32000,0);

  tr->Branch("correctionE",&correctionE,32000,0);
  //tr->BranchRef();
  
  TEnv *set = new TEnv(SetFile);
  double targetz = set->GetValue("TargetPosition",0.0);
  double targetthick = set->GetValue("TargetThickness",703.0);
  double targetdensity = set->GetValue("TargetDensity",1.85);
  targetdensity*=100; //to convert to mg/cm^2/mm
  double betabeam = set->GetValue("BetaBeam",0.6);
  char* targetprojrange = (char*)set->GetValue("TargetProjRange",(char*)"/home/wimmer/simulation/lineshape/range70Kr.dat");
  char* targetejecrange = (char*)set->GetValue("TargetEjecRange",(char*)"/home/wimmer/simulation/lineshape/range70Kr.dat");
  double egamma0 = set->GetValue("EGamma",1000.);
  double massprojinu = set->GetValue("MassProj",69.956040);
  double massejecinu = set->GetValue("MassEjec",69.956040);
  double lifetime = set->GetValue("LifeTime",0.010);
  double betareconstruct = set->GetValue("BetaReconstruct",0.010);
  if(clifetime>0)
    lifetime = clifetime;
  double resolutionE = set->GetValue("ResolutionE",0.01/2.355); 
  double rsphere = set->GetValue("RadiusSphere",150.);
  double resolutionP = set->GetValue("ResolutionP",2./2.355); 
  
  int plunger = set->GetValue("Plunger",0);
  double degraderthick = set->GetValue("DegraderThickness",200.0);
  double degraderdensity = set->GetValue("DegraderDensity",16.601);
  degraderdensity*=100; //to convert to mg/cm^2/mm
  char* degraderrange = (char*)set->GetValue("DegraderRange",(char*)"/home/wimmer/simulation/lineshape/70KronTa.dat");
  double degraderdistance = set->GetValue("DegraderDistance",2.0);
  if(plunger==1)
    targetz =0;
    
  TList* hlist = new TList();
  TH1F* egamdc = new TH1F("egamdc","egamdc",8000,0,8000);hlist->Add(egamdc);
  TH1F* egamdc_4pi = new TH1F("egamdc_4pi","egamdc_4pi",8000,0,8000);hlist->Add(egamdc_4pi);
  TH1F* egamdc_1pi = new TH1F("egamdc_1pi","egamdc_1pi",8000,0,8000);hlist->Add(egamdc_1pi);
  TH2F* egamdc_theta = new TH2F("egamdc_theta","egamdc_theta",180,0,180,8000,0,8000);hlist->Add(egamdc_theta);

  Range *rangeprojtarget = new Range(massprojinu,targetprojrange);
  Range *rangeejectarget = new Range(massejecinu,targetejecrange);
  Range *rangedegrader = new Range(massejecinu,degraderrange);
  cout << "simulating "<< NEvent << " events "<< endl;
  for(int i=0;i<NEvent;i++){
    if(signal_received){
      break;
    }
    //find the reaction point
    reactionZ = targetthick*myrandom->Uniform(0,1);
    //beta at reaction
    reactionB = rangeprojtarget->GetBetaAfter(betabeam,reactionZ);
    reactionZ/=targetdensity;
    reactionZ+=targetz;
    //time of gamma emission
    if(lifetime>0)
      emissionT = myrandom->Exp(lifetime);
    else
      emissionT = 0;
    
    //approximation: use velocity at creation time to get emission position
    emissionZ = reactionB*SPEEDOFLIGHT*emissionT+reactionZ;

    //set emission velocity (after target if ejectile has left the target)
    if(emissionZ>targetthick/targetdensity+targetz){
      //after target
      emissionB = rangeejectarget->GetBetaAfter(reactionB,targetthick);
      if(plunger==1){
	//after distance
	if(emissionZ > (targetthick/targetdensity+degraderdistance)){
	  //in degrader
	  if(emissionZ - (targetthick/targetdensity+degraderdistance) <degraderthick/degraderdensity)
	    emissionB = rangedegrader->GetBetaAfter(emissionB,(emissionZ - (targetthick/targetdensity+degraderdistance))*degraderdensity);
	  else
	    emissionB = rangedegrader->GetBetaAfter(emissionB,degraderthick);
	}
      }
    }
    else
      emissionB = rangeejectarget->GetBetaAfter(reactionB,emissionZ*targetdensity);

    //isotropic in cm system
    emissionD_cm = randdirection();
    emissionD = new TVector3(0,0,1);
    emissionD->SetPhi(emissionD_cm->Phi());
    //boost to labsystem
    emissionD->SetTheta(acos((cos(emissionD_cm->Theta())+emissionB)/(1+emissionB*cos(emissionD_cm->Theta()))));

    emissionE_cm = egamma0;
    emissionE = emissionE_cm/(rangeejectarget->beta2gamma(emissionB)*(1-emissionB*cos(emissionD->Theta())));
    detectionE=emissionE+myrandom->Gaus(0,detectionE*resolutionE);

    //find intersection of emission direction and sphere
    //https://en.wikipedia.org/wiki/Line–sphere_intersection
    TVector3 spherecenter(0,0,0);
    TVector3 emissionpoint(0,0,emissionZ);
    double a = (*emissionD) * (emissionpoint-spherecenter);
    double b = a*a-(emissionpoint-spherecenter)*(emissionpoint-spherecenter)+rsphere*rsphere;
    TVector3 f(sqrt(-1),sqrt(-1),sqrt(-1));
    if(b>0){
      double c = -a+sqrt(b);
      double d = -a-sqrt(b);
      if (d>0)
	cout << "d >0" << endl;
      if(c>0){
	f = emissionpoint+c*(*emissionD);
      }
    }
    detectionD->SetX(f.X()+myrandom->Gaus(0,resolutionP));
    detectionD->SetY(f.Y()+myrandom->Gaus(0,resolutionP));
    detectionD->SetZ(f.Z()+myrandom->Gaus(0,resolutionP));

    correctionE = doppler(detectionE,betareconstruct,detectionD);
    egamdc->Fill(correctionE);
    egamdc_theta->Fill(detectionD->Theta()*180/PI,correctionE);
    if(detectionD->Theta()>0.15 && detectionD->Theta() <3.00)
      egamdc_4pi->Fill(correctionE);
    if(detectionD->Theta()>0.69 && detectionD->Theta() <1.97)
      egamdc_1pi->Fill(correctionE);

    tr->Fill();
    if(i%10000 == 0){
      double time_end = get_time();
      cout << setw(5) << setiosflags(ios::fixed) << setprecision(1) << (100.*i)/NEvent <<
	" % done\t" << (Float_t)i/(time_end - time_start) << " events/s " << 
	(NEvent-i)*(time_end - time_start)/(Float_t)i << "s to go \r" << flush;
    }
  }
  
  myrandom->WriteRandom("random.dat");
  file->cd();
  tr->Write();
  hlist->Write();
  file->Close();
  double time_end = get_time();
  cout << "Program Run time: " << time_end - time_start << " s." << endl;
  timer.Stop();
  cout << "CPU time: " << timer.CpuTime() << "\tReal time: " << timer.RealTime() << endl;
}

double doppler(double lab, double beta, TVector3* pg){
  double gamma = 1/sqrt(1-beta*beta);
  return lab*gamma*(1-beta*cos(pg->Theta()));
}
TVector3* randdirection(){
  double costheta = myrandom->Uniform(-1,1);
  double phi = myrandom->Uniform(-PI,PI);
  TVector3 *v = new TVector3(0,0,1);
  v->SetTheta(acos(costheta));
  v->SetPhi(phi);
  return v;
}
void signalhandler(int sig){
  if (sig == SIGINT){
    signal_received = true;
  }
}

double get_time(){  
    struct timeval t;  
    gettimeofday(&t, NULL);  
    double d = t.tv_sec + (double) t.tv_usec/1000000;  
    return d;  
}  
