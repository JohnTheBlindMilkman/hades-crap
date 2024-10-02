#include "Includes.h" // A lot of header files from Root and Hydra
#include "GeomFunct.h" // Some fuctions to do DCA
#include "HParticleTool.h"
#include <vector>
#include <algorithm>
#include "TStyle.h"
#include "TSystem.h"
#include "henergylosscorrpar.h"

#include "hparcond.h"

#include "/lustre/nyx/hades/user/mgrunwal/eventChara/hparticleevtchara.h"

using namespace std; // cout - The most important debugging tool ;D

const double D2R = 1.74532925199432955e-02;
const double R2D = 57.2957795130823229; 


struct NParticle {
  Long_t evtId;
  TLorentzVector vec;
  Int_t nToFRPC;
  Int_t CENTRALITY;
  Float_t Xvertex;
  Float_t Yvertex;
  Float_t E;
  Float_t p;
  Float_t px;
  Float_t py;
  Float_t pz;
  Float_t mass;
  Float_t theta;
  Float_t phi;
  Float_t pT;
  Float_t rap;
  Float_t pseudorap;
  Float_t Mag;
  Int_t plate;
} particle_data;

int GetRapBin(double rap){
  if(rap <= 0.4) return 0;
  else if(rap >= 0.4 && rap <= 0.8) return 1;
  else if(rap >= 0.8 && rap <= 1.2) return 2;
  else if(rap >= 1.2 && rap <= 1.6) return 3;
  else if(rap >= 1.6 && rap <= 2.0) return 4;
  return 5;
}

int GetRapDiffClassBin(double deltarap){
  if(deltarap >= 0.0 && deltarap <= 0.1) return 0;
  else if(deltarap >= 0.1 && deltarap <= 0.2) return 1;
  else if(deltarap >= 0.2 && deltarap <= 0.3) return 2;
  else if(deltarap >= 0.3 && deltarap <= 0.4) return 3;
  else if(deltarap >= 0.4 && deltarap <= 0.5) return 4;
  return 5;
}

int GetkTDiffClassBin(double deltakT){
  if(deltakT >= 0.0 && deltakT <= 400.0) return 0;
  else if(deltakT >= 400.0 && deltakT <= 800.0) return 1;
  else if(deltakT >= 800.0 && deltakT <= 1200.0) return 2;
  else if(deltakT >= 1200.0 && deltakT <= 1600.0) return 3;
  return 4;
}

int SetTargeMomLlateNumber(double eventZvertex)
{
  int nVertex = 0;
  if(eventZvertex <= -63        && eventZvertex > -60)   nVertex = 1;
  else if(eventZvertex <= -60   && eventZvertex > -56.5) nVertex = 2;
  else if(eventZvertex <= -56.5 && eventZvertex > -53)   nVertex = 3;
  else if(eventZvertex <= -53   && eventZvertex > -49)   nVertex = 4;
  else if(eventZvertex <= -49   && eventZvertex > -45.5) nVertex = 5;
  else if(eventZvertex <= -45.5 && eventZvertex > -42)   nVertex = 6;
  else if(eventZvertex <= -42   && eventZvertex > -39)   nVertex = 7;
  else if(eventZvertex <= -39   && eventZvertex > -35.5) nVertex = 8;
  else if(eventZvertex <= -35.5 && eventZvertex > -31.5) nVertex = 9;
  else if(eventZvertex <= -31.5 && eventZvertex > -28)   nVertex = 10;
  else if(eventZvertex <= -28   && eventZvertex > -24.5) nVertex = 11;
  else if(eventZvertex <= -24.5 && eventZvertex > -21)   nVertex = 12;
  else if(eventZvertex <= -21   && eventZvertex > -16.5) nVertex = 13;
  else if(eventZvertex <= -16.5 && eventZvertex > -13.5) nVertex = 14;
  else nVertex = 15;
  return nVertex;
}
Int_t getTargeMomLlateNumber(Double_t eventZvertex){ //for mar19 data
  if(eventZvertex >= -63 && eventZvertex < -60) return 1;
  else if(eventZvertex >= -60 && eventZvertex < -56.5) return 2;
  else if(eventZvertex >= -56.5 && eventZvertex < -53) return 3;
  else if(eventZvertex >= -53 && eventZvertex < -49) return 4;
  else if(eventZvertex >= -49 && eventZvertex < -45.5) return 5;
  else if(eventZvertex >= -45.5 && eventZvertex < -42) return 6;
  else if(eventZvertex >= -42 && eventZvertex < -39) return 7;
  else if(eventZvertex >= -39 && eventZvertex < -35.5) return 8;
  else if(eventZvertex >= -35.5 && eventZvertex < -31.5) return 9;
  else if(eventZvertex >= -31.5 && eventZvertex < -28) return 10;
  else if(eventZvertex >= -28 && eventZvertex < -24.5) return 11;
  else if(eventZvertex >= -24.5 && eventZvertex < -21) return 12;
  else if(eventZvertex >= -21 && eventZvertex < -16.5) return 13;
  else if(eventZvertex >= -16.5 && eventZvertex < -13.5) return 14;
  else if(eventZvertex >= -13.5 && eventZvertex < -11) return 15;
  else return 1000;
}

double openingangle(const TLorentzVector& a, const TLorentzVector& b){
  return TMath::ACos( (a.Px()*b.Px() + a.Py()*b.Py() +  a.Pz()*b.Pz() ) / ( a.Vect().Mag() * b.Vect().Mag() ) );
}

double openingangle(const TVector3& a, const TVector3& b){
  return TMath::ACos( (a.Px()*b.Px() + a.Py()*b.Py() +  a.Pz()*b.Pz() ) / ( a.Mag() * b.Mag() ) );
}

Double_t correlateQinv(TLorentzVector v1, TLorentzVector v2){
  double dpx = v1.Px() - v2.Px();
  double dpy = v1.Py() - v2.Py();
  double dpz = v1.Pz() - v2.Pz();
  double dE = v1.E() - v2.E();
  return TMath::Sqrt(TMath::Abs(dpx*dpx + dpy*dpy + dpz*dpz - dE*dE));
}


Double_t correlatePairPRF(TLorentzVector v1, TLorentzVector v2){
  Double_t fpx1 = v1.Px();   Double_t fpx2 = v2.Px();
  Double_t fpy1 = v1.Py();   Double_t fpy2 = v2.Py();
  Double_t fpz1 = v1.Pz();   Double_t fpz2 = v2.Pz();
  Double_t fe1  = v1.E();   Double_t fe2  = v2.E();

  Double_t MomLx = fpx1 + fpx2;
  Double_t MomLy = fpy1 + fpy2;
  Double_t MomLz = fpz1 + fpz2;
  Double_t tE  = fe1 + fe2;
  Double_t MomLt = MomLx * MomLx + MomLy * MomLy;
  Double_t tMt = tE * tE - MomLz * MomLz;  // mCVK;
  Double_t tM  = TMath::Sqrt(tMt - MomLt);
  tMt          = TMath::Sqrt(tMt);
  MomLt          = TMath::Sqrt(MomLt);

  Double_t tBeta  = MomLz / tE;
  Double_t tGamma = tE / tMt;
  Double_t fZ     = tGamma * (fpz1 - tBeta * fe1);
  Double_t tE1L   = tGamma * (fe1 - tBeta * fpz1);


  // Rotate in transverse plane
  Double_t fX = (fpx1 * MomLx + fpy1 * MomLy) / MomLt;
  Double_t fY = (-fpx1 * MomLy + fpy1 * MomLx) / MomLt;

  // Boost to pair cms
  fX = tMt / tM * (fX - MomLt / tMt * tE1L);

  Double_t fT = fX > 0. ? 1. : -1.;
  return TMath::Sqrt(fY * fY + fX * fX + fZ * fZ);
}

Double_t correlatePairPRFGKine(HGeantKine *v1, HGeantKine *v2){
  float vx1, vy1, vz1;
  v1->getMomentum(vx1,vy1,vz1);
  float vx2, vy2, vz2;
  v2->getMomentum(vx2,vy2,vz2);
  Double_t fpx1 = vx1;   Double_t fpx2 = vx2;
  Double_t fpy1 = vy1;   Double_t fpy2 = vy2;
  Double_t fpz1 = vz1;   Double_t fpz2 = vz2;
  Double_t fe1  = v1->getE();   Double_t fe2  = v2->getE();

  Double_t MomLx = fpx1 + fpx2;
  Double_t MomLy = fpy1 + fpy2;
  Double_t MomLz = fpz1 + fpz2;
  Double_t tE  = fe1 + fe2;
  Double_t MomLt = MomLx * MomLx + MomLy * MomLy;
  Double_t tMt = tE * tE - MomLz * MomLz;  // mCVK;
  Double_t tM  = TMath::Sqrt(tMt - MomLt);
  tMt          = TMath::Sqrt(tMt);
  MomLt          = TMath::Sqrt(MomLt);

  Double_t tBeta  = MomLz / tE;
  Double_t tGamma = tE / tMt;
  Double_t fZ     = tGamma * (fpz1 - tBeta * fe1);
  Double_t tE1L   = tGamma * (fe1 - tBeta * fpz1);


  // Rotate in transverse plane
  Double_t fX = (fpx1 * MomLx + fpy1 * MomLy) / MomLt;
  Double_t fY = (-fpx1 * MomLy + fpy1 * MomLx) / MomLt;

  // Boost to pair cms
  fX = tMt / tM * (fX - MomLt / tMt * tE1L);

  Double_t fT = fX > 0. ? 1. : -1.;
  //Double_t fT *= TMath::Sqrt(fY * fY + fX * fX + fZ * fZ);
  return TMath::Sqrt(fY * fY + fX * fX + fZ * fZ);
}

void FillParticleStruc(NParticle &particle, Long64_t event1, TLorentzVector TLv, int nTOFRPC, double eventZvertex, double X_vertex, double Y_vertex, Int_t centrality_Class){
  particle.evtId      = event1;
  particle.nToFRPC    = nTOFRPC;
  particle.CENTRALITY = centrality_Class;
  particle.Xvertex    = X_vertex;
  particle.Yvertex    = Y_vertex;
  particle.p          = TLv.P();
  particle.E          = TLv.E();
  particle.px         = TLv.Px();
  particle.py         = TLv.Py();
  particle.pz         = TLv.Pz();
  particle.pT         = TLv.Pt();
  particle.mass       = TLv.M();
  particle.phi        = TLv.Phi()*TMath::DegToRad();
  particle.theta      = TLv.Theta()*TMath::DegToRad();
  particle.rap        = TLv.Rapidity();
  particle.pseudorap  = TLv.PseudoRapidity();
  particle.vec        = TLv;
  particle.Mag        = TLv.Vect().Mag();
  particle.plate      = SetTargeMomLlateNumber(eventZvertex);
}

// event mixing selection
bool EvMixSelection(NParticle p1, NParticle p2){
  if(p1.CENTRALITY == p2.CENTRALITY && p1.plate == p2.plate &&
  (abs(p1.Xvertex-p2.Xvertex) < 0.2) && (abs(p1.Yvertex-p2.Yvertex) < 0.2) )
    return 1;
  else return 0;
}

void MakeSignal(vector <NParticle> particle1, vector <NParticle> particle2, TH1D* hsignal){
  for(int ipart = 0; ipart < particle1.size(); ipart++){
    NParticle p1 = particle1[ipart];
    for(int jpart = 0; jpart < particle2.size(); jpart++){
      NParticle p2 = particle2[jpart];

      double phi1, phi2, theta1, theta2;

      phi1=p1.phi; phi2=p2.phi; theta1=p1.theta; theta2=p2.theta;

      double dTHETA = theta1 - theta2;
      double dPHI = phi1 - phi2;
      double dPt  = p1.pT - p2.pT;                          //---------- delta-Pt ------------------
      double Pt2  = 2*(p1.pT*p2.pT);                        //---------- Pt2 -----------------------
      double dphi  = TMath::ACos(1 - ((dPt*dPt)/(Pt2)));    //---------- delta-Phi -----------------
      double dpRap = p1.pseudorap - p2.pseudorap;           //---------- delta-Pseudo Rapidity -----
      double PairRap = (p1.rap + p2.rap)/2.0;                     //---------- Total Pair Rapidity -------
      double PairPt  = (p1.pT + p2.pT)/2.0;
      double kstar = fabs(p1.p - p2.p)/2.0;                 //---------- |p1 - p2|/2 ---------------

      if(p2.evtId == p1.evtId) hsignal->Fill(correlatePairPRF(p1.vec, p2.vec));
    }
  }
} 

void MakeBackground(vector <NParticle> particle1, vector <NParticle> particle2, TH1D* hback, int maxMixEvents){
  for(int ipart = 0; ipart < particle1.size(); ipart++){
    int mixEvCounter = 0;
    NParticle p1 = particle1[ipart];
    int tempEv = p1.evtId;
    for(int jpart = 0; jpart < particle2.size(); jpart++){
    NParticle p2 = particle2[jpart];

    double phi1, phi2, theta1, theta2;
    phi1=p1.phi; phi2=p2.phi; theta1=p1.theta; theta2=p2.theta; 

    double dTHETA = theta1 - theta2;
    double dPHI = phi1 - phi2;
    double dPt  = p1.pT - p2.pT;                          //---------- delta-Pt -----------------
    double Pt2  = 2*(p1.pT*p2.pT);                        //---------- Pt2 ----------------------
    double dphi  = TMath::ACos(1 - ((dPt*dPt)/(Pt2)));    //---------- delta-Phi ----------------
    double dpRap = p1.pseudorap - p2.pseudorap;           //---------- delta-Pseudo Rapidity ----

    if(p2.evtId != p1.evtId && EvMixSelection(p1,p2)!=0){
      if(mixEvCounter == maxMixEvents) {
        cout << "Break mixing!" << endl;
        break;
      }

      if(mixEvCounter >= maxMixEvents) cout << "mixEvCounter overflow!!!! " << mixEvCounter << " : " << maxMixEvents << endl;

      if(tempEv!=p2.evtId){
        tempEv=p2.evtId;
        mixEvCounter++;
      }   
      hback->Fill(correlatePairPRF(p1.vec, p2.vec));
      }
    }
  }
}

void MakeSignal_Purity(vector <NParticle> particle1, vector <NParticle> particle2, TH1D* hsignal){
  for(int ipart = 0; ipart < particle1.size(); ipart++){
    NParticle p1 = particle1[ipart];
    for(int jpart = 0; jpart < particle2.size(); jpart++){
      NParticle p2 = particle2[jpart];
                                
      if(p2.evtId == p1.evtId) hsignal->Fill(correlatePairPRF(p1.vec, p2.vec));
    }
  }
}
void MakeBackground_Purity(vector <NParticle> particle1, vector <NParticle> particle2, TH1D* hback, int maxMixEvents){
  for(int ipart = 0; ipart < particle1.size(); ipart++){
  int mixEvCounter = 0;
  NParticle p1 = particle1[ipart];
  int tempEv = p1.evtId;
  for(int jpart = 0; jpart < particle2.size(); jpart++){
    NParticle p2 = particle2[jpart];
    if(tempEv != p2.evtId) mixEvCounter++;
                              
    if(p2.evtId != p1.evtId && EvMixSelection(p1,p2)!=0){
      if(mixEvCounter == maxMixEvents) break;

      if(tempEv!=p2.evtId){
        tempEv==p2.evtId;
        mixEvCounter++;	
      }   
      hback->Fill(correlatePairPRF(p1.vec, p2.vec));
    }
    }
  }
}

void MakeSignal_cent(vector <NParticle> particle1,vector <NParticle> particle2,TH1D* hsignal[]){
  for(int ipart = 0; ipart < particle1.size(); ipart++){
    NParticle p1 = particle1[ipart];
    for(int jpart = 0; jpart < particle2.size(); jpart++){
      NParticle p2 = particle2[jpart];

      double phi1, phi2, theta1, theta2;
      phi1=p1.phi; phi2=p2.phi; theta1=p1.theta; theta2=p2.theta;
        
      double dTHETA = theta1 - theta2;
      double dPHI = phi1 - phi2;
      double dPt  = p1.pT - p2.pT;                          //---------- delta-Pt -----------------
      double Pt2  = 2*(p1.pT*p2.pT);                        //---------- Pt2 ----------------------
      double dphi  = TMath::ACos(1 - ((dPt*dPt)/(Pt2)));    //---------- delta-Phi ----------------
      double dpRap = p1.pseudorap - p2.pseudorap;           //---------- delta-Pseudo Rapidity ----

      if(p2.evtId == p1.evtId) hsignal[p1.CENTRALITY-1]->Fill(correlatePairPRF(p1.vec, p2.vec));
      
    }
  }
} 

void MakeBackground_cent(vector <NParticle> particle1,vector <NParticle> particle2,TH1D* hback[],int maxMixEvents){
  for(int ipart = 0; ipart < particle1.size(); ipart++){
  int mixEvCounter = 0;
  NParticle p1 = particle1[ipart];
  int tempEv = p1.evtId;
    for(int jpart = 0; jpart < particle2.size(); jpart++){
      NParticle p2 = particle2[jpart];

      if(tempEv != p2.evtId) mixEvCounter++;

      double phi1, phi2, theta1, theta2;
      phi1=p1.phi; phi2=p2.phi; theta1=p1.theta; theta2=p2.theta; 

      double dTHETA = theta1 - theta2;
      double dPHI = phi1 - phi2;
      double dPt  = p1.pT - p2.pT;                          //---------- delta-Pt -----------------
      double Pt2  = 2*(p1.pT*p2.pT);                        //---------- Pt2 ----------------------
      double dphi  = TMath::ACos(1 - ((dPt*dPt)/(Pt2)));    //---------- delta-Phi ----------------
      double dpRap = p1.pseudorap - p2.pseudorap;           //---------- delta-Pseudo Rapidity ----

      if(p2.evtId != p1.evtId && EvMixSelection(p1,p2)!=0){
        if(mixEvCounter == maxMixEvents) break;

        if(tempEv!=p2.evtId){
          tempEv==p2.evtId;
          mixEvCounter++;
        }   

      hback[p1.CENTRALITY-1]->Fill(correlatePairPRF(p1.vec, p2.vec));
      }
    }
  }
}

void MakeSignal_cent_Purity(vector <NParticle> particle1, vector <NParticle> particle2, TH1D* hsignal[]){
  for(int ipart = 0; ipart < particle1.size(); ipart++){
    NParticle p1 = particle1[ipart];
    for(int jpart = 0; jpart < particle2.size(); jpart++){
      NParticle p2 = particle2[jpart];
                                
      if(p2.evtId == p1.evtId) hsignal[p1.CENTRALITY-1]->Fill(correlatePairPRF(p1.vec, p2.vec));
    }
  }
}
void MakeBackground_cent_Purity(vector <NParticle> particle1, vector <NParticle> particle2, TH1D* hback[], int maxMixEvents){
  for(int ipart = 0; ipart < particle1.size(); ipart++){
  int mixEvCounter = 0;
  NParticle p1 = particle1[ipart];
  int tempEv = p1.evtId;
  for(int jpart = 0; jpart < particle2.size(); jpart++){
    NParticle p2 = particle2[jpart];
    if(tempEv != p2.evtId) mixEvCounter++;
                              
    if(p2.evtId != p1.evtId && EvMixSelection(p1,p2)!=0){
      if(mixEvCounter == maxMixEvents) break;

      if(tempEv!=p2.evtId){
        tempEv==p2.evtId;
        mixEvCounter++; 
      }   
      hback[p1.CENTRALITY-1]->Fill(correlatePairPRF(p1.vec, p2.vec));
    }
    }
  }
}


void MakeSignal_kT(vector <NParticle> particle1,vector <NParticle> particle2,TH1D* hsignal_kT[])
{
  for(int ipart = 0; ipart < particle1.size(); ipart++){
    NParticle p1 = particle1[ipart];
    for(int jpart = 0; jpart < particle2.size(); jpart++){
      NParticle p2 = particle2[jpart];

      double kt   = (p1.pT + p2.pT)/2.;                        //---------- Average of Pt ----------------------

      if(p2.evtId == p1.evtId){
        int kTBin = GetkTDiffClassBin(kt);
        hsignal_kT[kTBin]->Fill(correlatePairPRF(p1.vec, p2.vec));
      }
    }
  }
} 


void MakeBackground_kT(vector <NParticle> particle1,vector <NParticle> particle2, TH1D* hback_kT[], int maxMixEvents){
  for(int ipart = 0; ipart < particle1.size(); ipart++){
    int mixEvCounter = 0;
    NParticle p1 = particle1[ipart];
    int tempEv = p1.evtId;
    for(int jpart = 0; jpart < particle2.size(); jpart++){
      NParticle p2 = particle2[jpart];

      if(tempEv != p2.evtId) mixEvCounter++;

      double kt   = (p1.pT + p2.pT)/2.;                        //---------- Pt2 ----------------------

      if(p2.evtId != p1.evtId && EvMixSelection(p1,p2)!=0){
        if(mixEvCounter == maxMixEvents) break;

        if(tempEv!=p2.evtId){
          tempEv==p2.evtId;
          mixEvCounter++;
        }
        int kTBin = GetkTDiffClassBin(kt);
        hback_kT[kTBin]->Fill(correlatePairPRF(p1.vec, p2.vec));
      }
    }
  }
}


void MakeSignal_kT_Purity(vector <NParticle> particle1, vector <NParticle> particle2, TH1D* hsignal_kT[]){
  for(int ipart = 0; ipart < particle1.size(); ipart++){
    NParticle p1 = particle1[ipart];
    for(int jpart = 0; jpart < particle2.size(); jpart++){
      NParticle p2 = particle2[jpart];

      double kt   = (p1.pT + p2.pT)/2.;   
                                
       if(p2.evtId == p1.evtId){
        int kTBin = GetkTDiffClassBin(kt);
        hsignal_kT[kTBin]->Fill(correlatePairPRF(p1.vec, p2.vec));
      }
    }
  }
}
void MakeBackground_kT_Purity(vector <NParticle> particle1, vector <NParticle> particle2, TH1D* hback_kT[], int maxMixEvents){
  for(int ipart = 0; ipart < particle1.size(); ipart++){
  int mixEvCounter = 0;
  NParticle p1 = particle1[ipart];
  int tempEv = p1.evtId;
  for(int jpart = 0; jpart < particle2.size(); jpart++){
    NParticle p2 = particle2[jpart];
    if(tempEv != p2.evtId) mixEvCounter++;

    double kt   = (p1.pT + p2.pT)/2.;  
                              
    if(p2.evtId != p1.evtId && EvMixSelection(p1,p2)!=0){
      if(mixEvCounter == maxMixEvents) break;

      if(tempEv!=p2.evtId){
        tempEv==p2.evtId;
        mixEvCounter++; 
      }   
      int kTBin = GetkTDiffClassBin(kt);
      hback_kT[kTBin]->Fill(correlatePairPRF(p1.vec, p2.vec));
    }
    }
  }
}


struct MomResParameters{
  double PhiA;
  double PhiB;
  double PhiC;
  double PhiD;
  double PhiE;
  double PhiF;
  double ThetaA;
  double ThetaB;
  double ThetaC;
  double ThetaD;
  double ThetaE;
  double ThetaF;
  double PResA;
  double PResB;
  double PResC;
  double PResD;
  double PResE;
  double PResF;
};

void setResParamsSigma(MomResParameters *lambda, MomResParameters *deuteron){
  lambda->PResA = 1826.4;
  lambda->PResB = -361.082;
  lambda->PResC = 1.11948;
  lambda->PResD = -0.000556701;
  lambda->ThetaA = 0.00155009;
  lambda->ThetaB = 2.27965;
  lambda->ThetaC = 204.219;
  lambda->PhiA = 0.00434149;
  lambda->PhiB = -1.57389;
  lambda->PhiC = 9648.34;
  lambda->PhiD = -7.07828e+06;

  deuteron->PResA = -348.439;
  deuteron->PResB = 71.6403;
  deuteron->PResC = -0.178639;
  deuteron->PResD = 7.31801e-05;
  deuteron->ThetaA = 0.000481703;
  deuteron->ThetaB = -0.16082;
  deuteron->ThetaC = 4793.77;
  deuteron->PhiA = 0.00406378;
  deuteron->PhiB = -9.47511;
  deuteron->PhiC = 21984.6;
  deuteron->PhiD = -8.68765e+06;
}

void setResParamsMean(MomResParameters *lambda, MomResParameters *deuteron){
  lambda->PResA = 78.403;
  lambda->PResB = -0.16415;
  lambda->PResC = 0.000126503;
  lambda->PResD = -3.23618e-08;
  lambda->ThetaA = -0.000494939;
  lambda->ThetaB = 3.69832e-07;
  lambda->PhiA = -0.000191174;
  lambda->PhiB = 9.76148e-08;

  deuteron->PResA = 59.3621;
  deuteron->PResB = -0.117937;
  deuteron->PResC = 7.58512e-05;
  deuteron->PResD = -1.29528e-08;
  deuteron->ThetaA = -0.000392047;
  deuteron->ThetaB = 1.68789e-07;
  deuteron->PhiA = -3.87536e-05;
  deuteron->PhiB = 1.45387e-08;
}

//================================================================================
//================================================================================
// Main Function of the Macro with 4 parameters
// 1. Depends on value of 4th parameter:
//   -1 (default): Comma separated DST file list to run the macro on the cluster
//   -2: Regular expression string describing DST file to be added
//   -3: Name / Path of a file containing one DST file per line to be added
//   Any other case: Unused
// 2. Name / Path of the ouMomLutfile to be created
// 3. Amount of events to process (-1 means all events)
// 4. If positive amount of random DST files to select, otherwise see 1st parameter.
//================================================================================
Int_t dLambda_SIM_momRes(TString input = "", TString outFile = "dLambda_SIM_v1.root", Long64_t nEventsDesired = -1, Int_t nInFiles = -1) {
  //--------------------------------------------------------------------------------
  // Initialization of the global ROOT object and the Hades Loop
  // The Hades Loop used as an interface to the DST data (Basically a container of a TChain).
  // kTRUE - The global HADES object is being created if not existing
  //--------------------------------------------------------------------------------

  TROOT trDSTAnalysis("trDSTAnalysis", "Off Vertex Particle Reconstruction Macro");
  HLoop *gloop = new HLoop(kTRUE); // Later accessed via the global gLoop variable!

  // For mom res
  MomResParameters par_sigma_lambda;
  MomResParameters par_mean_lambda;
  MomResParameters par_sigma_deuteron;
  MomResParameters par_mean_deuteron;

  setResParamsSigma(&par_sigma_lambda, &par_sigma_deuteron);
  setResParamsMean(&par_mean_lambda, &par_mean_deuteron);

  //--------------------------------------------------------------------------------
  // The following block finds / adds the input DST files to the HLoop
  //--------------------------------------------------------------------------------
  if (nInFiles == -1) gloop->addMultFiles(input);
  else {  
    Int_t nFiles = 0;

    TString inputFolder = "/lustre/nyx/hades/dstsim/mar19/ag1580ag/gen5/bmax8/no_enhancement_gcalor/root"; //path to data simulated
    TSystemDirectory* inputDir = new TSystemDirectory("inputDir",inputFolder);
    TList* files = inputDir->GetListOfFiles();

    for(Int_t i = 0; i <= files->LastIndex() && nFiles < nInFiles; i++){
      if(((TSystemFile*) files->At(i))->IsDirectory()) continue;

      gloop->addFile(inputFolder + "/" + ((TSystemFile*) files->At(i))->GetName());
      nFiles++;
    }
  }

  //--------------------------------------------------------------------------------
  // Booking the categories to be read from the DST files.
  // By default all categories are booked therefore -* (Unbook all) first and book the ones needed
  // All required categories have to be booked except the global Event Header which is always booked
  //--------------------------------------------------------------------------------
  if (!gloop->setInput("-*,+HParticleEvtInfo,+HParticleCandSim,+HGeantKine"))
  exit(1);  

  //--------------------------------------------------------------------------------
  // Setting the cache size of the HLoop internal TChain to read data from the DSTs
  // Improves performance of the lustre storage by decreasing load on lustre META servers
  //--------------------------------------------------------------------------------
  gloop->getChain()->SetCacheSize(8e6); // 8Mb
  gloop->getChain()->AddBranchToCache("*", kTRUE);
  gloop->getChain()->StopCacheLearningPhase();  
  gloop->printCategories(); // Just for informative purposes

  //--------------------------------------------------------------------------------
  // Creating the placeholder variables to read data from categories and getting categories
  //(They \ have to be booked!)
  //--------------------------------------------------------------------------------

  HParticleCandSim*     particle_cand;
  HGeantKine*               kine_cand;
  HEventHeader*         event_header;
  HParticleEvtInfo*     particle_info;

  HCategory* particle_info_cat  = (HCategory*) HCategoryManager::getCategory(catParticleEvtInfo);
  HCategory* particle_cand_cat  = (HCategory*) HCategoryManager::getCategory(catParticleCand);
  HCategory* kine_cand_cat      = (HCategory*) HCategoryManager::getCategory(catGeantKine);

  //--------------------------------------------------------------------------------

  //--------------------------------------------------------------------------------
  // Creating the placeholder variables to read data from categories and getting categories
  //(They \ have to be booked!)
  //--------------------------------------------------------------------------------
  if (!particle_cand_cat)  // If the category for the reconstructed trackes does not exist the macro makes no sense
  exit(1);    
  //--------------------------------------------------------------------------------
  // Lists to store the possible daugther tracks
  //--------------------------------------------------------------------------------
  vector<HParticleCandSim*> vDau1Tracks, vDau2Tracks, Deuteron;

  //=============================================================================================================================================
  // Put your object declarations here
  //=============================================================================================================================================

  //============================================================================================================================                                                                              
  //TFile* outFile = new TFile(outfile.Data(), "RECREATE");
  TFile *file = new TFile("/lustre/hades/user/mstefan/sub/loopDST/Mar19AgAg1580.root");
  TCutG *protonRPCCut = (TCutG*)file->Get("tcgPBetaProtonRPC2Sig");
  TCutG *protonToFCut = (TCutG*)file->Get("tcgPBetaProtonToF2Sig");
  TCutG *piNRPCCut = (TCutG*)file->Get("tcgPBetaPiMRPC2Sig");
  TCutG *piNToFCut = (TCutG*)file->Get("tcgPBetaPiMToF2Sig");
  TCutG *piPRPCCut = (TCutG*)file->Get("tcgPBetaPiPRPC2Sig");
  TCutG *piPToFCut = (TCutG*)file->Get("tcgPBetaPiPToF2Sig");

  TFile *file2 = new TFile("/lustre/hades/user/mstefan/sub/loopDST/Mar19AgAg1580_Gen5.root");
  TCutG *dRPCCut = (TCutG*)file2 -> Get("tcgPdEdxMDCDeuteronRPC2Sig");
  TCutG *dToFCut = (TCutG*)file2 -> Get("tcgPdEdxMDCDeuteronTOF2Sig");

  vector<NParticle> ProtonRPC_NPart;
  vector<NParticle> pionNRPC_NPart;
  vector<NParticle> HBMomLrotonRPC_NPart;
  vector<NParticle> HBMomLionNRPC_NPart;
  vector<NParticle> HBTLambdaRPC_NPart;
  vector<NParticle> HBTDeuteron_NPart;


  vector<TLorentzVector> usedProtonRPC;
  vector<TLorentzVector> usedPionNRPC;
  vector<TLorentzVector> HBTLambdaRPC;
  vector<TLorentzVector> HBMomLrotonRPC;

  TLorentzVector tlv_d;

  vector<NParticle> ProtonL;
  vector<NParticle> PionL;
  vector<NParticle> DeuteronP;

  vector<NParticle> Deuteron_NPart_Purity;
  vector<NParticle> HBTLambdaRPC_NPart_Purity;


  //--------------------------------------------------------------------------------   
  int maxMixEvents = 5;

  //--------------------------------------------------------------------------------
  // The following counter histogram is used to gather some basic information on the analysis
  //--------------------------------------------------------------------------------
  enum eCounters {
    cNumAllEvents      = 0,
    cNumSelectedEvents = 1,
    cNumAllTracks      = 2,
    cNumSelectedTracks = 3,
    cNumDau1Tracks     = 4,
    cNumDau2Tracks     = 5,
    cNumMotCan         = 6,
    cNumMotRec         = 7
  };

  TH1D* hCounter = new TH1D("hCounter", "", 8, 0., 8.);
  TH1D *hLambdamassRPC = new TH1D("hLambdamass_RPC_ToF","(RPC || ToF);Mass [MeV/c]" , 2000, 0.0, 2000.0);
  TH1D *hLambdamassRPC_checked = new TH1D("hLambdamass_RPC_ToF_checked","(RPC || ToF);Mass [MeV/c]" , 2000, 0.0, 2000.0);

  TH1D *hsignalRPC = new TH1D("hsignal_RPC_ToF","k* (RPC || ToF);k* [MeV/c]" , 1000, 0.0, 2000.0);
  TH1D *hbackgroundRPC = new TH1D("hbackground_RPC_ToF","k* (RPC || ToF);k* [MeV/c]" , 1000, 0.0, 2000.0);


  TH1D *hsignal_Purity = new TH1D("hsignal_RPC_ToF_Purity","k* (RPC || ToF);k* [MeV/c]" , 1000, 0.0, 2000.0);
  TH1D *hbackground_Purity = new TH1D("hbackground_RPC_ToF_Purity","k* (RPC || ToF);k* [MeV/c]" , 1000, 0.0, 2000.0);


  TH2D *hsigAlphaPt_NoMC = new TH2D("hsigAlphaPt_NoMC","Armenteros-Podolanski plot (no MC); #alpha; Pt [MeV/c]", 2000, -1, 1, 2000, 0.0, 2000.);
  TH2D *hsigAlphaPt_MC = new TH2D("hsigAlphaPt_MC","Armenteros-Podolanski plot (with MC); #alpha; Pt [MeV/c]", 2000, -1, 1, 2000, 0.0, 2000.);
  TH2D *hsigAlphaPt_MC2 = new TH2D("hsigAlphaPt_MC2","Armenteros-Podolanski plot (with MC tight); #alpha; Pt [MeV/c]",2000,-1,1,2000,0.0,2000.);

  TH1D *hMultiLambdaCounter = new TH1D("hMultiLambdaCounter","Number of Lambdas in same event", 4, 1, 5);
  TH1D *hCounterparticle = new TH1D("hCounterparticle","Number of particles(#lambda) in same event", 5, 0, 5);
  TH1D *LambdaProton = new TH1D("LambdaProton", "LambdaProton",10,0,10);

  //----------Mass
  TH1D *hProtonMass = new TH1D("hProtonMass", "hProtonMass; Mass [MeV/c2]; Counts", 2000, 0.0, 2000.0); 
  TH1D *hPionMass = new TH1D("hPionMass", "hPionMass; Mass [MeV/c2]; Counts", 2000, 0.0, 2000.0);
  TH1D *hLambdaMass = new TH1D("hLambdaMass", "hLambdaMass; Mass [MeV/c2]; Counts", 2000, 0.0, 2000.0);
  //----------Rapidity
  TH1D *hProtonRap = new TH1D("hProtonRap", "hProtonRap; Rapidity; Counts", 2000, -4, 4); 
  TH1D *hPionRap =   new TH1D("hPionRap", "hPionRap; Rapidity; Counts", 2000, -4, 4);
  TH1D *hLambdaRap = new TH1D("hLambdaRap", "hLambdaRap; Rapidity; Counts", 2000, -4, 4);
  //----------Transverse Momentum
  TH1D *hProtonPt = new TH1D("hProtonPt", "hProtonPt; Pt [MeV/c]; Counts", 2000, 0.0, 2000.0); 
  TH1D *hPionPt =   new TH1D("hPionPt", "hPionPt; Pt [MeV/c]; Counts", 2000, 0.0, 2000.0);
  TH1D *hLambdaPt = new TH1D("hLambdaPt", "hLambdaPt; Pt [MeV/c]; Counts", 2000, 0.0, 2000.0);
  //----------Momentum
  TH1D *hProtonP = new TH1D("hProtonP", "hProtonP; P [MeV/c]; Counts", 2000, 0.0, 2000.0); 
  TH1D *hPionP =   new TH1D("hPionP", "hPionP; P [MeV/c]; Counts", 2000, 0.0, 2000.0);
  TH1D *hLambdaP = new TH1D("hLambdaP", "hLambdaP; P [MeV/c]; Counts", 2000, 0.0, 2000.0);
  //----------Pion-Proton opening angle
  TH1D *hPPiOA  = new TH1D("hPPiOA", "hPPiOA; OA; Counts", 100, 0.0, 50.);


  hCounter->GetXaxis()->SetBinLabel(1, "All Events");
  hCounter->GetXaxis()->SetBinLabel(2, "Selected Events");
  hCounter->GetXaxis()->SetBinLabel(3, "All Tracks");
  hCounter->GetXaxis()->SetBinLabel(4, "Selected Tracks");
  hCounter->GetXaxis()->SetBinLabel(5, "Daughter 1 Tracks");
  hCounter->GetXaxis()->SetBinLabel(6, "Daughter 2 Tracks");
  hCounter->GetXaxis()->SetBinLabel(7, "Daughter combinations");
  hCounter->GetXaxis()->SetBinLabel(8, "Reconstructed Mothers");

  hCounter->Fill(-1., 0.); // For whatever reason this is required to merge the histogram correctly!

  TH2D *hbackPairRapvsPairPt = new TH2D("hbackPairRapvsPairPt",";kT [MeV/c]; Rapidity",1200,-3000, 9000,900,-3.,6.);

  // Event monitor //
  TH1D *hMultiplicity = new TH1D("hMultiplicity","Multiplicity;[Nch];", 400, 0, 401);
  // Track monitors //
  TH1D *hZVertex      = new TH1D("hZVertex","Z Vertex;[mm];",1000, -70, 0);
  TH2D *hBetaMomToF   = new TH2D("hBetaMomToF","#beta vs. p (ToF);p [MeV/c^{2}];#beta",6000,-3000, 3000, 1100, 0, 1.1);
  TH2D *hBetaMomRPC   = new TH2D("hBetaMomRPC","#beta vs. p (RPC);p [MeV/c^{2}];#beta",6000,-3000, 3000, 1100, 0, 1.1); 

  TH1D *hMomentum = new TH1D("hMomentum deuteron", "hMomentum deuteron", 500, 0, 2000);
  TH1D *hTransMomentum = new TH1D("hTransMomentum deuteron", "hTransMomentum deuteron", 500, 0, 2000);
  TH1D *hRapidity = new TH1D("hRapidity deuteron", "hRapidity deuteron", 500, -1, 4);

  //----------BetaMom
  TH2D *hProtonBetaMom = new TH2D("hProtonBetaMom", "hProtonBetaMom; p [MeV/c]; #beta", 2000, 0.0, 2000.0, 600, 0.0, 1.2);
  TH2D *hPionNBetaMom = new TH2D("hPionNBetaMom", "hPionNBetaMom; p [MeV/c]; #beta", 2000, 0.0, 2000.0, 600, 0.0, 1.2);

  TH2D *hDeuteronBetaMom = new TH2D("hDeuteronBetaMom", "hDeuteronBetaMom; p [MeV/c]; #beta", 2000, 0.0, 2000.0, 600, 0.0, 1.2);
  TH2D *hDeuterondEdx = new TH2D("hDeuterondEdx", "hDeuterondEdx; p [MeV/c]; #beta", 2000, 0.0, 2000.0, 600, 0.0, 100);

  TH2D *hDeuteronBetaMom_checked = new TH2D("hDeuteronBetaMom_checked", "hDeuteronBetaMom; p [MeV/c]; #beta", 2000, 0.0, 2000.0, 600, 0.0, 1.2);
  TH2D *hDeuterondEdx_checked = new TH2D("hDeuterondEdx_checked", "hDeuterondEdx; p [MeV/c]; #beta", 2000, 0.0, 2000.0, 600, 0.0, 100);

  int CentN = 4;

  TH1D *hSIGNALRPC[CentN];
  TH1D *hBACKGROUNDRPC[CentN];
  TH1D *hSIGNALRPC_Purity[CentN];
  TH1D *hBACKGROUNDRPC_Purity[CentN];

  for( int icent = 0; icent < CentN; icent++){
    hSIGNALRPC[icent] = new TH1D(Form("Signal_function_proton-lambda_c%d", icent), Form("Signal_function_proton-lambda_c%d; k*[MeV/c]; C(k*)",icent), 1000,0,2000);
    hBACKGROUNDRPC[icent] = new TH1D(Form("Background_function_proton-lambda_c%d", icent), Form("Background_function_proton-lambda_c%d; k*[MeV/c]; C(k*)",icent), 1000, 0, 2000);
    hSIGNALRPC_Purity[icent] = new TH1D(Form("Signal_function_proton-lambda_c%d_Purity", icent), Form("Signal_function_proton-lambda_c%d_Purity; k*[MeV/c]; C(k*)",icent), 1000,0,2000);
    hBACKGROUNDRPC_Purity[icent] = new TH1D(Form("Background_function_proton-lambda_c%d_Purity", icent), Form("Background_function_proton-lambda_c%d_Purity; k*[MeV/c]; C(k*)",icent), 1000, 0, 2000);
  }  

  int kTN = 5;

  TH1D *hSIGNALkT[kTN];
  TH1D *hBACKGROUNDkT[kTN];
  TH1D *hSIGNALkT_Purity[kTN];
  TH1D *hBACKGROUNDkT_Purity[kTN];

  for( int ikT = 0; ikT < kTN; ikT++){
    hSIGNALkT[ikT] = new TH1D(Form("Signal_function_proton-lambda_kT_%d", ikT), Form("Signal_function_proton-lambda_kT_%d; k*[MeV/c]; C(k*)",ikT), 1000,0,2000);
    hBACKGROUNDkT[ikT] = new TH1D(Form("Background_function_proton-lambda_kT_%d", ikT), Form("Background_function_proton-lambda_kT_%d; k*[MeV/c]; C(k*)",ikT), 1000, 0, 2000);
    hSIGNALkT_Purity[ikT] = new TH1D(Form("Signal_function_proton-lambda_kT_%d_Purity", ikT), Form("Signal_function_proton-lambda_kT_%d_Purity; k*[MeV/c]; C(k*)",ikT), 1000,0,2000);
    hBACKGROUNDkT_Purity[ikT] = new TH1D(Form("Background_function_proton-lambda_kT_%d_Purity", ikT), Form("Background_function_proton-lambda_kT_%d_Purity; k*[MeV/c]; C(k*)",ikT), 1000, 0, 2000);
  }  

  TH1D *hSIGNAL_Cent_kT[CentN][kTN];
  TH1D *hBACKGROUND_Cent_kT[CentN][kTN];

  for( int iC = 0; iC < CentN; iC++){
    for( int iCkT = 0; iCkT < kTN; iCkT++){
      hSIGNAL_Cent_kT[iC][iCkT] = new TH1D(Form("Signal_function_proton-lambda_c%d_kT%d", iC, iCkT), Form("Signal_function_proton-lambda_c%d_kT%d; k*[MeV/c]; C(k*)", iC, iCkT), 1000,0,2000);
      hBACKGROUND_Cent_kT[iC][iCkT] = new TH1D(Form("Background_function_proton-lambda_c%d_kT%d", iC, iCkT), Form("Background_function_proton-lambda_c%d_kT%d; k*[MeV/c]; C(k*)", iC, iCkT), 1000, 0, 2000);
    }  
  }

  //--------------------------------------------------------------------------------
  // Creating and initializing the track sorter and a simple stopwatch object
  //--------------------------------------------------------------------------------
  HParticleTrackSorter sorter;
  sorter.init();

  TStopwatch timer;
  timer.Reset();
  timer.Start();

  //--------------------------------------------------------------------------------
  // Accounting energy Loss Correction
  //--------------------------------------------------------------------------------
  HEnergyLossCorrPar enLossCorr;
  enLossCorr.setDefaultPar("mar19");

  //--------------------------------------------------------------------------------
  //  Event chara and centrality
  //--------------------------------------------------------------------------------
  HParticleEvtChara evtChara;
  TString ParameterfileCVMFS = "/lustre/nyx/hades/user/bkardan/param/centrality_epcorr_mar19_ag158ag_3200A_glauber_gen5_pass3_2021_08.root";
  if(!evtChara.setParameterFile(ParameterfileCVMFS)){
    cout << "Parameterfile not found !!! " << endl;
    return kFALSE;
  }
  if(!evtChara.init()) {
    cout << "HParticleEvtChara not init!!! " << endl;
    return kFALSE;
  }
  
  Int_t eCentEst   = HParticleEvtChara::kTOFRPC;
  Int_t eCentClass = HParticleEvtChara::k10;
  Int_t eEPcorr    = HParticleEvtChara::kDefault;
  cout << "\t selected EPcorrection method is:  "  << evtChara.getStringEventPlaneCorrection(eEPcorr) << endl;
  evtChara.printCentralityClass(eCentEst, eCentClass);

  int    nRPC,nTOF;
  //--------------------------------------------------------------------------------
  // The amount of events to be processed
  //--------------------------------------------------------------------------------
  Long64_t nEvents = gloop->getEntries();
  if (nEventsDesired < nEvents && nEventsDesired >= 0) nEvents = nEventsDesired;

  //--------------------------------------------------------------------------------
  // The global event loop which loops over all events in the DST files added to HLoop
  // The loop breaks if the end is reached
  //--------------------------------------------------------------------------------
  for (Long64_t event = 0; event < nEvents; event++) {
    if (gloop->nextEvent(event) <= 0) {
        cout << " Last event processed " << endl;
        break;
    }
    //--------------------------------------------------------------------------------
    // Just the progress of the analysis
    //--------------------------------------------------------------------------------
    HTool::printProgress(event, nEvents, 1, "Analyzed events: ");

    //--------------------------------------------------------------------------------
    // The Event Header Object containing general information
    //-------------------------------------------------------------------------------
    HEventHeader* event_header = gloop->getEventHeader();

    //--------------------------------------------------------------------------------
    // Getting the amount of tracks (Particle Candidates), the Particle event Info object and the reconstructed global event vertex
    //--------------------------------------------------------------------------------
    Int_t nTracks                   = particle_cand_cat->getEntries();
    HGeomVector EventVertex         = event_header->getVertexReco().getPos();
    particle_info = HCategoryManager::getObject(particle_info, particle_info_cat, 0); 

    double EventVertexZ = event_header->getVertexReco().getZ();
    //===================================================================================================================
    // --- EVENT LEVEL ANALYSIS ---
    //===================================================================================================================
    nRPC          = particle_info->getSumRpcMultHitCut();
    nTOF          = particle_info->getSumTofMultCut();  
    //===================================================================================================================  
    //--------------------------------------------------------------------------------
    // Discarding bad events with multiple criteria and counting amount of all / good events
    //--------------------------------------------------------------------------------
    hCounter->AddBinContent(cNumAllEvents);

    hCounter->AddBinContent(cNumSelectedEvents);

    int lambdaCounter = 0;

    Float_t  event_weight = evtChara.getEventWeight();   // event_weight dependent if PT2(down-scaled) or PT3
    Int_t centrality = evtChara.getCentralityClass(eCentEst, eCentClass); //1(0-10%) - 7(60-70%) ... 0:Overflow max:Underflow
    Int_t plate = getTargeMomLlateNumber(EventVertex.getZ()); //number between 1-15, -1 means outside the target
    if(centrality<1 || centrality>4 || plate<1) continue;

    // 10% Centrality-Classes:       1(0-10%) - 5(40-50%) ... 0:Overflow max:Underflow
    Float_t CentralityTOFRPC      = evtChara.getCentralityPercentile(eCentClass);

    hMultiplicity  ->  Fill(nRPC+nTOF);
    hZVertex       ->  Fill(EventVertex.Z());                                                                                                

    //--------------------------------------------------------------------------------
    // Resetting the track sorter and selecting hadrons ranked by Chi2 Runge Kutta
    //--------------------------------------------------------------------------------
    sorter.cleanUp();
    sorter.resetFlags(kTRUE, kTRUE, kTRUE, kTRUE);
    sorter.fill(HParticleTrackSorter::selectHadrons);

    sorter.selectBest(Particle::kIsBestRKSorter, Particle::kIsHadronSorter);
    //--------------------------------------------------------------------------------
    // Removing tracks from previous events from the lists
    //--------------------------------------------------------------------------------
    vDau1Tracks.clear();
    vDau2Tracks.clear();
    Deuteron.clear();

    int nTOFRPC;
    int lambdaCOUNTER = 0;

    TRandom2* mRand = new TRandom2;
    //--------------------------------------------------------------------------------
    // The loop over all tracks (Particle Candidates) in the current event
    //--------------------------------------------------------------------------------
    for (Int_t track = 0; track < nTracks; track++) {
      //--------******HParticleCand* particle_cand = (HParticleCand*) particle_cand_cat->getObject(track);

      HParticleCandSim* particle_cand = HCategoryManager::getObject(particle_cand, particle_cand_cat, track);
      nTOFRPC = particle_info->getSumRpcMultHitCut() + particle_info->getSumTofMultCut();

      //--------------------------------------------------------------------------------
      // Getting information on the current track (Not all of them might be used - commented to supress warnings)
      //--------------------------------------------------------------------------------

      Float_t dEdx = particle_cand -> getMdcdEdx();
      Short_t chargeCand = particle_cand -> getCharge();
      Float_t momCand = particle_cand -> getMomentum();
      Float_t betaCand = particle_cand ->getBeta();

      Float_t chi2Rk   = particle_cand->getChi2();
      Float_t mass     = particle_cand->getMass();
      Float_t theta = particle_cand->getTheta();
      Float_t phi = particle_cand->getPhi();
                               
      Short_t sysCand    = particle_cand -> getSystemUsed(); // detector: 0 RPC, 1 ToF
      Int_t   sector     = particle_cand -> getSector();
      Float_t RichInd    = particle_cand -> getRichInd();
      Float_t dTheta     = particle_cand -> getDeltaTheta();
      Float_t dPhi       = particle_cand -> getDeltaPhi();
      Float_t richPhi    = particle_cand -> getRichPhi();
      Float_t richTheta  = particle_cand -> getRichTheta();
      Float_t rapCand    = particle_cand -> Rapidity();   

      if(sysCand == 0)hBetaMomRPC->Fill(chargeCand*momCand, betaCand);
      if(sysCand == 1)hBetaMomToF->Fill(chargeCand*momCand, betaCand);

      //--------------------------------------------------------------------------------
      // Discarding all tracks that have been discarded by the track sorter and counting all / good tracks
      //--------------------------------------------------------------------------------
      if (!particle_cand->isFlagBit(Particle::kIsUsed))
      continue;

      hCounter->AddBinContent(cNumSelectedTracks);

      //================================================================================================================================================================
      // Selecting daugther 1 and 2 tracks, setting mass to nominal value and adding to the lists
      // Notice that the selection has to be exclusive
      //================================================================================================================================================================
      if (protonRPCCut->IsInside(chargeCand*momCand, betaCand) || protonToFCut->IsInside(chargeCand*momCand, betaCand)){ // RPC+ToF
        hCounter->AddBinContent(cNumDau1Tracks);

        particle_cand->calc4vectorProperties(HPhysicsConstants::mass(14));
        Double_t momCandproton = enLossCorr.getCorrMom(14,momCand,theta);
        particle_cand->setMomentum(momCandproton);

        vDau1Tracks.push_back(particle_cand);
      }

      else if (piNRPCCut->IsInside(chargeCand*momCand, betaCand) || piNToFCut->IsInside(chargeCand*momCand, betaCand)) { // RPC+ToF
        hCounter->AddBinContent(cNumDau2Tracks);

        particle_cand->calc4vectorProperties(HPhysicsConstants::mass(9));
        Double_t momCandpion = enLossCorr.getCorrMom(9,momCand,theta);
        particle_cand->setMomentum(momCandpion);

        vDau2Tracks.push_back(particle_cand);
      }

      // if(betaCand > 0.45 || momCand > 900) continue;

      if( ( mass > 1600 && mass < 2150) && (
      (sysCand == 0 && dRPCCut->IsInside(chargeCand*momCand, dEdx)) ||
      (sysCand == 1 && dToFCut->IsInside(chargeCand*momCand, dEdx)))){

        particle_cand->setMomentum(particle_cand->getCorrectedMomentumPID(45));
        particle_cand -> calc4vectorProperties(HPhysicsConstants::mass(45)); // d

        Float_t momCand = particle_cand -> P();
        Float_t pt = particle_cand -> Pt();
        Float_t rap    = particle_cand -> Rapidity();
        Deuteron.push_back(particle_cand);

        hTransMomentum->Fill(pt);
        hMomentum->Fill(momCand);
        hRapidity->Fill(rap);
        hDeuteronBetaMom->Fill(particle_cand->getMomentum(), particle_cand->getBeta());
        hDeuterondEdx->Fill(particle_cand->getMomentum(), particle_cand->getMdcdEdx());
      }
    } // End of track loop
    
    LambdaProton->Fill(lambdaCOUNTER);
    //--------------------------------------------------------------------------------
    // The analysis makes no sense if there is not at least one candidate track for both dauthers
    // Furthermore counting the amount of possible combinations (mother candidates)
    //--------------------------------------------------------------------------------
    if (vDau1Tracks.size() == 0 || vDau2Tracks.size() == 0) continue;

    hCounter->AddBinContent(cNumMotCan, vDau1Tracks.size() * vDau2Tracks.size());

    //--------------------------------------------------------------------------------
    // Looping over candiates for first daugther track and calculating base and direction vector
    //--------------------------------------------------------------------------------
    for (vector<HParticleCandSim*>::iterator itDau1 = vDau1Tracks.begin(); itDau1 != vDau1Tracks.end(); itDau1++) {
      HParticleCandSim* hpcDau1 = *itDau1;
      TLorentzVector P11 = *hpcDau1;	    
      HGeomVector hgvBaseDau1, hgvDirDau1;
      HParticleTool::calcSegVector(hpcDau1->getZ(), hpcDau1->getR(), (TMath::DegToRad()*hpcDau1->getPhi()), (TMath::DegToRad()*hpcDau1->getTheta()), hgvBaseDau1, hgvDirDau1);

      //=====================================================================================================================================
      // Distance of closest approach of first daugther track to the global event vertex
      //=====================================================================================================================================
      Double_t VerDistA = HParticleTool::calculateMinimumDistanceStraightToPoint(hgvBaseDau1, hgvDirDau1, EventVertex);

      //--------------------------------------------------------------------------------
      // Looping over candiates for second daugther track and calculating base and direction vector
      //--------------------------------------------------------------------------------
      for (vector<HParticleCandSim*>::iterator itDau2 = vDau2Tracks.begin(); itDau2 != vDau2Tracks.end(); itDau2++) {
        HParticleCandSim* hpcDau2 = *itDau2;
        TLorentzVector P22 = *hpcDau2;
        HGeomVector hgvBaseDau2, hgvDirDau2;
        HParticleTool::calcSegVector(hpcDau2->getZ(), hpcDau2->getR(), (TMath::DegToRad()*hpcDau2->getPhi()), (TMath::DegToRad()*hpcDau2->getTheta()), hgvBaseDau2, hgvDirDau2);

        //====================================================================================================================================
        // Distance of closest approach of second daugther track to the global event vertex
        //====================================================================================================================================
        Double_t VerDistB = HParticleTool::calculateMinimumDistanceStraightToPoint(hgvBaseDau2, hgvDirDau2, EventVertex);

        //--------------------------------------------------------------------------------
        // Calculating the median point of closest approach of the two daugther tracks - Estimated decay vertex
        // Can also be used as the base vector of the mother particle
        //--------------------------------------------------------------------------------
        HGeomVector hgvDecayVertex = HParticleTool::calcVertexAnalytical(hgvBaseDau1, hgvDirDau1, hgvBaseDau2, hgvDirDau2);

        //====================================================================================================================================
        // Distance between secondary vertex and global event vertex - Decay length of mother particle
        //====================================================================================================================================
        Double_t VerDistX = (hgvDecayVertex - EventVertex).length();

        //--------------------------------------------------------------------------------
        // Opening angle between the two daugther tracks - Necessary for mixed event background estimation
        //--------------------------------------------------------------------------------
        Double_t Alpha = TMath::RadToDeg()*TMath::ACos(hgvDirDau1.scalarProduct(hgvDirDau2) / (hgvDirDau1.length() * hgvDirDau2.length()));

        //--------------------------------------------------------------------------------
        // Lorentz Vector of mother particle and corresponding direction vector
        //--------------------------------------------------------------------------------
        TLorentzVector tlvMot = *hpcDau1 + *hpcDau2;
        TLorentzVector tlvMot1 = *hpcDau1;
        TLorentzVector tlvMot2 = *hpcDau2;
        HGeomVector hgvDirMot(tlvMot.Px(), tlvMot.Py(), tlvMot.Pz());

        //================================================================================================================================================================
        // Distance of closest approach of reconstructed mother track to the global event vertex
        //================================================================================================================================================================
        Double_t VerDistC = HParticleTool::calculateMinimumDistanceStraightToPoint(hgvDecayVertex, hgvDirMot, EventVertex);

        //====================================================================================================================================
        // Distance of closest approach between the two daugther tracks
        //====================================================================================================================================
        Double_t MinTrackDist = HParticleTool::calculateMinimumDistance(hgvBaseDau1, hgvDirDau1, hgvBaseDau2, hgvDirDau2);

        NParticle proton_DCA  =  NParticle();
        FillParticleStruc(proton_DCA, event, P11, nTOFRPC, EventVertexZ, EventVertex.X(), EventVertex.Y(), centrality);

        NParticle pion_DCA  =  NParticle();
        FillParticleStruc(pion_DCA, event, P11, nTOFRPC, EventVertexZ, EventVertex.X(), EventVertex.Y(), centrality);

        //====================================================================================================================================
        // Mass of the reconstructed mother track
        //====================================================================================================================================
        Double_t MMot = tlvMot.M();

        //----------Mass
        hProtonBetaMom->Fill(hpcDau1->getMomentum(), hpcDau1->getBeta());
        hPionNBetaMom->Fill(hpcDau2->getMomentum(), hpcDau2->getBeta());
        //----------Rapidity
        hProtonRap->Fill(tlvMot1.Rapidity()); 
        hPionRap->Fill(tlvMot2.Rapidity());
        //----------Transverse Momentum
        hProtonPt->Fill(tlvMot1.Pt()); 
        hPionPt->Fill(tlvMot2.Pt());
        //----------Momentum
        hProtonP->Fill(tlvMot1.P()); 
        hPionP->Fill(tlvMot2.P());


        //====================================================================================================================================
        // Armenteros – Podolanski plot parametrs
        //====================================================================================================================================
        TVector3 beta (0., 0., 0.99);
        TLorentzVector tlvMot1_AP = tlvMot1;
        TLorentzVector tlvMot2_AP = tlvMot2;
        tlvMot1_AP.Boost(beta);
        tlvMot2_AP.Boost(beta);

        TLorentzVector tlvMot_AP = tlvMot1_AP + tlvMot2_AP;

        TVector3 vecN = tlvMot2_AP.Vect();          // Negetive particles kinematics 
        TVector3 vecP = tlvMot1_AP.Vect();          // Positive particles kinematics 
        TVector3 vecM = tlvMot_AP.Vect();           // (Positive+Negetive) mother particles kinematics 

        Double_t thetaP = TMath::ACos( ( vecP*vecM ) / ( vecP.Mag()*vecM.Mag() ) ) ;
        Double_t thetaN = TMath::ACos( ( vecN*vecM ) / ( vecN.Mag()*vecM.Mag() ) ) ;

        Double_t alfa = ( ( vecP.Mag() ) * TMath::Cos(thetaP) - ( vecN.Mag() ) * TMath::Cos(thetaN) ) / ( ( vecP.Mag() ) * TMath::Cos(thetaP) + ( vecN.Mag() ) * TMath::Cos(thetaN) );
        Double_t qt = vecP.Mag() * TMath::Sin(thetaP);
        //cout << "alfa : " <<alfa<<"   "<<" qt : "<<qt<<endl;

        //================================================================================================================================================================
        // Put your analyses of the mother particles here
        //================================================================================================================================================================
        if(VerDistX > 65 && VerDistA > 8 && VerDistB > 24 && VerDistC < 5 && MinTrackDist < 6 && Alpha > 15){//---------Dr. Severus Snape
          if(find(usedProtonRPC.begin(), usedProtonRPC.end(), tlvMot1) != usedProtonRPC.end()) continue; //do not use pip twice
          if(find(usedPionNRPC.begin(), usedPionNRPC.end(), tlvMot2) != usedPionNRPC.end()) continue;

          double OA = R2D*openingangle(tlvMot1,tlvMot);
          hPPiOA->Fill(OA);
          NParticle selectedPart_PionsP  =  NParticle();
          FillParticleStruc(selectedPart_PionsP, event, tlvMot1, nTOFRPC, EventVertexZ, EventVertex.X(), EventVertex.Y(), centrality);
          ProtonL.push_back(selectedPart_PionsP);

          NParticle selectedPart_PionsN  =  NParticle();
          FillParticleStruc(selectedPart_PionsN, event, tlvMot2, nTOFRPC, EventVertexZ, EventVertex.X(), EventVertex.Y(), centrality);
          PionL.push_back(selectedPart_PionsN);

          NParticle p1 = selectedPart_PionsP;
          NParticle p2 = selectedPart_PionsN;
          double LambdaMassB = p1.mass + p2.mass;  

          hLambdamassRPC->Fill(MMot);
          hsigAlphaPt_NoMC->Fill(alfa, qt);


          if(MMot > 1100 && MMot < 1130) hsigAlphaPt_MC->Fill(alfa, qt);

          if(MMot > 1110 && MMot < 1120){
            usedProtonRPC.push_back(tlvMot1);
            usedPionNRPC.push_back(tlvMot2);


            //LAMBDA RESOLUTION
            Float_t MassL = tlvMot.M();
            Float_t MomL = tlvMot.P();
            Float_t ThetaL = tlvMot.Theta()*TMath::DegToRad();
            Float_t PhiL = tlvMot.Phi()*TMath::DegToRad();

            double per = par_sigma_lambda.PResA + par_sigma_lambda.PResB * TMath::Log(MomL) + par_sigma_lambda.PResC * MomL + par_sigma_lambda.PResD * MomL * MomL + par_sigma_lambda.PResE * MomL * MomL * MomL;
            double thetaer = par_sigma_lambda.ThetaA + par_sigma_lambda.ThetaB / ThetaL + par_sigma_lambda.ThetaC / ThetaL / ThetaL + par_sigma_lambda.ThetaD / ThetaL / ThetaL / ThetaL + par_sigma_lambda.ThetaE  / ThetaL / ThetaL / ThetaL / ThetaL;
            double phier = par_sigma_lambda.PhiA + par_sigma_lambda.PhiB / PhiL + par_sigma_lambda.PhiC / PhiL / PhiL + par_sigma_lambda.PhiD * PhiL;

            double per_mean = par_mean_lambda.PResA + par_mean_lambda.PResB * MomL + par_mean_lambda.PResC * MomL;
            double thetaer_mean = par_mean_lambda.ThetaA + par_mean_lambda.ThetaB * ThetaL;
            double phier_mean = par_mean_lambda.PhiA + par_mean_lambda.PhiB * PhiL;

            MomL = MomL + mRand->Gaus(per_mean,per);
            ThetaL = ThetaL + mRand->Gaus(thetaer_mean,thetaer);
            PhiL = PhiL + mRand->Gaus(phier_mean,phier);

            double px_L = MomL * TMath::Cos(PhiL) * TMath::Sin(ThetaL);
            double py_L = MomL * TMath::Sin(PhiL) * TMath::Cos(ThetaL);
            double pz_L = MomL * TMath::Cos(ThetaL);
            double e_L = TMath::Sqrt(px_L*px_L + py_L*py_L + pz_L*pz_L + MassL*MassL);

            tlvMot.SetPx(px_L); 
            tlvMot.SetPy(py_L); 
            tlvMot.SetPz(pz_L); 

            // LAMBDA RESOLUTION END

            NParticle selectedPart_Lambda  =  NParticle();
            FillParticleStruc(selectedPart_Lambda, event, tlvMot, nTOFRPC, EventVertexZ, EventVertex.X(), EventVertex.Y(), centrality);
            HBTLambdaRPC_NPart.push_back(selectedPart_Lambda);
            hCounter->AddBinContent(cNumMotRec);

            

            //----------Mass
            hLambdaMass->Fill(tlvMot.M());
            //----------Rapidity
            hLambdaRap->Fill(tlvMot.Rapidity());
            //----------Transverse Momentum
            hLambdaPt->Fill(tlvMot.Pt());
            //----------Momentum
            hLambdaP->Fill(tlvMot.P());

            if(hpcDau1->getGeantPID()==14 && hpcDau2->getGeantPID()==9){ //checking if we have true proton & pion

              if(hpcDau1->getGeantParentPID()==18 && hpcDau2->getGeantParentPID()==18 && (hpcDau1->getGeantParentTrackNum()==hpcDau2->getGeantParentTrackNum())){
                Int_t parentTrack = hpcDau1->getGeantParentTrackNum();
                // HGeantKine* kineCand = HCategoryManager::getObject(kineCand, kine_cand_cat, parentTrack);
                // if(kineCand->isPrimary()){
                  NParticle selectedPart_Lambda_PURITY  =  NParticle();
                  FillParticleStruc(selectedPart_Lambda_PURITY, event, tlvMot, nTOFRPC, EventVertexZ, EventVertex.X(), EventVertex.Y(), centrality);
                  HBTLambdaRPC_NPart_Purity.push_back(selectedPart_Lambda_PURITY); 

                  hLambdamassRPC_checked->Fill(MMot);
                // }			
              }
            }

            hsigAlphaPt_MC2->Fill(alfa, qt);

            lambdaCounter++;
          } //Lambda invariant mass cut loop
        } // DCA cuts
      } // End of second dauther track loop
    } // first daugther track loop


    hMultiLambdaCounter->Fill(lambdaCounter);
    hCounterparticle->Fill(lambdaCounter);

    for(vector<HParticleCandSim*>::iterator itDeut = Deuteron.begin(); itDeut != Deuteron.end(); itDeut++) {

      HParticleCandSim* hpcDeut = *itDeut;

      //DEUTERON RESOLUTION
      Float_t MassD = hpcDeut->getMass();
      Float_t MomD = hpcDeut->getMomentum();
      Float_t ThetaD = hpcDeut->getTheta()*TMath::DegToRad();
      Float_t PhiD = hpcDeut->getPhi()*TMath::DegToRad();

      double per = par_sigma_deuteron.PResA + par_sigma_deuteron.PResB * TMath::Log(MomD) + par_sigma_deuteron.PResC * MomD + par_sigma_deuteron.PResD * MomD * MomD + par_sigma_deuteron.PResE * MomD * MomD * MomD;
      double thetaer = par_sigma_deuteron.ThetaA + par_sigma_deuteron.ThetaB / ThetaD + par_sigma_deuteron.ThetaC / ThetaD / ThetaD + par_sigma_deuteron.ThetaD / ThetaD / ThetaD / ThetaD + par_sigma_deuteron.ThetaE  / ThetaD / ThetaD / ThetaD / ThetaD;
      double phier = par_sigma_deuteron.PhiA + par_sigma_deuteron.PhiB/ PhiD + par_sigma_deuteron.PhiC / PhiD / PhiD + par_sigma_deuteron.PhiD * PhiD;

      double per_mean = par_mean_deuteron.PResA + par_mean_deuteron.PResB * MomD + par_mean_deuteron.PResC * MomD;
      double thetaer_mean = par_mean_deuteron.ThetaA + par_mean_deuteron.ThetaB * ThetaD;
      double phier_mean = par_mean_deuteron.PhiA + par_mean_deuteron.PhiB * PhiD;

      MomD = MomD + mRand->Gaus(per_mean,per);
      ThetaD = ThetaD + mRand->Gaus(thetaer_mean,thetaer);
      PhiD = PhiD + mRand->Gaus(phier_mean,phier);

      double px_D = MomD * TMath::Cos(PhiD) * TMath::Sin(ThetaD);
      double py_D = MomD * TMath::Sin(PhiD) * TMath::Cos(ThetaD);
      double pz_D = MomD * TMath::Cos(ThetaD);
      double e_D = TMath::Sqrt(px_D*px_D + py_D*py_D + pz_D*pz_D + MassD*MassD);
      double mom_D = TMath::Sqrt(px_D*px_D + py_D*py_D + pz_D*pz_D);

      hpcDeut->setMomentum(mom_D); 

      // DEUTERON RESOLUTION END

      NParticle selectedPart_Deuterons = NParticle();
      FillParticleStruc(selectedPart_Deuterons, event, *hpcDeut, nTOFRPC, EventVertexZ, EventVertex.X(), EventVertex.Y(), centrality);
      HBTDeuteron_NPart.push_back(selectedPart_Deuterons);

      if(hpcDeut->getGeantPID()==45){
        NParticle selectedPart_deuteron_PURITY  =  NParticle();
        FillParticleStruc(selectedPart_deuteron_PURITY, event, *hpcDeut, nTOFRPC, EventVertexZ, EventVertex.X(), EventVertex.Y(), centrality);
        Deuteron_NPart_Purity.push_back(selectedPart_deuteron_PURITY);

        hDeuteronBetaMom_checked->Fill(hpcDeut->getMomentum(), hpcDeut->getBeta());
        hDeuterondEdx_checked->Fill(hpcDeut->getMomentum(), hpcDeut->getMdcdEdx());
      }
    } // end for for Deuterons


    usedPionNRPC.clear();
    usedProtonRPC.clear();

  } // End of event loop

  // static ProcInfo_t info;
  // constexpr float toGB = 1.f/1024.f/1024.f;

  // gSystem->GeMomLrocInfo(&info);
  // std::cout << "\n---=== Memory Usage ===---\n";
  // std::cout << "resident memory used: " << info.fMemResident*toGB << " GB\t virtual memory used: " << info.fMemVirtual*toGB << " GB\n\n";


  MakeSignal(HBTLambdaRPC_NPart,HBTDeuteron_NPart,hsignalRPC);
  MakeBackground(HBTLambdaRPC_NPart,HBTDeuteron_NPart,hbackgroundRPC,maxMixEvents);

  // MakeSignal_Purity(HBTLambdaRPC_NPart_Purity,Deuteron_NPart_Purity,hsignal_Purity);
  // MakeBackground_Purity(HBTLambdaRPC_NPart_Purity,Deuteron_NPart_Purity,hbackground_Purity,maxMixEvents);

  // MakeSignal_cent(HBTLambdaRPC_NPart,HBTDeuteron_NPart,hSIGNALRPC);
  // MakeBackground_cent(HBTLambdaRPC_NPart,HBTDeuteron_NPart,hBACKGROUNDRPC,maxMixEvents);

  // MakeSignal_kT(HBTLambdaRPC_NPart,HBTDeuteron_NPart,hSIGNALkT);
  // MakeBackground_kT(HBTLambdaRPC_NPart,HBTDeuteron_NPart,hBACKGROUNDkT,maxMixEvents);

  // MakeSignal_cent_Purity(HBTLambdaRPC_NPart_Purity,Deuteron_NPart_Purity,hSIGNALRPC_Purity);
  // MakeBackground_cent_Purity(HBTLambdaRPC_NPart_Purity,Deuteron_NPart_Purity,hBACKGROUNDRPC_Purity,maxMixEvents);

  // MakeSignal_kT_Purity(HBTLambdaRPC_NPart_Purity,Deuteron_NPart_Purity,hSIGNALkT_Purity);
  // MakeBackground_kT_Purity(HBTLambdaRPC_NPart_Purity,Deuteron_NPart_Purity,hBACKGROUNDkT_Purity,maxMixEvents);

  sorter.finalize();

  HBTLambdaRPC_NPart.clear();
  HBMomLrotonRPC_NPart.clear();

  HBTLambdaRPC_NPart_Purity.clear();
  Deuteron_NPart_Purity.clear();

  timer.Stop();
  cout << "Finished DST processing" << endl;

  TFile* out = new TFile(outFile.Data(), "RECREATE");
  out->cd();


  hCounter->Write();
  hLambdamassRPC->Write();

  hCounterparticle->Write();
  hMultiLambdaCounter->Write();

  hsignalRPC->Write();
  hbackgroundRPC->Write();

  hsignal_Purity->Write();
  hbackground_Purity->Write();

  hsigAlphaPt_NoMC->Write();
  hsigAlphaPt_MC->Write();
  hsigAlphaPt_MC2->Write();


  // for(int icount = 0; icount < 3; icount++){
  //   hSIGNALRPC[icount]->Write();
  //   hBACKGROUNDRPC[icount]->Write();
  //   hSIGNALRPC_Purity[icount]->Write();
  //   hBACKGROUNDRPC_Purity[icount]->Write();
  // }

  // for( int ikT = 0; ikT < 5; ikT++){
  //   hSIGNALkT[ikT]->Write();
  //   hBACKGROUNDkT[ikT]->Write();
  //   hSIGNALkT_Purity[ikT]->Write();
  //   hBACKGROUNDkT_Purity[ikT]->Write();
  // }  


  hProtonBetaMom->Write();
  hPionNBetaMom->Write();

  //----------Mass
  hProtonMass->Write(); 
  hPionMass->Write();
  hLambdaMass->Write();
  //----------Rapidity
  hProtonRap->Write(); 
  hPionRap->Write();
  hLambdaRap->Write();
  //----------Transverse Momentum
  hProtonPt->Write(); 
  hPionPt->Write();
  hLambdaPt->Write();
  //----------Momentum
  hProtonP->Write(); 
  hPionP->Write();
  hLambdaP->Write();

  hBetaMomToF->Write();
  hBetaMomRPC->Write();
  hMomentum->Write();
  hTransMomentum->Write();
  hRapidity->Write();

  hMultiplicity  ->  Write();
  hZVertex       ->  Write();

  hDeuteronBetaMom->Write();
  hDeuterondEdx->Write();
  hLambdamassRPC_checked->Write();

  hDeuteronBetaMom_checked->Write();
  hDeuterondEdx_checked->Write();

  out->Save();
  out->Close();

  return 0 ;
}