#ifndef garutils_H
#define garutils_H

//#include "TTreeStream.h"
#include "TStopwatch.h"
#include "TMatrix.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "TChain.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TPad.h"
#include "TCanvas.h"
#include "AliPID.h"
#include "fastSimulation.h"
#include "TDecompLU.h"
#include <iostream>
#include <string>
#include <fstream>

Int_t BuildParticlePoints(fastParticle &particle, double Center[3], fastGeometry geom, 
                   std::vector<TVector3> trajxyz, long PDGcode, double displaceX, double displaceY)
{
          uint fMaxLayer = 0;
          TParticlePDG *p = TDatabasePDG::Instance()->GetParticle(PDGcode);
          if (p == nullptr) {
            ::Error("fastParticle::simulateParticle", "Invalid pdgCode %ld", PDGcode);
            return -1;
          }
          //Short_t sign = 1 * p->Charge() / 3.;
          Float_t mass = p->Mass();
          particle.fMassMC=mass;
          size_t count = 0;

          for(size_t k=0;k<trajxyz.size();k++) 
          {             
              //Bool_t invert = kFALSE;
              Double_t xyz_conv[3]= {(trajxyz.at(k).Z()-(Center[2]-displaceX)),
                                (trajxyz.at(k).Y()-(Center[1]-displaceY)),
                                trajxyz.at(k).X()-Center[0]};


              Double_t alpha=TMath::ATan2(xyz_conv[1],xyz_conv[0]);
              Double_t radius=sqrt((xyz_conv[1]-displaceY)*(xyz_conv[1]-displaceY)+(xyz_conv[0]-displaceX)*(xyz_conv[0]-displaceX));
              Double_t X_loc =  xyz_conv[0]*cos(alpha) + xyz_conv[1]*sin(alpha);
              Double_t Y_loc = -xyz_conv[0]*sin(alpha) + xyz_conv[1]*cos(alpha);
              //double covar[21]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
              double ptemp[5] = {Y_loc,xyz_conv[2],0,0,0};
              particle.fDirection.resize(count+1);
              particle.fParamMC.resize(count+1);
              particle.fLayerIndex.resize(count+1);
              particle.fParamMC[count].SetParamOnly(X_loc,alpha,ptemp);
              uint indexR = uint(std::upper_bound (geom.fLayerRadius.begin(),geom.fLayerRadius.end(), radius)-geom.fLayerRadius.begin());
              particle.fLayerIndex[count] = indexR;
              if(k!=0)
              {
                double x_now = particle.fParamMC[count].GetX();
                double x_prev = particle.fParamMC[count-1].GetX();
                if((x_now/x_prev)>1) particle.fDirection[count] = +1;
                else particle.fDirection[count] = -1;
              }
              else particle.fDirection[count] = +1;

              if (indexR>fMaxLayer) fMaxLayer=indexR;

              count++;
          }
          if(particle.fDirection.size()>1) particle.fDirection[0]=particle.fDirection[1];

          return 1;
}

void ExpandMatrix(float mA[15], float mB[25]){
  
  TMatrixF m(5,5);
  size_t counter = 0;
  for(size_t i = 0; i<5; i++){
    for(size_t j = 0; j<=i; j++){
        m(i,j) = mA[counter];
        m(j,i) = mA[counter];
        counter++;
    }
  }

  //m.Print();

  counter = 0;
  for(size_t i = 0; i<5; i++){
    for(size_t j = 0; j<5; j++){
      mB[counter] = m(i,j);
      counter++;
    }
  }

}

#endif