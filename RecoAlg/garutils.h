#ifndef garutils_H
#define garutils_H

#include "TTreeStream.h"
#include "TStopwatch.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "TChain.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TPad.h"
#include "TCanvas.h"
#include "AliPID.h"
#include "TDecompLU.h"
#include <iostream>
#include <string>
#include <fstream>



void idx2ijsort2(size_t n, size_t k, size_t &i, size_t &j)
{
  if (n<2)
    {
      i=0;
      j=0;
      return;
    }
  if (k > n*(n-1)/2)
    {
      throw std::runtime_error(Form("idx2ijsort2: k too big: %lu %lu",k,n)) ;
    }
  i = n - 2 - TMath::Floor(TMath::Sqrt( (double) (-8*k + 4*n*(n-1)-7) )/2.0 - 0.5);
  j = k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2;
}

size_t ij2idxsort2(size_t n, size_t i_in, size_t j_in)
{
  if (i_in >= n)
    {
      throw std::runtime_error("ij2idxsort2 i_in >= n");
    }
  if (j_in >= n)
    {
      throw std::runtime_error("ij2idxsort2 j_in >= n");
    }
  size_t i = TMath::Min(i_in,j_in);
  size_t j = TMath::Max(i_in,j_in);
  size_t k = (n*(n-1)/2) - (n-i)*((n-i)-1)/2 + j - i - 1;
  return k;
}


bool sort2_check_cyclic(std::vector<int> &link1, std::vector<int> &link2, int i, int j)
{
  if (i == j) return true;

  int prev = i;
  int next = (link1.at(i) == -1) ? link2.at(i) : link1.at(i);
  while (next != -1)
    {
      if (next == j) return true;
      int ptmp = next;
      next = (link1.at(next) == prev) ? link2.at(next) : link1.at(next);
      prev = ptmp;
    }
  return false;
}


void sort_TPCClusters_along_track2(const std::vector<TVector3>  &TPCClusters,
                                        std::vector<int> &hlf,
                                        std::vector<int> &hlb,
                                        int /* printlevel */,
                                        float &lengthforwards,
                                        float &lengthbackwards,
					     float dcut)
{
  size_t nclus = TPCClusters.size();
  hlf.clear();
  hlb.clear();
  lengthforwards = 0;
  lengthbackwards = 0;
  if (nclus == 0) return;
  if (nclus == 1)
    {
      hlf.push_back(0);
      hlb.push_back(0);
      return;
    }

  size_t ndists = nclus*(nclus-1)/2;
  std::vector<float> dists(ndists);
  std::vector<int> dsi(ndists);
  std::vector<int> link1(nclus,-1);
  std::vector<int> link2(nclus,-1);
  for (size_t i=0; i<nclus-1; ++i)
    {
      TVector3 iclus(TPCClusters.at(i));
      for (size_t j=i+1; j<nclus; ++j)
        {
          TVector3 jclus(TPCClusters.at(j));
          float d = (iclus-jclus).Mag();
          dists.at(ij2idxsort2(nclus,i,j)) = d;
        }
    }
  TMath::Sort( (int) ndists, dists.data(), dsi.data(), false );
  //std::cout << "Dists sorting: " << std::endl;
  //for (size_t i=0; i<ndists; ++i)
  //  {
  //    std::cout << i << " " << dists[i] << " " << dists[dsi[i]] << std::endl;
  //  }

  // greedy pairwise association -- may end up leaving some hits stranded

  for (size_t k=0; k<ndists; ++k)
    {
      size_t i = 0;
      size_t j = 0;
      size_t k2 = dsi.at(k);
      if (dists.at(k2) > dcut) break;
      idx2ijsort2(nclus,k2,i,j);
      TVector3 iclus(TPCClusters.at(i));
      TVector3 jclus(TPCClusters.at(j));

      if (link1.at(i) == -1)
        {
          if (link1.at(j) == -1)
            {
              link1.at(i) = j;
              link1.at(j) = i;
            }
          else  // already have a link on j'th cluster, add a second one only if it's on the other side
            {
              if (link2.at(j) == -1)
                {
                  TVector3 j2clus(TPCClusters.at(link1.at(j)));

		  // old test on angle instead of cyclicness
                  //if ( (j2clus-jclus).Angle(iclus-jclus) > TMath::TMath::Pi()/2)

		  if (!sort2_check_cyclic(link1,link2,j,i))
                    {
                      link1.at(i) = j;
                      link2.at(j) = i;
                    }
                }
            }
        }
      else
        {
          if (link2.at(i) == -1)
            {
              if (link1.at(j) == -1)
                {
                  link2.at(i) = j;
                  link1.at(j) = i;
                }
              else  // already have a link on j'th cluster, add a second one only if it's on the other side
                {
                  if (link2.at(j) == -1)
                    {
                      TVector3 j2clus(TPCClusters.at(link1.at(j)));

		      // old test on angle instead of cyclicness
                      //if ( (j2clus-jclus).Angle(iclus-jclus) > TMath::TMath::Pi()/2)

		      if (!sort2_check_cyclic(link1,link2,j,i))
                        {
                          link2.at(i) = j;
                          link2.at(j) = i;
                        }
                    }
                }
            }
        }
    }

  //std::cout << "sorting clusters and links" << std::endl;
  //for (size_t i=0; i<nclus; ++i)
  //  {
  //    std::cout << i << " " << link1.at(i) << " " << link2.at(i) << " " << 
  //	TPCClusters.at(i).Position()[0] << " " <<
  // 	TPCClusters.at(i).Position()[1] << " " <<
  //	TPCClusters.at(i).Position()[2] << std::endl;
  //}

  std::vector<int> used(nclus,-1);
  std::vector<std::vector<int> > clistv;

  // find the first tpc cluster with just one link

  int ifirst = -1;
  for (size_t i=0; i<nclus; ++i)
    {
      int il1 = link1.at(i);
      int il2 = link2.at(i);

      if ( (il1 == -1 && il2 != -1) || (il1 != -1 && il2 == -1) ) 
	{
	  ifirst = i;
	  break;
	}
    }
  // this can happen if all TPC clusters are singletons -- no links anywhere.  Happens if
  // there are no TPC clusters at all or just one.
  if (ifirst == -1)
    {
      return;
    }

  // follow the chains

  for (size_t i=ifirst; i<nclus; ++i)
    {
      if (used.at(i) == 1) continue;
      std::vector<int> clist;
      clist.push_back(i);
      used.at(i) = 1;

      int idx = i;
      int il1 = link1.at(idx);
      int il2 = link2.at(idx);
      
      int u1 = -2;
      if (il1 != -1) u1 = used.at(il1);
      int u2 = -2;
      if (il2 != -1) u2 = used.at(il2);

      while (u1 == -1 || u2 == -1)
	{
	  idx =  (u1 == -1) ?  il1 : il2;
	  clist.push_back(idx);
	  used.at(idx) = 1;
          il1 = link1.at(idx);
          il2 = link2.at(idx);
          u1 = -2;
          if (il1 != -1) u1 = used.at(il1);
          u2 = -2;
          if (il2 != -1) u2 = used.at(il2);
	}
      clistv.push_back(clist);
    } 

  // temporary solution -- think about what to do.  Report the longest chain

  size_t msize=0;
  int ibest=-1;
  for (size_t ichain=0; ichain<clistv.size(); ++ichain)
    {
      size_t cursize = clistv.at(ichain).size();
      if (cursize > msize)
	{
	  ibest = ichain;
	  msize = cursize;  
	}
    }
  if (ibest == -1) return;

  for (size_t i=0; i<clistv.at(ibest).size(); ++i)
    {
      hlf.push_back(clistv[ibest].at(i));
      if (i>0)
	{
	  TVector3 lastpoint(TPCClusters.at(hlf.at(i-1)));
	  TVector3 curpoint(TPCClusters.at(hlf.at(i)));
	  lengthforwards += (curpoint-lastpoint).Mag();
	}
    } 
  for (size_t i=0; i< hlf.size(); ++i)
    {
      hlb.push_back(hlf[hlf.size()-1-i]);  // just invert the order for the backward sort
    }
  lengthbackwards = lengthforwards;
 
}
/*
bool GetLArm(fastParticle part, double cutoff=75){
  double yend = part.fParamMC[part.fParamMC.size()-1].GetX()*sin(part.fParamMC[part.fParamMC.size()-1].GetAlpha())+part.fParamMC[part.fParamMC.size()-1].GetParameter()[0]*cos(part.fParamMC[part.fParamMC.size()-1].GetAlpha());
  double ystart = part.fParamMC[0].GetX()*sin(part.fParamMC[0].GetAlpha())+part.fParamMC[0].GetParameter()[0]*cos(part.fParamMC[0].GetAlpha());
  double xend = part.fParamMC[part.fParamMC.size()-1].GetX()*cos(part.fParamMC[part.fParamMC.size()-1].GetAlpha())-part.fParamMC[part.fParamMC.size()-1].GetParameter()[0]*sin(part.fParamMC[part.fParamMC.size()-1].GetAlpha());
  double xstart = part.fParamMC[0].GetX()*cos(part.fParamMC[0].GetAlpha())-part.fParamMC[0].GetParameter()[0]*sin(part.fParamMC[0].GetAlpha());
  double lArmMC = sqrt((xend-xstart)*(xend-xstart)+(yend-ystart)*(yend-ystart));
  return lArmMC;
}

Double_t GetCenterR(Double_t x1,Double_t y1, Double_t x2,Double_t y2, Double_t x3,Double_t y3, Double_t &x0, Double_t &y0, Double_t &r){
  //-----------------------------------------------------------------
  // Initial approximation of the track curvature
  //-----------------------------------------------------------------
  x3 -=x1;
  x2 -=x1;
  y3 -=y1;
  y2 -=y1;
  //  
  Double_t det = x3*y2-x2*y3;
  if (TMath::Abs(det)<1e-10){
    return 100;
  }
  //
  Double_t u = 0.5* (x2*(x2-x3)+y2*(y2-y3))/det;
  x0 = (x3*0.5-y3*u);
  y0 = (y3*0.5+x3*u);
  r = sqrt(x0*x0+y0*y0);
  x0 +=x1;
  y0 +=y1;
  return 1;
}

int longestStreak(const std::vector<float>& numbers, bool positive) {
    int longestStreak = 0;
    int currentStreak = 0;
    
    for (int num : numbers) {
        if ((positive && num > 0) || (!positive && num < 0)) {
            // Current number matches the desired sign
            currentStreak++;
        } else {
            // Current number does not match the desired sign
            longestStreak = std::max(longestStreak, currentStreak);
            currentStreak = 0;
        }
    }
    
    // Check if the streak continues till the end of the vector
    longestStreak = std::max(longestStreak, currentStreak);
    
    return longestStreak;
}

bool CheckLooper(fastParticle particle){
  uint layer1 = uint(particle.fParamMC.size()-1);
  Double_t xyzS[3][3];
  Int_t step=layer1/3;
  if (step>20) step=20;

  for (int dLayer=0; dLayer<3  && layer1-step*dLayer>0; dLayer++) {
          Int_t index=layer1-step*dLayer-1;
          Int_t index0=layer1-1;
          //fParamMC[index1-step*dLayer-1].GetXYZ(xyzS[dLayer]);
          xyzS[dLayer][0]=particle.fParamMC[index].GetX();
          xyzS[dLayer][1]=particle.fParamMC[index].GetY();
          xyzS[dLayer][2]=particle.fParamMC[index].GetZ();
    }

  Double_t x0,y0,r = 0; 
  GetCenterR(xyzS[0][0],xyzS[0][1],xyzS[1][0],xyzS[1][1],xyzS[2][0],xyzS[2][1],x0,y0,r);
  double yrel0 = particle.fParamMC[0].GetY()-y0;
  double xrel0 = particle.fParamMC[0].GetX()-x0;
  double angle = atan2(yrel0,xrel0);
  std::vector<float> locangles;
  for(auto t: particle.fParamMC){
      double yrel = t.GetY()-y0;
      double xrel = t.GetX()-x0;
      Double_t X_loc =  xrel*cos(angle) + yrel*sin(angle);
      Double_t Y_loc = -xrel*sin(angle) + yrel*cos(angle);
      Double_t locangle = atan2(Y_loc,X_loc);
      locangles.push_back(locangle);
  }
  int pos = longestStreak(locangles,true);
  int neg = longestStreak(locangles,false);
  if(pos>10 && neg>10) return 1;
  else return 0;
}

double CheckLooperMC(fastParticle particle){
   double Length = particle.fLengthInRot;
   double Radius = 1/abs(particle.fParamMC[0].GetParameter()[4]*5*0.0003);
   double circ_frac = Length/(2*TMath::Pi()*Radius);
   return circ_frac;
}


void BuiltParticleUnsmeared(fastParticle &particle, double Center[3], std::vector<TVector3> trajxyz, std::vector<TVector3> trajxyz_uns,double displaceX, double displaceY)
{
  for(size_t k=0;k<trajxyz.size();k++) 
  {
    Double_t xyz_conv[3]= {trajxyz.at(k).Z()-(Center[2]-displaceX),
                              trajxyz.at(k).Y()-(Center[1]-displaceY),
                              trajxyz.at(k).X()-Center[0]};
    Double_t xyz_conv_uns[3]= {trajxyz_uns.at(k).Z()-(Center[2]-displaceX),
                              trajxyz_uns.at(k).Y()-(Center[1]-displaceY),
                              trajxyz_uns.at(k).X()-Center[0]};
    Double_t alpha=TMath::ATan2(xyz_conv[1],xyz_conv[0]);
    Double_t X_loc =  xyz_conv_uns[0]*cos(alpha) + xyz_conv_uns[1]*sin(alpha);
    Double_t Y_loc = -xyz_conv_uns[0]*sin(alpha) + xyz_conv_uns[1]*cos(alpha);
    Double_t Z_loc =  xyz_conv_uns[2];
    double ptemp[] = {Y_loc,Z_loc,0,0,0};
    particle.fParamMC.resize(k+1);
    particle.fParamMC[k].SetParamOnly(X_loc,alpha,ptemp);  
  }
}

void BuiltParticleUnsmearedNoRot(fastParticle &particle, double Center[3], std::vector<TVector3> trajxyz, std::vector<TVector3> trajxyz_uns, double alpha=0)
{
  for(size_t k=0;k<trajxyz.size();k++) 
  {
    Double_t xyz_conv[3]= {(trajxyz.at(k).Z()-Center[2]),
                              trajxyz.at(k).Y()-Center[1],
                              trajxyz.at(k).X()-Center[0]};
    Double_t xyz_conv_uns[3]= {(trajxyz_uns.at(k).Z()-Center[2]),
                              trajxyz_uns.at(k).Y()-Center[1],
                              trajxyz_uns.at(k).X()-Center[0]};
    //Double_t alpha=TMath::ATan2(xyz_conv[1],xyz_conv[0]);
    Double_t X_loc =  xyz_conv_uns[0]*cos(alpha) + xyz_conv_uns[1]*sin(alpha);
    Double_t Y_loc = -xyz_conv_uns[0]*sin(alpha) + xyz_conv_uns[1]*cos(alpha);
    Double_t Z_loc =  xyz_conv_uns[2];
    double ptemp[] = {Y_loc,Z_loc,0,0,0};
    particle.fParamMC.resize(k+1);
    particle.fParamMC[k].SetParamOnly(X_loc,alpha,ptemp);  
  }
}

size_t ClosestPoint(std::vector<TVector3> trajxyz, TVector3 ClusterXYZ)
{
  Double_t dist = 10000;
  size_t index=0;
  for(size_t i=0; i<trajxyz.size(); i++)
  {
    double cl[] = {ClusterXYZ.X(),ClusterXYZ.Y(),ClusterXYZ.Z()};
    double traj[] = {trajxyz.at(i).X(),trajxyz.at(i).Y(),trajxyz.at(i).Z()};
    TVector3 diff = trajxyz.at(i)-ClusterXYZ;
    Double_t checkdist = sqrt(diff.Mag2());
    //checkdist = sqrt(pow(ClusterXYZ.X()-trajxyz.at(i).X(),2)+pow(ClusterXYZ.Y()-trajxyz.at(i).Y(),2));
    if(checkdist<dist)
    {
      dist=checkdist;
      index=i;
    }

  }
  return index;
}

size_t ClosestPoint2D(std::vector<TVector3> trajxyz, TVector3 ClusterXYZ)
{
  Double_t dist = 10000;
  size_t index=0;
  for(size_t i=0; i<trajxyz.size(); i++)
  {
    double cl[] = {0,ClusterXYZ.Y(),ClusterXYZ.Z()};
    double traj[] = {0,trajxyz.at(i).Y(),trajxyz.at(i).Z()};
    TVector3 diff = trajxyz.at(i)-ClusterXYZ;
    Double_t checkdist = sqrt(diff.Mag2());
    //checkdist = sqrt(pow(ClusterXYZ.X()-trajxyz.at(i).X(),2)+pow(ClusterXYZ.Y()-trajxyz.at(i).Y(),2));
    if(checkdist<dist)
    {
      dist=checkdist;
      index=i;
    }

  }
  return index;
}

TVector3 FindPoca(TVector3 p, TVector3 xyz, TVector3 TPCxyz, Double_t extdist, size_t steps){
  TVector3 unitp = p.Unit();
  Double_t dist = sqrt((xyz-TPCxyz).Mag2());
  Double_t angle = TPCxyz.Angle(unitp);
  Double_t angleini = angle;
  Double_t distini= dist;
  TVector3 poca = xyz;
  for(size_t s=0; s<steps; s++){
    TVector3 nextpoint = xyz + unitp*s*(extdist/steps);
    Double_t distnext = sqrt((nextpoint-TPCxyz).Mag2());
    TVector3 distvectnext =  nextpoint-TPCxyz;
    Double_t anglenext = distvectnext.Angle(unitp);
    if(abs(anglenext-TMath::Pi()/2)<abs(angle-TMath::Pi()/2)){
      poca = nextpoint;
      dist = distnext;
      angle = anglenext;
    }
    TVector3 prevpoint = xyz - unitp*s*(extdist/steps);
    Double_t distprev = sqrt((prevpoint-TPCxyz).Mag2());
    TVector3 distvectprev =  prevpoint-TPCxyz;
    Double_t angleprev = distvectprev.Angle(unitp);
    if(abs(angleprev-TMath::Pi()/2)<abs(angle-TMath::Pi()/2)){
      poca = prevpoint;
      dist = distprev;
      angle = angleprev;
    }
  }
  return poca;
}


void ClosestTrajectory(std::vector<TVector3> trajxyz, std::vector<TVector3> trajpxyz, std::vector<TVector3> ClusterXYZ, 
                        std::vector<TVector3> &trajxyz_closest, std::vector<TVector3> &trajpxyz_closest,Bool_t &IsForward)
{
  size_t FirstTrajIndex,LastTrajIndex;
  for(size_t t=0;t<ClusterXYZ.size();t++)
  {
    size_t cl_index= ClosestPoint(trajxyz,ClusterXYZ[t]);
    Double_t distmax=1;
    if(cl_index>0) distmax = TMath::Max(sqrt((trajxyz[cl_index]-trajxyz[cl_index-1]).Mag2()),sqrt((trajxyz[cl_index]-trajxyz[cl_index+1]).Mag2()));
    else distmax = sqrt((trajxyz[cl_index]-trajxyz[cl_index+1]).Mag2());
    TVector3 poca = FindPoca(trajpxyz[cl_index],trajxyz[cl_index],ClusterXYZ[t],distmax,500);
    trajpxyz_closest.push_back(trajpxyz[cl_index]);
    trajxyz_closest.push_back(poca);
    if(t==0) FirstTrajIndex=ClosestPoint(trajxyz,ClusterXYZ[t]);
    if(t==(ClusterXYZ.size()-1)) LastTrajIndex=ClosestPoint(trajxyz,ClusterXYZ[t]);
  }
  if(LastTrajIndex>FirstTrajIndex) IsForward=kTRUE;
  else IsForward=kFALSE;
}

Int_t BuildParticleCov(fastParticle &particle, double Center[3], fastGeometry geom, 
                   std::vector<TVector3> trajxyz, std::vector<TVector3> trajpxyz, long PDGcode, std::vector<float> CovXX,std::vector<float> CovYY)
{
          uint fMaxLayer = 0;
          TParticlePDG *p = TDatabasePDG::Instance()->GetParticle(PDGcode);
          if (p == nullptr) {
            ::Error("fastParticle::simulateParticle", "Invalid pdgCode %ld", PDGcode);
            return -1;
          }
          Short_t sign = 1 * p->Charge() / 3.;
          Float_t mass = p->Mass();
          particle.fMassMC=mass;

          for(size_t k=0;k<trajxyz.size();k++) 
          {
              Bool_t invert = kFALSE;
              Double_t xyz_conv[3]= {(trajxyz.at(k).Z()-(Center[2]-278)),
                                trajxyz.at(k).Y()-Center[1],
                                trajxyz.at(k).X()-Center[0]};

              Double_t pxyz_conv[3]= {(trajpxyz.at(k).Z()),
                                      trajpxyz.at(k).Y(),
                                      trajpxyz.at(k).X()};

              Double_t alpha=TMath::ATan2(xyz_conv[1],xyz_conv[0]);
              Double_t radius=sqrt(xyz_conv[1]*xyz_conv[1]+(xyz_conv[0]-278)*(xyz_conv[0]-278));
              Double_t X_loc =  xyz_conv[0]*cos(alpha) + xyz_conv[1]*sin(alpha);
              Double_t Y_loc = -xyz_conv[0]*sin(alpha) + xyz_conv[1]*cos(alpha);
              //Double_t param[5]={Y_loc,xyz_conv.Z(),0,0,0};              
              //particle.fParamMC[k].SetParamOnly(X_loc,alpha,param);
              double covar[21]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
              AliExternalTrackParam param(xyz_conv,pxyz_conv,covar,sign);
              AliExternalTrackParam4D param4D(param,mass,1);
              Double_t pxyz_conv_test[3];
              param4D.GetPxPyPz(pxyz_conv_test);
              Bool_t status = param4D.Rotate(alpha);
              Double_t pxyz_inv[3] = {-pxyz_conv[0],-pxyz_conv[1],-pxyz_conv[2]};

              if(abs(pxyz_conv_test[0]-pxyz_conv[0])>0.00001)
              {
                AliExternalTrackParam param_inv(xyz_conv,pxyz_inv,covar,sign);
                AliExternalTrackParam4D param4D_inv(param_inv,mass,1);
                status = param4D_inv.Rotate(alpha);
                // Double_t pxyz_inv_test[3];
                // param4D_inv.GetPxPyPz(pxyz_inv_test);
                // param4D = param4D_inv;
                // Double_t pxyz_inv_test2[3];
                // param4D.GetPxPyPz(pxyz_inv_test2);
                invert=kTRUE;
              }

              particle.fDirection.resize(k+1);
              particle.fParamMC.resize(k+1);
              // particle.fResolRPhi.resize(k+1);
              // particle.fResolZ.resize(k+1);
              particle.fLayerIndex.resize(k+1);
              if(status)
              {
                if(!invert)
                {
                  particle.fParamMC[k].Set(xyz_conv,pxyz_conv,covar,sign);
                  particle.fParamMC[k].Rotate(alpha);
                  // particle.fResolRPhi[k]=CovYY[k];
                  // particle.fResolZ[k]=CovXX[k];
                }
                else
                {
                  particle.fParamMC[k].Set(xyz_conv,pxyz_inv,covar,-sign);
                  particle.fParamMC[k].Rotate(alpha);
                  // particle.fResolRPhi[k]=CovYY[k];
                  // particle.fResolZ[k]=CovXX[k];
                }
              }
              else
              {
                double ptemp[] = {Y_loc,xyz_conv[2],0,0,0};
                particle.fParamMC[k].SetParamOnly(X_loc,alpha,ptemp);
                // particle.fResolRPhi[k]=CovYY[k];
                // particle.fResolZ[k]=CovXX[k];
              }
              //particle.fParamMC[k].Rotate(alpha);
              uint indexR = uint(std::upper_bound (geom.fLayerRadius.begin(),geom.fLayerRadius.end(), radius)-geom.fLayerRadius.begin());
              particle.fLayerIndex[k] = indexR;
              if(k!=0)
              {
                // TVector3 xyz_prev(trajxyz.at(k-1).Z()-Center[2],
                //                   trajxyz.at(k-1).Y()-Center[1],
                //                   trajxyz.at(k-1).X()-Center[0]);
                // Double_t r_prev = sqrt(xyz_prev.Y()*xyz_prev.Y()+xyz_prev.X()*xyz_prev.X());
                double x_now = particle.fParamMC[k].GetX();
                double x_prev = particle.fParamMC[k-1].GetX();
                if((x_now/x_prev)>1) particle.fDirection[k] = +1;
                else particle.fDirection[k] = -1;
              }
              else particle.fDirection[k] = +1;

              if (indexR>fMaxLayer) fMaxLayer=indexR;
          }
          if(particle.fDirection.size()>1) particle.fDirection[0]=particle.fDirection[1];

          return 1;
}

Int_t BuildParticle(fastParticle &particle, double Center[3], fastGeometry geom, 
                   std::vector<TVector3> trajxyz, std::vector<TVector3> trajpxyz, long PDGcode, double displaceX, double displaceY, int Sampling = 1)
{
          uint fMaxLayer = 0;
          TParticlePDG *p = TDatabasePDG::Instance()->GetParticle(PDGcode);
          if (p == nullptr) {
            ::Error("fastParticle::simulateParticle", "Invalid pdgCode %ld", PDGcode);
            return -1;
          }
          Short_t sign = 1 * p->Charge() / 3.;
          Float_t mass = p->Mass();
          particle.fMassMC=mass;
          double outofTPC = 0;
          size_t count = 0;

          for(size_t k=0;k<trajxyz.size();k++) 
          {

              if(k % Sampling != 0) continue;

              
              Bool_t invert = kFALSE;
              Double_t xyz_conv[3]= {(trajxyz.at(k).Z()-(Center[2]-displaceX)),
                                (trajxyz.at(k).Y()-(Center[1]-displaceY)),
                                trajxyz.at(k).X()-Center[0]};

              Double_t pxyz_conv[3]= {(trajpxyz.at(k).Z()),
                                      trajpxyz.at(k).Y(),
                                      trajpxyz.at(k).X()};

              Double_t alpha=TMath::ATan2(xyz_conv[1],xyz_conv[0]);
              Double_t radius=sqrt((xyz_conv[1]-displaceY)*(xyz_conv[1]-displaceY)+(xyz_conv[0]-displaceX)*(xyz_conv[0]-displaceX));
              Double_t X_loc =  xyz_conv[0]*cos(alpha) + xyz_conv[1]*sin(alpha);
              Double_t Y_loc = -xyz_conv[0]*sin(alpha) + xyz_conv[1]*cos(alpha);
              if(radius>250 || abs(xyz_conv[2])>250) break;
              //if(radius>250) outofTPC+=1;
              //Double_t param[5]={Y_loc,xyz_conv.Z(),0,0,0};              
              //particle.fParamMC[k].SetParamOnly(X_loc,alpha,param);
              double covar[21]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
              AliExternalTrackParam param(xyz_conv,pxyz_conv,covar,sign);
              AliExternalTrackParam4D param4D(param,mass,1);
              Double_t pxyz_conv_test[3];
              param4D.GetPxPyPz(pxyz_conv_test);
              Bool_t status = param4D.Rotate(alpha);
              Double_t pxyz_inv[3] = {-pxyz_conv[0],-pxyz_conv[1],-pxyz_conv[2]};

              if(abs(pxyz_conv_test[0]-pxyz_conv[0])>0.00001)
              {
                AliExternalTrackParam param_inv(xyz_conv,pxyz_inv,covar,sign);
                AliExternalTrackParam4D param4D_inv(param_inv,mass,1);
                status = param4D_inv.Rotate(alpha);
                // Double_t pxyz_inv_test[3];
                // param4D_inv.GetPxPyPz(pxyz_inv_test);
                // param4D = param4D_inv;
                // Double_t pxyz_inv_test2[3];
                // param4D.GetPxPyPz(pxyz_inv_test2);
                invert=kTRUE;
              }

              particle.fDirection.resize(count+1);
              particle.fParamMC.resize(count+1);
              // particle.fResolRPhi.resize(count+1);
              // particle.fResolZ.resize(count+1);
              particle.fLayerIndex.resize(count+1);
              if(status)
              {
                if(!invert)
                {
                  particle.fParamMC[count].Set(xyz_conv,pxyz_conv,covar,sign);
                  particle.fParamMC[count].Rotate(alpha);
                }
                else
                {
                  particle.fParamMC[count].Set(xyz_conv,pxyz_inv,covar,-sign);
                  particle.fParamMC[count].Rotate(alpha);
                }
              }
              else
              {
                double ptemp[] = {Y_loc,xyz_conv[2],0,0,0};
                particle.fParamMC[count].SetParamOnly(X_loc,alpha,ptemp);
              }
              //particle.fParamMC[k].Rotate(alpha);
              uint indexR = uint(std::upper_bound (geom.fLayerRadius.begin(),geom.fLayerRadius.end(), radius)-geom.fLayerRadius.begin());
              particle.fLayerIndex[count] = indexR;
              if(k!=0)
              {
                // TVector3 xyz_prev(trajxyz.at(k-1).Z()-Center[2],
                //                   trajxyz.at(k-1).Y()-Center[1],
                //                   trajxyz.at(k-1).X()-Center[0]);
                // Double_t r_prev = sqrt(xyz_prev.Y()*xyz_prev.Y()+xyz_prev.X()*xyz_prev.X());
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
          // double tot = trajxyz.size();
          // outofTPC*=1./tot;
          // std::cout<<"Percentage out of bounds: "<<outofTPC*100<<" %"<<std::endl;

          return 1;
}
*/
Int_t BuildParticlePoints(fastParticle &particle, double Center[3], fastGeometry geom, 
                   std::vector<TVector3> trajxyz, long PDGcode, double displaceX, double displaceY)
{
          uint fMaxLayer = 0;
          TParticlePDG *p = TDatabasePDG::Instance()->GetParticle(PDGcode);
          if (p == nullptr) {
            ::Error("fastParticle::simulateParticle", "Invalid pdgCode %ld", PDGcode);
            return -1;
          }
          Short_t sign = 1 * p->Charge() / 3.;
          Float_t mass = p->Mass();
          particle.fMassMC=mass;
          double outofTPC = 0;
          size_t count = 0;

          for(size_t k=0;k<trajxyz.size();k++) 
          {             
              Bool_t invert = kFALSE;
              Double_t xyz_conv[3]= {(trajxyz.at(k).Z()-(Center[2]-displaceX)),
                                (trajxyz.at(k).Y()-(Center[1]-displaceY)),
                                trajxyz.at(k).X()-Center[0]};


              Double_t alpha=TMath::ATan2(xyz_conv[1],xyz_conv[0]);
              Double_t radius=sqrt((xyz_conv[1]-displaceY)*(xyz_conv[1]-displaceY)+(xyz_conv[0]-displaceX)*(xyz_conv[0]-displaceX));
              Double_t X_loc =  xyz_conv[0]*cos(alpha) + xyz_conv[1]*sin(alpha);
              Double_t Y_loc = -xyz_conv[0]*sin(alpha) + xyz_conv[1]*cos(alpha);
              double covar[21]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
              double ptemp[5] = {Y_loc,xyz_conv[2],0,0,0};
              particle.fDirection.resize(count+1);
              particle.fParamMC.resize(count+1);
              // particle.fResolRPhi.resize(count+1);
              // particle.fResolZ.resize(count+1);
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

/*
Int_t BuildParticleNoRotation(fastParticle &particle, double Center[3], fastGeometry geom, 
                   std::vector<TVector3> trajxyz, std::vector<TVector3> trajpxyz, long PDGcode, float angle = 0.)
{
          uint fMaxLayer = 0;
          TParticlePDG *p = TDatabasePDG::Instance()->GetParticle(PDGcode);
          if (p == nullptr) {
            ::Error("fastParticle::simulateParticle", "Invalid pdgCode %ld", PDGcode);
            return -1;
          }
          Short_t sign = 1 * p->Charge() / 3.;
          Float_t mass = p->Mass();
          particle.fMassMC=mass;

          for(size_t k=0;k<trajxyz.size();k++) 
          {
              Bool_t invert = kFALSE;
              Double_t xyz_conv[3]= {(trajxyz.at(k).Z()-(Center[2]-278)),
                                trajxyz.at(k).Y()-Center[1],
                                trajxyz.at(k).X()-Center[0]};

              Double_t pxyz_conv[3]= {(trajpxyz.at(k).Z()),
                                      trajpxyz.at(k).Y(),
                                      trajpxyz.at(k).X()};
              //Double_t totalp = sqrt(pxyz_conv[0]*pxyz_conv[0]+pxyz_conv[0]*pxyz_conv[0]+pxyz_conv[2]*pxyz_conv[2]);
              //if(totalp<0.01) break;

              Double_t alpha=TMath::ATan2(xyz_conv[1],xyz_conv[0]);
              Double_t radius=sqrt((xyz_conv[1])*(xyz_conv[1])+(xyz_conv[0]-278)*(xyz_conv[0]-278));
              Double_t X_loc =  xyz_conv[0]*cos(angle) + xyz_conv[1]*sin(angle);
              Double_t Y_loc = -xyz_conv[0]*sin(angle) + xyz_conv[1]*cos(angle);
              double covar[21]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

              particle.fDirection.resize(k+1);
              particle.fParamMC.resize(k+1);
              particle.fLayerIndex.resize(k+1);
              particle.fLoop.resize(k+1);

              particle.fLoop[k]=0;
              Float_t dirX=(trajxyz.at(trajxyz.size()-1).Z())-(trajxyz.at(0).Z()); 
              dirX>0? dirX=1:dirX=-1;
              if(dirX<0){
                pxyz_conv[0]*=-1;
                pxyz_conv[1]*=-1;
                pxyz_conv[2]*=-1;
              }
              AliExternalTrackParam param(xyz_conv,pxyz_conv,covar,sign*dirX);
              bool statuscheck = param.Rotate(angle);
              if(statuscheck){ 
                particle.fParamMC[k].Set(xyz_conv,pxyz_conv,covar,sign*dirX);
                particle.fParamMC[k].Rotate(angle);
                bool x = 0;
              }
              else {
                double ptemp[] = {Y_loc,xyz_conv[2],0,0,0};
                particle.fParamMC[k].SetParamOnly(X_loc,angle,ptemp);
                double x = 0;
              }
              //bool statuscheck = particle.fParamMC[k].Rotate(0);
              
              //particle.fParamMC[k].Rotate(alphca);
              uint indexR = uint(std::upper_bound (geom.fLayerRadius.begin(),geom.fLayerRadius.end(), radius)-geom.fLayerRadius.begin());
              particle.fLayerIndex[k] = indexR;
              if(k!=0)
              {
                TVector3 xyz_prev(trajxyz.at(k-1).Z()-Center[2],
                                  trajxyz.at(k-1).Y()-Center[1],
                                  trajxyz.at(k-1).X()-Center[0]);
                Double_t r_prev = sqrt(xyz_prev.Y()*xyz_prev.Y()+xyz_prev.X()*xyz_prev.X());
                if((radius/r_prev)>1) particle.fDirection[k] = +1;
                else particle.fDirection[k] = -1;
              }
              else particle.fDirection[k] = +1;

              if (indexR>fMaxLayer) fMaxLayer=indexR;
          }
          if(particle.fDirection.size()>1) particle.fDirection[0]=particle.fDirection[1];

          return 1;
}

Int_t BuildInRotfromIn(fastParticle &particle){
  particle.fLengthInRot=0;
  particle.fParamInRot.resize(particle.fParamMC.size());
  for(int i=0; i<particle.fParamMC.size(); i++){
    particle.fParamInRot[i]=particle.fParamIn[i];
    if(i==0) particle.fLengthInRot=0;
    else{
      double xyzp[3];
      double xyz[3];
      particle.fParamInRot[i-1].GetXYZ(xyzp);
      particle.fParamInRot[i].GetXYZ(xyz);
      if(xyz[0]!=0&&xyz[1]!=0&&xyz[2]!=0&&xyzp[0]!=0&&xyzp[1]!=0&&xyzp[2]!=0) particle.fLengthInRot+=sqrt(pow(xyz[0]-xyzp[0],2)+pow(xyz[1]-xyzp[1],2)+pow(xyz[1]-xyzp[1],2));
    }

  }
  return 1;
}
Int_t SetParticleLoop(fastParticle &particle){

 Int_t sign = particle.fParamMC[0].GetParameter()[4]>0? 1:-1;
 Int_t loopcounter = 0;
  for(size_t k=0;k<particle.fParamMC.size();k++) {
     Int_t signtest = particle.fParamMC[k].GetParameter()[4];
     if(signtest>0) signtest=1;
     else if(signtest<0) signtest=-1;
     else if(signtest==0) signtest=0;

     signtest*=sign;
     if(signtest>0){
      particle.fLoop.resize(k+1);
      particle.fLoop[k]=loopcounter;
     }
     if(signtest<0){
      sign*=-1;
      loopcounter++;
      particle.fLoop.resize(k+1);
      particle.fLoop[k]=loopcounter;
     }
     if(signtest==0){
      int zerocounter=0;
      float checkparam=0;
      while(checkparam==0){
        zerocounter++;
        if((k+zerocounter)==particle.fParamMC.size()) return 0;
        checkparam = particle.fParamMC[k+zerocounter].GetParameter()[4];
      }
      for(int i=0;i<=zerocounter;i++){
        if(i<=int(zerocounter/2)){
          particle.fLoop.resize(k+i+1);
          particle.fLoop[k+i]=loopcounter+1;
        }else{
          particle.fLoop.resize(k+i+1);
          particle.fLoop[k+i]=loopcounter+1;
        }
      }
      loopcounter++;
      sign*=-1;
      k+=(zerocounter);
     }

  }
  return 1;
}


Int_t BuildParamNDGArReco(AliExternalTrackParam4D &paramNDGAr4D, double Center[3], 
                          TVector3 trajxyz, TVector3 trajpxyz, long PDGcode,
                          Short_t signr=0, double displaceX=0, double displaceY=0)
{
        Int_t inv = 1;
        TParticlePDG *p = TDatabasePDG::Instance()->GetParticle(PDGcode);
        if (p == nullptr) {
          ::Error("fastParticle::simulateParticle", "Invalid pdgCode %ld", PDGcode);
          return -1;
        }
        Float_t mass = p->Mass();
        Short_t sign;
        if(signr==0) sign = 1 * p->Charge() / 3.;
        else sign=signr;

        Bool_t invert = kFALSE;
        Double_t xyz_conv[3]= {(trajxyz.Z()-(Center[2]-displaceX)),
                                (trajxyz.Y()-(Center[1]-displaceY)),
                          trajxyz.X()-Center[0]};

        Double_t pxyz_conv[3]= {(trajpxyz.Z()),
                                  trajpxyz.Y(),
                                  trajpxyz.X()};

        Double_t alpha=TMath::ATan2(xyz_conv[1],xyz_conv[0]);
        Double_t radius=sqrt(xyz_conv[1]*xyz_conv[1]+xyz_conv[0]*xyz_conv[0]);
        Double_t X_loc =  xyz_conv[0]*cos(alpha) + xyz_conv[1]*sin(alpha);
        Double_t Y_loc = -xyz_conv[0]*sin(alpha) + xyz_conv[1]*cos(alpha);

        double covar[21]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        AliExternalTrackParam param(xyz_conv,pxyz_conv,covar,sign);
        AliExternalTrackParam4D param4D(param,mass,1);
        Double_t pxyz_conv_test[3];
        param4D.GetPxPyPz(pxyz_conv_test);
        Bool_t status = param4D.Rotate(alpha);
        Double_t pxyz_inv[3] = {-pxyz_conv[0],-pxyz_conv[1],-pxyz_conv[2]};

        if((abs(pxyz_conv_test[0]-pxyz_conv[0])>0.00001) || (!status))
        {
          AliExternalTrackParam param_inv(xyz_conv,pxyz_inv,covar,sign);
          AliExternalTrackParam4D param4D_inv(param_inv,mass,1);
          status = param4D_inv.Rotate(alpha);
          invert=kTRUE;
        }

        if(status)
        {
          if(!invert)
          {
            paramNDGAr4D.Set(xyz_conv,pxyz_conv,covar,sign);
            paramNDGAr4D.Rotate(alpha);
          }
          else
          {
            paramNDGAr4D.Set(xyz_conv,pxyz_inv,covar,-sign);
            paramNDGAr4D.Rotate(alpha);
            inv = -1;
          }
        }
        else
        {
          double ptemp[] = {Y_loc,xyz_conv[2],0,0,0};
          paramNDGAr4D.SetParamOnly(X_loc,alpha,ptemp);
        }


    return inv;
}

void MoveParamToPoint(AliExternalTrackParam4D &paramNDGAr4D, double xyz[3], float xx0, float xrho, float ResolRPhi, float ResolZ, float Bz, double chi2Cut, double mass)
{
      AliExternalTrackParam4D param_prev = paramNDGAr4D;
      AliExternalTrackParam4D param = paramNDGAr4D;
      double alpha=TMath::ATan2(xyz[1],xyz[0]);
      double radius = TMath::Sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);

      bool status = kTRUE;

      status = param.Rotate(alpha);
      if(!status){       
        ::Error("fastParticle::MoveParamToPoint:", "Rotation failed");
        return;       
      }

      status = param.PropagateTo(radius,Bz,1);
      if(!status){       
        ::Error("fastParticle::MoveParamToPoint:", "Propagation failed");
        return;       
      }

      double pos[2]={0,xyz[2]};
      double cov[3]={ResolRPhi*ResolRPhi,0, ResolZ*ResolZ};
      float chi2 =  param.GetPredictedChi2(pos, cov);
      if(chi2>chi2Cut){
        ::Error("fastParticle::MoveParamToPoint:", "Too big chi2 %f", chi2);
        return;
      }

      status = param.Update(pos, cov);
      if(!status){
        ::Error("fastParticle::MoveParamToPoint:", "Update failed");
        return;
      }

      float tanPhi2 = param.GetParameter()[2]*param.GetParameter()[2];
      tanPhi2/=(1-tanPhi2);
      float deltaR = abs(param_prev.GetX()-param.GetX());
      float crossLength=deltaR*TMath::Sqrt((1.+tanPhi2)*(1+param.GetParameter()[3]*param.GetParameter()[3]));
      for (Int_t ic=0;ic<5; ic++) {
          status*= param.CorrectForMeanMaterial(crossLength * xx0/5., crossLength * xrho/5., mass, 0.01, 0+0x2*kTRUE);
      }
      if(!status){
        ::Error("fastParticle::MoveParamToPoint:", "Correct for material failed");
        return;
      }

      if(status) paramNDGAr4D = param;

}

void Vertexing3D(std::vector<TLorentzVector> p4_ev, std::vector<TVector3> xyz_ev, Double_t * xyz){
  
  std::vector<TVector3> exyz;
  TVector3 ex(1., 0.,0.);
  TVector3 ey(0., 1.,0.);
  TVector3 ez(0., 0.,1.);
  exyz.push_back(ex);
  exyz.push_back(ey);
  exyz.push_back(ez);

  std::vector<TVector3> epxyz;
  for(int t=0; t<p4_ev.size(); t++){
     TVector3 pxyz(p4_ev[t].Px(),p4_ev[t].Py(),p4_ev[t].Pz());
     epxyz.push_back(pxyz.Unit());
  }
  //std::cout<<epxyz.size()<<std::endl;
  //for(int t=0; t<epxyz.size(); t++) std::cout<<"epxyz: "<<epxyz[t].Px()<<" "<<epxyz[t].Py()<<" "<<epxyz[t].Pz()<<" "<<epxyz[t].Mag()<<"\n";
  TMatrixD a(3,3);
  for(int i=0;i<3;i++){
    for(int t=0;t<3;t++){
      Double32_t l=0;
      for(int g=0;g<3;g++) {
        TVector3 cross1 = exyz[i].Cross(epxyz[g]);
        TVector3 cross2 = exyz[t].Cross(epxyz[g]);
        l += cross1.Dot(cross2);
      }
      a(i,t)=l;
    }
  }

   
  Double_t det2;
  //a.Print();

  a.Invert(&det2);
  //a.Print();
   
  TMatrixD b(3,1);
  for(int i=0;i<3;i++){
      double l=0;
      for(int g=0;g<3;g++) {
        //std::cout<<"xyzev "<<xyz_ev[g].X()<<" "<<xyz_ev[g].Y()<<" "<<xyz_ev[g].Z()<<"\n";
        //std::cout<<"epxyz "<<-epxyz[g].X()<<" "<<-epxyz[g].Y()<<" "<<-epxyz[g].Z()<<"\n";
        TVector3 cross1 = xyz_ev[g].Cross(epxyz[g]);
        TVector3 cross2 = exyz[i].Cross(epxyz[g]);
        l+=cross1.Dot(cross2);
      }
      b(i,0)=l;
    }
  
  //b.Print();

  TMatrix vert = a * b;
  //vert.Print();

  xyz[0] = vert(0,0);
  xyz[1] = vert(1,0);
  xyz[2] = vert(2,0);
}


void findIntersection(const TVector3& position1, const TLorentzVector& direction1,
                               const TVector3& position2, const TLorentzVector& direction2,
                               const TVector3& position3, const TLorentzVector& direction3, double *intersection) {
    // Check if the directions are parallel
    TVector3 crossProduct = direction1.Vect().Cross(direction2.Vect());
    if (crossProduct.Mag() < 1e-6) {
        std::cerr << "Lines are parallel; no unique intersection point exists." << std::endl;
        return ; // Return a default TLorentzVector
    }

    // Use the equations of three lines to find the intersection point
    double t1 = ((position3.x() - position1.x()) * direction3.Py() * direction2.Pz() +
                  (position2.y() - position1.y()) * direction3.Py() * direction2.Pz() -
                  (position2.z() - position1.z()) * direction3.Py() * direction2.Px() -
                  (position3.y() - position1.y()) * direction3.Px() * direction2.Pz() +
                  (position3.z() - position1.z()) * direction3.Px() * direction2.Py()) /
                 (direction1.Px() * direction3.Py() * direction2.Pz() -
                  direction1.Py() * direction3.Px() * direction2.Pz() +
                  direction1.Pz() * direction3.Px() * direction2.Py() -
                  direction1.Px() * direction3.Pz() * direction2.Py() +
                  direction1.Py() * direction3.Pz() * direction2.Px() -
                  direction1.Pz() * direction3.Py() * direction2.Px());

    double t2 = ((position2.x() - position1.x()) * direction1.Py() * direction3.Pz() -
                  (position2.x() - position1.x()) * direction1.Py() * direction2.Pz() +
                  (position3.x() - position1.x()) * direction1.Py() * direction2.Pz() -
                  (position2.y() - position1.y()) * direction1.Px() * direction3.Pz() +
                  (position2.y() - position1.y()) * direction1.Px() * direction2.Pz() -
                  (position3.y() - position1.y()) * direction1.Px() * direction2.Pz() +
                  (position2.z() - position1.z()) * direction1.Px() * direction3.Py() -
                  (position2.z() - position1.z()) * direction1.Px() * direction2.Py() +
                  (position3.z() - position1.z()) * direction1.Px() * direction2.Py()) /
                 (direction1.Px() * direction3.Py() * direction2.Pz() -
                  direction1.Py() * direction3.Px() * direction2.Pz() +
                  direction1.Pz() * direction3.Px() * direction2.Py() -
                  direction1.Px() * direction3.Pz() * direction2.Py() +
                  direction1.Py() * direction3.Pz() * direction2.Px() -
                  direction1.Pz() * direction3.Py() * direction2.Px());

    // Calculate the intersection point
    
    intersection[0]=(position1.x() + t1 * direction1.Px());
    intersection[1]=(position1.y() + t1 * direction1.Py());
    intersection[2]=(position1.z() + t1 * direction1.Pz());
}





Int_t BuildParamNDGArRecoNoRotation(AliExternalTrackParam4D &paramNDGAr4D, double Center[3], 
                          TVector3 trajxyz, TVector3 trajpxyz, long PDGcode,
                          Short_t signr=0, double angle=0)
{
    TParticlePDG *p = TDatabasePDG::Instance()->GetParticle(PDGcode);
    if (p == nullptr) {
      ::Error("fastParticle::simulateParticle", "Invalid pdgCode %ld", PDGcode);
      return -1;
    }
    Float_t mass = p->Mass();
    Short_t sign;
    if(signr==0) sign = 1 * p->Charge() / 3.;
    else sign=signr;

    Bool_t invert = kFALSE;
    Double_t xyz_conv[3]= {(trajxyz.Z()-(Center[2]-278)),
                      trajxyz.Y()-Center[1],
                      trajxyz.X()-Center[0]};

    Double_t pxyz_conv[3]= {(trajpxyz.Z()),
                              trajpxyz.Y(),
                              trajpxyz.X()};
    Double_t pxyz_inv[3] = {-pxyz_conv[0],-pxyz_conv[1],-pxyz_conv[2]};

    //Double_t alpha=TMath::ATan2(xyz_conv[1],xyz_conv[0]);
    Double_t radius=sqrt(xyz_conv[1]*xyz_conv[1]+(xyz_conv[0]-278)*(xyz_conv[0]-278));
    //Double_t X_loc =  xyz_conv[0]*cos(alpha) + xyz_conv[1]*sin(alpha);
    //Double_t Y_loc = -xyz_conv[0]*sin(alpha) + xyz_conv[1]*cos(alpha);

    double covar[21]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    AliExternalTrackParam param(xyz_conv,pxyz_conv,covar,sign);
    AliExternalTrackParam4D param4D(param,mass,1);
    bool status = param4D.Rotate(angle);

    Double_t pxyz_conv_test[3];
    param4D.GetPxPyPz(pxyz_conv_test);
    AliExternalTrackParam4D param4D_inv_perm;

    if((abs(pxyz_conv_test[0]-pxyz_conv[0])>0.00001) || (!status))
    {
          AliExternalTrackParam param_inv(xyz_conv,pxyz_inv,covar,sign);
          AliExternalTrackParam4D param4D_inv(param_inv,mass,1);
          status = param4D_inv.Rotate(angle);
          param4D_inv_perm=param4D_inv;
          Double_t pxyz_conv_test_inv[3];
          param4D_inv.GetPxPyPz(pxyz_conv_test_inv);
          invert=kTRUE;
    }

        if(status)
        {
          if(!invert)
          {
            //paramNDGAr4D.Set(xyz_conv,pxyz_conv,covar,sign);
            //paramNDGAr4D.Rotate(0);
            paramNDGAr4D=param4D;
          }
          else
          {
            //paramNDGAr4D.Set(xyz_conv,pxyz_inv,covar,-sign);
            //paramNDGAr4D.Rotate(0);
            paramNDGAr4D=param4D_inv_perm;
            Double_t pxyz_conv_test_inv_again[3];
            paramNDGAr4D.GetPxPyPz(pxyz_conv_test_inv_again);
            int x=0;
          }
        }

    //paramNDGAr4D.Set(xyz_conv,pxyz_conv,covar,-sign);
    //bool statuscheck = paramNDGAr4D.Rotate(0);



    return 1;
}
*/
bool inFiducial(double x, double y, double z) {
	double const Rcut = 400.0 /2.0;		double const Zcut = 439.4 /2.0;
	double r = sqrt( z*z + y*y );
	if ( r > Rcut ) return false;
	if (fabs(x) > Zcut ) return false;
	return true;
};

bool inTPC(double x, double y, double z) {
	double const Rcut = 500 /2.0;		double const Zcut = 470 /2.0;
	double r = sqrt( z*z + y*y );
	if ( r > Rcut ) return false;
	if (fabs(x) > Zcut ) return false;
	return true;
};

  /// part.fStatusMaskIn[0]&8192>0&&part.fParamMC[0].GetP()<1000 conditions for drawing

  /*
  remake plots
  gSystem->AddIncludePath("-I\"$fastMCKalman/fastMCKalman/aliKalman/test/\"")
  gSystem->AddIncludePath("-I\"$fastMCKalman/fastMCKalman/MC/\"")
  gSystem->Load("$fastMCKalman/fastMCKalman/aliKalman/test/AliExternalTrackParam.so");
  .L $fastMCKalman/fastMCKalman/MC/fastSimulation.cxx++g
  .L $ConvertMC/testNDGAr.C++g
  AliPDG::AddParticlesToPdgDataBase();
  initTreeFast()
  treeFast->SetAlias("deltaP","(part.fParamMC[0].GetP()-part.fParamIn[0].GetP())/part.fParamMC[0].GetP()")
  treeFast->SetAlias("deltaPorig","(part.fParamMC[0].GetP()-paramSt.GetP())/part.fParamMC[0].GetP()")
  treeFast->Draw("deltaP:sqrt(part.fNPointsIn[0]/part.@fParamMC.size())>>hv(10,0,1,100,-0.4,0.4)","part.fStatusMaskIn[0]&8192>0&&part.fParamMC[0].GetP()<1000&&part.fNPointsIn[0]>50","colz")
  hv->FitSlicesY()
  treeFast->Draw("deltaPorig:sqrt(part.fNPointsIn[0]/part.@fParamMC.size())>>hk(10,0,1,100,-0.4,0.4)","part.fStatusMaskIn[0]&8192>0&&part.fParamMC[0].GetP()<1000&&part.fNPointsIn[0]>50","colz")
  hk->FitSlicesY()
  hv_2->Divide(hk_2)
  hv_2->SetTitle("Impact of skipping points on momentum resolution;#sqrt{N_{New}/N_{Old}};res_{New}/res_{Old}")
  hv_2->Draw()

  */

#endif