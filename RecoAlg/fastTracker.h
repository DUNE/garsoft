#ifndef fastTracker_H
#define fastTracker_H
/*
 .L $fastMCKalman/fastMCKalman/MC/fastTracker.h
 code based on the AliTPCtracker.
*/
// basic tracking functionality extracked
// https://github.com/miranov25/AliRoot/blob/master/TPC/TPCrec/AliTPCtracker.cxx#L1325
// https://github.com/miranov25/AliRoot/blob/2f249b4d1cff4e54e0b17f2b23ddb4d9b5e4bd9e/STEER/ESD/AliTrackerBase.cxx#L673
#include "AliExternalTrackParam.h"

class fastTracker{
public:
  static AliExternalTrackParam* makeSeed(double xyz0[3], double xyz1[3], double xyz2[3], double sy, double sz, float bz);
  static AliExternalTrackParam* makeSeedMB(double xyz0[3], double xyz1[3], double xyz2[3], double sy, double sz, float bz, float xx0, float xrho, float mass,int nSteps=5,float cross1=0, float cross2=0);
  static AliExternalTrackParam* makeSeedMB4(double xyz0[3], double xyz1[3], double xyz2[3], double sy, double sz, float bz, float xx0, float xrho, float mass,int nSteps=5, bool MSsmearing = kTRUE, float cross1=0, float cross2=0);
  static Double_t makeC(Double_t x1,Double_t y1, Double_t x2,Double_t y2, Double_t x3,Double_t y3);     // F1
  static Double_t makeSnp(Double_t x1,Double_t y1,Double_t x2,Double_t y2,Double_t x3,Double_t y3);     // F2
  static Double_t makeYC(Double_t x1,Double_t y1,Double_t x2,Double_t y2,Double_t x3,Double_t y3);     // YC
  static Double_t makeTgln(Double_t x1,Double_t y1, Double_t x2,Double_t y2, Double_t z1,Double_t z2,Double_t c);   // F3n
};


//_____________________________________________________________________________
Double_t fastTracker::makeC(Double_t x1,Double_t y1, Double_t x2,Double_t y2, Double_t x3,Double_t y3){
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
  Double_t x0 = x3*0.5-y3*u;
  Double_t y0 = y3*0.5+x3*u;
  Double_t c2 = 1/TMath::Sqrt(x0*x0+y0*y0);
  if (det<0) c2*=-1;
  return c2;
}


Double_t fastTracker::makeSnp(Double_t x1,Double_t y1, Double_t x2,Double_t y2, Double_t x3,Double_t y3){
  //-----------------------------------------------------------------
  // Initial approximation of the track snp at position x1
  //-----------------------------------------------------------------
  x3 -=x1;
  x2 -=x1;
  y3 -=y1;
  y2 -=y1;
  //  
  Double_t det = x3*y2-x2*y3;
  if (TMath::Abs(det)<1e-10) {
    return 100;
  }
  //
  Double_t u = 0.5* (x2*(x2-x3)+y2*(y2-y3))/det;
  Double_t x0 = x3*0.5-y3*u; 
  Double_t y0 = y3*0.5+x3*u;
  Double_t c2 = 1/TMath::Sqrt(x0*x0+y0*y0);
  if (det>0) c2*=-1;
  x0*=c2;  
  return x0;
}

Double_t fastTracker::makeYC(Double_t x1,Double_t y1, Double_t x2,Double_t y2, Double_t x3,Double_t y3){
  //-----------------------------------------------------------------
  // Initial approximation of the y coordinate of the center of the track circumference 
  // in the xy plane, with respects to the first point (x1,y1). Used to check consistency 
  // between points (i.e. if they are in the the same semiplane), not in the seeding itself. 
  // If the sign of yC is the same, the points are in the same semiplane.
  //-----------------------------------------------------------------
  x3 -=x1;
  x2 -=x1;
  y3 -=y1;
  y2 -=y1;
  //  
  Double_t det = x3*y2-x2*y3;
  if (TMath::Abs(det)<1e-10) {
    return 100;
  }
  //
  Double_t u = 0.5* (x2*(x2-x3)+y2*(y2-y3))/det;
  Double_t y0 = y3*0.5+x3*u;
  return y0;
}

//_____________________________________________________________________________
Double_t fastTracker::makeTgln(Double_t x1,Double_t y1, Double_t x2,Double_t y2,Double_t z1,Double_t z2,Double_t c){
  //-----------------------------------------------------------------
  // Initial approximation of the tangent of the track dip angle
  //-----------------------------------------------------------------
  Double_t d  =  TMath::Sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
  if (TMath::Abs(d*c*0.5)>1) return 0;
  Double_t   angle2    = asin(d*c*0.5);

  angle2  = (z1-z2)*c/(angle2*2.);    //dz /(R*dPhi)
  return angle2;
  //return (z1 - z2)/sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}

AliExternalTrackParam * fastTracker::makeSeed(double xyz0[3], double xyz1[3], double xyz2[3], double sy, double sz, float bz){
  Double_t sy2=sy*sy;
  Double_t sz2=sz*sz;
  Double_t param[5];
  Double_t c[15];
  // calculate initial param
  param[0]=xyz0[1];              
  param[1]=xyz0[2];
  param[2]=makeSnp(xyz0[0],xyz0[1],xyz1[0],xyz1[1],xyz2[0],xyz2[1]);
  param[4]=makeC(xyz0[0],xyz0[1],xyz1[0],xyz1[1],xyz2[0],xyz2[1]);
  param[3]=makeTgln(xyz0[0],xyz0[1],xyz1[0],xyz1[1],xyz0[2],xyz1[2],param[4]);
  // x[3]=F3n(x1,y1,x2,y2,z1,z2,x[4]);
  // x[4]=F1(x1,y1,x2,y2,x3,y3);
  //
  Double_t f40=(makeC(xyz0[0],xyz0[1]+sy,xyz1[0],xyz1[1],xyz2[0],xyz2[1])-param[4])/sy;
  Double_t f42=(makeC(xyz0[0],xyz0[1],xyz1[0],xyz1[1]+sy,xyz2[0],xyz2[1])-param[4])/sy;
  Double_t f43=(makeC(xyz0[0],xyz0[1],xyz1[0],xyz1[1],xyz2[0],xyz2[1]+sy)-param[4])/sy;
  //
  Double_t f20=(makeSnp(xyz0[0],xyz0[1]+sy,xyz1[0],xyz1[1],xyz2[0],xyz2[1])-param[2])/sy;
  Double_t f22=(makeSnp(xyz0[0],xyz0[1],xyz1[0],xyz1[1]+sy,xyz2[0],xyz2[1])-param[2])/sy;
  Double_t f23=(makeSnp(xyz0[0],xyz0[1],xyz1[0],xyz1[1],xyz2[0],xyz2[1]+sy)-param[2])/sy;
  //
  //  makeTgln(xyz0[0],xyz0[1],xyz1[0],xyz1[1],xyz0[2],xyz1[2],param[4]);
  Double_t f30=(makeTgln(xyz0[0],xyz0[1]+sy,xyz1[0],xyz1[1],xyz0[2],xyz1[2],param[4])-param[3])/sy;
  Double_t f31=(makeTgln(xyz0[0],xyz0[1],xyz1[0],xyz1[1],xyz0[2]+sz,xyz1[2],param[4])-param[3])/sz;
  Double_t f32=(makeTgln(xyz0[0],xyz0[1],xyz1[0],xyz1[1]+sy,xyz0[2],xyz1[2],param[4])-param[3])/sy;
  Double_t f34=(makeTgln(xyz0[0],xyz0[1],xyz1[0],xyz1[1],xyz0[2],xyz1[2]+sz,param[4])-param[3])/sz;
  c[0]=sy2;
  c[1]=0.;       c[2]=sz2;
  c[3]=f20*sy2;   c[4]=0.;       c[5]=f20*sy2*f20+f22*sy2*f22+f23*sy2*f23;
  c[6]=f30*sy2;  c[7]=f31*sz2;  c[8]=f30*sy2*f20+f32*sy2*f22;
  c[9]=f30*sy2*f30+f31*sz2*f31+f32*sy2*f32+f34*sz2*f34;
  c[10]=f40*sy2; c[11]=0.; c[12]=f40*sy2*f20+f42*sy2*f22+f43*sy2*f23;
  c[13]=f30*sy2*f40+f32*sy2*f42;
  c[14]=f40*sy2*f40+f42*sy2*f42+f43*sy2*f43;

  if (TMath::Abs(bz)>kAlmost0Field) {
    c[14]/=(bz*kB2C)*(bz*kB2C);
    c[13]/=(bz*kB2C);
    c[12]/=(bz*kB2C);
    c[10]/=(bz*kB2C);
    param[4]/=(bz*kB2C); // transform to 1/pt
  }
  else { // assign 0.6 GeV pT
    const double kq2pt = 1./0.6;
    param[4] = kq2pt;
    c[14] = (0.5*0.5)*kq2pt;
  }

  //
  //AliExternalTrackParam *param=new AliExternalTrackParam(
  //AliExternalTrackParam AliExternalTrackParam(Double_t x, Double_t alpha, const Double_t[5] param, const Double_t[15] covar)
  AliExternalTrackParam *extParam=new AliExternalTrackParam(xyz0[0],0,param,c);
  return extParam;
}

/// make seed correcting for the material budget  -assuming particle mass and mean x0 and density
/// \param xyz0
/// \param xyz1
/// \param xyz2
/// \param sy
/// \param sz
/// \param bz
/// \param xx0tocm
/// \param xrhotocm
/// \param mass
/// \return
AliExternalTrackParam* fastTracker::makeSeedMB(double xyz0[3], double xyz1[3], double xyz2[3], double sy, double sz, float bz, float xx0tocm, float xrhotocm, float mass, int nSteps, float cross1, float cross2){
  // calculate momentum loss between the seeding points
  // linear approximation   p0 -> p1 -> p2
  // in case seeding radius is homogenous - pSeed ~ (p0+p1+p2)/3    p0~ pSeed*3-p1-p2
  // in case different gaps - TODO later
  AliExternalTrackParam *extParam = makeSeed(xyz0,xyz1,xyz2,sy,sz,bz);   // first estimation
  AliExternalTrackParam paramFull=*extParam;
  paramFull.SetUseLogTermMS(kTRUE);
  AliExternalTrackParam paramFullnoMC=paramFull;
  Double_t *xyz[3]={xyz0,xyz1,xyz2};
  Double_t deltaCovar[15] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  Bool_t propStatus=kTRUE;
  Double_t p0[3]={paramFull.Pt()};
  for (int i=1; i<3; i++) {

    double crossLength = ((xyz[i][0] - xyz[i - 1][0]) * (xyz[i][0] - xyz[i - 1][0]) +
                         (xyz[i][1] - xyz[i - 1][1]) * (xyz[i][1] - xyz[i - 1][1]) +
                         (xyz[i][2] - xyz[i - 1][2]) * (xyz[i][2] - xyz[i - 1][2]));
    crossLength = TMath::Sqrt(crossLength);
    if(cross1!=0 && i==1) crossLength=cross1;
    if(cross2!=0 && i==2) crossLength=cross2;
    for (int j = 0; j < nSteps; j++) {
      propStatus &= paramFull.PropagateTo((xyz[i-1][0]+((j+1)*(xyz[i][0]-xyz[i-1][0])/nSteps)), bz);    
      propStatus &= paramFull.CorrectForMeanMaterial(crossLength * xx0tocm / nSteps, crossLength * xrhotocm / nSteps, mass, kFALSE);
      propStatus &= paramFullnoMC.PropagateTo((xyz[i-1][0]+((j+1)*(xyz[i][0]-xyz[i-1][0])/nSteps)), bz);
    }

    p0[i]=paramFull.Pt();
    if (i == 1) {
      for (int iCovar = 0; iCovar < 15; iCovar++) {
        deltaCovar[iCovar] = paramFull.GetCovariance()[iCovar] - paramFullnoMC.GetCovariance()[iCovar];
      }
    }
  }
  // Formula below approximation in case equal material distance of seeding layer
  Double_t ratio1= p0[1]/p0[0];
  Double_t ratio2= p0[2]/p0[0];
  //Double_t p0N=3.*p0[0]/(ratio2+ratio1+1.);
  Double_t p0NRatio=3./(ratio2+ratio1+1.);
  /// This hack - we should get proper curvature/sagita formula
  //
  ((double*)extParam->GetParameter())[4]/=  p0NRatio;
  ((double*)extParam->GetCovariance())[5] +=deltaCovar[5];
  ((double*)extParam->GetCovariance())[9] +=deltaCovar[9];
  ((double*)extParam->GetCovariance())[14]+=deltaCovar[14];
  return extParam;
}


AliExternalTrackParam* fastTracker::makeSeedMB4(double xyz0[3], double xyz1[3], double xyz2[3], double sy, double sz, float bz, float xx0tocm, float xrhotocm, float mass, int nSteps, bool MSsmearing, float cross1, float cross2){
  // calculate momentum loss between the seeding points
  // linear approximation   p0 -> p1 -> p2
  // in case seeding radius is homogenous - pSeed ~ (p0+p1+p2)/3    p0~ pSeed*3-p1-p2
  // in case different gaps - TODO later
  AliExternalTrackParam *extParam = makeSeed(xyz0,xyz1,xyz2,sy,sz,bz);   // first estimation
  AliExternalTrackParam paramFulls=*extParam;
  AliExternalTrackParam4D paramFull(paramFulls,mass,1);

  paramFull.SetUseLogTermMS(kTRUE);
  AliExternalTrackParam4D paramFullnoMC=paramFull;
  Double_t *xyz[3]={xyz0,xyz1,xyz2};
  Double_t deltaCovar[15];

  Bool_t propStatus=kTRUE;
  Double_t p0[3]={paramFull.Pt()};
  for (int i=1; i<3; i++) {

    double crossLength = ((xyz[i][0] - xyz[i - 1][0]) * (xyz[i][0] - xyz[i - 1][0]) +
                         (xyz[i][1] - xyz[i - 1][1]) * (xyz[i][1] - xyz[i - 1][1]) +
                         (xyz[i][2] - xyz[i - 1][2]) * (xyz[i][2] - xyz[i - 1][2]));
    crossLength = TMath::Sqrt(crossLength);
    if(cross1!=0 && i==1) crossLength=cross1;
    if(cross2!=0 && i==2) crossLength=cross2;
    
    for (int j = 0; j < nSteps; j++) {
      propStatus &= paramFull.PropagateTo((xyz[i-1][0]+((j+1)*(xyz[i][0]-xyz[i-1][0])/nSteps)), bz,1);    
      if(i==2) propStatus &= paramFull.CorrectForMeanMaterial(crossLength * xx0tocm/5., crossLength * xrhotocm/5., mass, 0.01, 0+0x2*MSsmearing);
      propStatus &= paramFullnoMC.PropagateTo((xyz[i-1][0]+((j+1)*(xyz[i][0]-xyz[i-1][0])/nSteps)), bz, 1);
    }

    p0[i]=paramFull.Pt();
    if (i == 2) {
      for (int iCovar = 0; iCovar < 15; iCovar++) {
        deltaCovar[iCovar] = paramFull.GetCovariance()[iCovar] - paramFullnoMC.GetCovariance()[iCovar];
      }
    }
  }
  // Formula below approximation in case equal material distance of seeding layer
  Double_t ratio1= p0[1]/p0[0];
  Double_t ratio2= p0[2]/p0[0];
  //Double_t p0N=3.*p0[0]/(ratio2+ratio1+1.);
  Double_t p0NRatio=3./(ratio2+ratio1+1.);
  /// This hack - we should get proper curvature/sagita formula
  //
  ((double*)extParam->GetParameter())[4]/=  p0NRatio;
  ((double*)extParam->GetCovariance())[5] +=deltaCovar[5];
  ((double*)extParam->GetCovariance())[9] +=deltaCovar[9];
  ((double*)extParam->GetCovariance())[14]+=deltaCovar[14];
  return extParam;
}



#endif