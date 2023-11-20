/*
    .L $fastMCKalman/fastMCKalman/MC/fastSimulation.cxx+
    geom= fastGeometry(201.)
    particle=fastParticle(200);
    geom.setLayerRadiusPower(0,10,1,20,1.2,0);
    geom.setLayerRadiusPower(11,200,20,150,1.2,0);
*/

#include "fastSimulation.h"
//#include "TTreeStream.h"
#include "TRandom.h"
#include "fastTracker.h"
#include <TMatrixD.h>
#include <TVector3.h>
#include <iostream>


TTreeSRedirector* fastParticle::fgStreamer = nullptr;

// ClassImp(fastGeometry)
// ClassImp(fastParticle)
// ClassImp(AliExternalTrackParam4D)


AliExternalTrackParam4D::AliExternalTrackParam4D():
  AliExternalTrackParam(),
  fZ(1),
  fMass(0),
  fLength(0),
  fTime(0){

}
AliExternalTrackParam4D::AliExternalTrackParam4D(const AliExternalTrackParam &t,Double_t mass,Int_t z,float length,float time):
  AliExternalTrackParam(t),
  fZ(z),
  fMass(mass),
  fLength(length),
  fTime(time){
}

AliExternalTrackParam4D::~AliExternalTrackParam4D(){}

// /// Estimate shape size for particle - shape is determined by the diffusion, angles and deposited charge
// /// \param sigma           -  diffusion sigma ~ induction gap - if as in the TPC
// /// \param width           -  pixel width
// /// \param threshold       -  threhsold in units of MIP
// /// \return                - mean number of pixels above threshold
// Double_t AliExternalTrackParam4D::GetOverThr(Float_t sigma, Float_t width, Float_t threshold){
//   return 0; /// TODO add implemetation
// }


/// propagate to radius and update Track length and time
/// \param xk
/// \param b
/// \param timeDir
/// \return
Bool_t AliExternalTrackParam4D::PropagateTo(Double_t xk, Double_t b, Int_t timeDir){
  static const Double_t kcc = 2.99792458e-2;
  Double_t xyzIn[3],xyzOut[3];
  GetXYZ(xyzIn);
  bool status = AliExternalTrackParam::PropagateTo(xk,b);
  if (status && (fX!=xk)){
    //::Error("AliExternalTrackParam4D::PropagateTo","Incosistent propagate %f\t%f",fX,xk);
    return false;
  }
  GetXYZ(xyzOut);
  Double_t length=0;
  for (Int_t i=0;i<3;i++) length+=(xyzIn[i]-xyzOut[i])*(xyzIn[i]-xyzOut[i]);
  length=TMath::Sqrt(length);
  //
  fLength += length;
  Double_t p2inv = fP[4]*fP[4]/(1+fP[3]*fP[3]);
  Double_t mBeta = TMath::Sqrt( 1. + fMass*fMass*p2inv ); // 1/beta
  Double_t time = length * mBeta / kcc;
  fTime += timeDir*time;
  return status;
}

Double_t AliExternalTrackParam4D::PropagateToMirrorX(Double_t b, Float_t dir, Double_t  sy, Double_t sz)
{
   //----------------------------------------------------------------
  // Impose a "flip" on the parameter vector by a rotation defined 
  // by a diagonal matrix with diagonal elements: R = {1,1,-1,-1,-1}
  //----------------------------------------------------------------

  Double_t &fP1=fP[1], &fP2=fP[2], &fP3=fP[3], &fP4=fP[4];

  Double_t 
  &fC00=fC[0],   &fC11=fC[2], 
  &fC20=fC[3],   &fC21=fC[4],  
  &fC30=fC[6],   &fC31=fC[7],    
  &fC40=fC[10],  &fC41=fC[11];

  Double_t xyz_Old[3];
  GetXYZ(xyz_Old);
  
  Double_t &fA=fAlpha;

  
  Double_t xc, yc, rc, x0, y0;
  Double_t cs= cosf(fAlpha); Double_t sn=sinf(fAlpha); // RS use float versions: factor 2 in CPU speed
  Double_t crv= GetC(b);    // Curvature 

  rc = 1/crv;  xc = fX-fP[2]*rc;

  //update alpha

  Double_t dummy = 1-(fX-xc)*(fX-xc)*crv*crv;
  if (dummy<0) dummy = 0;
  yc  =  fP[0]+TMath::Sqrt(dummy)/crv;

  x0 = xc*cs - yc*sn; y0 = xc*sn + yc*cs;
  
  float alphaC    = TMath::ATan2(y0,x0);
  float dAlpha    = fAlpha-alphaC;
  if (dAlpha>TMath::Pi()) dAlpha-=TMath::TwoPi();
  if (dAlpha<-TMath::Pi()) dAlpha+=TMath::TwoPi();
  fA = alphaC-dAlpha;

  ///Propagate z
     
  Double_t xyz_New[3];
  GetXYZ(xyz_New);
  double cross = sqrt((xyz_New[0]-xyz_Old[0])*(xyz_New[0]-xyz_Old[0])+(xyz_New[1]-xyz_Old[1])*(xyz_New[1]-xyz_Old[1]));
  Double_t sinphic = 0.5*cross / rc;
  if (abs(sinphic)>1) 
  {
    return 0;
  }
  Double_t dPhic = 2*asin(sinphic);
  Double_t darchxy = dir*abs(dPhic*rc);
  fP1 += darchxy*fP[3];
  fC00+=(sy*darchxy)*(sy*darchxy);
  fC11+=(sz*darchxy)*(sz*darchxy);


  //Flip parameters 2/4
  //fP=R(fP)
  fP2 *= -1;
  fP3 *= -1;
  fP4 *= -1;

  //C=RCR^T
  fC20*=-1; fC21*=-1;
  fC30*=-1; fC31*=-1;
  fC40*=-1; fC41*=-1;



  
 
  CheckCovariance();

  Double_t dArch = abs(darchxy*sqrt(1+fP[3]*fP[3]));

  return dArch; 

}


///  Clone of the original method removing one protection - here we disable  TMath::Abs(cosT)>kAlmost1 check (lAlmost1 was very restrictive)
///  // This method has 3 modes of behaviour
///  // 1) xyz[3] array is provided but alpSect pointer is 0: calculate the position of track intersection
///  //    with circle of radius xr and fill it in xyz array
///  // 2) alpSect pointer is provided: find alpha of the sector where the track reaches local coordinate xr
///  //    Note that in this case xr is NOT the radius but the local coordinate.
///  //    If the xyz array is provided, it will be filled by track lab coordinates at local X in this sector
///  // 3) Neither alpSect nor xyz pointers are provided: just check if the track reaches radius xr
/// \param xr
/// \param bz
/// \param xyz
/// \param alpSect
/// \return
Bool_t AliExternalTrackParam4D::GetXYZatR(Double_t xr,Double_t bz, Double_t *xyz, Double_t* alpSect) const {
  double crv = GetC(bz);
  if ((TMath::Abs(bz)) < kAlmost0Field) crv = 0.;
  const double &fy = fP[0];
  const double &fz = fP[1];
  const double &sn = fP[2];
  const double &tgl = fP[3];
  //
  // general circle parameterization:
  // x = (r0+tR)cos(phi0) - tR cos(t+phi0)
  // y = (r0+tR)sin(phi0) - tR sin(t+phi0)
  // where qb is the sign of the curvature, tR is the track's signed radius and r0
  // is the DCA of helix to origin
  //
  double tR = 1. / crv;            // track radius signed
  double cs = TMath::Sqrt((1 - sn) * (1 + sn));
  double x0 = fX - sn * tR;        // helix center coordinates
  double y0 = fy + cs * tR;
  double phi0 = TMath::ATan2(y0, x0);  // angle of PCA wrt to the origin
  if (tR < 0) phi0 += TMath::Pi();
  if (phi0 > TMath::Pi()) phi0 -= 2. * TMath::Pi();
  else if (phi0 < -TMath::Pi()) phi0 += 2. * TMath::Pi();
  double cs0 = TMath::Cos(phi0);
  double sn0 = TMath::Sin(phi0);
  double r0 = x0 * cs0 + y0 * sn0 - tR; // DCA to origin
  double r2R = 1. + r0 / tR;
  //
  //
  if (r2R < kAlmost0) return kFALSE;  // helix is centered at the origin, no specific intersection with other concetric circle
  if (!xyz && !alpSect) return kTRUE;
  double xr2R = xr / tR;
  double r2Ri = 1. / r2R;
  // the intersection cos(t) = [1 + (r0/tR+1)^2 - (r0/tR)^2]/[2(1+r0/tR)]
  double cosT = 0.5 * (r2R + (1 - xr2R * xr2R) * r2Ri);
  if ((cosT * cosT) >= 1.) {  /// TODO - fix protection or change calculation
    //    printf("Does not reach : %f %f\n",r0,tR);
    return kFALSE; // track does not reach the radius xr
  }
  //
  double t = TMath::ACos(cosT);
  if (tR < 0) t = -t;
  // intersection point
  double xyzi[3];
  xyzi[0] = x0 - tR * TMath::Cos(t + phi0);
  xyzi[1] = y0 - tR * TMath::Sin(t + phi0);
  if (xyz) { // if postition is requested, then z is needed:
    double t0 = TMath::ATan2(cs, -sn) - phi0;
    double z0 = fz - t0 * tR * tgl;
    xyzi[2] = z0 + tR * t * tgl;
  } else xyzi[2] = 0;
  //
  Local2GlobalPosition(xyzi, fAlpha);
  //
  if (xyz) {
    xyz[0] = xyzi[0];
    xyz[1] = xyzi[1];
    xyz[2] = xyzi[2];
  }
  //
  if (alpSect) {
    double &alp = *alpSect;
    // determine the sector of crossing
    double phiPos = TMath::Pi() + TMath::ATan2(-xyzi[1], -xyzi[0]);
    int sect = ((Int_t)(phiPos * TMath::RadToDeg())) / 20;
    alp = TMath::DegToRad() * (20 * sect + 10);
    double x2r, f1, f2, r1, r2, dx, dy2dx, yloc = 0, ylocMax = xr * TMath::Tan(TMath::Pi() / 18); // min max Y within sector at given X
    //
    while (1) {
      Double_t ca = TMath::Cos(alp - fAlpha), sa = TMath::Sin(alp - fAlpha);
      if ((cs * ca + sn * sa) < 0) {
        //AliDebug(1, Form("Rotation to target sector impossible: local cos(phi) would become %.2f", cs * ca + sn * sa));
        return kFALSE;
      }
      //
      f1 = sn * ca - cs * sa;
      if (TMath::Abs(f1) >= kAlmost1) {
        //AliDebug(1, Form("Rotation to target sector impossible: local sin(phi) would become %.2f", f1));
        return kFALSE;
      }
      //
      double tmpX = fX * ca + fy * sa;
      double tmpY = -fX * sa + fy * ca;
      //
      // estimate Y at X=xr
      dx = xr - tmpX;
      x2r = crv * dx;
      f2 = f1 + x2r;
      if (TMath::Abs(f2) >= kAlmost1) {
        //AliDebug(1, Form("Propagation in target sector failed ! %.10e", f2));
        return kFALSE;
      }
      r1 = TMath::Sqrt((1. - f1) * (1. + f1));
      r2 = TMath::Sqrt((1. - f2) * (1. + f2));
      dy2dx = (f1 + f2) / (r1 + r2);
      yloc = tmpY + dx * dy2dx;
      if (yloc > ylocMax) {
        alp += 2 * TMath::Pi() / 18;
        sect++;
      }
      else if (yloc < -ylocMax) {
        alp -= 2 * TMath::Pi() / 18;
        sect--;
      }
      else break;
      if (alp >= TMath::Pi()) alp -= 2 * TMath::Pi();
      else if (alp < -TMath::Pi()) alp += 2 * TMath::Pi();
      //      if (sect>=18) sect = 0;
      //      if (sect<=0) sect = 17;
    }
    //
    // if alpha was requested, then recalculate the position at intersection in sector
    if (xyz) {
      xyz[0] = xr;
      xyz[1] = yloc;
      if (TMath::Abs(x2r) < 0.001) xyz[2] = fz + dx * (r2 + f2 * dy2dx) * tgl;
      else {
        // for small dx/R the linear apporximation of the arc by the segment is OK,
        // but at large dx/R the error is very large and leads to incorrect Z propagation
        // angle traversed delta = 2*asin(dist_start_end / R / 2), hence the arc is: R*deltaPhi
        // The dist_start_end is obtained from sqrt(dx^2+dy^2) = x/(r1+r2)*sqrt(2+f1*f2+r1*r2)
        // Similarly, the rotation angle in linear in dx only for dx<<R
        double chord = dx * TMath::Sqrt(1 + dy2dx * dy2dx);   // distance from old position to new one
        double rot = 2 * TMath::ASin(0.5 * chord * crv); // angular difference seen from the circle center
        xyz[2] = fz + rot / crv * tgl;
      }
      Local2GlobalPosition(xyz, alp);
    }
  }
  return kTRUE;
  //
}
/// get radiial direction sign
Int_t AliExternalTrackParam4D::GetDirectionSign(){
  Float_t dir = GetX()*Px()+GetY()*Py();
  return (dir>0)? 1:-1;
}




/// Correct for mean material
/// WARNING - we use strange ALICE convention signing Z==2 particle with negative mass - TODO - replace it with explicit Q
/// \param xOverX0        - X/X0, the thickness in units of the radiation length.
/// \param xTimesRho      - is the product length*density (g/cm^2).
//                        - It should be passed as negative when propagating tracks
//                        - from the interaction point to the outside of the central barrel.
/// \param mass           - the mass of this particle (GeV/c^2). Negative mass means charge=2 particle
/// \param mcSwitch    - can be used for the simulation - bit 2 = addSmearing, bit 1 = stop on covariance
/// \param f              - dEdx formula
/// \param stepFraction   - step fraction  - above some limits RungeKuta instead of the Euler Method used
/// \return  CorrectForMeanMaterial status  (kFalse - Failed, kTrue - Success)
Bool_t AliExternalTrackParam4D::CorrectForMeanMaterial(Double_t xOverX0, Double_t xTimesRho, Double_t mass, Float_t stepFraction, int mcSwitch, Double_t (*f)(Double_t)){
  bool addMSSmearing = (mcSwitch&0x2)>0;
  bool isMC = (mcSwitch&0x1)>0;
  const Double_t kBGStop=0.000040;
  double pOld=GetP();
  Double_t p=GetP();
  Double_t q=(mass<0)?2.:1.;   // q=2 particle in ALICE convention
  mass=TMath::Abs(mass);
  Double_t mass2=mass*mass;
  p*=q;
  if ((p/mass)<kBGStop) {
    return kFALSE;
  }
  Double_t p2=p*p;

  Double_t Ein=TMath::Sqrt(p2+mass2);
  Double_t dP= dPdxEulerStep(p,mass,xTimesRho,stepFraction,f);
  if (dP==0) {
    return kFALSE;
  }
  if (dP>0 &&isMC){
    //::Error("aliExternalTrackParam4D", "Incorrect energy loss %f -> %f ", p,dP);
    //dPdxEulerStep(p,mass,xTimesRho,stepFraction,f); // THIS was debug symbol - TODO remove it later
  }

  Double_t pOut=p+dP;
  if ((pOut/mass)<kBGStop) {
    if (!isMC) {return kFALSE; }
  }
  Double_t Eout=TMath::Sqrt(pOut*pOut+mass2);
  //p=(p+pOut)*0.5; /// TODO -this is turtle - should be used the out value
  p=pOut;           /// use the momenta at the end of step  - approcimation dor in dPdxEulerStep
  // Use mean values for p2, p and beta2
  p2=p*p;
  Double_t beta2=p2/(p2+mass2);
  //
  Double_t &fP2=fP[2];
  Double_t &fP3=fP[3];
  Double_t &fP4=fP[4];
  Double_t &fC22=fC[5];
  Double_t &fC33=fC[9];
  Double_t &fC43=fC[13];
  Double_t &fC44=fC[14];
  //
  //Calculating the multiple scattering corrections******************
  Double_t cC22 = 0.;
  Double_t cC33 = 0.;
  Double_t cC43 = 0.;
  Double_t cC44 = 0.;
  Double_t theta2=0;
  Double_t sigmadPRel=0;
  if (addMSSmearing && xOverX0 != 0) {
    //Double_t theta2=1.0259e-6*14*14/28/(beta2*p2)*TMath::Abs(d)*9.36*2.33;
    theta2=0.0136*0.0136/(beta2*p2)*TMath::Abs(xOverX0);
    if (GetUseLogTermMS()) {
      double lt = 1+0.038*TMath::Log(TMath::Abs(xOverX0));
      // double Zeff=7.31201;
      // double lt = sqrt(1+0.105*TMath::Log(TMath::Abs(xOverX0))+0.0035*pow(TMath::Log(TMath::Abs(xOverX0)),2))*(1-0.24/(Zeff*(Zeff+1)));
      if (lt>0) theta2 *= lt*lt;
    }
    theta2 *= q*q;    // q=2 particle
    if(theta2>TMath::Pi()*TMath::Pi()) {
      return kFALSE;
    }
    // double Zeff=7.31201;
    // double thetaurban = 0.0136/(p*sqrt(beta2))*sqrt(TMath::Abs(xOverX0))*q;
    // double corr = sqrt(1+0.105*TMath::Log(TMath::Abs(xOverX0))+0.0035*pow(TMath::Log(TMath::Abs(xOverX0)),2))*(1-0.24/(Zeff*(Zeff+1)));
    // double lt = 1+0.038*TMath::Log(TMath::Abs(xOverX0));
    
    cC22 = theta2*((1.-fP2)*(1.+fP2))*(1. + fP3*fP3);
    cC33 = theta2*(1. + fP3*fP3)*(1. + fP3*fP3);
    cC43 = theta2*fP3*fP4*(1. + fP3*fP3);
    cC44 = theta2*fP3*fP4*fP3*fP4;
  }

  //Calculating the energy loss corrections************************
  Double_t cP4=1.;
  if ((xTimesRho != 0.) && (beta2 < 1.)) {
    Double_t dE=Eout-Ein;
    if ( (1.+ dE/p2*(dE + 2*Ein)) < 0. ) {
      return kFALSE;
    }
    cP4 = 1./TMath::Sqrt(1.+ dE/p2*(dE + 2*Ein));  //A precise formula by Ruben !  - //TODO -this formula is only first taylor approximation
    cP4 = pOld/pOut;                               /// TODO we use momentum loss not need to use "Ruben E loss approximation"
    //if (TMath::Abs(fP4*cP4)>100.) return kFALSE; //Do not track below 10 MeV/c -disable controlled by the BG cut
    // Approximate energy loss fluctuation (M.Ivanov)
    const Double_t knst=0.07; // TODO To be tuned.
    //Double_t sigmadE=knst*TMath::Abs(dE);                 /// TODO remove that part if momentm smearing working well
    //cC44 += ((sigmadE*Ein/p2*fP4)*(sigmadE*Ein/p2*fP4));
    //
    sigmadPRel=TMath::Abs(pOut-pOld)*knst/pOld;
    cC44+=sigmadPRel*sigmadPRel*fP[4]*fP[4];        // sigma due momentum loss fluctuation
  }

  //Applying the corrections*****************************
  if (addMSSmearing && isMC){
    const float kMaxP3=0.5;
    const float kMaxP4=0.3;
    if (TMath::Sqrt(cC44)>kMaxP4*TMath::Abs(fP4)) {
      if (!isMC) return kFALSE;
    }
    double phi = TMath::ASin(fP[2]);
    double lambda=TMath::ATan(fP[3]);
    float p2New=TMath::Sin(gRandom->Gaus(phi,sqrt(theta2)));
    float p3New=TMath::Tan(gRandom->Gaus(lambda,sqrt(theta2)));
    //Float_t p2New=fP[2]+gRandom->Gaus(0,TMath::Sqrt(cC22));
    //Float_t dp3New=gRandom->Gaus(0,TMath::Sqrt(cC33));
    if (TMath::Abs(p2New)>1.) {
      if (!isMC) return kFALSE;
    }
    if (TMath::Abs(p3New)>kMaxP3) {
      if (!isMC) return kFALSE;
    }
    Double_t cP4MS = sqrt((1+(p3New*p3New))/(1+(fP[3])*(fP[3]))); ////keep total momentum constant and modify q/pt accordingly (this factor is cos(lambda)/cos(lambda_new))
    double p4RelSmear=gRandom->Gaus(0,sigmadPRel);
    fP[2]=p2New;
    fP[3]=p3New;
    fP[4]*=cP4MS;
    fP[4]*=(1.+p4RelSmear);
    //fP[4]+=gRandom->Gaus(0,TMath::Sqrt(cC44)); //- TODO transform scattering in p2 and P3 to modification of qPt - can not be independent
  }
  fC22 += cC22;
  fC33 += cC33;
  fC43 += cC43;
  fC44 += cC44;
  fP4  *= cP4;
  if (isMC && (GetP()>pOld) ){
    //::Error("AliExternalTrackParam4D", "Incorrect energy loss %f -> %f ", pOld,GetP());
  }
  //CheckCovariance();
  return kTRUE;
}


/// First derivative of the dP/dx - for testing and visualization purposes
/// \param p       - particle momenta
/// \param mass    - particle mass
/// \param fdEdx   - dEdx function pointer
/// \return        - dP/dx
Double_t AliExternalTrackParam4D::dPdxEulerStep(double p, double mass,  Double_t xTimesRho, double step, Double_t (*fundEdx)(Double_t)){
    // const Double_t kBGStop=0.0040; // the position of non relybal BB  depends on the function ... not well defined BetheBlocAleph
    const Double_t kBGStop=0.0140;    // not well defined BetheBlocSolid
    Double_t bg=p/mass;
    double dEdxMin=fundEdx(3);
    float signCorr=(xTimesRho<0)? -1:1;
    double pOrig=p;
    if (bg<kBGStop) {
      return pOrig*signCorr;    /// if particle too low BG - we will let to loose full momenta
    }
    Double_t dPdx=TMath::Abs(fundEdx(bg))*TMath::Sqrt(1.+1./(bg*bg));

    Double_t dPdx2=TMath::Max(fundEdx(bg),dEdxMin)*TMath::Sqrt(1.+1./(bg*bg));
    //
    Int_t nSteps=1+(TMath::Abs(dPdx*xTimesRho)+TMath::Abs(dPdx2*xTimesRho))/step;
    if (nSteps==1) {
      float dP=0.5*(dPdx+dPdx2)*xTimesRho;
      if (abs(dP)>pOrig) return pOrig*signCorr;
      return dP; // can not loos more than full momenta -use linear approximation
    }
    Float_t xTimesRhoS=xTimesRho/nSteps;
    Float_t sumP=0;
    for (Int_t i=0; i<nSteps;i++){
      p+=dPdx*xTimesRhoS;
      sumP+=dPdx*xTimesRhoS;
      bg=p/mass;
      if (bg<kBGStop) return pOrig*signCorr;
      dPdx=TMath::Max(fundEdx(bg),dEdxMin)*TMath::Sqrt(1.+1./(bg*bg));
    }
    return TMath::Min(TMath::Abs(double(sumP)),pOrig)*signCorr;
}







void AliExternalTrackParam4D::UpdateTrack(AliExternalTrackParam4D &track1, const AliExternalTrackParam4D &track2){
  //
  // Update track 1 with track 2
  //
  //
  //
  TMatrixD vecXk(5,1);    // X vector
  TMatrixD covXk(5,5);    // X covariance 
  TMatrixD matHk(5,5);    // vector to mesurement
  TMatrixD measR(5,5);    // measurement error 
  TMatrixD vecZk(5,1);    // measurement
  //
  TMatrixD vecYk(5,1);    // Innovation or measurement residual
  TMatrixD matHkT(5,5);
  TMatrixD matSk(5,5);    // Innovation (or residual) covariance
  TMatrixD matKk(5,5);    // Optimal Kalman gain
  TMatrixD mat1(5,5);     // update covariance matrix
  TMatrixD covXk2(5,5);   // 
  TMatrixD covOut(5,5);
  //
  Double_t *param1=(Double_t*) track1.GetParameter();
  Double_t *covar1=(Double_t*) track1.GetCovariance();
  Double_t *param2=(Double_t*) track2.GetParameter();
  Double_t *covar2=(Double_t*) track2.GetCovariance();
  //
  // copy data to the matrix
  for (Int_t ipar=0; ipar<5; ipar++){
    for (Int_t jpar=0; jpar<5; jpar++){
      covXk(ipar,jpar) = covar1[track1.GetIndex(ipar, jpar)];
      measR(ipar,jpar) = covar2[track2.GetIndex(ipar, jpar)];
      matHk(ipar,jpar)=0;
      mat1(ipar,jpar)=0;
    }
    vecXk(ipar,0) = param1[ipar];
    vecZk(ipar,0) = param2[ipar];
    matHk(ipar,ipar)=1;
    mat1(ipar,ipar)=1;
  }
  //
  //
  //
  //
  //
  vecYk = vecZk-matHk*vecXk;                 // Innovation or measurement residual
  matHkT=matHk.T(); matHk.T();
  matSk = (matHk*(covXk*matHkT))+measR;      // Innovation (or residual) covariance
  matSk.Invert();
  matKk = (covXk*matHkT)*matSk;              //  Optimal Kalman gain
  vecXk += matKk*vecYk;                      //  updated vector 
  covXk2 = (mat1-(matKk*matHk));
  covOut =  covXk2*covXk; 
  //
  //
  //
  // copy from matrix to parameters
  if (0) {
    vecXk.Print();
    vecZk.Print();
    //
    measR.Print();
    covXk.Print();
    covOut.Print();
    //
    track1.Print();
    track2.Print();
  }


  for (Int_t ipar=0; ipar<5; ipar++){
    param1[ipar]= vecXk(ipar,0) ;
    for (Int_t jpar=0; jpar<5; jpar++){
      covar1[track1.GetIndex(ipar, jpar)]=covOut(ipar,jpar);
    }
  }
}

void getHelix(Double_t *fHelix,  AliExternalTrackParam t, float bz){
  // addapted code from AliHelix
  Double_t x,cs,sn;
  x=t.GetX();
  for (Int_t i=0;i<5; i++) fHelix[i]=t.GetParameter()[i];
  float alpha = t.GetAlpha();
  //
  //circle parameters
  //PH Sometimes fP4 and fHelix[4] are very big and the calculation
  //PH of the Sqrt cannot be done. To be investigated...
  fHelix[4]=fHelix[4]/(-1000/0.299792458/bz);    // C
  //  cs=TMath::Cos(alpha); sn=TMath::Sin(alpha);
  cs=cosf(alpha); sn=sinf(alpha); // RS use float versions: factor 2 in CPU speed
  Double_t xc, yc, rc;
  rc  =  1/fHelix[4];
  xc  =  x-fHelix[2]*rc;
  Double_t dummy = 1-(x-xc)*(x-xc)*fHelix[4]*fHelix[4];
  if (dummy<0) dummy = 0;
  yc  =  fHelix[0]+TMath::Sqrt(dummy)/fHelix[4];
  fHelix[6] = xc*cs - yc*sn;
  fHelix[7] = xc*sn + yc*cs;
  fHelix[8] =  TMath::Abs(rc);
  return;
}


/// make layer distribution in interval layer0(r0)->layerN(rN) following power low distance
/// \param layer0
/// \param layerN
/// \param r0
/// \param rN
/// \param power
/// \param X0
void fastGeometry::setLayerRadiusPower(int layer0, int layerN, float r0, float rN, float power, float X0,float rho, float resol[2]){
    float dLayerN=(layerN-layer0);
    for (Int_t iLayer=layer0; iLayer<=layerN; iLayer++){
        fLayerX0[iLayer]=X0;
        fLayerRho[iLayer]=rho;
        float radius=r0+TMath::Power(float(iLayer-layer0)/dLayerN,power)*(rN-r0);
        fLayerRadius[iLayer]=radius;
        fLayerResolRPhi[iLayer]=resol[0];
        fLayerResolZ[iLayer]=resol[1];
    }
    fLayerIndex = Argsort(fLayerRadius);
}

void fastGeometry::setLayer(int iLayer, float radius,  float X0,float rho, float resol[2]) {
  fLayerX0[iLayer] = X0;
  fLayerRho[iLayer] = rho;
  fLayerRadius[iLayer] = radius;
  fLayerResolRPhi[iLayer] = resol[0];
  fLayerResolZ[iLayer] = resol[1];
}


void fastParticle::refitParticle()
{
  fParamRefit = fParamIn;
  fParamRefit.resize(fParamIn.size());
  fStatusMaskRefit.resize(fParamIn.size());
  fNPointsRefit.resize(fParamIn.size());
  for(size_t i=0; i<fParamRefit.size();i++) fStatusMaskRefit[i]=0;
  for(size_t i=0; i<fParamRefit.size();i++)
  {
    Bool_t statusIn= kFALSE;
    Bool_t statusOut = kFALSE;

    std::uint16_t flagsIn=fStatusMaskIn[i];

    if((flagsIn & kTrackisOK)>0)
    {
      statusIn=kTRUE;
    }

    std::uint16_t flagsOut=fStatusMaskOut[i]; 

    if((flagsOut & kTrackisOK)>0)
    {
      statusOut=kTRUE;
    }

    if(statusIn && !statusOut)
    {
      fParamRefit[i]=fParamIn[i];
      fStatusMaskRefit[i]|=kTrackUsedIn;
      fStatusMaskRefit[i]|=kTrackisOK;
      fNPointsRefit[i]=fNPointsIn[i];
    }
    else if(!statusIn && statusOut)
    {
      fParamRefit[i]=fParamOut[i];
      fStatusMaskRefit[i]|=kTrackUsedOut;
      fStatusMaskRefit[i]|=kTrackisOK;
      fNPointsRefit[i]=fNPointsOut[i];
    }
    else if (statusIn && statusOut)
    {
      AliExternalTrackParam4D::UpdateTrack(fParamRefit[i],fParamOut[i]);
      fStatusMaskRefit[i]|=kTrackRefitted;
      fStatusMaskRefit[i]|=kTrackisOK;
      fNPointsRefit[i]=fNPointsIn[i]+fNPointsOut[i];
    }


  }
}



///
/// \param geom          - pointer to geometry to use
/// \param pdgCode       - pdgCode used in the reconstruction
/// \param layerStart    - starting layer to do tracking
/// \return   -  TODO  status flags to be decides
int fastParticle::reconstructParticle(fastGeometry  &geom, long pdgCode, uint indexStart){
  const Float_t chi2Cut=100;
  const float kMaxSnp=0.95;
  const float kMaxLoss=0.3;
  const float kCovarFactor=2;
  fLengthIn=0;
  float_t mass=0;
  fPdgCodeRec   =pdgCode;
  if (pdgCode==0){
    fMassRec=fMassMC;
    mass=fMassMC;
  }else {
    TParticlePDG *particle = TDatabasePDG::Instance()->GetParticle(pdgCode);
    if (particle == nullptr) {
      //::Error("fastParticle::reconstructParticle", "Invalid pdgCode %ld", pdgCode);
      return -1;
    }
    mass = particle->Mass();
    fMassRec=particle->Mass();
  }
  uint index1 = TMath::Min(indexStart,uint(fParamMC.size()-1));
  fMaxLayerRec=index1;
  for (;index1>0; index1--){
    if (fLoop[index1]>0) continue;                                   // loop number - only closet part of helix to the production vertex
    if (TMath::Abs(fParamMC[index1].GetSnp())>kMaxSnp) continue;     // track inclination anlge has to be smaller than kmaxSnp
    fMaxLayerRec=index1;
    if ((fParamMC[index1].P())>kMaxLoss*fParamMC[0].P()) break;      // if particle momneta is bigger than kMaxLoss production momenta - define starting index
  }

  if (index1<=3){
    ::Info("fastParticle::reconstructParticle","short track");
    return -1;
  }

  //AliExternalTrackParam4D param(fParamMC[index1],mass,1);
  Double_t xyzS[3][3];
  Int_t step=index1/3;
  if (step>10) step=10;
  float sign0=(fParamMC[index1-1].GetX()>fParamMC[index1-2].GetX())? 1.:-1.;
  Float_t alpha0=fParamMC[index1-1].GetAlpha();
  for (int dLayer=0; dLayer<3  && index1-step*dLayer>0; dLayer++) {
        Int_t index=index1-step*dLayer-1;

        int layer = fLayerIndex[index];
        //fParamMC[index1-step*dLayer-1].GetXYZ(xyzS[dLayer]);
        xyzS[dLayer][0]=fParamMC[index].GetX();
        xyzS[dLayer][1]=fParamMC[index].GetY()+gRandom->Gaus(0,geom.fLayerResolRPhi[layer]);;
        xyzS[dLayer][2]=fParamMC[index].GetZ()+gRandom->Gaus(0,geom.fLayerResolZ[layer]);;
        fParamMC[index].Local2GlobalPosition(xyzS[dLayer],fParamMC[index].GetAlpha()-alpha0);

  }
  /// seeds in alpha0 coordinate frame
  Int_t indexst = fLayerIndex[index1-1];
  AliExternalTrackParam * paramSeedI = fastTracker::makeSeed(xyzS[0],xyzS[1],xyzS[2],geom.fLayerResolRPhi[indexst],geom.fLayerResolZ[indexst],geom.fBz);
  AliExternalTrackParam * paramSeed = fastTracker::makeSeedMB(xyzS[0],xyzS[1],xyzS[2],geom.fLayerResolRPhi[indexst],geom.fLayerResolZ[indexst],
                                                              geom.fBz,geom.fLayerX0[indexst],geom.fLayerRho[indexst],fMassMC);
  AliExternalTrackParam   paramRot(paramSeed->GetX(),alpha0, paramSeed->GetParameter(),paramSeed->GetCovariance());
  AliExternalTrackParam4D param(paramRot,mass,1);
  if (sign0<0) {
    ((double*)param.GetParameter())[4]*=-1;
    ((double*)param.GetParameter())[3]*=-1;
    ((double*)param.GetParameter())[2]*=-1;
    
    ////Rotate Covariance accordingly RCR^T
    ((double*)param.GetCovariance())[3]*=-1;
    ((double*)param.GetCovariance())[4]*=-1;
    ((double*)param.GetCovariance())[6]*=-1;
    ((double*)param.GetCovariance())[7]*=-1;
    ((double*)param.GetCovariance())[10]*=-1;
    ((double*)param.GetCovariance())[11]*=-1;
  }

  //param.fMass=.fMass;

  //Double_t dEdx=AliExternalTrackParam::BetheBlochAleph(param.P()/mass);
  //Double_t dPdx=AliExternalTrackParam4D::dPdx(fParamMC[index1-1]);
  //  int version=0;
  // (*fgStreamer)<<"seedDump"<<   // seeding not ideal in case significant energy loss
  //   "version="<<version<<
  //   "gid="<<gid<<
  //   "sign0="<<sign0<<
  //   "fMassMC="<<fMassMC<<
  //   "dEdx="<<dEdx<<
  //   "step="<<step<<
  //   "seed.="<<&param<<
  //   "input.="<<&fParamMC[index1-1]<<
  //   "input1.="<<&fParamMC[index1-step-1]<<
  //   "input2.="<<&fParamMC[index1-2*step-1]<<
  //   "paramSeed.="<<paramSeed<<
  //   "paramSeedI.="<<paramSeedI<<
  //   "paramRot.="<<&paramRot<<
  //   "\n";
  delete paramSeed;
  delete paramSeedI;
  ///scale covaraince to remove double counting
  for (int i=0; i<15; i++) ((double*)param.GetCovariance())[i]*=kCovarFactor;
  //double *covar = (double*)param.GetCovariance();
  //for (int i=0; i<15; i++)covar[i]*=2;

//  Double_t resFactor=(geom.fLayerResolRPhi[0]*geom.fLayerResolRPhi[0]/0.01)+0.2; // this is hack - we will need to
//  double covar0[5]={0.1,0.1,(1.+2.*dEdx/param.P())/(1.+LArm),(1.+2.*dEdx/param.P())/(1.+LArm),resFactor*400*(1.+2.*dEdx/param.P())/(1.+LArm*LArm)};
//  covar[0]+=covar0[0]*covar0[0];
//  covar[2]+=covar0[1]*covar0[1];
//  covar[5]+=covar0[2]*covar0[2];
//  covar[9]+=covar0[3]*covar0[3];
//  covar[14]+=covar0[4]*covar0[4];
  //
  //float length=0, time=0;
  //float radius = sqrt(param.GetX()*param.GetX()+param.GetY()*param.GetY());
  fParamIn.resize(index1+1);
  fStatusMaskIn.resize(index1+1);
  fNPointsIn.resize(index1+1);
  fChi2.resize(index1+1);
  fParamIn[index1]=param;
  double xyz[3];
  int status=0;
  const double *par = param.GetParameter();
  for (int index=index1-1; index>=0; index--){   // dont propagate to vertex , will be done later ...

      Int_t layer = fLayerIndex[index];
      AliExternalTrackParam & p = fParamMC[index];
      p.GetXYZ(xyz);
      double alpha=TMath::ATan2(xyz[1],xyz[0]);
      double radius = TMath::Sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
      fStatusMaskIn[index]=0;
      fStatusMaskIn[index]|=kTrackEnter;
      if (radius>0) {
        status = param.Rotate(alpha);
      }
      if (status) {
        fStatusMaskIn[index]|=kTrackRotate;
      }else{
        //::Error("fastParticle::reconstructParticle:", "Rotation failed");
        break;
      }
      status = param.PropagateTo(radius,geom.fBz,1);
      if (status) {
        fStatusMaskIn[index]|=kTrackPropagate;
      }else{
        //::Error("fastParticle::reconstructParticle:", "Propagation failed");
        break;
      }
      float xrho  =geom.fLayerRho[layer];
      float xx0  =geom.fLayerX0[layer];
      float deltaY = gRandom->Gaus(0,geom.fLayerResolRPhi[layer]);
      float deltaZ = gRandom->Gaus(0,geom.fLayerResolZ[layer]);
      double pos[2]={0+deltaY,xyz[2]+deltaZ};
      double cov[3]={geom.fLayerResolRPhi[layer]*geom.fLayerResolRPhi[layer],0, geom.fLayerResolZ[layer]*geom.fLayerResolZ[layer]};
      fParamIn[index]=param;
      if (TMath::Abs(param.GetX()-fParamMC[index].GetX())>0.0001){ /// TODO  add float almost0
        ///problem
        //::Error("fastParticle::reconstructParticle:","ICONSISTENT X, should never happen %f\t%f",param.GetX(),fParamMC[index].GetX());
      }
      float chi2 =  param.GetPredictedChi2(pos, cov);
      fChi2[index]=chi2;
      if (chi2<chi2Cut) {
        fStatusMaskIn[index]|=kTrackChi2;
      }else{
        //::Error("fastParticle::reconstructParticle:", "Too big chi2 %f", chi2);
        break;
      }
      if (cov[0]>0) {
        status = param.Update(pos, cov);
        if (status) {
          fStatusMaskIn[index]|=kTrackUpdate;
        }else{
          //::Error("fastParticle::reconstructParticle:", "Update failed");
          break;
        }
      }
      else{
          //::Error("fastParticle::reconstructParticle:", "Update failed");
          break;
      }

      float tanPhi2 = par[2]*par[2];
      tanPhi2/=(1-tanPhi2);
      float crossLength=TMath::Sqrt(1.+tanPhi2+par[3]*par[3]);                /// geometrical path assuming crossing cylinder
      //status = param.AliExternalTrackParam::CorrectForMeanMaterial(crossLength*xx0,crossLength*xrho,mass);
      for (Int_t ic=0;ic<5; ic++) {
        status*= param.CorrectForMeanMaterial(crossLength * xx0/5., crossLength * xrho/5., mass, 0.01);
      }
      //status = param.CorrectForMeanMaterialT4(crossLength*xx0,crossLength*xrho,mass);
      //if (gRandom->Rndm() <fracUnitTest) param.UnitTestDumpCorrectForMaterial(fgStreamer,crossLength*xx0,crossLength*xrho,mass,20);
      if (status) {
        fStatusMaskIn[index]|=kTrackCorrectForMaterial;
      }else{
        //::Error("fastParticle::reconstructParticle:", "Correct for material failed");
        break;
      }
      fLengthIn++;
      fStatusMaskIn[index]|=kTrackisOK;
      fNPointsIn[index]=fLengthIn;
  }
  return 1;
}


///
/// \param geom          - pointer to geometry to use
/// \param pdgCode       - pdgCode used in the reconstruction
/// \param layerStart    - starting layer to do tracking
/// \return   -  TODO  status flags to be decides
int fastParticle::reconstructParticleFull(fastGeometry  &geom, long pdgCode, uint indexStart){
  const Float_t chi2Cut=10000/(geom.fLayerResolZ[0]);
  const float kMaxSnp=0.975;
  const float kMaxSnpStart=0.85;
  //const float kMaxLoss=0.3;
  const float kCovarFactor=1.5;
  const int   kMaxSkipped=500;
  const int   kSeedStep=500;
  fLengthIn=0;
  float_t mass=0;
  fPdgCodeRec   =pdgCode;
  if (pdgCode==0){
    fMassRec=fMassMC;
    mass=fMassMC;
  }else {
    TParticlePDG *particle = TDatabasePDG::Instance()->GetParticle(pdgCode);
    if (particle == nullptr) {
      //::Error("fastParticle::reconstructParticleFull", "Invalid pdgCode %ld", pdgCode);
      return -1;
    }
    mass = particle->Mass();
    fMassRec=particle->Mass();
  }
  uint index1 = TMath::Min(indexStart,uint(fParamMC.size()-1));
  fStatusMaskIn.resize(index1+1);
  fParamIn.resize(index1+1);
  fNPointsIn.resize(index1+1);
  
  for(uint i=0;i<=index1;i++) fStatusMaskIn[i]=0;
  /// skip layers with too big erregy loss - to smalle BG
  if(fUseMCInfo)
  {
    for (int i=index1; i>0; i--){  /// TODO - make query on fraction of the energy loss
      if (fParamMC[i].Beta()<0.05) {
        continue; /// TODO this is hack
      }
      index1=i;
      break;
    }
  }
  

  fMaxLayerRec=index1;

  if (index1<=3){
    ::Info("fastParticle::reconstructParticleFull","short track");
    return -1;
  }
  //Double_t LArm=getStat(0);
  //AliExternalTrackParam4D param(fParamMC[index1],mass,1);
  Double_t xyzS[3][3];
  Float_t alpha0;
  float sign0;
  float semiplane = -1;
  float sphi1 = 0.1;
  Int_t step=index1/3;
  if (step>kSeedStep) step=kSeedStep;

  while(semiplane<0 || sphi1>kMaxSnpStart)
  {
    sign0=(fParamMC[index1-1].GetX()>fParamMC[index1-1-step].GetX())? 1.:-1.;
    alpha0=fParamMC[index1-1].GetAlpha();

    for (int dLayer=0; dLayer<3  && index1-step*dLayer>0; dLayer++) {
          Int_t index=index1-step*dLayer-1;
          //Int_t index0=index1-1;
          int layer = fLayerIndex[index];
          //fParamMC[index1-step*dLayer-1].GetXYZ(xyzS[dLayer]);
          xyzS[dLayer][0]=fParamMC[index].GetX();
          xyzS[dLayer][1]=fParamMC[index].GetY();
          xyzS[dLayer][2]=fParamMC[index].GetZ();
          if(fAddPadsmearing)
          {
            xyzS[dLayer][1]+=gRandom->Gaus(0,geom.fLayerResolRPhi[layer]);
            xyzS[dLayer][2]+=gRandom->Gaus(0,geom.fLayerResolZ[layer]);
          }
          fParamMC[index].Local2GlobalPosition(xyzS[dLayer],fParamMC[index].GetAlpha()-alpha0);
    }

    sphi1 = TMath::Abs(fastTracker::makeSnp(xyzS[0][0],xyzS[0][1],xyzS[1][0],xyzS[1][1],xyzS[2][0],xyzS[2][1]));
    if(sphi1>kMaxSnpStart) {
      index1-=1;
      if((index1)<=3) break;
      step=index1/3;
      if (step>kSeedStep) step=kSeedStep;
      continue;
    }
    
    float sp0 = fastTracker::makeYC(xyzS[0][0],xyzS[0][1],xyzS[1][0],xyzS[1][1],xyzS[2][0],xyzS[2][1]);
    float sp2 = fastTracker::makeYC(xyzS[2][0],xyzS[2][1],xyzS[1][0],xyzS[1][1],xyzS[0][0],xyzS[0][1]);  
    semiplane = sp0/sp2; ///if <0 the two points are in different semiplanes

    if(semiplane<0) step-=1;
    if(step==0) 
    {
      index1-=1;
      if((index1)<=3) break;
      step=index1/3;
      if (step>kSeedStep) step=kSeedStep;
      continue;
    }
  }

  if (step==0 || (index1)<=3)
  {
    //::Error("fastParticle::reconstructParticleFull", "Too few consecutive points in same semiplane");
    return -1;
  }

  //float totLength=0;
  //float cross1=0; //clang complains about set but unused cross1 and cross2, length and totlength
  //float cross2=0;
  for(int l=index1-2*step; l<int(index1);l++)
  {
    Double_t xyznow[3];
    Double_t xyznext[3];
    fParamMC[l].GetXYZ(xyznow);
    fParamMC[l+1].GetXYZ(xyznext);
    //float Length=sqrt(pow((xyznow[0]-xyznext[0]),2)+pow((xyznow[1]-xyznext[1]),2)+pow((xyznow[2]-xyznext[2]),2));
    //totLength+=Length;
    //if(l<int(index1-step)) cross1+=Length;
    //else cross2+=Length;
  }

  /// seeds in alpha0 coordinate frame
  Int_t indexst = fLayerIndex[index1-1];
  // float crossing = 0;
  // std::vector<TVector3> conv_points_bkw;
  // for(int c=(index1-1); c>=0; c--){
  //    Double_t xyztest[3];
  //    fParamMC[c].GetXYZ(xyztest);
  //    TVector3 tr(xyztest[2], //z
  //   -xyztest[0]*sin(alpha0) + xyztest[1]*cos(alpha0),   //y
  //    xyztest[0]*cos(alpha0) + xyztest[1]*sin(alpha0));  //x
  //    conv_points_bkw.push_back(tr);
  // }
  //AliExternalTrackParam * partestbkw = fastTrackerGAR::Helix_Fit(conv_points_bkw,geom.fBz,geom.fLayerResolRPhi[indexst],geom.fLayerResolZ[indexst],-1,0);
  Double_t resolrphi = geom.fLayerResolRPhi[indexst];
  Double_t resolz = geom.fLayerResolZ[indexst];
  // if(resolrphi<0.001 || resolrphi>10 || isnan(resolrphi) || resolz<0.001 || resolz>10 || isnan(resolz)) {
  //   resolrphi=geom.fLayerResolRPhi[indexst];
  //   resolz=geom.fLayerResolZ[indexst];
  // }else{
  //   resolrphi = resolrphi;
  //   resolz =resolz;
  // }

  AliExternalTrackParam * paramSeedI = fastTracker::makeSeed(xyzS[0],xyzS[1],xyzS[2],resolrphi,resolz,geom.fBz);
  AliExternalTrackParam * paramSeed = fastTracker::makeSeedMB(xyzS[0],xyzS[1],xyzS[2],resolrphi,resolz,
                                                              geom.fBz,geom.fLayerX0[indexst],geom.fLayerRho[indexst],mass,5);
  /*
  std::cout<<"xyzS[0]: "<<xyzS[0][0]<<" "<<xyzS[0][1]<<" "<<xyzS[0][2]<<std::endl;
  std::cout<<"xyzS[1]: "<<xyzS[1][0]<<" "<<xyzS[1][1]<<" "<<xyzS[1][2]<<std::endl;
  std::cout<<"xyzS[2]: "<<xyzS[2][0]<<" "<<xyzS[2][1]<<" "<<xyzS[2][2]<<std::endl;
  std::cout<<"resolrphi: "<<resolrphi<<std::endl;
  std::cout<<"resolz: "<<resolz<<std::endl;
  std::cout<<"geom.fBz: "<<geom.fBz<<std::endl;
  std::cout<<"geom.fLayerX0: "<<geom.fLayerX0[indexst]<<std::endl;
  std::cout<<"geom.fLayerRho: "<<geom.fLayerRho[indexst]<<std::endl;
  std::cout<<"fMassMC: "<<fMassMC<<std::endl;
  */

  //AliExternalTrackParam * paramSeed4D = fastTracker::makeSeedMB4(xyzS[0],xyzS[1],xyzS[2],resolrphi,resolz,
  //                                                            geom.fBz,geom.fLayerX0[indexst],geom.fLayerRho[indexst],fMassMC,1,fAddMSsmearing);
  AliExternalTrackParam   paramRot(paramSeed->GetX(),alpha0, paramSeed->GetParameter(),paramSeed->GetCovariance());
  AliExternalTrackParam4D param(paramRot,mass,1);

  if (sign0<0) {
    ((double*)param.GetParameter())[4]*=-1;
    ((double*)param.GetParameter())[3]*=-1;
    ((double*)param.GetParameter())[2]*=-1;
    
    ////Rotate Covariance accordingly RCR^T
    ((double*)param.GetCovariance())[3]*=-1;
    ((double*)param.GetCovariance())[4]*=-1;
    ((double*)param.GetCovariance())[6]*=-1;
    ((double*)param.GetCovariance())[7]*=-1;
    ((double*)param.GetCovariance())[10]*=-1;
    ((double*)param.GetCovariance())[11]*=-1;
  }
  //param.fMass=.fMass;
  //Double_t dEdx=AliExternalTrackParam::BetheBlochAleph(param.P()/mass);
  //float LArmSeed = sqrt((xyzS[2][0]-xyzS[0][0])*(xyzS[2][0]-xyzS[0][0])+(xyzS[2][1]-xyzS[0][1])*(xyzS[2][1]-xyzS[0][1]));
  //float deltaY = abs(xyzS[2][1]);
  //Double_t dPdx=AliExternalTrackParam4D::dPdx(fParamMC[index1-1]);
  //int version=1;
  // (*fgStreamer)<<"seedDump"<<   // seeding not ideal in case significant energy loss
  //   "version="<<version<<       // version 1 leg
  //   "indexst="<<indexst<<
  //   "gid="<<gid<<
  //   "sign0="<<sign0<<
  //   "LArm="<<LArmSeed<<
  //   "totLength="<<totLength<<
  //   "fMassMC="<<fMassMC<<
  //   "dEdx="<<dEdx<<
  //   "step="<<step<<
  //   "seed.="<<&param<<
  //   "input.="<<&fParamMC[index1-1]<<
  //   "input1.="<<&fParamMC[index1-step-1]<<
  //   "input2.="<<&fParamMC[index1-2*step-1]<<
  //   "paramSeed.="<<paramSeed<<
  //   "paramSeedI.="<<paramSeedI<<
  //   "paramRot.="<<&paramRot<<
  //   "\n";
  delete paramSeed;
  delete paramSeedI;
  ///scale covaraine to avoid double counting
  for (int i=0; i<15; i++) ((double*)param.GetCovariance())[i]*=kCovarFactor;

  //float length=0, time=0;
  //float radius = sqrt(param.GetX()*param.GetX()+param.GetY()*param.GetY());
  fFirstIndexIn = index1;
  fParamIn.resize(index1+1);
  fStatusMaskIn.resize(index1+1);
  fNPointsIn.resize(index1+1);
  fChi2.resize(index1+1);
  fParamIn[index1]=param;
  double xyz[3];
  int status=0;
  //const double *par = param.GetParameter();
  int checkloop=0;
  int mirrored =0;
  //Bool_t restarted = kFALSE;
  //Bool_t Propagate_Failed = kFALSE;  ///Used to avoid to PropagateToMirrorX after Propagate failed twice consecutively
  for (int index=index1-1; index>=0; index--){   // dont propagate to vertex , will be done later ...
      //if(fParamMC.size()==111)std::cout<<"Parameter: "<<param.GetParameter()[0]<<" "<<param.GetParameter()[1]<<" "<<param.GetParameter()[2]<<" "<<param.GetParameter()[3]<<" "<<param.GetParameter()[4]<<std::endl;
      //Bool_t Propagate_First = kFALSE;
      Bool_t SkipUpdate = kFALSE;
      //double resol=0;
      float crossLength = 0;
      Int_t layer = fLayerIndex[index];
      AliExternalTrackParam & p = fParamMC[index];
      p.GetXYZ(xyz);
      double alpha=TMath::ATan2(xyz[1],xyz[0]);
      double radius = TMath::Sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
      fStatusMaskIn[index]|=kTrackEnter;

      if(checkloop==0) 
      {
        if(fUseMCInfo) checkloop = fLoop[index]-fLoop[index+1];  /////// PropagateToMirror triggered for now using flag from MC information: not realistic reconstruction
        else if (index>20)
        {
          if (mirrored==0)
          {
            checkloop = (TMath::Abs(param.GetParameter()[2])>(kMaxSnp))? 1:0;
          }
        }
        else if (TMath::Abs(param.GetParameter()[2])>0.999) break;
      }
      if(mirrored>0) mirrored--;
      
      if (checkloop!=0){   
          mirrored=fParamMC.size()/6;
          float dir = 0;
          for (size_t d=1;d<TMath::Min((fDirection.size()-index),size_t(20));d++) dir+=fDirection[index+d];
          //dir = - fDirection[index+1];
          dir = dir>0? -1:1;
          crossLength = param.PropagateToMirrorX(geom.fBz, dir, geom.fLayerResolRPhi[layer],geom.fLayerResolZ[layer]);  ////Using direction from MC, not realistic reconstruction
          if (crossLength<0)
          {
            //::Error("fastParticle::reconstructParticleFull:", "PropagateToMirrorX failed");
            break;
          }      
          ///Find closest point in X after PropagateToMirrorX    
          int Skip = TMath::Min(kMaxSkipped,index);
          float dxyz_min = 9999;
          float dx_min = 9999;
          int  new_index = 0;
          for(int n=0; n<Skip; n++)
          {
            //Int_t layer_m = fLayerIndex[index-n];
            float dx_m = TMath::Abs(fParamMC[index-n].GetX()-param.GetX());
            if(dx_m<dx_min)
            {
              dx_min=dx_m;
              if(fUseMCInfo) new_index=index-n;
            }
            
            Double_t xyzMC[3];
            fParamMC[index-n].GetXYZ(xyzMC);
            Double_t xyzparam[3];
            param.GetXYZ(xyzparam);
            float dxyz_m = sqrt((xyzMC[0]-xyzparam[0])*(xyzMC[0]-xyzparam[0])+(xyzMC[1]-xyzparam[1])*(xyzMC[1]-xyzparam[1])+(xyzMC[2]-xyzparam[2])*(xyzMC[2]-xyzparam[2]));
            if(dxyz_m<dxyz_min)
            {
              dxyz_min=dxyz_m;
              if(!fUseMCInfo) new_index=index-n;
            }
          }
          //fLengthIn+=1+TMath::Abs(new_index-index);
          for(Int_t k=index;k>new_index;k--) fStatusMaskIn[k]|=kTrackSkip;
          index=new_index;
          fParamIn[index]=param;
          fStatusMaskIn[index]|=kTrackPropagatetoMirrorX;
          fStatusMaskIn[index]|=kTrackisOK;
          fNPointsIn[index]=fLengthIn;
          checkloop=0;
          continue;
      }
      else{

          if (radius>0) {
            status = param.Rotate(alpha);
          }
          if (status) {
            fStatusMaskIn[index]|=kTrackRotate;
          }else if(mirrored==0 && index>15){
              ///If propagation fails go back a step -> PropagateToMirrorX(); If this has already been tried in the previous step, propagation has failed.          
              param=fParamIn[index+1];
              index++;
              checkloop=2;
              continue;
          }
          else {
            int Skip = TMath::Min(kMaxSkipped-450,index);
            int  new_index = 0;
            AliExternalTrackParam4D paramtemp = param;
            for(int n=0; n<Skip; n++)
            {
                AliExternalTrackParam & p_sk = fParamMC[index-n];
                double xyz_sk[3];
                p_sk.GetXYZ(xyz_sk);
                double alpha_sk=TMath::ATan2(xyz_sk[1],xyz_sk[0]);
                //double radius_sk = TMath::Sqrt(xyz_sk[0]*xyz_sk[0]+xyz_sk[1]*xyz_sk[1]);
                status = paramtemp.Rotate(alpha_sk);
                if(status)
                {
                  new_index=index+1-n;
                  break;
                }
            }
            if(status){
              index=new_index;
              continue;
            }
            else{
              ////::Error("fastParticle::reconstructParticleFull:", "Rotation failed");
              break;
            }
          }
          status = param.PropagateTo(radius,geom.fBz,1);
          if (status) {
            fStatusMaskIn[index]|=kTrackPropagate;
            //if(Propagate_Failed) Propagate_Failed=kFALSE;
          }
          else if(mirrored==0 && index>15){
              ///If propagation fails go back a step -> PropagateToMirrorX(); If this has already been tried in the previous step, propagation has failed.          
              param=fParamIn[index+1];
              index++;
              checkloop=2;
              continue;
          }
          else 
          {
            int Skip = TMath::Min(kMaxSkipped-450,index);
            int  new_index = 0;
            for(int n=0; n<Skip; n++)
            {
                AliExternalTrackParam4D paramtemp = param;
                AliExternalTrackParam & p_sk = fParamMC[index-n];
                double xyz_sk[3];
                p_sk.GetXYZ(xyz_sk);
                double alpha_sk=TMath::ATan2(xyz_sk[1],xyz_sk[0]);
                double radius_sk = TMath::Sqrt(xyz_sk[0]*xyz_sk[0]+xyz_sk[1]*xyz_sk[1]);
                bool status1 = paramtemp.Rotate(alpha_sk);
                bool status2 = paramtemp.PropagateTo(radius_sk,geom.fBz,1);
                if(status1 && status2)
                {
                  status=1;
                  new_index=index+1-n;
                  break;
                }
            }
            if(status){
              index=new_index;
              continue;
            }
            else{
              //::Error("fastParticle::reconstructParticleFull:", "Propagation failed");
              //std::cout<<"Parameter: "<<param.GetParameter()[0]<<" "<<param.GetParameter()[1]<<" "<<param.GetParameter()[2]<<" "<<param.GetParameter()[3]<<" "<<param.GetParameter()[4]<<std::endl;
              break;
            }
          }
        }
      // }
      float xrho  =geom.fLayerRho[layer];
      float xx0  =geom.fLayerX0[layer];
      float deltaY = gRandom->Gaus(0,geom.fLayerResolRPhi[layer]);
      float deltaZ = gRandom->Gaus(0,geom.fLayerResolZ[layer]);
      double pos[2]={0,xyz[2]};
      if(fAddPadsmearing)
      {
        pos[0]+=deltaY;
        pos[1]+=deltaZ;
      }
      double cov[3]={geom.fLayerResolRPhi[layer]*geom.fLayerResolRPhi[layer],0, geom.fLayerResolZ[layer]*geom.fLayerResolZ[layer]};
      //fParamIn[index]=param;
      float chi2 =  param.GetPredictedChi2(pos, cov);
      fChi2[index]=chi2;
      if (chi2<chi2Cut) {
        fStatusMaskIn[index]|=kTrackChi2;
      }else{
            //::Error("fastParticle::reconstructParticleFull:", "Too big chi2 %f", chi2);
            break;
      }
      
      if(chi2>200 && mirrored==0 && abs(param.GetParameter()[2])>0.8)
      {
        param=fParamIn[index+1];
        index++;
        checkloop=2;
        continue;
      }

      if (cov[0]>0) {
        status = param.Update(pos, cov);
        if (status) {
          fStatusMaskIn[index]|=kTrackUpdate;
        }else{
            ///skip the Update
            fStatusMaskIn[index]|=kTrackSkipUpdate;
            SkipUpdate = kTRUE;
        }
      }
      else{
            //::Error("fastParticle::reconstructParticleFull:", "Update failed");
            break;
      }
      

      //float tanPhi2 = par[2]*par[2];
      float tanPhi2 = param.GetParameter()[2]*param.GetParameter()[2];
      tanPhi2/=(1-tanPhi2);
      float deltaR = 1;
      if(fParamIn[index+1].GetX()!=0) deltaR=abs(fParamIn[index+1].GetX()-param.GetX());
      Double_t xyz_prev[3];
      Double_t xyz_now[3];
      fParamMC[index+1].GetXYZ(xyz_prev);
      fParamMC[index].GetXYZ(xyz_now);
      //if(crossLength==0) crossLength=deltaR*TMath::Sqrt(1.+tanPhi2+param.GetParameter()[3]*param.GetParameter()[3]);                /// geometrical path assuming crossing cylinder
      if(crossLength==0) crossLength=deltaR*TMath::Sqrt((1.+tanPhi2)*(1+param.GetParameter()[3]*param.GetParameter()[3]));
      //if(crossLength==0) crossLength=TMath::Sqrt(pow(xyz_prev[0]-xyz_now[0],2)+pow(xyz_prev[1]-xyz_now[1],2)+pow(xyz_prev[2]-xyz_now[2],2));
      //param.SetUseLogTermMS(kTRUE);
      //status = param.AliExternalTrackParam::CorrectForMeanMaterial(crossLength*xx0,crossLength*xrho,mass);
      if(!SkipUpdate && crossLength>0)
      {
         for (Int_t ic=0;ic<5; ic++) {
           status*= param.CorrectForMeanMaterial(crossLength * xx0/5., crossLength * xrho/5., mass, 0.01, 0+0x2*fAddMSsmearing);
         }
        //status*= param.CorrectForMeanMaterial(crossLength * xx0, crossLength * xrho, mass, 0.01, 0+0x2*fAddMSsmearing);
        //status = param.CorrectForMeanMaterialT4(crossLength*xx0,crossLength*xrho,mass);
        //if (gRandom->Rndm() <fracUnitTest) param.UnitTestDumpCorrectForMaterial(fgStreamer,crossLength*xx0,crossLength*xrho,mass,20);
        if (status) {
          fStatusMaskIn[index]|=kTrackCorrectForMaterial;
        }else{
          //::Error("fastParticle::reconstructParticleFull:", "Correct for material failed");
          break;
        }
      }

      fParamIn[index]=param;
      fLengthIn++;
      fNPointsIn[index]=fLengthIn;
      fStatusMaskIn[index]|=kTrackisOK;
  }
  if(fParamMC.size()!=fParamIn.size())
  {
    fParamIn.resize(fParamMC.size());
    fStatusMaskIn.resize(fParamMC.size());
    fNPointsIn.resize(fParamMC.size());
  }
  return 1;
  
}


///
/// \param geom          - pointer to geometry to use
/// \param pdgCode       - pdgCode used in the reconstruction
/// \param layerStart    - starting layer to do tracking
/// \return   -  TODO  status flags to be decides
int fastParticle::reconstructParticleFullOut(fastGeometry  &geom, long pdgCode, uint lastPoint){
  const Float_t chi2Cut=10000/(geom.fLayerResolZ[0]);
  const float kMaxSnp=0.975;
  const float kMaxSnpStart=0.85;
  //const float kMaxLoss=0.3;
  const float kCovarFactor=1.5;
  const int   kMaxSkipped=500;
  const int   kSeedStep=500;
  fLengthOut=0;
  float_t mass=0;
  fPdgCodeRec   =pdgCode;
  if (pdgCode==0){
    fMassRec=fMassMC;
    mass=fMassMC;
  }else {
    TParticlePDG *particle = TDatabasePDG::Instance()->GetParticle(pdgCode);
    if (particle == nullptr) {
      //::Error("fastParticle::reconstructParticleFullOut", "Invalid pdgCode %ld", pdgCode);
      return -1;
    }
    mass = particle->Mass();
    fMassRec=particle->Mass();
  }
  uint indexlast = TMath::Min(lastPoint,uint(fParamMC.size()-1));
  fStatusMaskOut.resize(indexlast+1);
  fParamOut.resize(indexlast+1);
  fNPointsOut.resize(indexlast+1);
  
  for(uint i=0;i<=indexlast;i++) fStatusMaskOut[i]=0;
  /// skip layers with too big erregy loss - to smalle BG
  if(fUseMCInfo)
  {
    for (int i=indexlast; i>0; i--){  /// TODO - make query on fraction of the energy loss
      if (fParamMC[i].Beta()<0.05) {
        continue; /// TODO this is hack
      }
      indexlast=i;
      break;
    }
  }
  

  fMaxLayerRec=indexlast;

  if (indexlast<=3){
    ::Info("fastParticle::reconstructParticleFullOut","short track");
    return -1;
  }
  //Double_t LArm=getStat(0);
  //AliExternalTrackParam4D param(fParamMC[index1],mass,1);
  Double_t xyzS[3][3];
  Float_t alpha0;
  float sign0;
  uint index1 = 0;
  float semiplane = -1;
  float sphi1 = 0.1;
  Int_t step=indexlast/3;
  if (step>kSeedStep) step=kSeedStep;

  while(semiplane<0 || sphi1>kMaxSnpStart)
  {
    sign0=(fParamMC[index1+1].GetX()>fParamMC[index1+1+step].GetX())? 1.:-1.;
    alpha0=fParamMC[index1+1].GetAlpha();

    for (int dLayer=0; dLayer<3 ; dLayer++) {
          Int_t index=step*dLayer+index1+1;
          //Int_t index0=index1+1;
          int layer = fLayerIndex[index];
          //fParamMC[index1-step*dLayer-1].GetXYZ(xyzS[dLayer]);
          xyzS[dLayer][0]=fParamMC[index].GetX();
          xyzS[dLayer][1]=fParamMC[index].GetY();
          xyzS[dLayer][2]=fParamMC[index].GetZ();
          if(fAddPadsmearing)
          {
            xyzS[dLayer][1]+=gRandom->Gaus(0,geom.fLayerResolRPhi[layer]);
            xyzS[dLayer][2]+=gRandom->Gaus(0,geom.fLayerResolZ[layer]);
          }
          fParamMC[index].Local2GlobalPosition(xyzS[dLayer],fParamMC[index].GetAlpha()-alpha0);
    }

    sphi1 = TMath::Abs(fastTracker::makeSnp(xyzS[0][0],xyzS[0][1],xyzS[1][0],xyzS[1][1],xyzS[2][0],xyzS[2][1]));
    if(sphi1>kMaxSnpStart) {
      index1+=1;
      if((indexlast-index1)<=3) break;
      step=(indexlast-index1)/3;
      if (step>kSeedStep) step=kSeedStep;
      continue;
    }
    
    float sp0 = fastTracker::makeYC(xyzS[0][0],xyzS[0][1],xyzS[1][0],xyzS[1][1],xyzS[2][0],xyzS[2][1]);
    float sp2 = fastTracker::makeYC(xyzS[2][0],xyzS[2][1],xyzS[1][0],xyzS[1][1],xyzS[0][0],xyzS[0][1]);  
    semiplane = sp0/sp2; ///if <0 the two points are in different semiplanes

    if(semiplane<0) step-=1;
    if(step==0) 
    {
      index1+=1;
      if((indexlast-index1)<=3) break;
      step=(indexlast-index1)/3;
      if (step>kSeedStep) step=kSeedStep;
      continue;
    }
  }

  if (step==0 || (indexlast-index1)<=3)
  {
    //::Error("fastParticle::reconstructParticleFullOut", "Too few consecutive points in same semiplane");
    return -1;
  }

  //float totLength=0;
  //float cross1=0; //clang complains about set but unused cross1 and cross2
  //float cross2=0;
  for(int l=index1; l<int(index1+2*step);l++)
  {
    Double_t xyznow[3];
    Double_t xyznext[3];
    fParamMC[l].GetXYZ(xyznow);
    fParamMC[l+1].GetXYZ(xyznext);
    //float Length=sqrt(pow((xyznow[0]-xyznext[0]),2)+pow((xyznow[1]-xyznext[1]),2)+pow((xyznow[2]-xyznext[2]),2));
    //totLength+=Length;
    //if(l<int(index1+step)) cross1+=Length;
    //else cross2+=Length;
  }

  /// seeds in alpha0 coordinate frame
  Int_t indexst = fLayerIndex[index1+1];
  Double_t resolrphi = geom.fLayerResolRPhi[indexst];
  Double_t resolz = geom.fLayerResolZ[indexst];

  AliExternalTrackParam * paramSeedI = fastTracker::makeSeed(xyzS[0],xyzS[1],xyzS[2],resolrphi,resolz,geom.fBz);
  AliExternalTrackParam * paramSeed = fastTracker::makeSeedMB(xyzS[0],xyzS[1],xyzS[2],resolrphi,resolz,
                                                              geom.fBz,geom.fLayerX0[indexst],geom.fLayerRho[indexst],mass,5);
  //AliExternalTrackParam * paramSeed4D = fastTracker::makeSeedMB4(xyzS[0],xyzS[1],xyzS[2],resolrphi,resolz,
  //                                                            geom.fBz,geom.fLayerX0[indexst],geom.fLayerRho[indexst],fMassMC,1,fAddMSsmearing);
  AliExternalTrackParam   paramRot(paramSeed->GetX(),alpha0, paramSeed->GetParameter(),paramSeed->GetCovariance());
  AliExternalTrackParam4D param(paramRot,mass,1);

  if (sign0<0) {
    ((double*)param.GetParameter())[4]*=-1;
    ((double*)param.GetParameter())[3]*=-1;
    ((double*)param.GetParameter())[2]*=-1;
    
    ////Rotate Covariance accordingly RCR^T
    ((double*)param.GetCovariance())[3]*=-1;
    ((double*)param.GetCovariance())[4]*=-1;
    ((double*)param.GetCovariance())[6]*=-1;
    ((double*)param.GetCovariance())[7]*=-1;
    ((double*)param.GetCovariance())[10]*=-1;
    ((double*)param.GetCovariance())[11]*=-1;
  }
  //param.fMass=.fMass;
  //Double_t dEdx=AliExternalTrackParam::BetheBlochAleph(param.P()/mass);
  //float LArmSeed = sqrt((xyzS[2][0]-xyzS[0][0])*(xyzS[2][0]-xyzS[0][0])+(xyzS[2][1]-xyzS[0][1])*(xyzS[2][1]-xyzS[0][1]));
  //float deltaY = abs(xyzS[2][1]);
  //Double_t dPdx=AliExternalTrackParam4D::dPdx(fParamMC[index1-1]);
  // int version=2;
  // (*fgStreamer)<<"seedDumpOut"<<   // seeding not ideal in case significant energy loss
  //   "version="<<version<<       // version 1 leg
  //   "indexst="<<indexst<<
  //   "gid="<<gid<<
  //   "sign0="<<sign0<<
  //   "LArm="<<LArmSeed<<
  //   "totLength="<<totLength<<
  //   "fMassMC="<<fMassMC<<
  //   "dEdx="<<dEdx<<
  //   "step="<<step<<
  //   "seed.="<<&param<<
  //   "input.="<<&fParamMC[index1+1]<<
  //   "input1.="<<&fParamMC[index1+step+1]<<
  //   "input2.="<<&fParamMC[index1+2*step+1]<<
  //   "paramSeed.="<<paramSeed<<
  //   "paramSeedI.="<<paramSeedI<<
  //   "paramRot.="<<&paramRot<<
  //   "\n";
  delete paramSeed;
  delete paramSeedI;
  ///scale covaraine to avoid double counting
  for (int i=0; i<15; i++) ((double*)param.GetCovariance())[i]*=kCovarFactor;

  //float length=0, time=0;
  //float radius = sqrt(param.GetX()*param.GetX()+param.GetY()*param.GetY());
  fFirstIndexOut = index1;
  fParamOut.resize(indexlast+1);
  fStatusMaskOut.resize(indexlast+1);
  fNPointsOut.resize(indexlast+1);
  fChi2Out.resize(indexlast+1);
  fParamOut[index1]=param;
  double xyz[3];
  int status=0;
  //const double *par = param.GetParameter();
  int checkloop=0;
  int mirrored =0;
  //Bool_t restarted = kFALSE;
  //Bool_t Propagate_Failed = kFALSE;  ///Used to avoid to PropagateToMirrorX after Propagate failed twice consecutively
  for (int index=index1+1; index<=int(indexlast); index++){   // dont propagate to vertex , will be done later ...
      //Bool_t Propagate_First = kFALSE;
      Bool_t SkipUpdate = kFALSE;
      //double resol=0;
      float crossLength = 0;
      Int_t layer = fLayerIndex[index];
      AliExternalTrackParam & p = fParamMC[index];
      p.GetXYZ(xyz);
      double alpha=TMath::ATan2(xyz[1],xyz[0]);
      double radius = TMath::Sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
      fStatusMaskOut[index]|=kTrackEnter;

      if(checkloop==0) 
      {
        if(fUseMCInfo) checkloop = fLoop[index]-fLoop[index+1];  /////// PropagateToMirror triggered for now using flag from MC information: not realistic reconstruction
        else if (index<(int(indexlast)-20))
        {
          if (mirrored==0)
          {
            checkloop = (TMath::Abs(param.GetParameter()[2])>(kMaxSnp))? 1:0;
          }
        }
        else if (TMath::Abs(param.GetParameter()[2])>0.999) break;
      }
      if(mirrored>0) mirrored--;
      
      if (checkloop!=0){   
          mirrored=fParamMC.size()/6;
          float dir = 0;
          for (size_t d=1;d<TMath::Min(size_t(index),size_t(20));d++) dir+=fDirection[index-d];
          //dir = - fDirection[index+1];
          dir = dir>0? 1:-1;
          crossLength = param.PropagateToMirrorX(geom.fBz, dir, geom.fLayerResolRPhi[layer],geom.fLayerResolZ[layer]);  ////Using direction from MC, not realistic reconstruction
          if (crossLength<0)
          {
            //::Error("fastParticle::reconstructParticleFullOut:", "PropagateToMirrorX failed");
            break;
          }      
          ///Find closest point in X after PropagateToMirrorX    
          int Skip = TMath::Min(kMaxSkipped,int(indexlast-index));
          float dxyz_min = 9999;
          float dx_min = 9999;
          int  new_index = 0;
          for(int n=0; n<Skip; n++)
          {
            //Int_t layer_m = fLayerIndex[index+n];
            float dx_m = TMath::Abs(fParamMC[index+n].GetParameter()[1]-param.GetX());
            if(dx_m<dx_min)
            {
              dx_min=dx_m;
              if(fUseMCInfo) new_index=index+n;
            }
            
            Double_t xyzMC[3];
            fParamMC[index+n].GetXYZ(xyzMC);
            Double_t xyzparam[3];
            param.GetXYZ(xyzparam);
            float dxyz_m = sqrt((xyzMC[0]-xyzparam[0])*(xyzMC[0]-xyzparam[0])+(xyzMC[1]-xyzparam[1])*(xyzMC[1]-xyzparam[1])+(xyzMC[2]-xyzparam[2])*(xyzMC[2]-xyzparam[2]));
            if(dxyz_m<dxyz_min)
            {
              dxyz_min=dxyz_m;
              if(!fUseMCInfo) new_index=index+n;
            }
          }
          //fLengthOut+=1+TMath::Abs(new_index-index);
          for(Int_t k=index;k<new_index;k++) fStatusMaskOut[k]|=kTrackSkip;
          index=new_index;
          fParamOut[index]=param;
          fStatusMaskOut[index]|=kTrackPropagatetoMirrorX;
          fStatusMaskOut[index]|=kTrackisOK;
          fNPointsOut[index]=fLengthOut;
          checkloop=0;
          continue;
      }
      else{

          if (radius>0) {
            status = param.Rotate(alpha);
          }
          if (status) {
            fStatusMaskOut[index]|=kTrackRotate;
          }else if(mirrored==0 && index<(int(indexlast)-15)){
              ///If propagation fails go back a step -> PropagateToMirrorX(); If this has already been tried in the previous step, propagation has failed.          
              param=fParamOut[index-1];
              index--;
              checkloop=2;
              continue;
          }
          else {
            int Skip = TMath::Min(kMaxSkipped-450,int(indexlast-index));
            int  new_index = index;
            AliExternalTrackParam4D paramtemp = param;
            for(int n=0; n<Skip; n++)
            {
                AliExternalTrackParam & p_sk = fParamMC[index+n];
                double xyz_sk[3];
                p_sk.GetXYZ(xyz_sk);
                double alpha_sk=TMath::ATan2(xyz_sk[1],xyz_sk[0]);
                //double radius_sk = TMath::Sqrt(xyz_sk[0]*xyz_sk[0]+xyz_sk[1]*xyz_sk[1]);
                status = paramtemp.Rotate(alpha_sk);
                if(status)
                {
                  new_index=index+n;
                  break;
                }
            }
            if(status){
              index=new_index;
              continue;
            }
            else{
              //::Error("fastParticle::reconstructParticleFullOut:", "Rotation failed");
              break;
            }
          }
          status = param.PropagateTo(radius,geom.fBz,1);
          if (status) {
            fStatusMaskOut[index]|=kTrackPropagate;
            //if(Propagate_Failed) Propagate_Failed=kFALSE;
          }
          else if(mirrored==0 && index<int(indexlast-15)){
              ///If propagation fails go back a step -> PropagateToMirrorX(); If this has already been tried in the previous step, propagation has failed.          
              param=fParamOut[index-1];
              index--;
              checkloop=2;
              continue;
          }
          else 
          {
            int Skip = TMath::Min(kMaxSkipped-450,int(indexlast-index));
            int  new_index = index;
            for(int n=0; n<Skip; n++)
            {
                AliExternalTrackParam4D paramtemp = param;
                AliExternalTrackParam & p_sk = fParamMC[index+n];
                double xyz_sk[3];
                p_sk.GetXYZ(xyz_sk);
                double alpha_sk=TMath::ATan2(xyz_sk[1],xyz_sk[0]);
                double radius_sk = TMath::Sqrt(xyz_sk[0]*xyz_sk[0]+xyz_sk[1]*xyz_sk[1]);
                bool status1 = paramtemp.Rotate(alpha_sk);
                bool status2 = paramtemp.PropagateTo(radius_sk,geom.fBz,1);
                if(status1 && status2)
                {
                  status=1;
                  new_index=index+n;
                  break;
                }
            }
            if(status){
              index=new_index;
              continue;
            }
            else{
              //::Error("fastParticle::reconstructParticleFullOut:", "Propagation failed");
              break;
            }
          }
        }
      // }
      float xrho  =geom.fLayerRho[layer];
      float xx0  =geom.fLayerX0[layer];
      float deltaY = gRandom->Gaus(0,geom.fLayerResolRPhi[layer]);
      float deltaZ = gRandom->Gaus(0,geom.fLayerResolZ[layer]);
      double pos[2]={0,xyz[2]};
      if(fAddPadsmearing)
      {
        pos[0]+=deltaY;
        pos[1]+=deltaZ;
      }
      double cov[3]={geom.fLayerResolRPhi[layer]*geom.fLayerResolRPhi[layer],0, geom.fLayerResolZ[layer]*geom.fLayerResolZ[layer]};
      //fParamIn[index]=param;
      float chi2 =  param.GetPredictedChi2(pos, cov);
      fChi2Out[index]=chi2;
      if (chi2<chi2Cut) {
        fStatusMaskOut[index]|=kTrackChi2;
      }else{
            //::Error("fastParticle::reconstructParticleFullOut:", "Too big chi2 %f", chi2);
            break;
      }
      
      if(chi2>200 && mirrored==0 && abs(param.GetParameter()[2])>0.8)
      {
        param=fParamOut[index-1];
        index--;
        checkloop=2;
        continue;
      }

      if (cov[0]>0) {
        status = param.Update(pos, cov);
        if (status) {
          fStatusMaskOut[index]|=kTrackUpdate;
        }else{
            ///skip the Update
            fStatusMaskOut[index]|=kTrackSkipUpdate;
            SkipUpdate = kTRUE;
        }
      }
      else{
            //::Error("fastParticle::reconstructParticleFullOut:", "Update failed");
            break;
      }
      

      //float tanPhi2 = par[2]*par[2];
      float tanPhi2 = param.GetParameter()[2]*param.GetParameter()[2];
      tanPhi2/=(1-tanPhi2);
      float deltaR = 1;
      if(fParamOut[index-1].GetX()!=0) deltaR=abs(fParamOut[index-1].GetX()-param.GetX());
      Double_t xyz_prev[3];
      Double_t xyz_now[3];
      fParamMC[index-1].GetXYZ(xyz_prev);
      fParamMC[index].GetXYZ(xyz_now);
      //if(crossLength==0) crossLength=deltaR*TMath::Sqrt(1.+tanPhi2+param.GetParameter()[3]*param.GetParameter()[3]);                /// geometrical path assuming crossing cylinder
      if(crossLength==0) crossLength=deltaR*TMath::Sqrt((1.+tanPhi2)*(1+param.GetParameter()[3]*param.GetParameter()[3]));
      //if(crossLength==0) crossLength=TMath::Sqrt(pow(xyz_prev[0]-xyz_now[0],2)+pow(xyz_prev[1]-xyz_now[1],2)+pow(xyz_prev[2]-xyz_now[2],2));
      //param.SetUseLogTermMS(kTRUE);
      //status = param.AliExternalTrackParam::CorrectForMeanMaterial(crossLength*xx0,crossLength*xrho,mass);
      if(!SkipUpdate && crossLength>0)
      {
         for (Int_t ic=0;ic<5; ic++) {
           status*= param.CorrectForMeanMaterial(crossLength * xx0/5., -crossLength * xrho/5., mass, 0.01, 0+0x2*fAddMSsmearing);
         }
        //status*= param.CorrectForMeanMaterial(crossLength * xx0, crossLength * xrho, mass, 0.01, 0+0x2*fAddMSsmearing);
        //status = param.CorrectForMeanMaterialT4(crossLength*xx0,crossLength*xrho,mass);
        //if (gRandom->Rndm() <fracUnitTest) param.UnitTestDumpCorrectForMaterial(fgStreamer,crossLength*xx0,crossLength*xrho,mass,20);
        if (status) {
          fStatusMaskOut[index]|=kTrackCorrectForMaterial;
        }else{
          //::Error("fastParticle::reconstructParticleFullOut:", "Correct for material failed");
          break;
        }
      }

      fParamOut[index]=param;
      fLengthOut++;
      fNPointsOut[index]=fLengthOut;
      fStatusMaskOut[index]|=kTrackisOK;
  }
  if(fParamMC.size()!=fParamOut.size())
  {
    fParamOut.resize(fParamMC.size());
    fStatusMaskOut.resize(fParamMC.size());
    fNPointsOut.resize(fParamOut.size());
  }
  return 1;
  
}


///
/// \param geom          - pointer to geometry to use
/// \param pdgCode       - pdgCode used in the reconstruction
/// \param layerStart    - starting layer to do tracking
/// \return   -  TODO  status flags to be decides
int fastParticle::reconstructParticleRotate0(fastGeometry  &geom, long pdgCode, uint layerStart){
  const Float_t chi2Cut=10000;
  const float kMaxSnp=0.95;
  //const float kMaxLoss=0.3;
  const float kCovarFactor=2;
  //double covar0[5]={0.1,0.1,0.001,0.001,0.01};
  fLengthInRot=0;
  fPdgCodeRec   =pdgCode;
     float_t mass=0;
  if (pdgCode==0){
    fMassRec=fMassMC;
    mass=fMassMC;
  }else {
      TParticlePDG *particle = TDatabasePDG::Instance()->GetParticle(pdgCode);
      if (particle == nullptr) {
        //::Error("fastParticle::simulateParticle", "Invalid pdgCode %ld", pdgCode);
        return -1;
      }
      mass = particle->Mass();
      fMassRec=mass;
    }
  uint layer1 = TMath::Min(layerStart,uint(fParamMC.size()-1));
  fMaxLayerRec=layer1;
  // for (;layer1>0; layer1--){
  //   if (fLoop[layer1]>0) continue;
  //   if (TMath::Abs(fParamMC[layer1].GetSnp())>kMaxSnp) continue;
  //   if ((fParamMC[layer1].P()/fParamMC[0].P())>kMaxLoss) break;
  //   fMaxLayerRec=layer1;
  // }
  //layer1=1500;
  if (layer1<=3){
    ::Info("fastParticle::reconstructParticleRotate0","short track");
    return -1;
  }
  // AliExternalTrackParam4D param(fParamMC[layer1],mass,1);
  // Double_t dEdx=AliExternalTrackParam::BetheBlochAleph(param.P()/mass);
  // double *covar = (double*)param.GetCovariance();
  // Double_t LArm=getStat(0);
  //  Double_t resFactor=(geom.fLayerResolRPhi[0]*geom.fLayerResolRPhi[0]/0.01)+0.2; // this is hack - we will need to
  // double covar0[5]={0.1,0.1,(1+2.*dEdx/param.P())/(1+LArm),(1+2.*dEdx/param.P())/(1+LArm),resFactor*20*20*(1+2.*dEdx/param.P())/(1+LArm*LArm)};
  // for (int i=0; i<15; i++)covar[i]*=5;
  // covar[0]+=covar0[0]*covar0[0];
  // covar[2]+=covar0[1]*covar0[1];
  // covar[5]+=covar0[2]*covar0[2];
  // covar[9]+=covar0[3]*covar0[3];
  // covar[14]+=covar0[4]*covar0[4];

  Double_t xyzS[3][3];
  Int_t step=layer1/3;
  if (step>10) step=10;

  for (int dLayer=0; dLayer<3  && layer1-step*dLayer>0; dLayer++) {
          Int_t index=layer1-step*dLayer-1;
  //        Int_t index0=layer1-1;
          int layer = fLayerIndex[index];
          //fParamMC[index1-step*dLayer-1].GetXYZ(xyzS[dLayer]);
          xyzS[dLayer][0]=fParamMC[index].GetX();
          xyzS[dLayer][1]=fParamMC[index].GetY();
          xyzS[dLayer][2]=fParamMC[index].GetZ();
          if(fAddPadsmearing)
          {
            xyzS[dLayer][1]+=gRandom->Gaus(0,geom.fLayerResolRPhi[layer]);
            xyzS[dLayer][2]+=gRandom->Gaus(0,geom.fLayerResolZ[layer]);
          }
    }

  Int_t indexst = fLayerIndex[layer1-1];
  Double_t alpha0=fParamMC[layer1-1].GetAlpha();
  float sign0=(fParamMC[layer1-1].GetX()>fParamMC[layer1-2].GetX())? 1.:-1.;
  std::vector<TVector3> points;
  for(int t = layer1-1; t>0;t--){
    double zyx[3];
    zyx[0]=fParamMC[t].GetZ();
    zyx[1]=fParamMC[t].GetY();
    zyx[2]=fParamMC[t].GetX();
    TVector3 ht(zyx[0],zyx[1],zyx[2]);
    points.push_back(ht);
  }
  AliExternalTrackParam * paramSeedI = fastTracker::makeSeed(xyzS[0],xyzS[1],xyzS[2],geom.fLayerResolRPhi[indexst],geom.fLayerResolZ[indexst],geom.fBz);
  AliExternalTrackParam * paramSeed = fastTracker::makeSeedMB(xyzS[0],xyzS[1],xyzS[2],geom.fLayerResolRPhi[indexst],geom.fLayerResolZ[indexst],
                                                              geom.fBz,geom.fLayerX0[indexst],geom.fLayerRho[indexst],fMassMC);
  // AliExternalTrackParam * paramSeedGAr = fastTrackerGAR::Helix_Fit(points,geom.fBz,geom.fLayerResolRPhi[indexst],geom.fLayerResolZ[indexst],1,0); 
  // Double_t c[15] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  // fastTrackerGAR::CalculateCovariance(c,points,paramSeedGAr,geom.fBz,geom.fLayerResolRPhi[indexst],geom.fLayerResolZ[indexst],1,0);
  // if(fParamMC.size()<150){
  //   AliExternalTrackParam paramRot(paramSeedGAr->GetX(),alpha0, paramSeedGAr->GetParameter(),paramSeed->GetCovariance());
  //   AliExternalTrackParam4D param(paramRot,mass,1);
  //   fParamInRot.resize(layer1+1);
  //   fStatusMaskInRot.resize(layer1+1);
  //   fParamInRot[0] = param;
  //   return 1;
  // }                                                      
  AliExternalTrackParam   paramRot;
  // if(paramSeedGAr!=NULL && fParamMC.size()<=75) paramRot = AliExternalTrackParam(paramSeedGAr->GetX(),alpha0, paramSeedGAr->GetParameter(),c);
  // else paramRot = AliExternalTrackParam(paramSeed->GetX(),alpha0, paramSeed->GetParameter(),paramSeed->GetCovariance());
  paramRot = AliExternalTrackParam(paramSeed->GetX(),alpha0, paramSeed->GetParameter(),paramSeed->GetCovariance());
  //AliExternalTrackParam   paramGAr(paramSeedGAr->GetX(),alpha0, paramSeedGAr->GetParameter(),paramSeedGAr->GetCovariance());
  AliExternalTrackParam4D param(paramRot,mass,1);

  if (sign0<0) {
    ((double*)param.GetParameter())[4]*=-1;
    ((double*)param.GetParameter())[3]*=-1;
    ((double*)param.GetParameter())[2]*=-1;
    
    ////Rotate Covariance accordingly RCR^T
    ((double*)param.GetCovariance())[3]*=-1;
    ((double*)param.GetCovariance())[4]*=-1;
    ((double*)param.GetCovariance())[6]*=-1;
    ((double*)param.GetCovariance())[7]*=-1;
    ((double*)param.GetCovariance())[10]*=-1;
    ((double*)param.GetCovariance())[11]*=-1;
  }


  //int version=4;
  //Double_t dEdx=AliExternalTrackParam::BetheBlochAleph(param.P()/mass);
  // (*fgStreamer)<<"seedDump"<<   // seeding not ideal in case significant energy loss
  //   "version="<<version<<
  //   "gid="<<gid<<
  //   "sign0="<<sign0<<
  //   "fMassMC="<<fMassMC<<
  //   "dEdx="<<dEdx<<
  //   "step="<<step<<
  //   "seed.="<<&param<<
  //   "input.="<<&fParamMC[layer1-1]<<
  //   "input1.="<<&fParamMC[layer1-step-1]<<
  //   "input2.="<<&fParamMC[layer1-2*step-1]<<
  //   "paramSeed.="<<paramSeed<<
  //   "paramSeedI.="<<paramSeedI<<
  //   "paramRot.="<<&paramRot<<
  //   "\n";
  delete paramSeed;
  delete paramSeedI;
  ///scale covaraince to remove double counting
  for (int i=0; i<15; i++) ((double*)param.GetCovariance())[i]*=kCovarFactor;
  //
  Double_t sinT=TMath::Sin(param.GetAlpha());
  Double_t cosT=TMath::Cos(param.GetAlpha());
  //
  //float length=0, time=0;
  //float radius = sqrt(param.GetX()*param.GetX()+param.GetY()*param.GetY());
  fParamInRot.resize(layer1+1);
    fStatusMaskInRot.resize(layer1+1);
  fParamInRot[layer1]=param;
  double xyz[3];
  int status=0;
  const double *par = param.GetParameter();
  for (int layer=layer1-1; layer>=0; layer--){   // dont propagate to vertex , will be done later ...
  //  double resol=0;
    Int_t index = fLayerIndex[layer];
    AliExternalTrackParam & p = fParamMC[layer];
    p.GetXYZ(xyz);
  //  double alpha=TMath::ATan2(xyz[1],xyz[0]);
  //  double radius = TMath::Sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
    Double_t localX, localY;
    localX=cosT*xyz[0]+sinT*xyz[1];
    localY=-sinT*xyz[0]+cosT*xyz[1];
    fStatusMaskInRot[layer]=0;
    status = param.PropagateTo(localX,geom.fBz,1);
    if (status) {
      fStatusMaskInRot[layer]|=kTrackPropagate;
    }else{
      //::Error("fastParticle::reconstructParticle:", "Propagation failed");
      param.PropagateTo(localX,geom.fBz,1);
      break;
    }
    //

    float xrho  =geom.fLayerRho[index];
    float xx0  =geom.fLayerX0[index];
    float deltaY = gRandom->Gaus(0, geom.fLayerResolRPhi[index]);
    float deltaZ = gRandom->Gaus(0, geom.fLayerResolZ[index]);
    double pos[2]={localY,xyz[2]};
    if(fAddPadsmearing)
      {
        pos[0]+=deltaY;
        pos[1]+=deltaZ;
      }
    double cov[3]={geom.fLayerResolRPhi[index]*geom.fLayerResolRPhi[index],0, geom.fLayerResolZ[index]*geom.fLayerResolZ[index]};
    fParamInRot[layer]=param;    ///
    float chi2 =  param.GetPredictedChi2(pos, cov);
    fChi2[layer]=chi2;
    if (chi2<chi2Cut) {
      fStatusMaskInRot[layer] |= kTrackChi2;
    }else {
      //::Error("fastParticle::reconstructParticle:", "Chi2 -  failed");
      break;
    }
    //
    if (TMath::Abs(param.GetSnp())<kMaxSnp &&cov[0]>0)  {
      status = param.Update(pos, cov);
      if (status) {
        fStatusMaskInRot[layer] |= kTrackUpdate;
      } else {
        //::Error("fastParticle::reconstructParticle:", "Update failed");
        break;
      }
    }
    float tanPhi2 = pow(fParamInRot[layer+1].GetParameter()[2],2);
    tanPhi2/=(1-tanPhi2);
    float deltaR = abs(fParamInRot[layer+1].GetX()-fParamInRot[layer].GetX());
    float deltax = fParamInRot[layer+1].GetX()-fParamInRot[layer].GetX();
    //float deltay = fParamInRot[layer+1].GetY()-fParamInRot[layer].GetY();
    //float deltaz = fParamInRot[layer+1].GetZ()-fParamInRot[layer].GetZ();
    float deltaxMC = fParamMC[layer+1].GetX()-fParamMC[layer].GetX();
    float deltayMC = fParamMC[layer+1].GetY()-fParamMC[layer].GetY();
    float deltazMC = fParamMC[layer+1].GetZ()-fParamMC[layer].GetZ();
    if(deltax!=0){
      //float crossLength=TMath::Sqrt(deltax*deltax+deltay*deltay+deltaz*deltaz);                /// geometrical path assuming crossing cylinder  = to be racalculated
      float crossLength=deltaR*TMath::Sqrt(1.+tanPhi2+par[3]*par[3]);
      //status = param.AliExternalTrackParam::CorrectForMeanMaterial(crossLength*xx0,crossLength*xrho,mass);
      for (Int_t ic=0;ic<5; ic++) {
        status *= param.CorrectForMeanMaterial(crossLength * xx0/5., crossLength * xrho/5., mass, 0.01);
      }
      //status = param.CorrectForMeanMaterialT4(crossLength*xx0,crossLength*xrho,mass);
      if (status) {
        fStatusMaskInRot[layer] |= kTrackCorrectForMaterial;
      } else {
        //::Error("fastParticle::reconstructParticle:", "Correct for material failed");
        break;
      }
      fLengthInRot+= TMath::Sqrt(deltaxMC*deltaxMC+deltayMC*deltayMC+deltazMC*deltazMC);
    }
      
  }
  return 1;
}

///
/// \param valueType      0 - <pt>,  1 - <1/pt> - 2 -<dEdxExp>, 3 - 1/dEdx, 4 - dEdx/p**2
/// \return
Float_t fastParticle::getMean(Int_t valueType,Float_t beginF, Float_t endF, Int_t powerType){
  Float_t valueMean=0;
  Float_t nPoints=0;
  if (beginF<0 || beginF>=1 || beginF>endF) return 0;
  if (endF<0 || endF>1) return 0;
  Int_t sizeRec = fMaxLayerRec;
  for (UInt_t i=beginF*sizeRec; i<endF*sizeRec;i++){
    Float_t value=0;
    if (valueType==0) value=fParamMC[i].Pt();
    if (valueType==1) value=1/fParamMC[i].Pt();
    if (valueType==2) value=TMath::Max(AliExternalTrackParam4D::BetheBlochAleph(fParamMC[i].GetP()/fMassMC),0.);
    if (valueType==3) value=TMath::Max(1/AliExternalTrackParam4D::BetheBlochAleph(fParamMC[i].GetP()/fMassMC),0.);
    if (valueType==4) value=TMath::Max(AliExternalTrackParam4D::BetheBlochAleph(fParamMC[i].GetP()/fMassMC),.0)/(fParamMC[i].P()*fParamMC[i].P());
    //
    if (powerType==2) value*=value;
    if (powerType==-1) value=1./value;
    valueMean+=value;
    nPoints++;
  }
  if (nPoints==0) return 0;
  return  valueMean/nPoints;
}
///
/// \param valueType      0 - LArmMC  1 - LArmIn 2-MC length, 3 - MC legnth for the last point
/// \param averageType    dummy
/// \return
Float_t fastParticle::getStat(Int_t valueType){
  //Float_t valueMean=0;
  //Float_t nPoints=0;
  Float_t rMin = -1;
  Float_t rMax = -1;
  Double_t vecMin[3],vecMax[3];
  if (valueType==2) return fParamMC.size();
  if (valueType==3) {return fParamMC[fParamMC.size()-1].fLength;}
  if (valueType==0) {
    for (UInt_t i = 0; i < fParamMC.size(); i++) {
      Float_t x=  fParamMC[i].GetX();
      if (rMin>x || rMin<0)  rMin=x;
      if (rMax<x || rMax==0) rMax=x;
    }
    return rMax-rMin;
  }
  if (valueType==1) {
    for (UInt_t i = 0; i < fParamIn.size(); i++) {
      Float_t x=  fParamIn[i].GetX();
      if (rMin>x || rMin<0)  rMin=x;
      if (rMax<x || rMax==0) rMax=x;
    }
    return rMax-rMin;
  }
  if (valueType==2) {
    for (UInt_t i = 0; i < fParamIn.size(); i++) {
      Float_t x=  fParamIn[i].GetX();
      if (rMin>x || rMin<0)  {
        rMin=x; fParamIn[i].GetXYZ(vecMin);
      }
      if (rMax<x || rMax==0) {
        rMax=x; fParamIn[i].GetXYZ(vecMax);
      }
    }
    return TMath::Sqrt((vecMin[0]-vecMax[0])*(vecMin[0]-vecMax[0])+ (vecMin[1]-vecMax[1])*(vecMin[1]-vecMax[1]));
  }

  return 0;
}




