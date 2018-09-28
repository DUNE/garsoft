////////////////////////////////////////////////////////////////////////
/// \file  ECALUtils.cxx
/// \brief Interface to algorithm class ECAL specific
///
/// \author  eldwan.brianne@desy.de
////////////////////////////////////////////////////////////////////////

#define SHIFT_CALOID 61
#define SHIFT_NEG_I 60
#define SHIFT_I 34
#define SHIFT_NEG_J 33
#define SHIFT_J 7
#define SHIFT_LAYER 0

#include "ReadoutSimulation/ECALUtils.h"
#include <cmath>

#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoBBox.h"
#include "TGeoTube.h"
#include "TGeoMatrix.h"

namespace gar {
  namespace rosim {

    //----------------------------------------------------------------------
    ECALUtils::ECALUtils()
    : fNeffPx(0),
    fRatio(0)
    {
    }

    //----------------------------------------------------------------------
    ECALUtils::ECALUtils(double NeffPx, double ratio)
    : fNeffPx(NeffPx),
    fRatio(ratio)
    {
    }

    //----------------------------------------------------------------------
    ECALUtils::~ECALUtils()
    {
    }

    //----------------------------------------------------------------------------
    double ECALUtils::Saturate(const double unsat_px)
    {
      /*  protect against negative/zero input*/
      if ( unsat_px <= 0 ) return unsat_px;

      /* calculate saturated signal using simple exponential formula */
      double saturatedSignal = fNeffPx * ( 1. - std::exp( - unsat_px / fNeffPx ) );

      return saturatedSignal;
    }

    //----------------------------------------------------------------------------
    double ECALUtils::DeSaturate(const double sat_px)
    {
      /*  protect against negative/zero input*/
      if ( sat_px <= 0 ) return sat_px;

      if(sat_px < fRatio * fNeffPx)
      {
        /* desaturate using inverse function */
        float unSaturatedSignal = - fNeffPx * std::log(1 - sat_px / fNeffPx);
        return unSaturatedSignal;
      }
      else
      {
        /* desaturate using linear continuation function for hits above linearisation threshold */
        float unSaturatedSignal = 1/( 1 - fRatio ) * (sat_px - fRatio * fNeffPx) - fNeffPx * std::log( 1 - fRatio );
        return unSaturatedSignal;
      }
    }

    //----------------------------------------------------------------------------
    int ECALUtils::PositionToBin(double position, double cellsize, double offset)
    {
      return int(std::floor((position + 0.5 * cellsize - offset) / cellsize));
    }

    //----------------------------------------------------------------------------
    double ECALUtils::BinToPosition(int bin, double cellsize, double offset)
    {
      return double(bin * cellsize + offset);
    }

    //----------------------------------------------------------------------------
    unsigned long long int ECALUtils::MakeCellID(int id, int binI, int binJ, unsigned int layer)
    {
      unsigned long int kEndcap = 0;
      unsigned long int kNegI = 0;
      unsigned long int kNegJ = 0;

      //id is related to Endcap or Barrel
      if(id == 3) kEndcap = 1;

      //Encode the cellID in 64 bits
      //1) take absolute of binI and binJ
      //2) make flag if binI / binJ negative
      //3) encode on 64 bits
      //  - 26 bits for binI + flag negative
      //  - 26 bits for binJ + flag negative
      //  - 7 bits for layer (up to 128 layers)
      //  - 1 bit for Endcap or Barrel (1 Endcap / 0 Barrel)

      if(binI < 0)
      kNegI = 1;

      if(binJ < 0)
      kNegJ = 1;

      unsigned long int absbinI = std::abs(binI);
      unsigned long int absbinJ = std::abs(binJ);

      unsigned long long int cellID = (
        ( kEndcap << SHIFT_CALOID ) |
        ( kNegI << SHIFT_NEG_I ) |
        ( absbinI << SHIFT_I ) |
        ( kNegJ << SHIFT_NEG_J ) |
        ( absbinJ << SHIFT_J ) |
        ( layer << SHIFT_LAYER )
      );

      return cellID;
    }

    //----------------------------------------------------------------------------
    void ECALUtils::GetLayerThickness(const gar::geo::GeometryCore *geo, double &layer_thickness)
    {
      TGeoManager *gGeoManager = geo->ROOTGeoManager();
      TGeoVolume *vol = gGeoManager->FindVolumeFast("IECLayer_vol");

      if(vol)
      layer_thickness = ((TGeoBBox*)vol->GetShape())->GetDZ() * 2;

      return;
    }

    //----------------------------------------------------------------------------
    void ECALUtils::GetRadius(const gar::geo::GeometryCore *geo, const std::string name, double &R)
    {
      TGeoManager *gGeoManager = geo->ROOTGeoManager();
      TGeoVolume *vol = gGeoManager->FindVolumeFast(name.c_str());

      if(vol)
      R = ((TGeoTube*)vol->GetShape())->GetRmin();

      return;
    }

    //----------------------------------------------------------------------------
    void ECALUtils::GetEndcapStartPosition(const gar::geo::GeometryCore *geo, double &pos)
    {
      TGeoManager *gGeoManager = geo->ROOTGeoManager();
      TGeoVolume *det_vol = gGeoManager->FindVolumeFast("volNDHPgTPC");
      TGeoVolume *endcap_vol = gGeoManager->FindVolumeFast("InnerEndcapECal_vol");

      double thickness = 0.;

      if(endcap_vol)
      thickness = ((TGeoTube*)endcap_vol->GetShape())->GetDZ();

      if(det_vol)
      {
        TGeoNode *endcap_node = det_vol->FindNode("InnerEndcapECal_vol_0");
        TGeoMatrix *mat = endcap_node->GetMatrix();
        const double *origin = mat->GetTranslation();
        pos = origin[0] - thickness;
      }

      return;
    }

  }
} // gar
