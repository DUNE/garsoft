////////////////////////////////////////////////////////////////////////
/// \file ColorDrawingOptions_service.cc
///
/// \version $Id: ColorDrawingOptions_service.cc,v 1.1 2010/11/11 18:11:22 p-novaart Exp $
/// \author  brebel@fnal.gov

// Framework includes

/// GArSoft includes
#include "EventDisplay/ColorDrawingOptions.h"
#include "EventDisplay/RawDrawingOptions.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Event.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace gar {
namespace evd{
  
  //......................................................................
  ColorDrawingOptions::ColorDrawingOptions(fhicl::ParameterSet   const& pset,
                                           art::ActivityRegistry      & /* reg */)
    
    : evdb::Reconfigurable{pset}
    , fColorOrGray  (pset.get< int    >("ColorOrGrayScale"))
    , fRawDiv       (pset.get< int    >("RawDiv")          )
    , fRecoDiv      (pset.get< int    >("RecoDiv")         )
    , fRawQLow      (pset.get< double >("RawQLow")         )
    , fRawQHigh     (pset.get< double >("RawQHigh")        )
    , fRecoQLow     (pset.get< double >("RecoQLow")        )
    , fRecoQHigh    (pset.get< double >("RecoQHigh")       )
    , fColorScaleRaw(evdb::ColorScale(fRawQLow, fRawQHigh,
				      evdb::kBlueToRedII, evdb::kLinear,
				      fRawDiv,
				      285.0, 135.0, // angle in the color wheel
				      0.65, 0.25))  // intensity from light to dark,
      // starting with low color wheel value
    , fGrayScaleRaw(evdb::ColorScale(fRawQLow, fRawQHigh,
				     evdb::kLinGray, evdb::kLinear,
				     fRawDiv,
				     270.0, 0.0, // angle in the color wheel
				     0.5, 0.5))  // intensity from light to dark,
      // starting with low color wheel value
    , fColorScaleReco(evdb::ColorScale(fRecoQLow, fRecoQHigh,
				       evdb::kBlueToRedII, evdb::kLinear,
				       fRecoDiv,
				       285.0, 135.0,
				       0.65, 0.25))
    , fGrayScaleReco(evdb::ColorScale(fRecoQLow, fRecoQHigh,
				      evdb::kLinGray, evdb::kLinear,
				      fRecoDiv,
				      270.0, 0.0,
				      0.5, 0.5))
    {
      return;
    }
  
  //......................................................................
  ColorDrawingOptions::~ColorDrawingOptions()
  {
  }

  //......................................................................
  void ColorDrawingOptions::reconfigure(fhicl::ParameterSet const& pset)
  {
    fColorOrGray = pset.get< int    >("ColorOrGrayScale");
    fRawDiv      = pset.get< int    >("RawDiv");
    fRecoDiv     = pset.get< int    >("RecoDiv");
    fRawQLow     = pset.get< double >("RawQLow");
    fRawQHigh    = pset.get< double >("RawQHigh");
    fRecoQLow    = pset.get< double >("RecoQLow");
    fRecoQHigh   = pset.get< double >("RecoQHigh");

    fColorScaleRaw = evdb::ColorScale(fRawQLow, fRawQHigh,
                                      evdb::kBlueToRedII, evdb::kLinear,
                                      fRawDiv,
                                      285.0, 135.0, // angle in the color wheel
                                      0.65, 0.25);  // intensity from light to dark,
                                                    // starting with low color wheel value
    
    fGrayScaleRaw = evdb::ColorScale(fRawQLow, fRawQHigh,
                                     evdb::kLinGray, evdb::kLinear,
                                     fRawDiv,
                                     270.0, 0.0, // angle in the color wheel
                                     0.5, 0.5);  // intensity from light to dark,
                                                 // starting with low color wheel value
    
    fColorScaleReco = evdb::ColorScale(fRecoQLow, fRecoQHigh,
                                       evdb::kBlueToRedII, evdb::kLinear,
                                       fRecoDiv,
                                       285.0, 135.0,
                                       0.65, 0.25);
    
    fGrayScaleReco = evdb::ColorScale(fRecoQLow, fRecoQHigh,
                                      evdb::kLinGray, evdb::kLinear,
                                      fRecoDiv,
                                      270.0, 0.0,
                                      0.5, 0.5);
    
    return;
  }

  //......................................................................
  const evdb::ColorScale& ColorDrawingOptions::RawQ() const
  {
    if(fColorOrGray > 0) return fGrayScaleRaw;
    
    return fColorScaleRaw;
  }

  //......................................................................
  const evdb::ColorScale& ColorDrawingOptions::CalQ() const
  {
    if(fColorOrGray > 0) return fGrayScaleReco;
    
    return fColorScaleReco;
  }

  //......................................................................
  const evdb::ColorScale& ColorDrawingOptions::RawT() const
  {
    if(fColorOrGray > 0) return fGrayScaleRaw;
    
    return fColorScaleRaw;
  }

  //......................................................................
  const evdb::ColorScale& ColorDrawingOptions::CalT() const
  {
    if(fColorOrGray > 0) return fGrayScaleReco;
    
    return fColorScaleReco;
  }
}// namespace

namespace evd {

  DEFINE_ART_SERVICE(ColorDrawingOptions)

} // namespace evd
}
////////////////////////////////////////////////////////////////////////
