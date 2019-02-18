/// \file  ColorDrawingOptions.h
/// \brief The color scales used by the event display
/// \author messier@indiana.edu
/// \version $Id: ColorScaleTable.h,v 1.1.1.1 2010/11/10 19:44:54 p-novaart Exp $
#ifndef EVD_COLORDRAWINGOPTIONS_H
#define EVD_COLORDRAWINGOPTIONS_H
#ifndef __CINT__
#include "nutools/EventDisplayBase/ColorScale.h"
#include "Geometry/Geometry.h"

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "nutools/EventDisplayBase/Reconfigurable.h"

namespace gar {
namespace evd {
  class ColorDrawingOptions :public evdb::Reconfigurable
  {
  public:

    explicit ColorDrawingOptions(fhicl::ParameterSet   const& pset,
				 art::ActivityRegistry      & reg);
    ~ColorDrawingOptions();

    void reconfigure(fhicl::ParameterSet const& pset);

    evdb::ColorScale const& RawQ() const;
    evdb::ColorScale const& CalQ() const;
    evdb::ColorScale const& RawT() const;
    evdb::ColorScale const& CalT() const;

    int    fColorOrGray; ///< 0 = color, 1 = gray
    int    fRawDiv;      ///< number of divisions in raw
    int    fRecoDiv;     ///< number of divisions in raw
    double fRawQLow;     ///< low  edge of ADC values for drawing raw digits
    double fRawQHigh;    ///< high edge of ADC values for drawing raw digits
    double fRecoQLow;    ///< low  edge of ADC values for drawing raw digits
    double fRecoQHigh;   ///< high edge of ADC values for drawing raw digits

  private:

    evdb::ColorScale fColorScaleRaw;
    evdb::ColorScale fGrayScaleRaw;
    evdb::ColorScale fColorScaleReco;
    evdb::ColorScale fGrayScaleReco;
  };
}
}
#endif // __CINT__
DECLARE_ART_SERVICE(gar::evd::ColorDrawingOptions, LEGACY)
#endif
////////////////////////////////////////////////////////////////////////
