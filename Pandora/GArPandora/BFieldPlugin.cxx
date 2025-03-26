#include "BFieldPlugin.h"

#include "Helpers/XmlHelper.h"

#include "CoreUtils/ServiceUtil.h"

namespace gar {
  namespace gar_pandora {

    BFieldPlugin::BFieldPlugin()
    {
      /* no op */
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    float BFieldPlugin::GetBField(const pandora::CartesianVector& positionVector) const
    {
      auto const* fieldService = gar::providerFrom<mag::MagneticFieldService>();
      G4ThreeVector PosVec(positionVector.GetX(), positionVector.GetY(), positionVector.GetZ());
      G4ThreeVector magfield = fieldService->FieldAtPoint(PosVec);
      return magfield[0];
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    pandora::StatusCode BFieldPlugin::Initialize()
    {
      /* nop */
      return pandora::STATUS_CODE_SUCCESS;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    pandora::StatusCode BFieldPlugin::ReadSettings(const pandora::TiXmlHandle /*xmlHandle*/)
    {
      /* nop */
      return pandora::STATUS_CODE_SUCCESS;
    }
  }
}
