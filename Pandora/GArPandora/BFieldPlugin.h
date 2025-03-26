#ifndef BFIELD_PLUGIN_H
#define BFIELD_PLUGIN_H 1

#include "Objects/CartesianVector.h"
#include "Plugins/BFieldPlugin.h"
#include "nug4/MagneticFieldServices/MagneticFieldService.h"

namespace gar {
  namespace gar_pandora {

    class BFieldPlugin : public pandora::BFieldPlugin {
    public:
      BFieldPlugin();
      float GetBField(const pandora::CartesianVector& positionVector) const;

    private:
      pandora::StatusCode Initialize();
      pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    };
  }
}

#endif // #ifndef BFIELD_PLUGIN_H
