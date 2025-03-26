#ifndef TRIGGERDATA_CXX
#define TRIGGERDATA_CXX

#include "RawDataProducts/Trigger.h"

namespace gar {
  namespace raw {

    //****************************************************
    bool Trigger::Triggered(const unsigned char bit) const
    {

      if (bit > 32) { throw std::invalid_argument("\n\nCannot access bit higher than 32!\n"); }

      return ((fTriggerBits >> bit) & 0x1);
    }

  }
} // gar
#endif
