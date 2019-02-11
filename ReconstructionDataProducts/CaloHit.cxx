//
//  Hit.cxx
//
//  Created by Brian Rebel on 10/6/16.
//
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "CoreUtils/ServiceUtil.h"

#include "ReconstructionDataProducts/CaloHit.h"

namespace gar {
    namespace rec {

        //--------------------------------------------------------------------------
        CaloHit::CaloHit()
        : fEnergy(0.),
        fTime(0.),
        fCellID(0.),
        fPositionVector(0., 0., 0.)
        {
            fPosition[0] = 0.;
            fPosition[1] = 0.;
            fPosition[2] = 0.;

            return;
        }

        //--------------------------------------------------------------------------
        CaloHit::CaloHit(float energy, float time, float *pos, long long int cellID)
        : fEnergy  (energy  )
        , fTime   (time   )
        , fCellID     (cellID  )
        {

            fPosition[0] = pos[0];
            fPosition[1] = pos[1];
            fPosition[2] = pos[2];

            fPositionVector = CLHEP::Hep3Vector(pos[0], pos[1], pos[2]);

            return;
        }

        //--------------------------------------------------------------------------
        std::ostream& operator<< (std::ostream& o, gar::rec::CaloHit const& h)
        {

            o << "CaloHit "
            << "\n\tposition = ("
            << h.Position()[0]
            << ", "
            << h.Position()[1]
            << ", "
            << h.Position()[2]
            << ")"
            << "\n\tenergy = "
            << h.Energy()
            << "\n\t time: "
            << h.Time()
            << " cellID: "
            << h.CellID();

            return o;
        }

        //--------------------------------------------------------------------------
        const unsigned int CaloHit::GetLayer() const
        {
            gar::geo::GeometryCore const* fGeo = gar::providerFrom<geo::Geometry>();
            return fGeo->getIDbyCellID(this->CellID(), "layer");
            delete fGeo;
        }

    } // rec
} // gar
