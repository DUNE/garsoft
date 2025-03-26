//
//  StripSplitterAlg.cxx
//
//  Created by Eldwan Brianne on 18.04.2019
//

#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "CoreUtils/ServiceUtil.h"
#include "Geometry/GeometryGAr.h"
#include "Geometry/LocalTransformation.h"
#include "RecoAlg/StripSplitterAlg.h"

#include "canvas/Persistency/Common/Ptr.h"
#include <algorithm>
#include <array>
#include <functional>
#include <limits>

namespace gar {
  namespace rec {
    namespace alg {

      //----------------------------------------------------------------------------
      StripSplitterAlg::StripSplitterAlg(fhicl::ParameterSet const& pset)
      {
        ClearLists();

        fGeo = gar::providerFrom<geo::GeometryGAr>();
        fInnerSymmetry = 0;
        fStripWidth = -1;  //cm
        fStripLength = -1; //cm
        fnVirtual = 1;

        this->reconfigure(pset);

        return;
      }

      //----------------------------------------------------------------------------
      StripSplitterAlg::~StripSplitterAlg()
      {
        if (fFieldDecoder) delete fFieldDecoder;
        return;
      }

      //----------------------------------------------------------------------------
      void StripSplitterAlg::reconfigure(fhicl::ParameterSet const& pset)
      {
        fSSAAlgName = pset.get<std::string>("SSAAlgName");
        fDet = pset.get<std::string>("DetectorSystem"); //ECAL or MuID
        fSaveStripEnds = pset.get<bool>("SaveStripEnds", false);

        if (fDet == "ECAL") {
          fEncoding = fGeo->GetECALCellIDEncoding();
          fFieldDecoder = new gar::geo::BitFieldCoder(fEncoding);
          fInnerSymmetry = fGeo->GetECALInnerSymmetry();
        }
        else if (fDet == "MuID") {
          fEncoding = fGeo->GetMuIDCellIDEncoding();
          fFieldDecoder = new gar::geo::BitFieldCoder(fEncoding);
          fInnerSymmetry = fGeo->GetMuIDInnerSymmetry();
        }
        else {
          MF_LOG_ERROR("StripSplitterAlg::reconfigure")
            << "Detector system parameter incorrectly set. Should be ECAL or MuID -- Exiting!";
          throw cet::exception("StripSplitterAlg::reconfigure");
        }

        return;
      }

      //----------------------------------------------------------------------------
      void StripSplitterAlg::ClearLists()
      {
        m_CaloHitVecOdd.clear();
        m_CaloHitVecEven.clear();

        fStripEndsHits.clear();
        fIntersectionHits.clear();

        unSplitStripHits.clear();
        splitStripHits.clear();

        return;
      }

      //----------------------------------------------------------------------------
      bool StripSplitterAlg::SortByLayer(const gar::rec::CaloHit* rha, const gar::rec::CaloHit* rhb)
      {
        return rha->Layer() < rhb->Layer();
      }

      //----------------------------------------------------------------------------
      void StripSplitterAlg::PrepareAlgo(const std::vector<art::Ptr<gar::rec::CaloHit>>& hitVector)
      {
        MF_LOG_DEBUG("StripSplitterAlg") << "StripSplitterAlg::PrepareAlgo()";

        //Clear the lists
        ClearLists();

        //Loop over all hits
        for (std::vector<art::Ptr<gar::rec::CaloHit>>::const_iterator iter = hitVector.begin(),
                                                                      iterEnd = hitVector.end();
             iter != iterEnd;
             ++iter) {
          art::Ptr<gar::rec::CaloHit> hitPtr = *iter;
          const gar::rec::CaloHit* hit = hitPtr.get();

          //check the layer
          unsigned int layer = hit->Layer();
          if (layer % 2 == 0)
            m_CaloHitVecEven.emplace_back(hit);
          else {
            m_CaloHitVecOdd.emplace_back(hit);
          }
        }

        //Sort the hits by layer
        std::sort(m_CaloHitVecEven.begin(), m_CaloHitVecEven.end(), StripSplitterAlg::SortByLayer);
        std::sort(m_CaloHitVecOdd.begin(), m_CaloHitVecOdd.end(), StripSplitterAlg::SortByLayer);

        return;
      }

      //----------------------------------------------------------------------------
      void StripSplitterAlg::DoStripSplitting()
      {
        MF_LOG_DEBUG("StripSplitterAlg") << "StripSplitterAlg::DoStripSplitting()";

        //Collection to split
        std::vector<const gar::rec::CaloHit*>* toSplit;
        int orientation;
        //Orientation
        //even layers are segmented in Y (local) -> Transverse
        //odd layers are segmented in X (local) -> Longitudinal

        // loop over strip collections in even and odd layers (assumed to have perpendicular orientations)
        for (int icol = 0; icol < 2; icol++) {
          MF_LOG_DEBUG("StripSplitterAlg") << "Start Loop " << icol;

          switch (icol) {
          case 0:
            orientation = TRANSVERSE; //even layers
            toSplit = &m_CaloHitVecEven;
            break;
          case 1:
            orientation = LONGITUDINAL; //odd layers
            toSplit = &m_CaloHitVecOdd;
            break;
          default:
            MF_LOG_ERROR("StripSplitterAlg::DoStripSplitting") << "Crazy stuff!";
            throw cet::exception("StripSplitterAlg::DoStripSplitting");
          }

          // loop over hits of this type (long/trans) collection
          for (uint i = 0; i < toSplit->size(); i++) {
            const gar::rec::CaloHit* hit = toSplit->at(i);
            //1 == Barrel, 2 == Endcap, 4 == MuID
            unsigned int det_id = fFieldDecoder->get(hit->CellID(), "system");

            if (det_id != 1 && det_id != 2 && det_id != 4) {
              MF_LOG_ERROR("StripSplitterAlg::DoStripSplitting")
                << " Check det it " << det_id << " Problem with Hit " << i << " pointing at " << hit
                << " at position ( " << hit->Position()[0] << " cm , " << hit->Position()[1]
                << " cm , " << hit->Position()[2] << " cm )"
                << " orientation " << (orientation == LONGITUDINAL ? "LONGITUDINAL" : "TRANSVERSE");
              throw cet::exception("StripSplitterAlg::DoStripSplitting");
            }

            MF_LOG_DEBUG("StripSplitterAlg::DoStripSplitting")
              << "recohit " << hit << " with cellID " << hit->CellID() << " has energy "
              << hit->Energy() * CLHEP::MeV / CLHEP::GeV << " pos (" << hit->Position()[0] << ", "
              << hit->Position()[1] << ", " << hit->Position()[2] << ")";

            // split the hits
            std::vector<const gar::rec::CaloHit*> virtualhits;
            // As of Sep 2022, only barrel of yoke is instrumented so MuID is all barrel
            bool isBarrel = (det_id == 1 || det_id == 4) ? true : false;
            getVirtualHits(hit, orientation, isBarrel, virtualhits);

            // add (new) hits to collections
            if (virtualhits.size() == 0) {

              MF_LOG_DEBUG("StripSplitterAlg")
                << " Adding unsplit hit1 " << i << " pointing at " << hit << " orientation "
                << (orientation == LONGITUDINAL ? "LONGITUDINAL" : "TRANSVERSE");

              // not split, add original hit
              unSplitStripHits.emplace_back(hit);
            }
            else {
              // split was split, add the virtual hits
              for (uint hh = 0; hh < virtualhits.size(); hh++) {

                MF_LOG_DEBUG("StripSplitterAlg") << "adding virtual hit " << hh;

                splitStripHits.emplace_back(virtualhits.at(hh));
              }
            }
          } // loop over hits
        }   // long/trans loop

        return;
      }

      //----------------------------------------------------------------------------
      void StripSplitterAlg::getVirtualHits(const gar::rec::CaloHit* hit,
                                            int, // orientation,
                                            bool isBarrel,
                                            std::vector<const gar::rec::CaloHit*>& virtualhits)
      {
        MF_LOG_DEBUG("StripSplitterAlg") << "StripSplitterAlg::getVirtualHits()";

        // this splits the strip into zero or more hits along its length
        // by looking at nearby hits with different orientation (trans/long or square)

        int detid = fFieldDecoder->get(hit->CellID(), "system");
        int layer = hit->Layer();
        int module = fFieldDecoder->get(hit->CellID(), "module");
        int stave = fFieldDecoder->get(hit->CellID(), "stave");

        const std::array<double, 3> pt = {
          hit->Position()[0], hit->Position()[1], hit->Position()[2]};
        fStripWidth = fGeo->getStripWidth(pt);
        fStripLength = fGeo->getStripLength(pt, hit->CellID());
        fnVirtual = 0;

        if (fStripWidth != 0) fnVirtual = int(fStripLength / fStripWidth);

        MF_LOG_DEBUG("StripSplitterAlg") << " StripSplitterAlg::getVirtualHits()"
                                         << " Strip splitted in " << fnVirtual << " virtual cells";

        // get the ends of this strip
        std::pair<TVector3, TVector3> stripEnds = fGeo->GetStripEnds(pt, hit->CellID());
        TVector3 stripDir = stripEnds.first - stripEnds.second;

        if (fSaveStripEnds) {
          float ppp[3];
          ppp[0] = stripEnds.first.X();
          ppp[1] = stripEnds.first.Y();
          ppp[2] = stripEnds.first.Z();

          const gar::rec::CaloHit* StripEndHit1 = new gar::rec::CaloHit(0.01, 0., ppp, 0., 0);
          fStripEndsHits.emplace_back(StripEndHit1);

          ppp[0] = stripEnds.second.X();
          ppp[1] = stripEnds.second.Y();
          ppp[2] = stripEnds.second.Z();

          const gar::rec::CaloHit* StripEndHit2 = new gar::rec::CaloHit(0.01, 0., ppp, 0., 0);
          fStripEndsHits.emplace_back(StripEndHit2);
        }

        // decide which collections to use to split the strip
        int splitterOrientation;
        std::vector<const gar::rec::CaloHit*>* splitterVec;

        if (layer % 2 == 0) {
          splitterVec = &m_CaloHitVecOdd;
          splitterOrientation = LONGITUDINAL;
        }
        else {
          splitterVec = &m_CaloHitVecEven;
          splitterOrientation = TRANSVERSE;
        }

        std::map<int, float> virtEnergy;
        int nSplitters(0);

        // loop over splitter cols, find nearby hits
        // strips, cells
        for (int jj = 0; jj < 2; jj++) {
          std::vector<const gar::rec::CaloHit*>* splitter = splitterVec;

          for (uint i = 0; i < splitter->size(); i++) {
            const gar::rec::CaloHit* hit2 = splitter->at(i);

            int detid2 = fFieldDecoder->get(hit2->CellID(), "system");
            int layer2 = hit2->Layer();
            int module2 = fFieldDecoder->get(hit2->CellID(), "module");
            int stave2 = fFieldDecoder->get(hit2->CellID(), "stave");

            int ddetid = std::abs(detid2 - detid);
            int dlayer = std::abs(layer2 - layer);
            int dstave = std::abs(stave2 - stave);
            int dmodule = std::abs(module2 - module);

            TVector3 point2(hit2->Position()[0], hit2->Position()[1], hit2->Position()[2]);

            //if they are in different part of the det, ignore
            if (ddetid != 0) continue;

            //TODO CHECK!
            // are the two hits close enough to look at further?
            // if hits in same module and same stave, require that only one layer difference
            if (dmodule == 0 && dstave == 0 && dlayer > 1) continue;

            if (isBarrel) {
              dstave = std::min(dstave, fInnerSymmetry - dstave);
              // Leo Asks: Maybe we should just split inside a single stave & module
              if (dstave == 0 && dmodule > 1) continue; // allow same stave and +- 1 module
              if (dmodule == 0 && dstave > 1) continue; // or same module +- 1 stave
              // Leo Asks: Maybe we should just require dlayer == 1 always?
              if (dstave == 0 && dlayer > 1) continue; // if in same stave, require dlayer==1
            }
            else {
              //TODO CHECK!
              //module == 0 or 6 for endcap in hexagonal ECALs.  But not in dodecagonal ones.
              //stave is from 1 to 4
              dstave = std::min(dstave, 4 - dstave);
              if (dmodule != 0) continue; // different endcap
              // Leo Asks: Maybe we should just split inside a single stave & module
              if (dstave > 1) continue; // more than 1 stave (=quarter endcap) apart
              if (dlayer > 1) continue; // more than 1 layer apart
            }

            MF_LOG_DEBUG("StripSplitterAlg::getVirtualHits()")
              << " Getting hit2 " << i << " pointing at " << hit2 << " orientation "
              << (splitterOrientation == LONGITUDINAL ? "LONGITUDINAL" : "TRANSVERSE") << " dtave "
              << dstave << " dmodule " << dmodule << " dlayer " << dlayer;

            // simple distance check for remaining hit pairs
            float dist = std::sqrt(std::pow(hit2->Position()[0] - hit->Position()[0], 2) +
                                   std::pow(hit2->Position()[1] - hit->Position()[1], 2) +
                                   std::pow(hit2->Position()[2] - hit->Position()[2], 2));

            if (dist > 2 * fStripLength) {
              MF_LOG_DEBUG("StripSplitterAlg::getVirtualHits()")
                << " Distance between hit1 and hit2 " << dist << " > 2*fStripLength "
                << 2 * fStripLength;
              continue;
            }

            // for remaining hits, check if they overlap
            TVector3 stripDir2(0, 0, 0);
            if (jj == 0) {
              //strip
              std::array<double, 3> pt2 = {
                hit2->Position()[0], hit2->Position()[1], hit2->Position()[2]};
              std::pair<TVector3, TVector3> stripEnds2 = fGeo->GetStripEnds(pt2, hit2->CellID());
              stripDir2 = stripEnds2.first - stripEnds2.second;
            }

            // check if strips intersect
            TVector3 intercept = stripIntersect(hit, stripDir, hit2, stripDir2);
            // intercept found, calculate in which virtual cell
            if (intercept.Mag() > 0) {
              nSplitters++;
              float frac(-1);
              for (int ii = 0; ii < 3; ii++) {
                float dx = stripEnds.second[ii] - stripEnds.first[ii];
                if (std::fabs(dx) > 0.1) {
                  frac = (intercept[ii] - stripEnds.first[ii]) / dx;
                  break;
                }
              }

              if (frac >= 0.0 && frac <= 1.0) {
                int segment = int(frac * fnVirtual);
                if (segment >= 0 && segment < fnVirtual) {
                  if (virtEnergy.find(segment) != virtEnergy.end()) {
                    virtEnergy[segment] += hit2->Energy();
                  }
                  else {
                    virtEnergy[segment] = hit2->Energy();
                  }

                  if (fSaveStripEnds) {
                    float pos[3];
                    pos[0] = intercept.X();
                    pos[1] = intercept.Y();
                    pos[2] = intercept.Z();

                    const gar::rec::CaloHit* interhit = new gar::rec::CaloHit(0.1, 0., pos, 0., 0);
                    fIntersectionHits.emplace_back(interhit);
                  }
                }
                else {
                  MF_LOG_DEBUG("StripSplitterAlg::getVirtualHits()")
                    << "strange segment " << segment << " frac = " << frac
                    << " nvirt = " << fnVirtual;
                }
              }
              else {
                MF_LOG_DEBUG("StripSplitterAlg::getVirtualHits()") << "strange frac " << frac;
              }
            }
          }
        }

        MF_LOG_DEBUG("StripSplitterAlg::getVirtualHits()") << " Number of splitters " << nSplitters;

        // now create the virtual cells, and assign energy
        float totenergy(0);
        for (std::map<int, float>::iterator it = virtEnergy.begin(); it != virtEnergy.end(); ++it) {
          totenergy += it->second;
        }

        for (std::map<int, float>::iterator it = virtEnergy.begin(); it != virtEnergy.end(); ++it) {
          // energy of hit
          float energy = hit->Energy() * it->second / totenergy;
          // position of hit
          TVector3 virtualCentre = stripEnds.second - stripEnds.first;
          virtualCentre *= (it->first + 0.5) / fnVirtual;
          virtualCentre += stripEnds.first;

          float pos[3];
          pos[0] = virtualCentre.X();
          pos[1] = virtualCentre.Y();
          pos[2] = virtualCentre.Z();

          // make the new hit
          const gar::rec::CaloHit* newhit =
            new gar::rec::CaloHit(energy, hit->Time(), pos, hit->CellID(), hit->Layer());

          MF_LOG_DEBUG("StripSplitterAlg::getVirtualHits()")
            << " Creating new virtual hit pointing at " << newhit << " Energy " << energy
            << " Position " << pos[0] << " " << pos[1] << " " << pos[2];

          virtualhits.emplace_back(newhit);
        }
      }

      //----------------------------------------------------------------------------
      TVector3 StripSplitterAlg::stripIntersect(const gar::rec::CaloHit* hit0,
                                                const TVector3& dir0,
                                                const gar::rec::CaloHit* hit1,
                                                const TVector3& dir1)
      {
        // find intercept of hit1 with hit0
        // hit0 must be a strip
        // hit1 can be an orthogonal strip, or a cell
        // dir0,1 are direction of strip

        // centre position of cell/strip
        TVector3 stripCentre[2];
        stripCentre[0].SetXYZ(hit0->Position()[0], hit0->Position()[1], hit0->Position()[2]);
        stripCentre[1].SetXYZ(hit1->Position()[0], hit1->Position()[1], hit1->Position()[2]);

        const std::array<double, 3> pt = {
          hit0->Position()[0], hit0->Position()[1], hit0->Position()[2]};
        fStripLength = fGeo->getStripLength(pt, hit0->CellID());

        // direction of strip long axis
        // 0,0,0 for square cell
        TVector3 stripDir[2];
        stripDir[0] = dir0;
        stripDir[1] = dir1;

        // deal with cell case
        // define it's direction as perpendicular to strip and vector cell centre to origin
        bool isStrip[2];
        for (int i = 0; i < 2; i++) {
          if (stripDir[i].Mag() > 1e-10) { isStrip[i] = true; }
          else {
            isStrip[i] = false;
            stripDir[i] = stripDir[1 - i].Cross(stripCentre[i]);
          }
          // ensure dir is normalised
          stripDir[i] *= 1. / stripDir[i].Mag();
        }

        if (not isStrip[0]) {
          MF_LOG_ERROR("StripSplitterAlg::stripIntersect") << "first hit should be a strip";
          throw cet::exception("StripSplitterAlg::stripIntersect");
        }

        TVector3 p[2][2]; // ends of strips
        for (int j = 0; j < 2; j++) {
          for (int i = 0; i < 2; i++) {
            float ll = isStrip[j] ? fStripLength : fStripWidth * 1.1; // inflate a little for cells
            p[j][i] = stripCentre[j] - std::pow(-1, i) * 0.5 * ll * stripDir[j];
          }
        }

        TVector3 pNorm[2][2];
        for (int j = 0; j < 2; j++) {
          for (int i = 0; i < 2; i++) {
            float mm = p[j][i].Mag();
            pNorm[j][i] = p[j][i];
            pNorm[j][i] *= 1. / mm;
          }
        }

        TVector3 inPlane[2]; // difference between points: line inside plane
        for (int j = 0; j < 2; j++) {
          inPlane[j] = p[j][0] - p[j][1];
        }

        // vector normal to both lines (this is normal to the plane we're interested in)
        TVector3 normal = inPlane[0].Cross(inPlane[1]);
        float mag = normal.Mag();
        normal *= 1. / mag;

        // point on line [0]
        TVector3 point(p[0][0] + p[0][1]);
        point *= 0.5;

        // calculate the projected positions of ends of [1] on
        // the plane which contains "point" with normal "normal"
        TVector3 qPrime[2];
        for (int i = 0; i < 2; i++) {
          float d = ((point - p[1][i]).Dot(normal)) / (pNorm[1][i].Dot(normal));
          qPrime[i] = p[1][i] + d * pNorm[1][i];
        }

        // find the intersection of lines qPrime and p (they are in same plane)
        TVector3 a(p[0][1] - p[0][0]);
        TVector3 b(qPrime[1] - qPrime[0]);
        TVector3 c(qPrime[0] - p[0][0]);

        float factor = (c.Cross(b)).Dot(a.Cross(b)) / ((a.Cross(b)).Mag2());

        TVector3 x = p[0][0] + a * factor;

        // check that two lines really intercept
        bool intersect = true;
        for (int ii = 0; ii < 2; ii++) {
          for (int j = 0; j < 3; j++) {
            float d0 = ii == 0 ? x[j] - p[0][0][j] : x[j] - qPrime[0][j];
            float d1 = ii == 0 ? x[j] - p[0][1][j] : x[j] - qPrime[1][j];
            if (d0 * d1 > 1e-10) { intersect = false; }
          }
        }

        MF_LOG_DEBUG("StripSplitterAlg::stripIntersect") << "Intersection found " << intersect;

        if (intersect)
          return x;
        else
          return TVector3(0, 0, 0);
      }

    } // end namespace alg
  }   // end namespace rec
} // end namespace gar
