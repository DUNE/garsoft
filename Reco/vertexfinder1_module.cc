////////////////////////////////////////////////////////////////////////
// Class:       vertexfinder1
// Plugin Type: producer (art v2_11_02)
// File:        vertexfinder1_module.cc
//
// Generated at Thu Jul 12 16:15:30 2018 by Thomas Junk using cetskelgen
// from cetlib version v3_03_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Persistency/Common/PtrMaker.h"

#include "ReconstructionDataProducts/Track.h"
#include "ReconstructionDataProducts/Vertex.h"
#include "Reco/TrackPar.h"

#include "TMath.h"
#include "TVectorF.h"
#include "TMatrixF.h"

#include <memory>

namespace gar {
  namespace rec {

    class vertexfinder1 : public art::EDProducer {
    public:
      explicit vertexfinder1(fhicl::ParameterSet const & p);
      // The compiler-generated destructor is fine for non-base
      // classes without bare pointers or other resource use.

      // Plugins should not be copied or assigned.
      vertexfinder1(vertexfinder1 const &) = delete;
      vertexfinder1(vertexfinder1 &&) = delete;
      vertexfinder1 & operator = (vertexfinder1 const &) = delete;
      vertexfinder1 & operator = (vertexfinder1 &&) = delete;

      // Required functions.
      void produce(art::Event & e) override;

    private:

      float fChisquaredCut;
      int fPrintLevel;
      float fRCut;
      std::string fTrackLabel;

      // Stub to check that the track is vertexable
      bool goodTrack(TrackPar trk, bool usebeg);


      // fit an arbitrary number of tracks to a single vertex
      // returns 0 if success, 1 if failure
      // covmat is 3x3

      int fitVertex(std::vector<TrackPar> &tracks, std::vector<float> &xyz, float &chisquared, std::vector<std::vector<float> > &covmat, ULong64_t &time,
            std::vector<bool> usebeg);

    };


    vertexfinder1::vertexfinder1(fhicl::ParameterSet const & p)
    // :
    {
      produces< std::vector<rec::Vertex> >();
      produces< art::Assns<rec::Track, rec::Vertex> >();

      fTrackLabel = p.get<std::string>("TrackLabel","track");
      fChisquaredCut = p.get<float>("ChisquaredPerDOFCut",5.0);
      fRCut = p.get<float>("RCut", 12.0); // in cm
      fPrintLevel = p.get<int>("PrintLevel",0);
    }

    void vertexfinder1::produce(art::Event & e) {
      std::unique_ptr< std::vector<rec::Vertex> > vtxCol(new std::vector<rec::Vertex>);
      std::unique_ptr< art::Assns<rec::Track,rec::Vertex> > trkVtxAssns(new ::art::Assns<rec::Track,rec::Vertex>);

      auto trackHandle = e.getValidHandle< std::vector<Track>  >(fTrackLabel);
      auto const& tracks = *trackHandle;

      auto const vtxPtrMaker = art::PtrMaker<rec::Vertex>(e);
      auto const trackPtrMaker = art::PtrMaker<rec::Track>(e, trackHandle.id());



      size_t nTrack = tracks.size();
      if (nTrack>=2) {
        // For each track, build a vector of ends that can be vertexed with it.  Include
        // the number of the track for later ouroboros test
        std::vector< std::tuple<TrackPar,bool,int> > trackParEnds;
        for (size_t iTrack=0; iTrack<nTrack-1; ++iTrack) {
          TrackPar thisTrack = tracks.at(iTrack);
          // Build that vector for each of the 2 ends of the track
          for (int fin=0; fin<2; ++fin) {
            bool thisUseBeg = fin==0 ? true : false;
            trackParEnds.clear();
            if (goodTrack(thisTrack,thisUseBeg)) {
              trackParEnds.emplace_back( std::make_tuple(thisTrack,thisUseBeg,iTrack) );
            }
            TVector3 thisEnd;
            if (thisUseBeg) {
              thisEnd = thisTrack.getXYZBeg();
            } else {
              thisEnd = thisTrack.getXYZEnd();
            }
            // Loop over other trackends to see what is close
            for (size_t jTrack=iTrack+1; jTrack<nTrack; ++jTrack) {
              TrackPar thatTrack = tracks.at(jTrack);
              TVector3 thatEnd;
              if (goodTrack(thatTrack,true)) {
                thatEnd = thatTrack.getXYZBeg();
                double dist = (thisEnd -thatEnd).Mag();
                if (dist<fRCut) {
                  trackParEnds.emplace_back( std::make_tuple(thatTrack,true,jTrack) );
                }
              }
              if (goodTrack(thatTrack,false)) {
                thatEnd = thatTrack.getXYZEnd();
                double dist = (thisEnd -thatEnd).Mag();
                if (dist<fRCut) {
                  trackParEnds.emplace_back( std::make_tuple(thatTrack,false,jTrack) );
                }
              }
            } // End loop over 2 trackends for thatTrack
 


            // Having found all the track ends that can be vertexed with thisTrack on this
            // particular end (fin), it is time for the Ouroboros test: be certain that we
            // are not vertexing one end of a looping track to its own other end
            bool Ouroboros = false;
            size_t nTrackEnds = trackParEnds.size();
            for (size_t iTrackEnd=0; iTrackEnd<nTrackEnds-1; ++iTrackEnd) {
              for (size_t jTrackEnd=iTrackEnd+1; jTrackEnd<nTrackEnds; ++jTrackEnd) {
                Ouroboros = (std::get<int>(trackParEnds[iTrackEnd]) == std::get<int>(trackParEnds[jTrackEnd]));
                // Double break!
                if (Ouroboros) break;
              }
              if (Ouroboros) break;
            }



            if (!Ouroboros && nTrackEnds>=2) {
              // Feed the vertexable ends into the vertex finder
              std::vector<TrackPar> vtracks;
              std::vector<bool> usebeg;
              for (size_t iTrackEnd=0; iTrackEnd<nTrackEnds; ++iTrackEnd) {
                vtracks.push_back( std::get<TrackPar>(trackParEnds[iTrackEnd]) );            
                usebeg.push_back ( std::get<bool>(trackParEnds[iTrackEnd]) );            
              }
              std::vector<float> xyz(3,0);
              std::vector< std::vector<float> > covmat;
              float chisquared=0;
              ULong64_t time;
              if (fitVertex(vtracks,xyz,chisquared,covmat,time,usebeg) == 0) {
                float cmv[9];
                size_t icounter=0;
                for (size_t i=0; i<3;++i) {
                  for (size_t j=0; j<3; ++j) {
                    cmv[icounter] = covmat.at(i).at(j);
                  }
                }

                vtxCol->emplace_back(xyz.data(),cmv,time);  
                auto const vtxpointer = vtxPtrMaker(vtxCol->size()-1);
                for (size_t i=0; i<trackParEnds.size(); ++i) {
                  auto const trackpointer = trackPtrMaker(i);
                  trkVtxAssns->addSingle(trackpointer,vtxpointer);
                }
              } // end if (fitVertex...)
            }
          } // end loop on fin
        } // end loop on iTrack
      }  // End test for nTrack>=2; might put empty vertex list & Assns into event
      e.put(std::move(vtxCol));
      e.put(std::move(trkVtxAssns));
    } // end produce method.

    //--------------------------------------------------------------------------------------------------

    // fitVertex returns 0 if success
    // It does not cut tracks from the vertex candidate set -- may need to do that in the caller.

    int vertexfinder1::fitVertex(std::vector<TrackPar> &tracks, 
                   std::vector<float> &xyz, 
                   float &chisquared, 
                   std::vector< std::vector<float> > &covmat, 
                   ULong64_t &time,
                   std::vector<bool> usebeg)
    {
      // find the ends of the tracks that are closest to the other tracks' ends
      // pick the end that minimizes the sums of the min(dist) to the other tracks' endpoints

      // todo -- iterate on this as the tracks are helices and not lines -- move to close point on track and re-do linear extrapolation

      if (tracks.size() < 2) return(1);  // need at least two tracks to make a vertex

      // find the average time, protecting against overflows of the ULong64_t by subtracting off the smallest.
      std::vector<ULong64_t> times;
      for (size_t i=0; i<tracks.size(); ++i)
    {
      times.push_back(tracks.at(i).getTime());
    }
      std::sort(times.begin(),times.end());
      ULong64_t diffsum=0;
      for (size_t i=1; i<times.size(); ++i)
    {
      diffsum += (times.at(i) - times.at(0)); 
    }
      time = times.at(0) + diffsum/times.size();

      xyz.resize(3);
      covmat.resize(9);

      std::vector<TVectorF> dir;
      std::vector<TVectorF> p;

      for (size_t itrack = 0; itrack < tracks.size(); ++itrack)
    {


      TVector3 trackbeg = tracks.at(itrack).getXYZBeg();
      TVector3 trackend = tracks.at(itrack).getXYZEnd();
      float phi=0;
      float si=0;
      if (usebeg.at(itrack))
        {
          float tmppos[3] = {(float) trackbeg.X(), (float) trackbeg.Y(), (float) trackbeg.Z()};
          p.emplace_back(3,tmppos);
          phi = tracks.at(itrack).getTrackParametersBegin()[3];
          si =   TMath::Tan(tracks.at(itrack).getTrackParametersBegin()[4]);
        }
      else
        {
          float tmppos[3] = {(float) trackend.X(), (float) trackend.Y(), (float) trackend.Z()};
          p.emplace_back(3,tmppos);
          phi = tracks.at(itrack).getTrackParametersEnd()[3];
          si =   TMath::Tan(tracks.at(itrack).getTrackParametersEnd()[4]);
        }
      TVectorF dtmp(3);
      dtmp[0] = si;
      dtmp[1] = TMath::Sin(phi);
      dtmp[2] = TMath::Cos(phi);
      float norminv = 1.0/TMath::Sqrt(dtmp.Norm2Sqr());
      dtmp *= norminv;
      dir.push_back(dtmp);
    }

      TMatrixF I(3,3);
      I[0][0] = 1;
      I[1][1] = 1;
      I[2][2] = 1;

      TMatrixF A(3,3);
      TMatrixF VVT(3,3);
      TMatrixF IMVVT(3,3);
      TVectorF vsum(3);
      vsum.Zero();
      A *= 0;
       
      for (size_t itrack=0; itrack < tracks.size(); ++itrack)
    {
      for (size_t i=0; i<3; ++i)
        {
          for (size_t j=0; j<3; ++j)
        {
          VVT[i][j] = dir.at(itrack)[i]*dir.at(itrack)[j];
        }
        }
      IMVVT = I - VVT;
      A += IMVVT;
      vsum += IMVVT*p.at(itrack);
    }
      double det;
      A.Invert(&det);
      if (det == 0) return(1);
      TVectorF xyzsol = A*vsum;
      xyz.at(0) = xyzsol[0];
      xyz.at(1) = xyzsol[1];
      xyz.at(2) = xyzsol[2];

      chisquared = 0;
      for (size_t itrack=0; itrack < tracks.size(); ++itrack)
    {
      
      TVector3 dir3(dir.at(itrack)[0],dir.at(itrack)[1],dir.at(itrack)[2]);
      TVectorF diff = xyzsol - p.at(itrack);
      TVector3 diff3(diff[0],diff[1],diff[2]);
      chisquared +=  (diff3.Cross(dir3)).Mag2();   // todo -- include track errors in chisquared calc
    }

      // to iterate -- extrapolate helical tracks to the closest point to the first found vertex and
      // run the fitter again

      // to do: figure out what the vertex covariance matrix is

      covmat.clear();
      for (size_t i=0;i<3;++i)
    {
      std::vector<float> cmr;
      for (size_t j=0; j<3; ++j)
        {
          cmr.push_back(0);
        }
      covmat.push_back(cmr);
    }
      return 0;
       
    }

    //--------------------------------------------------------------------------------------------------

    // goodTrack returns true if the track can be used in vertexing

    bool vertexfinder1::goodTrack(TrackPar trk, bool usebeg)
    {
      return true;     // just a stub for now
    }

    DEFINE_ART_MODULE(vertexfinder1)

  } // namespace rec
} // namespace gar
