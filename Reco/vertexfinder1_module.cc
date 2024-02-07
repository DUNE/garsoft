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
#include "canvas/Persistency/Common/Assns.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Persistency/Common/PtrMaker.h"

#include "ReconstructionDataProducts/Vertex.h"
#include "Reco/vertex1algs.h"

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
      float fDCut;
      std::string fTrackLabel;
      int fNTPCClusCut;

      // Stub to check that the track is vertexable
      bool goodTrack(TrackPar trk, TrackEnd usebeg);

    };


    vertexfinder1::vertexfinder1(fhicl::ParameterSet const & p) : EDProducer{p} 
    {
      fTrackLabel = p.get<std::string>("TrackLabel","track");
      fChisquaredCut = p.get<float>("ChisquaredPerDOFCut",5.0);
      fRCut = p.get<float>("RCut", 12.0); // endpoints must be within this distnace in cm
      fDCut = p.get<float>("DCut", 2.0); // doca must be smaller than this in cm
      fNTPCClusCut = p.get<int>("NTPCClusCut",20);  // dimensionless
      fPrintLevel = p.get<int>("PrintLevel",0);

      art::InputTag trackTag(fTrackLabel);
      consumes< std::vector<rec::Track> >(trackTag);

      // Produce vertex, associated tracks & ends.  TrackEnd is defined in Track.h
      produces< std::vector<rec::Vertex> >();
      produces< art::Assns<rec::Track, rec::Vertex, rec::TrackEnd> >();
    }

    void vertexfinder1::produce(art::Event & e) {
      std::unique_ptr< std::vector<rec::Vertex> >
        vtxCol(new std::vector<rec::Vertex>);
      std::unique_ptr< art::Assns<rec::Track, rec::Vertex, rec::TrackEnd> >
        trkVtxAssns(new art::Assns<rec::Track, rec::Vertex, rec::TrackEnd>);

      auto trackHandle = e.getValidHandle< std::vector<Track>  >(fTrackLabel);
      auto const& tracks = *trackHandle;

      auto const vtxPtrMaker   = art::PtrMaker<rec::Vertex>(e);
      auto const trackPtrMaker = art::PtrMaker<rec::Track>(e, trackHandle.id());

      size_t nTrack = tracks.size();
      std::vector<int> usedflag_beg(nTrack,0);  // to tell if we've used a track yet or not.
      std::vector<int> usedflag_end(nTrack,0);  // to tell if we've used a track yet or not.

      if (nTrack>=2) {

        // For each track, build a vector of ends that can be vertexed with it.  Include
        // the number of the track for later ouroboros test

        std::vector< std::tuple<TrackPar,TrackEnd,int> > trackParEnds;
        for (size_t iTrack=0; iTrack<nTrack-1; ++iTrack) {
          TrackPar thisTrack = tracks.at(iTrack);
          // Build that vector for each of the 2 ends of the track
          for (int fin=0; fin<2; ++fin) {
            trackParEnds.clear();
            TrackEnd thisTrackEnd(fin == 0 ? 1 : 0);
            if (! goodTrack(thisTrack,thisTrackEnd)) continue;
            trackParEnds.emplace_back( std::make_tuple(thisTrack,thisTrackEnd,iTrack) );
            TVector3 thisEndinSpace;
            if (thisTrackEnd == TrackEndBeg) 
              { 
                usedflag_beg.at(iTrack) = 1; 
                thisEndinSpace = thisTrack.getXYZBeg();
              }
            else
              { 
                usedflag_end.at(iTrack) = 1; 
                thisEndinSpace = thisTrack.getXYZEnd();
              }

            // Loop over other trackends to see what is close
            for (size_t jTrack=iTrack+1; jTrack<nTrack; ++jTrack) {
              TrackPar thatTrack = tracks.at(jTrack);
              for (int fin2=0; fin2<2; ++fin2)
                {
                  TrackEnd thatTrackEnd(fin2 == 0 ? 1 : 0);
                  if (thatTrackEnd == TrackEndBeg && usedflag_beg.at(jTrack) != 0) continue;
                  if (thatTrackEnd == TrackEndEnd && usedflag_end.at(jTrack) != 0) continue;
                  if (!goodTrack(thatTrack,thatTrackEnd)) continue;
                  TVector3 thatEndinSpace;
                  if (thatTrackEnd == TrackEndBeg) 
                    { 
                      thatEndinSpace = thatTrack.getXYZBeg();
                    }
                  else
                    { 
                      thatEndinSpace = thatTrack.getXYZEnd();
                    }

                  if ( (thisEndinSpace-thatEndinSpace).Mag() > fRCut) continue;


                  // call the fitter to see what the docas are and cut on the doca

                  std::vector<TrackPar> vtracks;
                  std::vector<TrackEnd> usebeg;
                  vtracks.push_back(thisTrack);
                  vtracks.push_back(thatTrack);
                  usebeg.push_back(thisTrackEnd);
                  usebeg.push_back(thatTrackEnd);
                  std::vector<float> xyz(3,0);
                  std::vector< std::vector<float> > covmat;
                  std::vector<float> doca;
                  float chisquared=0;
                  double time;

                  if (gar::rec::fitVertex(vtracks,xyz,chisquared,covmat,time,usebeg,doca) == 0) 
                    {
                      if (doca.at(0) < fDCut && doca.at(1) < fDCut)
                        {
                          if (thatTrackEnd == TrackEndBeg)
                            {
                              usedflag_beg.at(jTrack) = 1;
                            }
                          else
                            {
                              usedflag_end.at(jTrack) = 1;
                            }
                          trackParEnds.emplace_back( std::make_tuple(thatTrack,thatTrackEnd,jTrack) );
                          continue;  // don't add the other end of the same track
                        } // end if doca test
                    } // end if vertex fit succeeded
                } // end loop over jtrack end
            } // End loop over jTrack
 
            size_t nTrackEnds = trackParEnds.size();

            if ( nTrackEnds>=2) {
              // Feed the vertexable ends into the vertex finder
              std::vector<TrackPar> vtracks;
              std::vector<TrackEnd> usebeg;
              for (size_t iTrackEnd=0; iTrackEnd<nTrackEnds; ++iTrackEnd) {
                vtracks.push_back( std::get<TrackPar>(trackParEnds[iTrackEnd]) );            
                usebeg.push_back ( std::get<1>(trackParEnds[iTrackEnd]) );            
              }
              std::vector<float> xyz(3,0);
              std::vector< std::vector<float> > covmat;
              std::vector<float> doca;
              float chisquared=0;
              double time;
              if (gar::rec::fitVertex(vtracks,xyz,chisquared,covmat,time,usebeg,doca) == 0) {
                float cmv[9];
                size_t icounter=0;
                for (size_t i=0; i<3;++i) {
                  for (size_t j=0; j<3; ++j) {
                    cmv[icounter] = covmat.at(i).at(j);
                  }
                }

                // Save vertex & its tracks into vectors that will be written 
                // to event.
                
                size_t IDcheck = 0;
                size_t nIDs = 0;

                for (size_t i=0; i<trackParEnds.size(); ++i) {
                  //TrackEnd endInVertex = std::get<1>(trackParEnds[i]);
                  auto const trackpointer = trackPtrMaker( std::get<2>(trackParEnds[i]) );
                  if(IDcheck!=trackpointer->getIDNumber()){
                    IDcheck = trackpointer->getIDNumber();
                    nIDs++;
                  }
                }
                if(nIDs<=1) continue;

                vtxCol->emplace_back(xyz.data(),cmv,time);
                auto const vtxpointer = vtxPtrMaker(vtxCol->size()-1);
               
                for (size_t i=0; i<trackParEnds.size(); ++i) {
                  TrackEnd endInVertex = std::get<1>(trackParEnds[i]);
                  auto const trackpointer = trackPtrMaker( std::get<2>(trackParEnds[i]) );
                  trkVtxAssns->addSingle(trackpointer,vtxpointer,endInVertex);
                }
              } // end if (fitVertex...)
            }
          } // end loop on 2 ends of iTrack
        } // end loop on iTrack
      }  // End test for nTrack>=2; might put empty vertex list & Assns into event
      for(size_t i = 0; i<trkVtxAssns->size(); i++){
        //auto Track = std::get<0>(trkVtxAssns->at(i));
        //auto Track = std::get<0>(trkVtxAssns->at(i));  
        //std::cout<<"Track ID Number: "<<Track->getIDNumber()<<"\n ";
        //std::cout<<"Trk ID Number: "<<Track->getIDNumber()<<" \n ";
      }
      e.put(std::move(vtxCol));
      e.put(std::move(trkVtxAssns));
    } // end produce method.


      //--------------------------------------------------------------------------------------------------

      // goodTrack returns true if the track can be used in vertexing

    bool vertexfinder1::goodTrack(TrackPar trk, TrackEnd /* usebeg */)
    {
      int trknclus = trk.getNTPCClusters();
      if ( trknclus < fNTPCClusCut ) return false;
      return true;     // just a stub for now
    }

    DEFINE_ART_MODULE(vertexfinder1)

  } // namespace rec
} // namespace gar
