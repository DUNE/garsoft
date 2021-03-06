////////////////////////////////////////////////////////////////////////
// Class:       tpcpatrec2
// Plugin Type: producer (art v3_00_00)
// File:        tpcpatrec2_module.cc
//
// Generated at Tue Feb  5 08:57:00 2019 by Thomas Junk using cetskelgen
// from cetlib version v3_04_00.
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
#include "canvas/Persistency/Common/FindManyP.h"

#include <memory>
#include <set>

#include "TMath.h"
#include "TVector3.h"

#include "Geant4/G4ThreeVector.hh"
#include "nutools/MagneticField/MagneticField.h"

// GArSoft Includes
#include "ReconstructionDataProducts/TPCCluster.h"
#include "ReconstructionDataProducts/VecHit.h"
#include "ReconstructionDataProducts/Track.h"
#include "Reco/TrackPar.h"
#include "Reco/tracker2algs.h"

namespace gar {
  namespace rec {

    class tpcpatrec2 : public art::EDProducer {
    public:
      explicit tpcpatrec2(fhicl::ParameterSet const& p);
      // The compiler-generated destructor is fine for non-base
      // classes without bare pointers or other resource use.

      // Plugins should not be copied or assigned.
      tpcpatrec2(tpcpatrec2 const&) = delete;
      tpcpatrec2(tpcpatrec2&&) = delete;
      tpcpatrec2& operator=(tpcpatrec2 const&) = delete;
      tpcpatrec2& operator=(tpcpatrec2&&) = delete;

      // Required functions.
      void produce(art::Event& e) override;

    private:

      // Declare member data here.

      std::string fVecHitLabel;     ///< label of module to get the vector hits and associations with hits
      std::string fTPCClusterLabel; ///< label of module to get the TPC clusters
      int    fPrintLevel;           ///< debug printout:  0: none, 1: selected, 2: all
      float  fVecHitMatchCos;       ///< matching condition for pairs of vector hits cos angle between directions
      float  fVecHitMatchPos;       ///< matching condition for pairs of vector hits -- 3D distance (cm)
      float  fVecHitMatchPEX;       ///< matching condition for pairs of vector hits -- miss distance (cm)
      float  fVecHitMatchEta;       ///< matching condition for pairs of vector hits -- eta match (cm)
      float  fVecHitMatchLambda;    ///< matching condition for pairs of vector hits -- dLambda (radians)
      unsigned int fInitialTPNTPCClusters; ///< number of hits to use for initial trackpar estimate, if present
      size_t fMinNumTPCClusters;           ///< minimum number of hits for a patrec track
      float  fSortTransWeight;      ///< for use in the hit sorting algorithm -- transverse distance weight factor
      float  fSortDistBack;         ///< for use in the hit sorting algorithm -- how far to go back before raising the distance figure of merit
      float  fCloseEtaUnmatch;      ///< distance to look for vector hits that don't match in eta. 

      bool fXBasedSecondPass;       ///< switch to tell us to use the X-based second-pass of patrec
      int fXBasedMinTPCClusters;    ///< min number of TPC clusters to require on new tracks found with the second pass
      size_t fPatRecLookBack1;      ///< n hits to look backwards to make a linear extrapolation
      size_t fPatRecLookBack2;      ///< extrapolate from lookback1 to lookback2 and see how close the new hit is to the line
      float  fHitResolYZ;           ///< resolution in cm of a hit in YZ (pad size)
      float  fHitResolX;            ///< resolution in cm of a hit in X (drift direction)
      float  fSigmaRoad;            ///< how many sigma away from a track a hit can be and still add it during patrec

      // criteria for associating vector hits together to form clusters
      bool vhclusmatch(const std::vector<gar::rec::VecHit> &vechits, std::vector<size_t> &cluster, size_t vh);

      // rough estimate of track parameters
      int makepatrectrack(std::vector<gar::rec::TPCCluster> &tpcclusters, gar::rec::TrackPar &trackpar);

      float calceta2d(gar::rec::VecHit &vhtest, gar::rec::VecHit &vh);

      // test to see if a cluster looks like a conversion and split it in two if it does.  Returns true
      // if the cluster is identified as a conversion.

      bool conversion_test_split(const std::vector<gar::rec::VecHit> &vechits,
                                 std::vector<size_t> &cluster,
                                 std::vector<size_t> &splitclus1,
                                 std::vector<size_t> &splitclus2);

    };


    tpcpatrec2::tpcpatrec2(fhicl::ParameterSet const& p) : EDProducer{p}  
      {
        fVecHitLabel       = p.get<std::string>("VecHitLabel","vechit");
        fTPCClusterLabel   = p.get<std::string>("TPCClusterLabel","tpcclusterpass1");
        fPrintLevel        = p.get<int>("PrintLevel",0);
        fVecHitMatchCos    = p.get<float>("VecHitMatchCos",0.9);
        fVecHitMatchPos    = p.get<float>("VecHitMatchPos",20.0);
        fVecHitMatchPEX    = p.get<float>("VecHitMatchPEX",5.0);
        fVecHitMatchEta    = p.get<float>("VecHitMatchEta",1.0);
        fVecHitMatchLambda = p.get<float>("VecHitMatchLambda",0.1);
        fInitialTPNTPCClusters    = p.get<unsigned int>("InitialTPNTPCClusters",100);
        fMinNumTPCClusters        = p.get<size_t>("MinNumTPCClusters",20);
        fSortTransWeight   = p.get<float>("SortTransWeight",0.1);
        fSortDistBack      = p.get<float>("SortDistBack",2.0);
        fCloseEtaUnmatch   = p.get<float>("CloseEtaUnmatch",20.0);
        fXBasedSecondPass  = p.get<bool>("XBasedSecondPass",true);
        fXBasedMinTPCClusters = p.get<int>("XBasedMinTPCClusters",10);
        fPatRecLookBack1     = p.get<size_t>("PatRecLookBack1",5);
        fPatRecLookBack2     = p.get<size_t>("PatRecLookBack2",10);
        if (fPatRecLookBack1 == fPatRecLookBack2)
          {
            throw cet::exception("tpcpatrec2_module: PatRecLookBack1 and PatRecLookBack2 are the same");
          }
        fHitResolYZ        = p.get<float>("HitResolYZ",1.0); // TODO -- think about what this value is
        fHitResolX         = p.get<float>("HitResolX",0.5);  // this is probably much better
        fSigmaRoad         = p.get<float>("SigmaRoad",5.0);

        art::InputTag vechitTag(fVecHitLabel);
        consumes< std::vector<gar::rec::VecHit> >(vechitTag);
        consumes< art::Assns<gar::rec::TPCCluster, gar::rec::VecHit> >(vechitTag);
        produces< std::vector<gar::rec::Track> >();
        produces< art::Assns<gar::rec::TPCCluster, gar::rec::Track> >();
        produces< art::Assns<gar::rec::VecHit, gar::rec::Track> >();
      }

    void tpcpatrec2::produce(art::Event& e)
    {
      std::unique_ptr< std::vector<gar::rec::Track> > trkCol(new std::vector<gar::rec::Track>);
      std::unique_ptr< art::Assns<gar::rec::VecHit,gar::rec::Track> > vhTrkAssns(new ::art::Assns<gar::rec::VecHit,gar::rec::Track>);
      std::unique_ptr< art::Assns<gar::rec::TPCCluster,gar::rec::Track> > TPCClusterTrkAssns(new ::art::Assns<gar::rec::TPCCluster,gar::rec::Track>);

      auto vechitHandle = e.getValidHandle< std::vector<gar::rec::VecHit> >(fVecHitLabel);
      auto const& vechits = *vechitHandle;

      auto TPCClusterHandle = e.getValidHandle< std::vector<gar::rec::TPCCluster> >(fTPCClusterLabel);
      auto const& e_tpcclusters = *TPCClusterHandle;  // TPC Clusters from the event

      auto const trackPtrMaker = art::PtrMaker<gar::rec::Track>(e);
      auto const vhPtrMaker = art::PtrMaker<gar::rec::VecHit>(e, vechitHandle.id());
      auto const clusPtrMaker = art::PtrMaker<gar::rec::TPCCluster>(e,TPCClusterHandle.id());

      const art::FindManyP<gar::rec::TPCCluster> TPCClustersFromVecHits(vechitHandle,e,fVecHitLabel);

      art::ServiceHandle<mag::MagneticField> magFieldService;
      G4ThreeVector zerovec(0,0,0);
      G4ThreeVector magfield = magFieldService->FieldAtPoint(zerovec);

      std::set<gar::rec::IDNumber> TPCClusNumerosUsedInFirstPass;

      // stitch together vector hits into tracks
      // question -- do we need to iterate this, first looking for small-angle matches, and then
      // loosening up?  Also need to address tracks that get stitched across the primary vertex -- may need to split
      // these after other tracks have been found.

      std::vector< std::vector< size_t > > vhclusters;  // change this to use just indices of vh's

      for (size_t ivh = 0; ivh< vechits.size(); ++ ivh)
        {
          if (fPrintLevel>1)
            {
              std::cout << " vhprint " << vechits[ivh].Position()[0] << " " <<  vechits[ivh].Position()[1] << " " <<  vechits[ivh].Position()[2] << " " <<
                vechits[ivh].Direction()[0] << " " <<  vechits[ivh].Direction()[1] << " " <<  vechits[ivh].Direction()[2] << " "  << std::endl;
            }

          std::vector<size_t> clusmatchlist;
          for (size_t iclus=0; iclus<vhclusters.size(); ++iclus)
            {
              if (vhclusmatch(vechits,vhclusters[iclus],ivh))
                {
                  clusmatchlist.push_back(iclus);
                }
            }
          if (clusmatchlist.size() == 0)
            {
              std::vector<size_t> newclus;
              newclus.push_back(ivh);
              vhclusters.push_back(newclus);
            }
          else if (clusmatchlist.size() == 1)
            {
              vhclusters[clusmatchlist[0]].push_back(ivh);
            }
          else   // multiple matches -- merge clusters togetehr
            {
              for (size_t icm=1; icm<clusmatchlist.size(); ++icm)
                {
                  for (size_t ivh2=0; ivh2<vhclusters[clusmatchlist[icm]].size(); ++ivh2)
                    {
                      vhclusters[clusmatchlist[0]].push_back(vhclusters[clusmatchlist[icm]][ivh2]);
                    }
                }
              // remove the merged vh clusters, using the new indexes after removing earlier ones
              for (size_t icm=1; icm<clusmatchlist.size(); ++icm)
                {
                  vhclusters.erase(vhclusters.begin() + (clusmatchlist[icm]-icm+1));
                }
            }
        }

      //std::cout << "Before conversion check, clusters: " << std::endl;
      //for (size_t iclus=0; iclus<vhclusters.size(); ++iclus)
      //    {
      // std::cout << "Cluster " << iclus << std::endl;
      // for (size_t ivh = 0; ivh< vhclusters.at(iclus).size(); ++ ivh)
      //   {
      //     size_t i=vhclusters.at(iclus).at(ivh);
      //     std::cout << "   " << ivh << " " <<
      //     vechits[i].Position()[0] << " " <<  vechits[i].Position()[1] << " " <<  vechits[i].Position()[2] << " " <<
      //     vechits[i].Direction()[0] << " " <<  vechits[i].Direction()[1] << " " <<  vechits[i].Direction()[2] << " "  << std::endl;
      //   }
      //    }

      // look for conversions with two tracks that have been joined together and split them in two.
      
      std::vector<size_t> identified_conversions;       // index into vhclusters of identified conversions
      std::vector< std::vector< size_t > > splitclus1;  // vhcluster for one leg after split 
      std::vector< std::vector< size_t > > splitclus2;  // vhcluster for the other leg after split

      for (size_t iclus=0; iclus<vhclusters.size(); ++iclus)
        {
          std::vector<size_t> scc1;
          std::vector<size_t> scc2;
          if (conversion_test_split(vechits, vhclusters.at(iclus), scc1, scc2))
            {
              identified_conversions.push_back(iclus);
              splitclus1.push_back(scc1);
              splitclus2.push_back(scc2);
            }
        }

      // rearrange the cluster list -- put the two legs separately in the list, replacing the incorrectly-merged one

      for (size_t icc=0; icc<identified_conversions.size(); ++icc)
        {
          vhclusters.at(identified_conversions.at(icc)) = splitclus1.at(icc);
          vhclusters.push_back(splitclus2.at(icc));
        }


      // make a local list of TPCClusters for each track and find initial track parameters
      for (size_t iclus=0; iclus < vhclusters.size(); ++iclus)
        {
          std::vector<gar::rec::TPCCluster> TPCClusters;
          std::vector<art::Ptr<gar::rec::TPCCluster> > TPCClusterptrs;
          for (size_t ivh=0; ivh<vhclusters[iclus].size(); ++ivh)
            {
              for (size_t iTPCCluster=0; iTPCCluster<TPCClustersFromVecHits.at(vhclusters[iclus][ivh]).size(); ++iTPCCluster)
                {
                  TPCClusterptrs.push_back(TPCClustersFromVecHits.at(vhclusters[iclus][ivh]).at(iTPCCluster));
                  TPCClusters.push_back(*TPCClusterptrs.back());
                }
            }

          if (TPCClusters.size() >= fMinNumTPCClusters)
            {
              gar::rec::TrackPar trackpar;
              if ( makepatrectrack(TPCClusters,trackpar) == 0 )
                {
                  trkCol->push_back(trackpar.CreateTrack());
                  auto const trackpointer = trackPtrMaker(trkCol->size()-1);
                  for (size_t iTPCCluster=0; iTPCCluster<TPCClusters.size(); ++iTPCCluster)
                    {
                      TPCClusNumerosUsedInFirstPass.insert(TPCClusters.at(iTPCCluster).getIDNumber());
                      TPCClusterTrkAssns->addSingle(TPCClusterptrs.at(iTPCCluster),trackpointer);
                    }
                  for (size_t ivh=0; ivh<vhclusters[iclus].size(); ++ivh)
                    {
                      auto const vhpointer = vhPtrMaker(ivh);
                      vhTrkAssns->addSingle(vhpointer,trackpointer);
                    }
                }
            }
        }

      // to think about -- remove dodgy tracks so that the second pass may have a better chance at picking them up.
      // any remaining misclassified conversions, for instance.  Or just don't put them on the list in the
      // first place.

      // Second pass patrec: 
      // go back and use the X-based pattern recognition on TPC clusters not included in tracks found above

      if (fXBasedSecondPass)
        {
          // first find unassociated TPC clusters and sort them by their X positions
          // keep them in place, just make a vector of sorted indices

          std::vector<gar::rec::TPCCluster> tcc2pass;
          std::vector<float> clusx;
          std::vector<size_t> tcc2passidx;

          float roadsq = fSigmaRoad*fSigmaRoad;
          float resolSq = fHitResolYZ*fHitResolYZ;

          for (size_t iclus=0; iclus < e_tpcclusters.size(); ++iclus)
            {
              if (TPCClusNumerosUsedInFirstPass.count(e_tpcclusters.at(iclus).getIDNumber()) == 0)
                {
                  tcc2pass.push_back(e_tpcclusters.at(iclus));
                  tcc2passidx.push_back(iclus);  // for use in making associations
                  clusx.push_back(tcc2pass.back().Position()[0]);
                }
            }

          if (tcc2pass.size()>0)
            {
              std::vector<int> csi(clusx.size());
              TMath::Sort((int) clusx.size(),clusx.data(),csi.data());

              // do this twice, once going forwards through the tcc2pass, once backwards, and then collect hits
              // in groups and split them when either the forwards or the backwards list says to split them.

              std::vector< std::vector<int> > cluslistf;
              std::vector<int> trackf(tcc2pass.size());
              for (size_t iclus=0; iclus<tcc2pass.size(); ++iclus)
                {
                  const float *cpos = tcc2pass[csi[iclus]].Position();
                  TVector3 hpvec(cpos);

                  float bestsignifs = -1;
                  int ibest = -1;
                  for (size_t itcand = 0; itcand < cluslistf.size(); ++itcand)
                    {
                      float signifs = 1E9;
                      //size_t clsiz = cluslistf[itcand].size();
                      //if (clsiz > fPatRecLookBack1 && clsiz > fPatRecLookBack2)
                      //  {
                      //  TVector3 pt1( tcc2pass[csi[cluslistf[itcand][clsiz-fPatRecLookBack1]]].Position() );
                      //  TVector3 pt2( tcc2pass[csi[cluslistf[itcand][clsiz-fPatRecLookBack2]]].Position() );
                      //  TVector3 uv = pt1-pt2;
                      //  uv *= 1.0/uv.Mag();
                      //  signifs = ((hpvec-pt1).Cross(uv)).Mag2()/resolSq;
                      //  }
                      //else // not enough tcc2pass, just look how close we are to the last one
                      //  {
                      const float *hpos = tcc2pass[csi[cluslistf[itcand].back()]].Position();
                      signifs = (TMath::Sq( (hpos[1]-cpos[1]) ) +
                                 TMath::Sq( (hpos[2]-cpos[2]) ))/resolSq;
                      //  }
                      if (bestsignifs < 0 || signifs < bestsignifs)
                        {
                          bestsignifs = signifs;
                          ibest = itcand;
                        }
                    }
                  if (ibest == -1 || bestsignifs > roadsq)  // start a new track if we're not on the road, or if we had no tracks to begin with
                    {
                      ibest = cluslistf.size();
                      std::vector<int> vtmp;
                      cluslistf.push_back(vtmp);
                    }
                  cluslistf[ibest].push_back(iclus);
                  trackf[iclus] = ibest;
                }

              std::vector< std::vector<int> > cluslistb;
              std::vector<int> trackb(tcc2pass.size());
              for (int iclus=tcc2pass.size()-1; iclus >= 0; --iclus)
                {
                  const float *cpos = tcc2pass[csi[iclus]].Position();
                  TVector3 hpvec(cpos);

                  float bestsignifs = -1;
                  int ibest = -1;
                  for (size_t itcand = 0; itcand < cluslistb.size(); ++itcand)
                    {
                      float signifs = 1E9;
                      //size_t clsiz = cluslistb[itcand].size();
                      //if (clsiz > fPatRecLookBack1 && clsiz > fPatRecLookBack2)
                      //{
                      //    TVector3 pt1( tcc2pass[csi[cluslistb[itcand][clsiz-fPatRecLookBack1]]].Position() );
                      //    TVector3 pt2( tcc2pass[csi[cluslistb[itcand][clsiz-fPatRecLookBack2]]].Position() );
                      //    TVector3 uv = pt1-pt2;
                      //    uv *= 1.0/uv.Mag();
                      //    signifs = ((hpvec-pt1).Cross(uv)).Mag2()/resolSq;
                      //  }
                      //else // not enough tcc2pass, just look how close we are to the last one
                      //  {
                      const float *hpos = tcc2pass[csi[cluslistb[itcand].back()]].Position();
                      signifs = (TMath::Sq( (hpos[1]-cpos[1]) ) +
                                 TMath::Sq( (hpos[2]-cpos[2]) ))/resolSq;
                      //  }

                      if (bestsignifs < 0 || signifs < bestsignifs)
                        {
                          bestsignifs = signifs;
                          ibest = itcand;
                        }
                    }
                  if (ibest == -1 || bestsignifs > roadsq)  // start a new track if we're not on the road, or if we had no tracks to begin with
                    {
                      ibest = cluslistb.size();
                      std::vector<int> vtmp;
                      cluslistb.push_back(vtmp);
                    }
                  cluslistb[ibest].push_back(iclus);
                  trackb[iclus] = ibest;
                }

              // make a list of tracks that is the set of disjoint subsets

              std::vector< std::vector<int> > cluslist;
              for (size_t itrack=0; itrack<cluslistf.size(); ++itrack)
                {
                  int itrackl = 0;
                  for (size_t iclus=0; iclus<cluslistf[itrack].size(); ++iclus)
                    {
                      int ihif = cluslistf[itrack][iclus];
                      if (iclus == 0 || (trackb[ihif] != itrackl))
                        {
                          std::vector<int> vtmp;
                          cluslist.push_back(vtmp);
                          itrackl = trackb[ihif];
                        }
                      cluslist.back().push_back(ihif);
                    }
                }

              // now make patrec tracks out of these grouped TPC clusters

              for (size_t itrack=0; itrack<cluslist.size(); ++itrack)
                {
                  if (cluslist.at(itrack).size() < (size_t) fXBasedMinTPCClusters) continue;
                  std::vector<gar::rec::TPCCluster> tcl;
                  for (size_t iclus=0; iclus<cluslist.at(itrack).size(); ++iclus)
                    {
                      tcl.push_back(tcc2pass.at(csi[cluslist[itrack].at(iclus)]));
                    }
                  gar::rec::TrackPar trackpar;
                  if (makepatrectrack(tcl,trackpar) == 0)
                    {
                      trkCol->push_back(trackpar.CreateTrack());
                      auto const trackptr = trackPtrMaker(trkCol->size()-1);
                      for (size_t iclus=0; iclus<tcl.size(); ++iclus)
                        {
                          auto const clusptr = clusPtrMaker(tcc2passidx.at(csi[cluslist[itrack][iclus]]));
                          TPCClusterTrkAssns->addSingle(clusptr,trackptr);
                          if (fPrintLevel>1)
                            {
                              std::cout << " second-pass track " << itrack << " clus: " << iclus << " :  ";
                              TVector3 cp(tcc2pass[csi[cluslist[itrack][iclus]]].Position());
                              std::cout << cp.X() << " " << cp.Y() << " " << cp.Z() << std::endl;
                            }
                        }
                    }
                }

            }  // end check if we have any TPC clusters for pass 2
        } // end check if we want to run pass 2 patrec

      e.put(std::move(trkCol));
      e.put(std::move(vhTrkAssns));
      e.put(std::move(TPCClusterTrkAssns));

    }


    bool tpcpatrec2::vhclusmatch(const std::vector<gar::rec::VecHit> &vechits, std::vector<size_t> &cluster, size_t vhindex)
    {
      gar::rec::VecHit vh = vechits[vhindex];
      TVector3 vhdir(vh.Direction());
      TVector3 vhpos(vh.Position());

      bool foundmatch = false;
      for (size_t ivh=0; ivh<cluster.size(); ++ivh)
        {
          gar::rec::VecHit vhtest = vechits[cluster[ivh]];
          TVector3 vhtestdir(vhtest.Direction());
          TVector3 vhtestpos(vhtest.Position());

          if (fPrintLevel > 1)
            {
              std::cout << "Testing vh " << ivh << " in a cluster of size: " << cluster.size() << std::endl;
            }

          // require the two VH's directions to point along each other -- use dot product

          if (TMath::Abs((vhdir).Dot(vhtestdir)) < fVecHitMatchCos)
            {
              if (fPrintLevel > 1)
                {
                  std::cout << " Dot failure: " << TMath::Abs((vhdir).Dot(vhtestdir)) << std::endl;
                }
              continue;
            }

          // require the positions to be within fVecHitMatchPos of each other

          if ((vhpos-vhtestpos).Mag() > fVecHitMatchPos)
            {
              if (fPrintLevel > 1)
                {
                  std::cout << " Pos failure: " << (vhpos-vhtestpos).Mag() << std::endl;
                }
              continue;
            }

          // require the extrapolation of one VH's line to another VH's center to match up.  Do for
          // both VH's.

          if ( ((vhpos-vhtestpos).Cross(vhdir)).Mag() > fVecHitMatchPEX )
            {
              if (fPrintLevel > 1)
                {
                  std::cout << "PEX failure: " << ((vhpos-vhtestpos).Cross(vhdir)).Mag() << std::endl;
                }
              continue;
            }
          if ( ((vhpos-vhtestpos).Cross(vhtestdir)).Mag() > fVecHitMatchPEX )
            {
              if (fPrintLevel > 1)
                {
                  std::cout << "PEX failure: " << ((vhpos-vhtestpos).Cross(vhtestdir)).Mag() << std::endl;
                }
              continue;
            }

          float eta = calceta2d(vhtest,vh);

          if ( eta > fVecHitMatchEta )
            {
              if (fPrintLevel > 1)
                {
                  std::cout << "Eta failure: " << eta << std::endl;
                }
              continue;
            }


          //----------------
          // lambda requirement

          float vhpd = TMath::Sqrt( TMath::Sq(vhdir.Y()) + TMath::Sq(vhdir.Z()) );
          float vhxd = TMath::Abs( vhdir.X() );
          float vhlambda = TMath::Pi()/2.0;
          if (vhpd >0) vhlambda = TMath::ATan(vhxd/vhpd);

          float cvhpd = TMath::Sqrt( TMath::Sq(vhtestdir.Y()) + TMath::Sq(vhtestdir.Z()) );
          float cvhxd = TMath::Abs( vhtestdir.X() );
          float cvhlambda = TMath::Pi()/2.0;
          if (cvhpd >0) cvhlambda = TMath::ATan(cvhxd/cvhpd);

          if ( TMath::Abs(vhlambda - cvhlambda) > fVecHitMatchLambda )
            {
              if (fPrintLevel > 1)
                {
                  std::cout << "dlambda  failure: " << vhlambda << " " << cvhlambda << std::endl;
                }
              continue;
            }

          if ( vhdir.Dot(vhtestdir) * vhdir.X() * vhtestdir.X() < 0 &&
               TMath::Abs(vhdir.X()) > 0.01 && TMath::Abs(vhtestdir.X()) > 0.01)
            {
              if (fPrintLevel > 1)
                {
                  std::cout << "lambda sign failure" << std::endl;
                }
              continue;
            }

          if (fPrintLevel > 1)
            {
              std::cout << " vh cluster match " << std::endl;
            }
          foundmatch = true;
        }
      if (!foundmatch)
        {
          return false;
        }
      else
        {
          // we have a match, but let's check it to see if we should discard it because there's a
          // close-by mismatch, which happens when a photon converts
          // the above loop stops when we find a match, but this time we want to go through
          // all the VH's in the cluster and check them, so new loop.

          // for (size_t ivh=0; ivh<cluster.size(); ++ivh)
          //   {
          //     gar::rec::VecHit vhtest = vechits[cluster[ivh]];
          //     TVector3 vhtestpos(vhtest.Position());

          //     // look for close-by VH's with an eta mismatch
          //     if ((vhpos-vhtestpos).Mag() < fCloseEtaUnmatch)
          //     {
          //       float eta = calceta2d(vhtest,vh);
          //       if ( eta > fVecHitMatchEta )
          //         {
          //           return false;
          //         }
          //     }
          //   } 

        }
      return true;
    }

    // rough estimate of track parameters -- both ends

    int tpcpatrec2::makepatrectrack(std::vector<gar::rec::TPCCluster> &trackTPCClusters, gar::rec::TrackPar &trackpar)
    {
      // track parameters:  x is the independent variable
      // 0: y
      // 1: z
      // 2: curvature
      // 3: phi
      // 4: lambda = angle from the cathode plane
      // 5: x   /// added on to the end


      float lengthforwards = 0;
      std::vector<int> hlf;
      float lengthbackwards = 0;
      std::vector<int> hlb;

      gar::rec::sort_TPCClusters_along_track(trackTPCClusters,hlf,hlb,fPrintLevel,lengthforwards,lengthbackwards,fSortTransWeight,fSortDistBack);

      std::vector<float> tparbeg(6,0);
      float xother = 0;
      if ( gar::rec::initial_trackpar_estimate(trackTPCClusters, hlf, tparbeg[2], tparbeg[4], 
                                               tparbeg[3], tparbeg[5], tparbeg[0], tparbeg[1], xother, fInitialTPNTPCClusters, fPrintLevel) != 0) 
        {
          return 1;
        }

      std::vector<float> tparend(6,0);
      if ( gar::rec::initial_trackpar_estimate(trackTPCClusters, hlb, tparend[2], tparend[4], 
                                               tparend[3], tparend[5], tparend[0], tparend[1], xother, fInitialTPNTPCClusters, fPrintLevel) != 0)
        {
          return 1;
        }

      // no chisquare or covariance in patrec tracks
      float covmatbeg[25];
      float covmatend[25];
      for (size_t i=0; i<25; ++i) // no covmat in patrec tracks
        {
          covmatend[i] = 0;
          covmatbeg[i] = 0;
        }

      trackpar.setNTPCClusters(trackTPCClusters.size());
      trackpar.setTime(0);
      trackpar.setChisqForwards(0);
      trackpar.setChisqBackwards(0);
      trackpar.setLengthForwards(lengthforwards);
      trackpar.setLengthBackwards(lengthbackwards);
      trackpar.setCovMatBeg(covmatbeg);
      trackpar.setCovMatEnd(covmatend);
      trackpar.setTrackParametersBegin(tparbeg.data());
      trackpar.setXBeg(tparbeg[5]);
      trackpar.setTrackParametersEnd(tparend.data());
      trackpar.setXEnd(tparend[5]);

      return 0;
    }


    // compute a 2D eta

    float tpcpatrec2::calceta2d(gar::rec::VecHit &vhtest, gar::rec::VecHit &vh)
    {
      // normalized direction vector for the VH under test, just the components
      // perpendicular to X

      TVector3 vhtestdir(vhtest.Direction());
      TVector3 vhtestpos(vhtest.Position());
      TVector3 vhdir(vh.Direction());
      TVector3 vhpos(vh.Position());

      TVector3 vhdp(vhdir);
      vhdp.SetX(0);
      float norm = vhdp.Mag();
      if (norm > 0) vhdp *= (1.0/norm);

      // same for the VH in the cluster under test

      TVector3 vhcp(vhtestdir);
      vhcp.SetX(0);
      norm = vhcp.Mag();
      if (norm > 0) vhcp *= (1.0/norm);

      float relsign = 1.0;
      if (vhdp.Dot(vhcp) < 0) relsign = -1;

      TVector3 dcent = vhpos-vhtestpos;
      dcent.SetX(0);

      TVector3 avgdir1 = 0.5*(vhdp + relsign*vhcp);
      float amag = avgdir1.Mag();
      if (amag != 0) avgdir1 *= 1.0/amag;
      float eta = (dcent.Cross(avgdir1)).Mag();
      return eta;

    }

    // test to see if a cluster looks like a conversion and split it in two if it does.  Returns true
    // if the cluster is identified as a conversion.

    bool tpcpatrec2::conversion_test_split(const std::vector<gar::rec::VecHit> &vechits,
                                           std::vector<size_t> &cluster,
                                           std::vector<size_t> &splitclus1,
                                           std::vector<size_t> &splitclus2)
    {

      splitclus1.clear();
      splitclus2.clear();

      // find the VH the farthest away from all the others as an approximation of the end of one of the legs of the
      // conversion candidate.

      double dsummax=0;
      size_t ivhend=0;
      for (size_t ivh=0; ivh<cluster.size(); ++ivh)
        {
          TVector3 vhpos(vechits.at(cluster.at(ivh)).Position());
          double dsum = 0;
          for (size_t jvh=0; jvh<cluster.size(); ++jvh)
            {
              if (ivh == jvh) continue;
              TVector3 vhpos2(vechits.at(cluster.at(jvh)).Position());
              gar::rec::VecHit vh1 = vechits.at(cluster.at(jvh));
              gar::rec::VecHit vh2 = vechits.at(cluster.at(ivh));
              dsum += calceta2d(vh1,vh2); 
              dsum += (vhpos-vhpos2).Mag();
            }
          //std::cout << "Looking for end: " << cluster.at(ivh) << " " << vhpos.X() << " " << vhpos.Y() << " " << vhpos.Z() << " " << dsum << std::endl;
          if (dsum > dsummax)
            {
              ivhend = ivh;  // this is an index into cluster
              dsummax = dsum;
              //std::cout << "This is the new max" << std::endl;
            }
        }
      if (dsummax == 0) return(false);

      // step along the cluster, looking for the closest vh not on the list yet, and identify places where we turn around
      // and backtrack.

      std::vector<size_t> already;  // already looked at vh list
      std::set<size_t> yet;      // set of vector hits yet to look at.

      size_t ivhlast = 0; // the last one looked at -- index into vechits

      for(size_t ivh=0; ivh<cluster.size(); ++ivh)
        {
          if (ivh == ivhend) 
            {
              ivhlast = cluster.at(ivh);
              already.push_back(ivhlast);
              //std::cout << " pushed " << ivhlast << " to already" << std::endl;
            }
          else
            {
              yet.insert(cluster.at(ivh));
            }
        }

      TVector3 lastpdir(0,0,0);  // direction in which we're going.  Don't know yet, and the VH's
      //have a two-fold ambiguity as to what that means

      TVector3 lastpos(vechits.at(ivhlast).Position());  // position of ivhlast
      gar::rec::VecHit lastvhdp = vechits.at(ivhlast);   // save for calculating eta
      
      while(yet.size() > 0)
        {

          // find which VH in the yet-to-be-looked-at pile that is most likely the "next" one.
          // choose the closest vechit that satisfies all the matching criteria. 

          float dmin = 0;
          bool foundmin = false;
          size_t inext=0;
          for (auto iyet : yet)
            {
              TVector3 testpos(vechits.at(iyet).Position());
              TVector3 testdir(vechits.at(iyet).Direction());
              gar::rec::VecHit testvhdp = vechits.at(iyet);
              float dtest = (lastpos - testpos).Mag() + calceta2d(lastvhdp,testvhdp);
              if (dtest < dmin || !foundmin)
                {
                  foundmin = true;
                  dmin = dtest;
                  inext = iyet;
                }
            }
          if (!foundmin) return(false);
          TVector3 nextpos(vechits.at(inext).Position());
          TVector3 nextdir(vechits.at(inext).Direction());
          if (lastpdir.Mag() > 1E-3)
            {
              TVector3 testpdir = nextpos - lastpos;
              if (testpdir.Angle(lastpdir) > 3)
                {
                  // we turned around -- found a conversion.
                  splitclus1 = already;
                  for (auto iyet : yet)
                    {
                      splitclus2.push_back(iyet);
                    }
                  //std::cout << "Found a conversion" << std::endl;
                  return(true);
                }
            }
          lastpdir = nextpos - lastpos;
          already.push_back(inext);
          //std::cout << " pushed " << inext << " to already" << std::endl;
          yet.erase(inext);
          lastpos = nextpos;
          lastvhdp = vechits.at(inext);
          //std::cout << "conv hunting: " << lastpos.X() << " " << lastpos.Y() << " " << lastpos.Z() << std::endl;
        }
      return(false);
    }

    DEFINE_ART_MODULE(tpcpatrec2)

    
  } // namespace rec
} // namespace gar
