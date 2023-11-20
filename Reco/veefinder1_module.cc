////////////////////////////////////////////////////////////////////////
// Class:       veefinder1
// Plugin Type: producer (art v2_11_02)
// File:        veefinder1_module.cc
//
// created 3 May 2020 Thomas Junk based on the vertex finder module
// use Susan Born's analysis cuts
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

#include "ReconstructionDataProducts/Track.h"
#include "ReconstructionDataProducts/Vee.h"
#include "Reco/TrackPar.h"

#include "TMath.h"
#include "TVectorF.h"
#include "TMatrixF.h"

#include <memory>

namespace gar {
  namespace rec {

    class veefinder1 : public art::EDProducer {
    public:
      explicit veefinder1(fhicl::ParameterSet const & p);
      // The compiler-generated destructor is fine for non-base
      // classes without bare pointers or other resource use.

      // Plugins should not be copied or assigned.
      veefinder1(veefinder1 const &) = delete;
      veefinder1(veefinder1 &&) = delete;
      veefinder1 & operator = (veefinder1 const &) = delete;
      veefinder1 & operator = (veefinder1 &&) = delete;

      // Required functions.
      void produce(art::Event & e) override;

    private:

      size_t fNTPCClusPerTrack;     ///< min number of TPC clusters for a track to be identified as part of a vee
      float fPCut;                  ///< min momentum for a track to be identified as part of a vee
      float fRCut;                  ///< DCA cut from track to vee vertex
      float fDCut;                  ///< distance cut from endpoints of tracks
      size_t fRequireNearTrack;     ///< number of nearby tracks to require (primary vertex)
      size_t fNearTrackNTPCClus;    ///< number of nearby tracks' TPC cluster count requirement
      float fNearTrackDMinCut;      ///< minimum distance of nearby track to the vee vertex
      float fNearTrackDMaxCut;      ///< maximum distance of nearby track to the vee vertex
      float fNearTrackAngCut;       ///< pointing requirement of the vee at the nearby track endpoint
      int fPrintLevel;              ///< module print level
      float fOpenAngMinCut;         ///< minimum opening angle cut in 3D
      float fOpenAngMaxCut;         ///< opening angle cut in 3D
      std::string fTrackLabel;      ///< track module label
      bool fRequireOppositeCharge;  ///< whether or not to require opposite-sign tracks.

      // check that the track passes requirements for association with a vee

      bool goodTrack(Track &trk, TrackEnd usebeg);

      // check that a track passes requirements for being another nearby track (primary vertex candidate)

      bool goodOtherTrack(Track &trk, TrackEnd usebeg, TVector3 &veepos, TVector3 &veep);

      // fit two tracks to a single vertex

      int fitVertex(std::vector<TrackPar> &tracks, std::vector<float> &xyz, float &chisquared, std::vector<std::vector<float> > &covmat, double &time,
                    std::vector<TrackEnd> &usebeg, std::vector<float> &dca);

    };


    veefinder1::veefinder1(fhicl::ParameterSet const & p) : EDProducer{p} 
    {
      fTrackLabel            = p.get<std::string>("TrackLabel","track");
      fPCut                  = p.get<float>("PCut",0.1);  // in GeV
      fRCut                  = p.get<float>("RCut", 0.1); // in cm
      fDCut                  = p.get<float>("DCut",20.); // in cm
      fOpenAngMinCut         = p.get<float>("OpenAngMinCut", 0.05); // in radians
      fOpenAngMaxCut         = p.get<float>("OpenAngMaxCut", 3.09); // in radians
      fNTPCClusPerTrack      = p.get<size_t>("NTPCClusPerTrack",100);
      fRequireNearTrack      = p.get<size_t>("RequireNearTrack",0);      
      fNearTrackNTPCClus     = p.get<size_t>("NearTrackNTPCClus",100);     
      fNearTrackDMinCut      = p.get<float>("NearTrackDMinCut",2.0);    
      fNearTrackDMaxCut      = p.get<float>("NearTrackDMaxCut",20.0);    
      fNearTrackAngCut       = p.get<float>("NearTrackAngCut",0.2);
      fRequireOppositeCharge = p.get<bool>("ReqireOppositeCharge",true);
      fPrintLevel = p.get<int>("PrintLevel",0);

      art::InputTag trackTag(fTrackLabel);
      consumes< std::vector<rec::Track> >(trackTag);

      // Produce vertex, associated tracks & ends.  TrackEnd is defined in Track.h
      produces< std::vector<rec::Vee> >();
      produces< art::Assns<rec::Track, rec::Vee, rec::TrackEnd> >();
    }

    void veefinder1::produce(art::Event & e) {
      std::unique_ptr< std::vector<rec::Vee> >
        veeCol(new std::vector<rec::Vee>);

      auto trkVeeAssns = std::make_unique<  art::Assns<rec::Track, rec::Vee, rec::TrackEnd> >();
      //std::unique_ptr< art::Assns<rec::Track, rec::Vee, rec::TrackEnd> >
      // trkVeeAssns(new art::Assns<rec::Track, rec::Vee, rec::TrackEnd>);

      auto trackHandle = e.getValidHandle< std::vector<Track>  >(fTrackLabel);
      auto const& tracks = *trackHandle;

      auto const veePtrMaker   = art::PtrMaker<rec::Vee>(e);
      auto const trackPtrMaker = art::PtrMaker<rec::Track>(e, trackHandle.id());

      size_t nTrack = tracks.size();

      const float m_prot = 0.9382723; // in GeV
      const float m_pi = 0.13957018;  // in GeV

      // For each end of each track, look for candidate vees.  Only two tracks per vee, i and j, or this and that

      if (nTrack>=2) 
        {
          std::vector<TrackPar> tplist;
          std::vector<TrackEnd> usebeg;
          for (size_t iTrack=0; iTrack<nTrack-1; ++iTrack)
            {
              Track thisTrack = tracks.at(iTrack);
              for (int iEnd=0; iEnd<2; ++iEnd)
                {
	          int thisCharge = 0;
                  TrackEnd iTrackEnd(iEnd);   
                  if (!goodTrack(thisTrack,iEnd)) continue;
                  TVector3 thisEndinSpace;
                  TVector3 thisP;
                  if (iTrackEnd==TrackEndBeg) 
                    {
                      TVector3 tmpv(thisTrack.Vertex());
                      thisEndinSpace = tmpv;
                      TVector3 tmpv2(thisTrack.VtxDir());
                      tmpv2 *= thisTrack.Momentum_beg();
                      thisP = tmpv2;
		      thisCharge = thisTrack.ChargeBeg();
                    } 
                  else 
                    {
                      TVector3 tmpv(thisTrack.End());
                      thisEndinSpace = tmpv;
                      TVector3 tmpv2(thisTrack.EndDir());
                      tmpv2 *= thisTrack.Momentum_end();
                      thisP = tmpv2;
		      thisCharge = thisTrack.ChargeEnd();
                    }
                  tplist.clear();
                  tplist.emplace_back(thisTrack);
                  usebeg.clear();
                  usebeg.push_back(iTrackEnd);

                  for (size_t jTrack=iTrack+1; jTrack<nTrack; ++jTrack)
                    {
                      Track thatTrack = tracks.at(jTrack);
                      for (int jEnd = 0; jEnd<2; ++jEnd)
                        {
			  int thatCharge = 0;
                          TrackEnd jTrackEnd(jEnd);
                          if (!goodTrack(thatTrack,jEnd)) continue;
                          TVector3 thatEndinSpace;
                          TVector3 thatP;
                          if (jTrackEnd==TrackEndBeg) 
                            {
                              TVector3 tmpv(thatTrack.Vertex());
                              thatEndinSpace = tmpv;
                              TVector3 tmpv2(thatTrack.VtxDir());
                              tmpv2 *= thatTrack.Momentum_beg();
                              thatP = tmpv2;
			      thatCharge = thatTrack.ChargeBeg();
                            } 
                          else 
                            {
                              TVector3 tmpv(thatTrack.End());
                              thatEndinSpace = tmpv;
                              TVector3 tmpv2(thatTrack.EndDir());
                              tmpv2 *= thatTrack.Momentum_end();
                              thatP = tmpv2;
			      thatCharge = thatTrack.ChargeEnd();
                            }

                          // distance and opening angle cuts

                          if (fPrintLevel>1) std::cout << "veefinder: got a pair of tracks" << std::endl;
                          if ( (thisEndinSpace-thatEndinSpace).Mag() > fDCut ) continue;
                          if (fPrintLevel>1) std::cout << "veefinder: passed the endpoint distance cut" << std::endl;
                          if (fPrintLevel>1) std::cout << "veefinder: opening angle: " << thisP.Angle(thatP) << std::endl;
                          float openang = thisP.Angle(thatP); 
                          if ( openang < fOpenAngMinCut) continue;
                          if ( openang > fOpenAngMaxCut) continue;
                          if (fPrintLevel>1) std::cout << "veefinder: passed the opening angle cut" << std::endl;
			  if (fRequireOppositeCharge)
			    {
			      if (thisCharge == 0 || thatCharge == 0)
				{
				  if (fPrintLevel>1)
				    {
				      std::cout << "veefinder: One of the track charges is zero: " << thisCharge << " " << thatCharge << std::endl;
				    }
				}
			      if (thisCharge == thatCharge) continue;
			      if (fPrintLevel>1) std::cout << "veefinder: passed the opposite-charge cut" << std::endl;
			    }

                          tplist.resize(1);
                          tplist.emplace_back(thatTrack);
                          usebeg.resize(1);
                          usebeg.push_back(jTrackEnd);

                          std::vector<float> dcas;
                          std::vector<float> vpos;
                          std::vector<std::vector<float>> covmat;
                          double time;
                          float chisquared;
                          int iret = fitVertex(tplist,vpos,chisquared,covmat,time,usebeg,dcas);
                          if (fPrintLevel>1) std::cout << "veefinder: vertex fitted" << std::endl;
                          if (iret != 0) continue;
                          if (fPrintLevel>1) std::cout << "veefinder: vertex retcode okay" << std::endl;
                          if (dcas[0] > fRCut || dcas[1] > fRCut) continue;
                          if (fPrintLevel>1) std::cout << "veefinder: tracks passed rcut" << std::endl;

                          TVector3 sumP = thisP + thatP;
                          TVector3 vposv(vpos.data());

                          // see if there are nearby tracks that might be from the primary vertex

                          if (fRequireNearTrack > 0)
                            {
                              size_t otherTrackCount = 0;
                              for (size_t kTrack=0; kTrack<nTrack; ++kTrack)
                                {
                                  if (kTrack == iTrack || kTrack == jTrack) continue;
                                  Track otherTrack = tracks.at(kTrack);
                                  for (int kEnd = 0; kEnd<2; ++kEnd)
                                    {
                                      TrackEnd kTrackEnd(kEnd);
                                      if (!goodOtherTrack(otherTrack,kTrackEnd,vposv,sumP)) continue;
                                      ++otherTrackCount;
                                      break;
                                    }                              
                                }
                              if (otherTrackCount < fRequireNearTrack) continue;
                            }
                          // calculate invariant masses for the three hypotheses

                          std::vector<TLorentzVector> fourmomentumvec;
                          float energy_kshort_hyp = 
                            TMath::Sqrt(thisP.Mag2() + m_pi*m_pi) +
                            TMath::Sqrt(thatP.Mag2() + m_pi*m_pi);

			  // sort the lambda hypotheses so that the first one has a positive proton (lambda)
			  // and the second one has a negative proton (lambdabar) if we are requiring opposite charged tracks
			  // if we have two same-charged tracks, the sorting is arbitrary

                          float energy_lambda1_hyp =
                            TMath::Sqrt(thisP.Mag2() + m_prot*m_prot) +
                            TMath::Sqrt(thatP.Mag2() + m_pi*m_pi);
                          float energy_lambda2_hyp =
                            TMath::Sqrt(thisP.Mag2() + m_pi*m_pi) +
                            TMath::Sqrt(thatP.Mag2() + m_prot*m_prot);
			  if (thisCharge == -1)
			    {
			      float tmp = energy_lambda1_hyp;
			      energy_lambda1_hyp = energy_lambda2_hyp;
			      energy_lambda2_hyp = tmp;
			    }
                          fourmomentumvec.emplace_back(sumP,energy_kshort_hyp);
                          fourmomentumvec.emplace_back(sumP,energy_lambda1_hyp);
                          fourmomentumvec.emplace_back(sumP,energy_lambda2_hyp);
                          float covmatvec[9];
                          size_t icounter=0;
                          for (size_t icm=0; icm<3; ++icm)
                            {
                              for (size_t jcm=0;jcm<3; ++jcm)
                                {
                                  covmatvec[icounter++] = covmat.at(icm).at(jcm);
                                }
                            }
			  if (fPrintLevel > 1) std::cout << "veefinder: Creating a vee object" << std::endl;
                          veeCol->emplace_back(vpos.data(),covmatvec,chisquared,time,fourmomentumvec.data());
                          auto const veeptr = veePtrMaker(veeCol->size()-1);
                          auto const itrackptr = trackPtrMaker( iTrack );
                          auto const jtrackptr = trackPtrMaker( jTrack );
                          trkVeeAssns->addSingle(itrackptr,veeptr,iTrackEnd);
                          trkVeeAssns->addSingle(jtrackptr,veeptr,jTrackEnd);
                        } // end loop over jEnd
                    } // end loop over jTrack
                } // end loop over iEnd
            } // end loop over iTrack

        }  // End test for nTrack>=2; might put empty vertex list & Assns into event
      e.put(std::move(veeCol));
      e.put(std::move(trkVeeAssns));
    } // end produce method.

    //--------------------------------------------------------------------------------------------------

    // fitVertex returns 0 if success
    // It does not cut tracks from the vertex candidate set -- may need to do that in the caller.

    int veefinder1::fitVertex(std::vector<TrackPar> &tracks, 
                              std::vector<float> &xyz, 
                              float &chisquared, 
                              std::vector< std::vector<float> > &covmat, 
                              double &time,
                              std::vector<TrackEnd> &usebeg,
                              std::vector<float> &dca
                              )
    {
      // find the ends of the tracks that are closest to the other tracks' ends
      // pick the end that minimizes the sums of the min(dist) to the other tracks' endpoints

      // todo -- iterate on this as the tracks are helices and not lines -- move to close point on track and re-do linear extrapolation

      if (tracks.size() < 2) return(1);  // need at least two tracks to make a vertex

      // find the average time of the tracks

      time = 0;
      for (size_t i=0; i<tracks.size(); ++i)
        {
          time += tracks.at(i).getTime();
        }
      if (tracks.size() > 0)
        {
          time /= tracks.size();
        }

      xyz.resize(3);
      dca.clear();

      std::vector<TVectorF> dir;
      std::vector<TVectorF> p;

      for (size_t itrack = 0; itrack < tracks.size(); ++itrack)
        {
          TVector3 trackbeg = tracks.at(itrack).getXYZBeg();
          TVector3 trackend = tracks.at(itrack).getXYZEnd();
          float phi=0;
          float si=0;
          if (usebeg.at(itrack)==TrackEndBeg)
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
      bool checkdiagonal = true;
      for(int e=0;e<3;e++) if(A(e,e)==0) checkdiagonal = false;
      if(!checkdiagonal) return (1);
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
          float dcatrack = (diff3.Cross(dir3)).Mag();
          dca.push_back(dcatrack);
          chisquared += dcatrack*dcatrack;   // todo -- include track errors in chisquared calc
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

    bool veefinder1::goodTrack(Track &trk, TrackEnd usebeg)
    {
      if (trk.NHits() < fNTPCClusPerTrack)   // historically-named NHits gets the number of TPC clusters
        { return false; }
      if (usebeg == TrackEndBeg)
        {
          if (trk.Momentum_beg() < fPCut) 
            { return false; }
        }
      else
        {
          if (trk.Momentum_end() < fPCut) 
            { return false; }
        }

      return true;     // if we got here, then the track passed all the cuts
    }

    //--------------------------------------------------------------------------------------------------
    // goodOtherTrack returns true if this track is consistent with being a primary vertex track
    // that is, not part of the vee but close ot it.  There is no requirement on the primary vertex
    // as we might not even really have one -- it may consist of just one track.

    bool veefinder1::goodOtherTrack(Track &trk, TrackEnd usebeg, TVector3 &veepos, TVector3 &veep)
    {
      if (trk.NHits() < fNearTrackNTPCClus)   // historically-named NHits gets the number of TPC clusters
        { return false; }
      TVector3 tep(0,0,0);
      if (usebeg == TrackEndBeg)
        {
          tep.SetXYZ(trk.Vertex()[0],trk.Vertex()[1],trk.Vertex()[2]);
        }
      else
        {
          tep.SetXYZ(trk.End()[0],trk.End()[1],trk.End()[2]);
        }
      TVector3 dv = tep - veepos;
      float dve = dv.Mag();
      if (dve < fNearTrackDMinCut || dve > fNearTrackDMaxCut)
        { return false; }
      if (dv.Angle(veep) > fNearTrackAngCut)
        { return false; }
      return true;
    }


    DEFINE_ART_MODULE(veefinder1)

  } // namespace rec
} // namespace gar
