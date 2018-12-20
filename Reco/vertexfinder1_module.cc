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

      // fit an arbitrary number of tracks to a single vertex
      // returns 0 if success, 1 if failure
      // covmat is 3x3

      int fitVertex(std::vector<TrackPar> &tracks, std::vector<float> &xyz, float &chisquared, std::vector<std::vector<float> > &covmat, ULong64_t &time);

    };


    vertexfinder1::vertexfinder1(fhicl::ParameterSet const & p)
    // :
    {
      produces< std::vector<rec::Vertex> >();
      produces< art::Assns<rec::Track, rec::Vertex> >();

      fTrackLabel = p.get<std::string>("TrackLabel","track");
      fChisquaredCut = p.get<float>("ChisquaredPerDOFCut",5.0);
      fRCut = p.get<float>("RCut",10.0); // in cm
      fPrintLevel = p.get<int>("PrintLevel",0);
    }

    void vertexfinder1::produce(art::Event & e)
    {
      std::unique_ptr< std::vector<rec::Vertex> > vtxCol(new std::vector<rec::Vertex>);
      std::unique_ptr< art::Assns<rec::Track,rec::Vertex> > trkVtxAssns(new ::art::Assns<rec::Track,rec::Vertex>);

      auto trackHandle = e.getValidHandle< std::vector<Track>  >(fTrackLabel);
      auto const& tracks = *trackHandle;

      auto const vtxPtrMaker = art::PtrMaker<rec::Vertex>(e, *this);
      auto const trackPtrMaker = art::PtrMaker<rec::Track>(e, trackHandle.id());

      // copy track info into TrackPar objects

      std::vector<TrackPar> tpars;
      for (size_t i=0; i<tracks.size(); ++i)
	{
	  tpars.emplace_back(tracks.at(i));
	}

      // find groups of track endpoints within rcut of each other
      std::vector< std::vector<float> > pairdist[4];  // bb, be, eb, ee
      for (size_t itrack=0; itrack < tpars.size(); ++itrack)
	{
	  TVector3 trackbeg = tpars.at(itrack).getXYZBeg();
	  TVector3 trackend = tpars.at(itrack).getXYZEnd();
	  std::vector<float> dvec[4];

	  for (size_t jtrack=0; jtrack < tpars.size(); ++jtrack)
	    {
	      if (jtrack == itrack)
		{
		  for (size_t i=0; i<4; ++i) dvec[i].push_back(0);
		}
	      else
		{
	          TVector3 otbeg = tpars.at(jtrack).getXYZBeg();
	          TVector3 otend = tpars.at(jtrack).getXYZEnd();
		  dvec[0].push_back((trackbeg-otbeg).Mag());
		  dvec[1].push_back((trackbeg-otend).Mag());
		  dvec[2].push_back((trackend-otbeg).Mag());
		  dvec[3].push_back((trackend-otend).Mag());
		}
	    }
	  for (size_t i=0; i<4; ++i)
	    {
	      pairdist[i].push_back(dvec[i]);
	    }
	}

      for (size_t itrack=0; itrack < tpars.size(); ++itrack)
	{
	  size_t sizelast = 0;
	  do
	    {
	      std::vector<size_t> vtlist;
	      vtlist.push_back(itrack);
	      for (size_t jtrack=itrack+1; jtrack < tpars.size(); ++jtrack)
		{
		  for (size_t i=0; i<4; ++i)
		    {
		      float d = pairdist[i][itrack][jtrack];
		      if (d > 0 && d < fRCut)
			{
			  vtlist.push_back(jtrack);
			  pairdist[i][itrack][jtrack] *= -1;
			  break;
			}
		    }
		}
	      sizelast = vtlist.size();
	      if (sizelast>1)
		{
		  std::vector<TrackPar> vtracks;
		  for (size_t i=0; i<vtlist.size(); ++i)
		    {
		      vtracks.push_back(tpars.at(vtlist.at(i)));
		    }
		  std::vector<float> xyz(3,0);
		  std::vector< std::vector<float> > covmat;
		  float chisquared=0;
		  ULong64_t time;
		  if (fitVertex(vtracks,xyz,chisquared,covmat,time) == 0)
		    {
		      float cmv[9];
		      size_t icounter=0;
		      for (size_t i=0; i<3;++i)
			{
			  for (size_t j=0; j<3; ++j)
			    {
			      cmv[icounter] = covmat.at(i).at(j);
			    }
			}
		      
		      vtxCol->emplace_back(xyz.data(),cmv,time);  
		      auto const vtxpointer = vtxPtrMaker(vtxCol->size()-1);
		      for (size_t i=0; i<vtlist.size(); ++i)
			{
			  auto const trackpointer = trackPtrMaker(vtlist.at(i));
			  trkVtxAssns->addSingle(trackpointer,vtxpointer);
			}
		    }
		}
	    }
	  while (sizelast>1); // if we found a vertex candidate for this track, we may find another for the other end
	}
      e.put(std::move(vtxCol));
      e.put(std::move(trkVtxAssns));
    }

    //--------------------------------------------------------------------------------------------------

    // fitVertex returns 0 if success
    // It does not cut tracks from the vertex candidate set -- may need to do that in the caller.

      int vertexfinder1::fitVertex(std::vector<TrackPar> &tracks, 
				   std::vector<float> &xyz, 
				   float &chisquared, 
				   std::vector< std::vector<float> > &covmat, 
				   ULong64_t &time)
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

      std::vector<bool> usebeg;
      std::vector<TVectorF> dir;
      std::vector<TVectorF> p;

      for (size_t itrack = 0; itrack < tracks.size(); ++itrack)
	{
	  float dsum_beg = 0;
	  float dsum_end = 0;
	  TVector3 trackbeg = tracks.at(itrack).getXYZBeg();
	  TVector3 trackend = tracks.at(itrack).getXYZEnd();
	  for (size_t jtrack = 0; jtrack < tracks.size(); ++jtrack)
	    {
	      if (jtrack == itrack) continue;
	      TVector3 otbeg = tracks.at(jtrack).getXYZBeg();
	      TVector3 otend = tracks.at(jtrack).getXYZEnd();
	      dsum_beg += TMath::Min( (trackbeg-otbeg).Mag(), (trackbeg-otend).Mag() );
	      dsum_end += TMath::Min( (trackend-otbeg).Mag(), (trackend-otend).Mag() );
	    }
	  float phi=0;
	  float s=0;
	  if (dsum_beg <= dsum_end)
	    {
	      usebeg.push_back(true);
	      float tmppos[3] = {(float) trackbeg.X(), (float) trackbeg.Y(), (float) trackbeg.Z()};
	      p.emplace_back(3,tmppos);
	      phi = tracks.at(itrack).getTrackParametersBegin()[3];
	      s =   TMath::Tan(tracks.at(itrack).getTrackParametersBegin()[4]);
	      if (s == 0)
		{
		  s = 999.;
		}
	      else
		{
		  s = 1/s; 
		}
	    }
	  else
	    {
	      usebeg.push_back(false);
	      float tmppos[3] = {(float) trackend.X(), (float) trackend.Y(), (float) trackend.Z()};
	      p.emplace_back(3,tmppos);
	      phi = tracks.at(itrack).getTrackParametersEnd()[3];
	      s =   TMath::Tan(tracks.at(itrack).getTrackParametersEnd()[4]);
	      if (s == 0)
		{
		  s = 999.;
		}
	      else
		{
		  s = 1/s; 
		}
	    }
	  if (s == 0) s=1E-6;
	  TVectorF dtmp(3);
	  dtmp[0] = 1.0/s;
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
	  chisquared +=  (diff3.Cross(dir3)).Mag2();
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

    DEFINE_ART_MODULE(vertexfinder1)

  } // namespace rec
} // namespace gar
