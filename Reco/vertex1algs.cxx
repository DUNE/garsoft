#include<Reco/vertex1algs.h>

//--------------------------------------------------------------------------------------------------

// fitVertex returns 0 if success
// It does not cut tracks from the vertex candidate set -- may need to do that in the caller.

int gar::rec::fitVertex(std::vector<TrackPar> &tracks, 
                        std::vector<float> &xyz, 
                        float &chisquared, 
                        std::vector< std::vector<float> > &covmat, 
                        double &time,
                        std::vector<TrackEnd> usebeg,
                        std::vector<float> &doca
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
  covmat.resize(9);

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
  try {
      A.Invert(&det);
  } catch (...) {
      // Uninvertible matrix some other how
      return (1);
  }
  if (det == 0) return(1);
  TVectorF xyzsol = A*vsum;
  xyz.at(0) = xyzsol[0];
  xyz.at(1) = xyzsol[1];
  xyz.at(2) = xyzsol[2];

  chisquared = 0;
  //std::cout << "in vertexfit ntracks: " << tracks.size() << std::endl;
  for (size_t itrack=0; itrack < tracks.size(); ++itrack)
    {
      
      TVector3 dir3(dir.at(itrack)[0],dir.at(itrack)[1],dir.at(itrack)[2]);
      TVectorF diff = xyzsol - p.at(itrack);
      TVector3 diff3(diff[0],diff[1],diff[2]);
      double dctmp = (diff3.Cross(dir3)).Mag2();
      doca.push_back(dctmp);
      chisquared +=  dctmp*dctmp;   // todo -- include track errors in chisquared calc
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
