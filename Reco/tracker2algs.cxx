#include "Reco/tracker2algs.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// these two methods are also duplicated in patrec.  Maybe we can use the patrec tracks as the
// initial guess of the track parameters.  Make sure we get them in the right order

int gar::rec::initial_trackpar_estimate(const std::vector<gar::rec::TPCCluster>& TPCClusters,
                                        std::vector<int>& TPCClusterlist,
                                        float& curvature_init,
                                        float& lambda_init,
                                        float& phi_init,
                                        float& xpos,
                                        float& ypos,
                                        float& zpos,
                                        float& x_other_end,
                                        unsigned int initialtpnTPCClusters,
                                        int printlevel)
{

  size_t nTPCClusters = TPCClusterlist.size();
  size_t firstTPCCluster = 0;
  size_t farTPCCluster = TMath::Min(nTPCClusters - 1, (size_t)initialtpnTPCClusters);
  ;
  size_t intTPCCluster = farTPCCluster / 2;
  size_t lastTPCCluster = nTPCClusters - 1;

  float trackbeg[3] = {TPCClusters[TPCClusterlist[firstTPCCluster]].Position()[0],
                       TPCClusters[TPCClusterlist[firstTPCCluster]].Position()[1],
                       TPCClusters[TPCClusterlist[firstTPCCluster]].Position()[2]};

  float tp1[3] = {TPCClusters[TPCClusterlist[intTPCCluster]].Position()[0],
                  TPCClusters[TPCClusterlist[intTPCCluster]].Position()[1],
                  TPCClusters[TPCClusterlist[intTPCCluster]].Position()[2]};

  float tp2[3] = {TPCClusters[TPCClusterlist[farTPCCluster]].Position()[0],
                  TPCClusters[TPCClusterlist[farTPCCluster]].Position()[1],
                  TPCClusters[TPCClusterlist[farTPCCluster]].Position()[2]};

  if (printlevel > 1) {
    std::cout << "TPCCluster Dump in initial_trackpar_estimate: " << std::endl;
    for (size_t i = 0; i < nTPCClusters; ++i) {
      size_t ihf = i;
      std::cout << i << " : " << TPCClusters[TPCClusterlist[ihf]].Position()[0] << " "
                << TPCClusters[TPCClusterlist[ihf]].Position()[1] << " "
                << TPCClusters[TPCClusterlist[ihf]].Position()[2] << std::endl;
    }
  }
  if (printlevel > 0) {
    std::cout << "first TPCCluster: " << firstTPCCluster << ", inter TPCCluster: " << intTPCCluster
              << " "
              << " far TPCCluster: " << farTPCCluster << std::endl;
    std::cout << "in the TPCCluster list: " << TPCClusterlist[firstTPCCluster] << " "
              << TPCClusterlist[intTPCCluster] << " " << TPCClusterlist[farTPCCluster] << std::endl;
    std::cout << "First TPCCluster x, y, z: " << trackbeg[0] << " " << trackbeg[1] << " "
              << trackbeg[2] << std::endl;
    std::cout << "Inter TPCCluster x, y, z: " << tp1[0] << " " << tp1[1] << " " << tp1[2]
              << std::endl;
    std::cout << "Far   TPCCluster x, y, z: " << tp2[0] << " " << tp2[1] << " " << tp2[2]
              << std::endl;
  }

  xpos = trackbeg[0];
  ypos = trackbeg[1];
  zpos = trackbeg[2];
  x_other_end = TPCClusters[TPCClusterlist[lastTPCCluster]].Position()[0];

  float ycc = 0;
  float zcc = 0;
  curvature_init =
    gar::rec::capprox(trackbeg[1], trackbeg[2], tp1[1], tp1[2], tp2[1], tp2[2], ycc, zcc);
  //std::cout << " inputs to trackpar circ fit (y,z): " << trackbeg[1] << " " << trackbeg[2] << " : "
  //        << tp1[1] << " " << tp1[2] << " : " << tp2[1] << " " << tp2[2] << std::endl;
  //std::cout << "curvature output: " << curvature_init << std::endl;

  // redo calc with a circle fit to all points, not just three

  double dycc = 0;
  double dzcc = 0;
  double drcc = 0;
  double dchisq = 0;
  int ifail = 0;
  std::vector<double> dtpcclusy;
  std::vector<double> dtpcclusz;
  for (size_t i = 0; i < nTPCClusters; ++i) {
    dtpcclusy.push_back(TPCClusters[TPCClusterlist[i]].Position()[1]);
    dtpcclusz.push_back(TPCClusters[TPCClusterlist[i]].Position()[2]);
  }
  int ict = nTPCClusters;
  ouchef(dtpcclusy.data(), dtpcclusz.data(), ict, dycc, dzcc, drcc, dchisq, ifail);
  if (ifail == 0 && drcc != 0) {
    ycc = dycc;
    zcc = dzcc;
    if (curvature_init < 0) drcc = -std::abs(drcc);
    //if (ict > 3)
    //  {
    //    std::cout << "updating curvature " << ict << " " << curvature_init << " with " << 1.0/drcc << std::endl;
    //  }
    curvature_init = 1.0 / drcc;
  }

  phi_init = TMath::ATan2(trackbeg[2] - zcc, ycc - trackbeg[1]);
  float phi2 = phi_init;
  if (curvature_init < 0) phi_init += TMath::Pi();
  float radius_init = 10000;
  if (curvature_init != 0) radius_init = 1.0 / curvature_init;

  float dx1 = tp2[0] - xpos;
  if (dx1 != 0) {
    float dphi2 = TMath::ATan2(tp2[2] - zcc, ycc - tp2[1]) - phi2;
    if (dphi2 > TMath::Pi()) dphi2 -= 2.0 * TMath::Pi();
    if (dphi2 < -TMath::Pi()) dphi2 += 2.0 * TMath::Pi();
    lambda_init = TMath::ATan(1.0 / ((radius_init / dx1) * dphi2));
  }
  else {
    //std::cout << "initial track par estimate failure" << std::endl;
    lambda_init = 0;
    return 1;
  } // got fMinNumTPCClusters all at exactly the same value of x (they were sorted).  Reject track.

  if (printlevel > 0) {
    std::cout << "phi calc: dz, dy " << tp2[2] - trackbeg[2] << " " << tp2[1] - trackbeg[1]
              << std::endl;
    std::cout << "initial curvature, phi, lambda: " << curvature_init << " " << phi_init << " "
              << lambda_init << std::endl;
  }
  return 0;
}

float gar::rec::capprox(float x1,
                        float y1,
                        float x2,
                        float y2,
                        float x3,
                        float y3,
                        float& xc,
                        float& yc)
{
  //-----------------------------------------------------------------
  // Initial approximation of the track curvature -- copied from ALICE
  // here x is y and y is z for us
  //-----------------------------------------------------------------
  x3 -= x1;
  x2 -= x1;
  y3 -= y1;
  y2 -= y1;
  //
  float det = x3 * y2 - x2 * y3;
  if (TMath::Abs(det) < 1e-10) { return 100; }
  //
  float u = 0.5 * (x2 * (x2 - x3) + y2 * (y2 - y3)) / det;
  float x0 = x3 * 0.5 - y3 * u;
  float y0 = y3 * 0.5 + x3 * u;
  float c2 = 1 / TMath::Sqrt(x0 * x0 + y0 * y0);
  xc = x0 + x1;
  yc = y0 + y1;
  if (det < 0) c2 *= -1;
  return c2;
}

// given a set of TPC Clusters, find sort orders hlf and hlb, which are vectors of indices into TPC Clusters,
// that run along the track as best as we can.  Which way is "forwards" and which "backwards" is arbitrary

// sorttransweight is the fraction of additional distance given to the transverse displacement from the
// local line candidate in addtion to the longitudinal component of the distance used to pick the closest next hit
// sortdisback is the turnover point for negative longitudinal distances so we can go back and pick up mis-sorted hits.

void gar::rec::sort_TPCClusters_along_track(const std::vector<gar::rec::TPCCluster>& TPCClusters,
                                            std::vector<int>& hlf,
                                            std::vector<int>& hlb,
                                            int printlevel,
                                            float& lengthforwards,
                                            float& lengthbackwards,
                                            float /* sorttransweight */,
                                            float /* sortdistback */)
{

  // this sorting code appears in the patrec stage too, if only to make tracks that can be drawn
  // on the event display

  // find candidate endpoints and sort TPCClusters

  float cmin[3]; // min x, y, and z coordinates over all TPCClusters
  float cmax[3]; // max x, y, and z coordinates over all TPCClusters
  size_t ihex
    [6]; // index of TPCCluster which gave the min or max ("extreme") 0-2: (min xyz)  3-5 (max xyz)

  for (size_t iTPCCluster = 0; iTPCCluster < TPCClusters.size(); ++iTPCCluster) {
    for (int i = 0; i < 3; ++i) {
      float c = TPCClusters[iTPCCluster].Position()[i];
      if (iTPCCluster == 0) {
        cmin[i] = c;
        cmax[i] = c;
        ihex[i] = 0;
        ihex[i + 3] = 0;
      }
      else {
        if (c < cmin[i]) {
          cmin[i] = c;
          ihex[i] = iTPCCluster;
        }
        if (c > cmax[i]) {
          cmax[i] = c;
          ihex[i + 3] = iTPCCluster;
        }
      }
    }
  }
  // now we have six TPCClusters that have the min and max x, y, and z values.  Find out which of these six
  // TPCClusters has the biggest sum of distances to all the other TPCClusters (the most extreme)
  float sumdmax = 0;
  size_t imax = 0;
  for (size_t i = 0; i < 6; ++i) {
    float sumd = 0;
    TVector3 poshc(TPCClusters[ihex[i]].Position());
    for (size_t iTPCCluster = 0; iTPCCluster < TPCClusters.size(); ++iTPCCluster) {
      TVector3 hp(TPCClusters[iTPCCluster].Position());
      sumd += (poshc - hp).Mag();
    }
    if (sumd > sumdmax) {
      sumdmax = sumd;
      imax = i;
    }
  }

  //  Use this TPCCluster, indexed by ihex[imax] as a starting point
  //  For cases with less than dmindir length, sort TPCClusters in order of how
  //  far they are from the first TPCCluster.  If longer, calculate a curdir vector
  //  and add hits that are a minimum step along curdir

  hlf.clear();
  lengthforwards = 0;
  hlf.push_back(ihex[imax]);
  TVector3 lpos(TPCClusters[hlf[0]].Position());
  TVector3 lastpos = lpos;
  TVector3 cdpos = lpos; // position of the beginning of the direction vector
  size_t cdposindex = 0; // and its index in the hlf vector
  TVector3 hpadd;
  TVector3 curdir(0, 0, 0);
  bool havecurdir = false;
  float dmindir = 10; // promote to a fcl parameter someday.  Distance over which to compute curdir

  for (size_t inh = 1; inh < TPCClusters.size(); ++inh) {
    float dmin = 0;
    int jmin = -1;
    for (size_t jh = 0; jh < TPCClusters.size(); ++jh) {
      bool found = false;
      for (size_t kh = 0; kh < hlf.size(); ++kh) {
        if (hlf[kh] == (int)jh) {
          found = true;
          break;
        }
      }
      if (found) continue; // skip if we've already assigned this TPCCluster on this track
      TVector3 hpos(TPCClusters[jh].Position());

      float d = 0;
      if (havecurdir) {
        float dtmp = (hpos - lastpos).Dot(curdir);
        if (dtmp < -2) { dtmp = -4 - dtmp; }
        d = dtmp + 0.1 * ((hpos - lastpos).Cross(curdir)).Mag();
      }
      else {
        d = (hpos - lpos).Mag();
      }
      if (jmin == -1) {
        jmin = jh;
        dmin = d;
        hpadd = hpos;
      }
      else {
        if (d < dmin) {
          jmin = jh;
          dmin = d;
          hpadd = hpos;
        }
      }
    }
    //  std::cout << "dmin: " << dmin << std::endl;
    hlf.push_back(jmin);
    lengthforwards += (hpadd - lastpos).Mag(); // add up track length
    lastpos = hpadd;
    if ((hpadd - cdpos).Mag() > dmindir) {
      TVector3 cdcand(TPCClusters[hlf[cdposindex]].Position());
      size_t itctmp = cdposindex;
      for (size_t itc = cdposindex + 1; itc < hlf.size(); ++itc) {
        TVector3 cdcandt(TPCClusters[hlf[itc]].Position());
        if ((hpadd - cdcandt).Mag() < dmindir) { break; }
        else {
          itctmp = itc;
          cdcand = cdcandt;
        }
      }
      //std::cout << "Updated curdir: " << cdposindex << " " << itctmp << " " << curdir.X() << " " << curdir.Y() << " " << curdir.Z() << " ";
      cdposindex = itctmp;
      cdpos = cdcand;
      curdir = hpadd - cdcand;
      curdir *= (1.0 / curdir.Mag());
      //std::cout << "To: " << curdir.X() << " " << curdir.Y() << " " << curdir.Z() << std::endl;
      havecurdir = true;
    }
  }

  if (printlevel > 2) {
    for (size_t iTPCCluster = 0; iTPCCluster < TPCClusters.size(); ++iTPCCluster) {
      printf("Sort compare: %5d %10.3f %10.3f %10.3f  %5d %10.3f %10.3f %10.3f\n",
             (int)iTPCCluster,
             TPCClusters[iTPCCluster].Position()[0],
             TPCClusters[iTPCCluster].Position()[1],
             TPCClusters[iTPCCluster].Position()[2],
             hlf[iTPCCluster],
             TPCClusters[hlf[iTPCCluster]].Position()[0],
             TPCClusters[hlf[iTPCCluster]].Position()[1],
             TPCClusters[hlf[iTPCCluster]].Position()[2]);
    }
  }

  // now go backwards -- just invert the sort order

  hlb.clear();
  for (size_t i = 0; i < hlf.size(); ++i) {
    hlb.push_back(hlf[hlf.size() - 1 - i]); // just invert the order for the backward sort
  }
  lengthbackwards = lengthforwards;

  if (printlevel > 2) {
    for (size_t iTPCCluster = 0; iTPCCluster < TPCClusters.size(); ++iTPCCluster) {
      printf("Backward Sort compare: %5d %10.3f %10.3f %10.3f  %5d %10.3f %10.3f %10.3f\n",
             (int)iTPCCluster,
             TPCClusters[iTPCCluster].Position()[0],
             TPCClusters[iTPCCluster].Position()[1],
             TPCClusters[iTPCCluster].Position()[2],
             hlb[iTPCCluster],
             TPCClusters[hlb[iTPCCluster]].Position()[0],
             TPCClusters[hlb[iTPCCluster]].Position()[1],
             TPCClusters[hlb[iTPCCluster]].Position()[2]);
    }
  }
}

// check if we add a link between i and j, do we induce a cycle?

bool gar::rec::sort2_check_cyclic(std::vector<int>& link1, std::vector<int>& link2, int i, int j)
{
  if (i == j) return true;

  int prev = i;
  int next = (link1.at(i) == -1) ? link2.at(i) : link1.at(i);
  while (next != -1) {
    if (next == j) return true;
    int ptmp = next;
    next = (link1.at(next) == prev) ? link2.at(next) : link1.at(next);
    prev = ptmp;
  }
  return false;
}

// try another method to sort TPC clusters -- a greedy algorithm that associates nearby TPC clusters first
// stranded hits may be not included if the best path skirts around them.  hlf and hlb might have a shorter
// length than TPCClusters.

void gar::rec::sort_TPCClusters_along_track2(const std::vector<gar::rec::TPCCluster>& TPCClusters,
                                             std::vector<int>& hlf,
                                             std::vector<int>& hlb,
                                             int /* printlevel */,
                                             float& lengthforwards,
                                             float& lengthbackwards,
                                             float dcut)
{
  size_t nclus = TPCClusters.size();
  hlf.clear();
  hlb.clear();
  lengthforwards = 0;
  lengthbackwards = 0;
  if (nclus == 0) return;
  if (nclus == 1) {
    hlf.push_back(0);
    hlb.push_back(0);
    return;
  }

  size_t ndists = nclus * (nclus - 1) / 2;
  std::vector<float> dists(ndists);
  std::vector<int> dsi(ndists);
  std::vector<int> link1(nclus, -1);
  std::vector<int> link2(nclus, -1);
  for (size_t i = 0; i < nclus - 1; ++i) {
    TVector3 iclus(TPCClusters.at(i).Position());
    for (size_t j = i + 1; j < nclus; ++j) {
      TVector3 jclus(TPCClusters.at(j).Position());
      float d = (iclus - jclus).Mag();
      dists.at(ij2idxsort2(nclus, i, j)) = d;
    }
  }
  TMath::Sort((int)ndists, dists.data(), dsi.data(), false);
  //std::cout << "Dists sorting: " << std::endl;
  //for (size_t i=0; i<ndists; ++i)
  //  {
  //    std::cout << i << " " << dists[i] << " " << dists[dsi[i]] << std::endl;
  //  }

  // greedy pairwise association -- may end up leaving some hits stranded

  for (size_t k = 0; k < ndists; ++k) {
    size_t i = 0;
    size_t j = 0;
    size_t k2 = dsi.at(k);
    if (dists.at(k2) > dcut) break;
    idx2ijsort2(nclus, k2, i, j);
    TVector3 iclus(TPCClusters.at(i).Position());
    TVector3 jclus(TPCClusters.at(j).Position());

    if (link1.at(i) == -1) {
      if (link1.at(j) == -1) {
        link1.at(i) = j;
        link1.at(j) = i;
      }
      else // already have a link on j'th cluster, add a second one only if it's on the other side
      {
        if (link2.at(j) == -1) {
          TVector3 j2clus(TPCClusters.at(link1.at(j)).Position());

          // old test on angle instead of cyclicness
          //if ( (j2clus-jclus).Angle(iclus-jclus) > TMath::Pi()/2)

          if (!sort2_check_cyclic(link1, link2, j, i)) {
            link1.at(i) = j;
            link2.at(j) = i;
          }
        }
      }
    }
    else {
      if (link2.at(i) == -1) {
        if (link1.at(j) == -1) {
          link2.at(i) = j;
          link1.at(j) = i;
        }
        else // already have a link on j'th cluster, add a second one only if it's on the other side
        {
          if (link2.at(j) == -1) {
            TVector3 j2clus(TPCClusters.at(link1.at(j)).Position());

            // old test on angle instead of cyclicness
            //if ( (j2clus-jclus).Angle(iclus-jclus) > TMath::Pi()/2)

            if (!sort2_check_cyclic(link1, link2, j, i)) {
              link2.at(i) = j;
              link2.at(j) = i;
            }
          }
        }
      }
    }
  }

  //std::cout << "sorting clusters and links" << std::endl;
  //for (size_t i=0; i<nclus; ++i)
  //  {
  //    std::cout << i << " " << link1.at(i) << " " << link2.at(i) << " " <<
  //	TPCClusters.at(i).Position()[0] << " " <<
  // 	TPCClusters.at(i).Position()[1] << " " <<
  //	TPCClusters.at(i).Position()[2] << std::endl;
  //}

  std::vector<int> used(nclus, -1);
  std::vector<std::vector<int>> clistv;

  // find the first tpc cluster with just one link

  int ifirst = -1;
  for (size_t i = 0; i < nclus; ++i) {
    int il1 = link1.at(i);
    int il2 = link2.at(i);

    if ((il1 == -1 && il2 != -1) || (il1 != -1 && il2 == -1)) {
      ifirst = i;
      break;
    }
  }
  // this can happen if all TPC clusters are singletons -- no links anywhere.  Happens if
  // there are no TPC clusters at all or just one.
  if (ifirst == -1) { return; }

  // follow the chains

  for (size_t i = ifirst; i < nclus; ++i) {
    if (used.at(i) == 1) continue;
    std::vector<int> clist;
    clist.push_back(i);
    used.at(i) = 1;

    int idx = i;
    int il1 = link1.at(idx);
    int il2 = link2.at(idx);

    int u1 = -2;
    if (il1 != -1) u1 = used.at(il1);
    int u2 = -2;
    if (il2 != -1) u2 = used.at(il2);

    while (u1 == -1 || u2 == -1) {
      idx = (u1 == -1) ? il1 : il2;
      clist.push_back(idx);
      used.at(idx) = 1;
      il1 = link1.at(idx);
      il2 = link2.at(idx);
      u1 = -2;
      if (il1 != -1) u1 = used.at(il1);
      u2 = -2;
      if (il2 != -1) u2 = used.at(il2);
    }
    clistv.push_back(clist);
  }

  // temporary solution -- think about what to do.  Report the longest chain

  size_t msize = 0;
  int ibest = -1;
  for (size_t ichain = 0; ichain < clistv.size(); ++ichain) {
    size_t cursize = clistv.at(ichain).size();
    if (cursize > msize) {
      ibest = ichain;
      msize = cursize;
    }
  }
  if (ibest == -1) return;

  for (size_t i = 0; i < clistv.at(ibest).size(); ++i) {
    hlf.push_back(clistv[ibest].at(i));
    if (i > 0) {
      TVector3 lastpoint(TPCClusters.at(hlf.at(i - 1)).Position());
      TVector3 curpoint(TPCClusters.at(hlf.at(i)).Position());
      lengthforwards += (curpoint - lastpoint).Mag();
    }
  }
  for (size_t i = 0; i < hlf.size(); ++i) {
    hlb.push_back(hlf[hlf.size() - 1 - i]); // just invert the order for the backward sort
  }
  lengthbackwards = lengthforwards;
}

// this method and the next come from formulas on
//  https://stackoverflow.com/questions/27086195/linear-index-upper-triangular-matrix

size_t gar::rec::ij2idxsort2(size_t n, size_t i_in, size_t j_in)
{
  if (i_in >= n) { throw cet::exception("gar::rec::ij2idxsort2 i_in >= n"); }
  if (j_in >= n) { throw cet::exception("gar::rec::ij2idxsort2 j_in >= n"); }
  size_t i = TMath::Min(i_in, j_in);
  size_t j = TMath::Max(i_in, j_in);
  size_t k = (n * (n - 1) / 2) - (n - i) * ((n - i) - 1) / 2 + j - i - 1;
  return k;
}

// j > i

void gar::rec::idx2ijsort2(size_t n, size_t k, size_t& i, size_t& j)
{
  if (n < 2) {
    i = 0;
    j = 0;
    return;
  }
  if (k > n * (n - 1) / 2) {
    throw cet::exception("gar::rec::idx2ijsort2: k too big. ") << k << " " << n;
  }
  i = n - 2 - TMath::Floor(TMath::Sqrt((double)(-8 * k + 4 * n * (n - 1) - 7)) / 2.0 - 0.5);
  j = k + i + 1 - n * (n - 1) / 2 + (n - i) * ((n - i) - 1) / 2;
}

void gar::rec::ouchef(double* x,
                      double* y,
                      int np,
                      double& xcirccent,
                      double& ycirccent,
                      double& rcirc,
                      double& chisq,
                      int& ifail)
{

  // converted from the original Fortran by T. Junk and de-OPALified -- apologies for the gotos

  //      SUBROUTINE OUCHEF(X,Y,NP,XCIRCCENT,YCIRCCENT,RCIRC,CHISQ,IFAIL)
  // *
  // *...OUCHEF  circle fit by linear regression method
  // *.
  // *.            OUCHEF(X,Y,NP,XCIRCCENT,YCIRCCENT,RCIRC,CHISQ,IFAIL)
  // *.
  // *.  INPUT
  // *.       X(.)      x coordinates of points
  // *.       Y(.)      y coordinates of points
  // *.       NP        number of points
  // *.
  // *.  OUTPUT
  // *.       XMED, YMED   center of circle
  // *.       R            radius of circle
  // *.       CHISQ     residual sum of squares : to be a chi**2 ,it should
  // *.                 be divided by (NP-3)*(mean residual)**2
  // *.
  // *.       IFAIL     error flag : 0   OK
  // *.                              1   less than 3 points
  // *.                              2   no convergence
  // *.                              5   number of points exceeding
  // *.                                  arrays dimensions -> fit performed
  // *.                                  with allowed maximum (NPTMAX)
  // *.
  // *. AUTHOR :  N.Chernov, G.Ososkov
  // *.
  // *.     written by N. Chernov and G. Ososkov, Dubna,
  // *.     reference:  computer physics communications vol 33, p329
  // *.     adapted slightly by F. James, March, 1986
  // *.     further adapted slightly by M.Hansroul May,1989
  // *.
  // *. CALLED   : any
  // *.
  // *. CREATED  : 31 May 1989
  // *. LAST MOD : 11 January 1996
  // *.
  // *. Modification Log.:
  // *. 11-jan-96  M.Hansroul  preset chisq to 9999.
  // *. 12-oct-95  M.Hansroul  clean up
  // *. 12-jul-93  M.Hansroul  protection against huge d0
  // *.                        put OUATG2 in line
  // *. 31-May-89  M.Hansroul  use OPAL output parameters
  // *.*********************************************************************

  const double eps = 1.0e-10;
  const int itmax = 20;
  const int nptmax = 1000;
  //const double twopi=8.0*atan(1.0);
  //const double piby2=2.0*atan(1.0);

  int nptlim = 0;
  int npf = 0;
  int n = 0;
  int ihalf = 0;
  int i = 0;
  int iter = 0;
  int iswap = 0;
  //double sfi0=0;
  //double cfi0=0;
  double xmid = 0;
  double ymid = 0;
  double* xf;
  double* yf;
  double vlf[2];
  double vmf[2];
  //double vlm[2];

  double alf = 0;
  double alm = 0;
  double a0 = 0;
  double a1 = 0;
  double a2 = 0;
  double a22 = 0;
  double bem = 0;
  double bet = 0;
  double cur = 0;
  double cu2 = 0;
  double dd = 0;
  double den = 0;
  double det = 0;
  double dy = 0;
  double d2 = 0;
  double f = 0;
  double fact = 0;
  double fg = 0;
  double f1 = 0;
  double g = 0;
  double gam = 0;
  double gam0 = 0;
  double gmm = 0;
  double g1 = 0;
  double h = 0;
  double h2 = 0;
  double p2 = 0;
  double q2 = 0;
  double rm = 0;
  double rn = 0;
  double xa = 0;
  double xb = 0;
  double xd = 0;
  double xi = 0;
  double xm = 0;
  double xx = 0;
  double xy = 0;
  double x1 = 0;
  double x2 = 0;
  double ya = 0;
  double yb = 0;
  double yd = 0;
  double yi = 0;
  double ym = 0;
  double yy = 0;
  double y1 = 0;
  double y2 = 0;
  double xcd = 0;
  double ycd = 0;
  //double xcycsq=0;
  //double rcd;

  //*     -----------------------------------------------------------------
  //*     xmid,ymid       = centre of fitted circle of radius radfit
  //*     chisq           = residual sum of squares per point
  //*     alf,bet,gam,cur = internal parameters of the circle

  ifail = 0;
  chisq = 9999.;

  nptlim = np;
  if (nptlim > nptmax) {
    nptlim = nptmax;
    ifail = 5;
    return;
  }

  npf = nptlim + 1;
label_12:
  npf = npf - 1;
  if (npf < 2) {
    ifail = 1;
    return;
  }
  vlf[0] = x[npf - 1] - x[0];
  vlf[1] = y[npf - 1] - y[0];

  //*                get a point near the middle

  ihalf = (npf + 1) / 2;

  vmf[0] = x[ihalf - 1] - x[0];
  vmf[1] = y[ihalf - 1] - y[0];
  //vlm[0] = x[npf-1] - x[ihalf-1];
  //vlm[1] = y[npf-1] - y[ihalf-1];

  //*            check for more than half a cicle of points

  if (vlf[0] * vmf[0] + vlf[1] * vmf[1] < 0.) goto label_12;

  //*            if the track is nearer to the y axis
  //*            interchange x and y

  if (fabs(y[npf - 1] - y[0]) > fabs(x[npf - 1] - x[0])) {
    yf = x;
    xf = y;
    iswap = 1;
  }
  else {
    xf = x;
    yf = y;
    iswap = 0;
  }

  n = npf;
  rn = 1. / ((float)n);
  xm = 0.;
  ym = 0.;
  for (i = 0; i < n; ++i) {
    xm = xm + xf[i];
    ym = ym + yf[i];
  }

  xm = xm * rn;
  ym = ym * rn;
  x2 = 0.;
  y2 = 0.;
  xy = 0.;
  xd = 0.;
  yd = 0.;
  d2 = 0.;
  for (i = 0; i < n; ++i) {
    xi = xf[i] - xm;
    yi = yf[i] - ym;
    xx = xi * xi;
    yy = yi * yi;
    x2 = x2 + xx;
    y2 = y2 + yy;
    xy = xy + xi * yi;
    dd = xx + yy;
    xd = xd + xi * dd;
    yd = yd + yi * dd;
    d2 = d2 + dd * dd;
  }

  x2 = x2 * rn;
  y2 = y2 * rn;
  xy = xy * rn;
  d2 = d2 * rn;
  xd = xd * rn;
  yd = yd * rn;
  f = 3. * x2 + y2;
  g = 3. * y2 + x2;
  fg = f * g;
  h = xy + xy;
  h2 = h * h;
  p2 = xd * xd;
  q2 = yd * yd;
  gam0 = x2 + y2;
  fact = gam0 * gam0;
  a2 = (fg - h2 - d2) / fact;
  fact = fact * gam0;
  a1 = (d2 * (f + g) - 2. * (p2 + q2)) / fact;
  fact = fact * gam0;
  a0 = (d2 * (h2 - fg) + 2. * (p2 * g + q2 * f) - 4. * xd * yd * h) / fact;
  a22 = a2 + a2;
  yb = 1.0e30;
  iter = 0;
  xa = 1.;
  //*                main iteration
label_3:
  ya = a0 + xa * (a1 + xa * (a2 + xa * (xa - 4.0)));
  if (fabs(ya) > fabs(yb)) goto label_4;
  if (iter >= itmax) goto label_4;
  dy = a1 + xa * (a22 + xa * (4. * xa - 12.));
  xb = xa - ya / dy;
  if (fabs(xa - xb) <= eps) goto label_5;
  xa = xb;
  yb = ya;
  iter = iter + 1;
  goto label_3;

label_4:
  ifail = 2;
  //*      print 99, iter
  //*   99 format ('  circle - no convergence after',i4,' iterations.')
  xb = 1.;

label_5:
  gam = gam0 * xb;
  f1 = f - gam;
  g1 = g - gam;
  x1 = xd * g1 - yd * h;
  y1 = yd * f1 - xd * h;
  det = f1 * g1 - h2;
  den = 1. / sqrt(x1 * x1 + y1 * y1 + gam * det * det);
  cur = det * den;
  alf = -(xm * det + x1) * den;
  bet = -(ym * det + y1) * den;
  rm = xm * xm + ym * ym;
  gam = ((rm - gam) * det + 2. * (xm * x1 + ym * y1)) * den * 0.5;
  alm = alf + cur * xm;
  bem = bet + cur * ym;

  gmm = gam + alf * xm + bet * ym + 0.5 * cur * rm;
  chisq = ((0.5 * cur) * (0.5 * cur) * d2 + cur * (alm * xd + bem * yd + gmm * gam0) +
           alm * alm * x2 + bem * bem * y2 + 2. * alm * bem * xy + gmm * gmm) /
          rn;
  cu2 = 0.;
  if (cur != 0.) {
    cu2 = 1. / cur;
    xcd = -alf * cu2;
    ycd = -bet * cu2;
    //rcd    =  fabs(cu2);
    xmid = xcd;
    ymid = ycd;
    xcirccent = xmid;
    ycirccent = ymid;
    if (iswap == 1) {
      xcirccent = ymid;
      ycirccent = xmid;
    }
    rcirc = cu2;
  }
  else {
    rcirc = 1e6;
  }

  return;
}
