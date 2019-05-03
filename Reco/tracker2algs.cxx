#include "Reco/tracker2algs.h"


// these two methods are also duplicated in patrec.  Maybe we can use the patrec tracks as the
// initial guess of the track parameters.  Make sure we get them in the right order

int gar::rec::initial_trackpar_estimate(const std::vector<gar::rec::TPCCluster>  &TPCClusters,
std::vector<int> &TPCClusterlist,
float &curvature_init,
float &lambda_init,
float &phi_init,
float &xpos,
float &ypos,
float &zpos,
float &x_other_end,
unsigned int initialtpnTPCClusters,
int printlevel)
{

    size_t nTPCClusters = TPCClusterlist.size();
    size_t firstTPCCluster = 0;
    size_t farTPCCluster = TMath::Min(nTPCClusters-1, (size_t) initialtpnTPCClusters);;
    size_t intTPCCluster = farTPCCluster/2;
    size_t lastTPCCluster = nTPCClusters-1;

    float trackbeg[3] = {TPCClusters[TPCClusterlist[firstTPCCluster]].Position()[0],
    TPCClusters[TPCClusterlist[firstTPCCluster]].Position()[1],
    TPCClusters[TPCClusterlist[firstTPCCluster]].Position()[2]};

    float tp1[3] = {TPCClusters[TPCClusterlist[intTPCCluster]].Position()[0],
    TPCClusters[TPCClusterlist[intTPCCluster]].Position()[1],
    TPCClusters[TPCClusterlist[intTPCCluster]].Position()[2]};

    float tp2[3] = {TPCClusters[TPCClusterlist[farTPCCluster]].Position()[0],
    TPCClusters[TPCClusterlist[farTPCCluster]].Position()[1],
    TPCClusters[TPCClusterlist[farTPCCluster]].Position()[2]};
    
    if (printlevel>1)
    {
        std::cout << "TPCCluster Dump in initial_trackpar_estimate: " << std::endl;
        for (size_t i=0;i<nTPCClusters;++i)
        {
            size_t ihf = i;
            std::cout << i << " : " <<
            TPCClusters[TPCClusterlist[ihf]].Position()[0] << " " <<
            TPCClusters[TPCClusterlist[ihf]].Position()[1] << " " <<
            TPCClusters[TPCClusterlist[ihf]].Position()[2] << std::endl;
        }
    }
    if (printlevel>0)
    {
        std::cout << "first TPCCluster: " << firstTPCCluster << ", inter TPCCluster: " << intTPCCluster << " " << " far TPCCluster: " << farTPCCluster << std::endl;
        std::cout << "in the TPCCluster list: " << TPCClusterlist[firstTPCCluster] << " " << TPCClusterlist[intTPCCluster] << " " << TPCClusterlist[farTPCCluster] << std::endl;
        std::cout << "First TPCCluster x, y, z: " << trackbeg[0] << " " << trackbeg[1] << " " << trackbeg[2] << std::endl;
        std::cout << "Inter TPCCluster x, y, z: " << tp1[0] << " " << tp1[1] << " " << tp1[2] << std::endl;
        std::cout << "Far   TPCCluster x, y, z: " << tp2[0] << " " << tp2[1] << " " << tp2[2] << std::endl;
    }

    xpos = trackbeg[0];
    ypos = trackbeg[1];
    zpos = trackbeg[2];
    x_other_end = TPCClusters[TPCClusterlist[lastTPCCluster]].Position()[0];


    float ycc=0;
    float zcc=0;
    curvature_init = gar::rec::capprox(trackbeg[1],trackbeg[2],tp1[1],tp1[2],tp2[1],tp2[2],ycc,zcc);
    //std::cout << " inputs to trackpar circ fit (y,z): " << trackbeg[1] << " " << trackbeg[2] << " : "
    //	    << tp1[1] << " " << tp1[2] << " : " << tp2[1] << " " << tp2[2] << std::endl;
    //std::cout << "curvature output: " << curvature_init << std::endl;

    phi_init = TMath::ATan2( trackbeg[2] - zcc, ycc - trackbeg[1] );
    float phi2 = phi_init;
    if (curvature_init<0) phi_init += TMath::Pi();
    float radius_init = 10000;
    if (curvature_init != 0) radius_init = 1.0/curvature_init;

    float dx1 = tp2[0] - xpos;
    if (dx1 != 0)
    {
        float dphi2 = TMath::ATan2(tp2[2]-zcc,ycc-tp2[1])-phi2;
        if (dphi2 > TMath::Pi()) dphi2 -= 2.0*TMath::Pi();
        if (dphi2 < -TMath::Pi()) dphi2 += 2.0*TMath::Pi();
        lambda_init = TMath::ATan(1.0/((radius_init/dx1)*dphi2));
    }
    else
    {
        //std::cout << "initial track par estimate failure" << std::endl;
        lambda_init = 0;
        return 1;
    } // got fMinNumTPCClusters all at exactly the same value of x (they were sorted).  Reject track.

    if (printlevel>0)
    {
        std::cout << "phi calc: dz, dy " << tp2[2]-trackbeg[2] << " " <<  tp2[1]-trackbeg[1] << std::endl;
        std::cout << "initial curvature, phi, lambda: " << curvature_init << " " << phi_init << " " << lambda_init << std::endl;
    }
    return 0;

}


float gar::rec::capprox(float x1,float y1,
float x2,float y2,
float x3,float y3,
float &xc, float &yc)
{
    //-----------------------------------------------------------------
    // Initial approximation of the track curvature -- copied from ALICE
    // here x is y and y is z for us
    //-----------------------------------------------------------------
    x3 -=x1;
    x2 -=x1;
    y3 -=y1;
    y2 -=y1;
    //
    float det = x3*y2-x2*y3;
    if (TMath::Abs(det)<1e-10){
        return 100;
    }
    //
    float u = 0.5* (x2*(x2-x3)+y2*(y2-y3))/det;
    float x0 = x3*0.5-y3*u;
    float y0 = y3*0.5+x3*u;
    float c2 = 1/TMath::Sqrt(x0*x0+y0*y0);
    xc = x0 + x1;
    yc = y0 + y1;
    if (det<0) c2*=-1;
    return c2;
}


// given a set of TPC Clusters, find sort orders hlf and hlb, which are vectors of indices into TPC Clusters,
// that run along the track as best as we can.  Which way is "forwards" and which "backwards" is arbitrary

// sorttransweight is the fraction of additional distance given to the transverse displacement from the
// local line candidate in addtion to the longitudinal component of the distance used to pick the closest next hit
// sortdisback is the turnover point for negative longitudinal distances so we can go back and pick up mis-sorted hits.

void gar::rec::sort_TPCClusters_along_track(const std::vector<gar::rec::TPCCluster>  &TPCClusters,
std::vector<int> &hlf,
std::vector<int> &hlb,
int printlevel,
float &lengthforwards,
float &lengthbackwards,
float sorttransweight,
float sortdistback)
{

    // this sorting code appears in the patrec stage too, if only to make tracks that can be drawn
    // on the event display

    // find candidate endpoints and sort TPCClusters

    float cmin[3];  // min x, y, and z coordinates over all TPCClusters
    float cmax[3];  // max x, y, and z coordinates over all TPCClusters
    size_t ihex[6];  // index of TPCCluster which gave the min or max ("extreme") 0-2: (min xyz)  3-5 (max xyz)


    for (size_t iTPCCluster=0; iTPCCluster < TPCClusters.size(); ++iTPCCluster)
    {
        for (int i=0; i<3; ++i)
        {
            float c = TPCClusters[iTPCCluster].Position()[i];
            if (iTPCCluster==0)
            {
                cmin[i] = c;
                cmax[i] = c;
                ihex[i] = 0;
                ihex[i+3] = 0;
            }
            else
            {
                if (c<cmin[i])
                {
                    cmin[i] = c;
                    ihex[i] = iTPCCluster;
                }
                if (c>cmax[i])
                {
                    cmax[i] = c;
                    ihex[i+3] = iTPCCluster;
                }
            }
        }
    }
    // now we have six TPCClusters that have the min and max x, y, and z values.  Find out which of these six
    // TPCClusters has the biggest sum of distances to all the other TPCClusters (the most extreme)
    float sumdmax = 0;
    size_t imax = 0;
    for (size_t i=0; i<6; ++i)
    {
        float sumd = 0;
        TVector3 poshc(TPCClusters[ihex[i]].Position());
        for (size_t iTPCCluster=0; iTPCCluster<TPCClusters.size(); ++iTPCCluster)
        {
            TVector3 hp(TPCClusters[iTPCCluster].Position());
            sumd += (poshc - hp).Mag();
        }
        if (sumd > sumdmax)
        {
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
    TVector3 lastpos=lpos;
    TVector3 cdpos=lpos;  // position of the beginning of the direction vector
    size_t cdposindex=0;  // and its index in the hlf vector
    TVector3 hpadd;
    TVector3 curdir(0,0,0);
    bool havecurdir = false;
    float dmindir=10;       // promote to a fcl parameter someday.  Distance over which to compute curdir

    for (size_t inh=1;inh<TPCClusters.size();++inh)
    {
        float dmin=0;
        int jmin=-1;
        for (size_t jh=0;jh<TPCClusters.size();++jh)
        {
            bool found = false;
            for (size_t kh=0;kh<hlf.size();++kh)
            {
                if (hlf[kh] == (int) jh)
                {
                    found = true;
                    break;
                }
            }
            if (found) continue;   // skip if we've already assigned this TPCCluster on this track
            TVector3 hpos(TPCClusters[jh].Position());

            float d=0;
            if (havecurdir)
            {
                float dtmp = (hpos-lastpos).Dot(curdir);
                if (dtmp<-2)
                {
                    dtmp = -4 - dtmp;
                }
                d= dtmp + 0.1 * ((hpos-lastpos).Cross(curdir)).Mag();
            }
            else
            {
                d=(hpos-lpos).Mag();
            }
            if (jmin == -1)
            {
                jmin = jh;
                dmin = d;
                hpadd = hpos;
            }
            else
            {
                if (d<dmin)
                {
                    jmin = jh;
                    dmin = d;
                    hpadd = hpos;
                }
            }
        }
        //  std::cout << "dmin: " << dmin << std::endl;
        hlf.push_back(jmin);
        lengthforwards += (hpadd-lastpos).Mag();   // add up track length
        lastpos = hpadd;
        if ( (hpadd-cdpos).Mag() > dmindir)
        {
            TVector3 cdcand(TPCClusters[hlf[cdposindex]].Position());
            size_t itctmp = cdposindex;
            for (size_t itc=cdposindex+1; itc<hlf.size(); ++itc)
            {
                TVector3 cdcandt(TPCClusters[hlf[itc]].Position());
                if (  (hpadd-cdcandt).Mag()  < dmindir  )
                {
                    break;
                }
                else
                {
                    itctmp = itc;
                    cdcand = cdcandt;
                }
            }
            //std::cout << "Updated curdir: " << cdposindex << " " << itctmp << " " << curdir.X() << " " << curdir.Y() << " " << curdir.Z() << " ";
            cdposindex = itctmp;
            cdpos = cdcand;
            curdir = hpadd - cdcand;
            curdir *= (1.0/curdir.Mag());
            //std::cout << "To: " << curdir.X() << " " << curdir.Y() << " " << curdir.Z() << std::endl;
            havecurdir = true;
        }
    }

    if (printlevel>2)
    {
        for (size_t iTPCCluster=0; iTPCCluster<TPCClusters.size(); ++iTPCCluster)
        {
            printf("Sort compare: %5d %10.3f %10.3f %10.3f  %5d %10.3f %10.3f %10.3f\n",
            (int) iTPCCluster,
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
    for (size_t i=0; i< hlf.size(); ++i)
    {
        hlb.push_back(hlf[hlf.size()-1-i]);  // just invert the order for the backward sort
        TVector3 curpos(TPCClusters[hlf[i]].Position());
    }
    lengthbackwards = lengthforwards;

    if (printlevel>2)
    {
        for (size_t iTPCCluster=0; iTPCCluster<TPCClusters.size(); ++iTPCCluster)
        {
            printf("Backward Sort compare: %5d %10.3f %10.3f %10.3f  %5d %10.3f %10.3f %10.3f\n",
            (int) iTPCCluster,
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
