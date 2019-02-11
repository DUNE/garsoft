#include "Reco/tracker2algs.h"


// these two methods are also duplicated in patrec.  Maybe we can use the patrec tracks as the
// initial guess of the track parameters.  Make sure we get them in the right order

int gar::rec::initial_trackpar_estimate(const std::vector<Hit>  &hits,
					std::vector<int> &hitlist,
					float &curvature_init,
					float &lambda_init,
					float &phi_init,
					float &xpos,
					float &ypos,
					float &zpos,
					float &x_other_end,
					unsigned int initialtpnhits,
					int printlevel)
{

  size_t nhits = hitlist.size();
  size_t firsthit = 0;
  size_t farhit = TMath::Min(nhits-1, (size_t) initialtpnhits);;
  size_t inthit = farhit/2;
  size_t lasthit = nhits-1;

  float trackbeg[3] = {hits[hitlist[firsthit]].Position()[0],
		       hits[hitlist[firsthit]].Position()[1],
		       hits[hitlist[firsthit]].Position()[2]};

  float tp1[3] = {hits[hitlist[inthit]].Position()[0],
		  hits[hitlist[inthit]].Position()[1],
		  hits[hitlist[inthit]].Position()[2]};

  float tp2[3] = {hits[hitlist[farhit]].Position()[0],
		  hits[hitlist[farhit]].Position()[1],
		  hits[hitlist[farhit]].Position()[2]};

  if (printlevel>1)
    {
      std::cout << "Hit Dump in initial_trackpar_estimate: " << std::endl;
      for (size_t i=0;i<nhits;++i)
	{
	  size_t ihf = i;
	  std::cout << i << " : " <<
	    hits[hitlist[ihf]].Position()[0] << " " <<
	    hits[hitlist[ihf]].Position()[1] << " " <<
	    hits[hitlist[ihf]].Position()[2] << std::endl;
	}
    }
  if (printlevel>0)
    {
      std::cout << "first hit: " << firsthit << ", inter hit: " << inthit << " " << " far hit: " << farhit << std::endl;
      std::cout << "in the hit list: " << hitlist[firsthit] << " " << hitlist[inthit] << " " << hitlist[farhit] << std::endl;
      std::cout << "First hit x, y, z: " << trackbeg[0] << " " << trackbeg[1] << " " << trackbeg[2] << std::endl;
      std::cout << "Inter hit x, y, z: " << tp1[0] << " " << tp1[1] << " " << tp1[2] << std::endl;
      std::cout << "Far   hit x, y, z: " << tp2[0] << " " << tp2[1] << " " << tp2[2] << std::endl;
    }

  xpos = trackbeg[0];
  ypos = trackbeg[1];
  zpos = trackbeg[2];
  x_other_end = hits[hitlist[lasthit]].Position()[0];


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
    } // got fMinNumHits all at exactly the same value of x (they were sorted).  Reject track.

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


void gar::rec::sort_hits_along_track(const std::vector<Hit>  &hits,
	 	                     std::vector<int> &hlf,
				     std::vector<int> &hlb,
				     int printlevel,
				     float &lengthforwards,
				     float &lengthbackwards)
{

      // this sorting code appears in the patrec stage too, if only to make tracks that can be drawn
      // on the event display

      // find candidate endpoints and sort hits

      float cmin[3];  // min x, y, and z coordinates over all hits
      float cmax[3];  // max x, y, and z coordinates over all hits
      size_t ihex[6];  // index of hit which gave the min or max ("extreme") 0-2: (min xyz)  3-5 (max xyz)

	 
      for (size_t ihit=0; ihit < hits.size(); ++ihit)
	{
	  for (int i=0; i<3; ++i)
	    {
	      float c = hits[ihit].Position()[i];
	      if (ihit==0)
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
		      ihex[i] = ihit;
		    }
		  if (c>cmax[i])
		    {
		      cmax[i] = c;
		      ihex[i+3] = ihit;
		    }
		}
	    }
	}
      // now we have six hits that have the min and max x, y, and z values.  Find out which of these six
      // hits has the biggest sum of distances to all the other hits (the most extreme)
      float sumdmax = 0;
      size_t imax = 0;
      for (size_t i=0; i<6; ++i)
	{
	  float sumd = 0;
	  TVector3 poshc(hits[ihex[i]].Position());
	  for (size_t ihit=0; ihit<hits.size(); ++ihit)
	    {
	      TVector3 hp(hits[ihit].Position());
	      sumd += (poshc - hp).Mag();
	    }
	  if (sumd > sumdmax)
	    {
	      sumdmax = sumd;
	      imax = i;
	    }
	}

      //  Use this hit as a starting point -- find the closest hit to the last
      //  and add it to the newly sorted list hlf.  Change -- sort hits in order of how
      //  far they are from the first hit.  Prevents oscillations in position on sort order.

      hlf.clear();
      lengthforwards = 0;
      hlf.push_back(ihex[imax]);
      TVector3 lpos(hits[hlf[0]].Position());
      TVector3 lastpos=lpos;
      TVector3 hpadd;

      for (size_t inh=1;inh<hits.size();++inh)
	{
	  float dmin=0;
	  float jmin=-1;
	  for (size_t jh=0;jh<hits.size();++jh)
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
	      if (found) continue;   // skip if we've already assigned this hit on this track
	      TVector3 hpos(hits[jh].Position());
	      float d=(hpos-lpos).Mag();
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
	}

      // replace our hit list with our newly sorted hit list.

      if (printlevel>2)
	{
	  for (size_t ihit=0; ihit<hits.size(); ++ihit)
	    {
	      printf("Sort compare: %5d %10.3f %10.3f %10.3f  %5d %10.3f %10.3f %10.3f\n",
		     (int) ihit,
		     hits[ihit].Position()[0],
		     hits[ihit].Position()[1],
		     hits[ihit].Position()[2],
		     hlf[ihit],
		     hits[hlf[ihit]].Position()[0],
		     hits[hlf[ihit]].Position()[1],
		     hits[hlf[ihit]].Position()[2]);
	    }
	}

      // now go backwards -- start at the end hit and use that as a starting point

      hlb.clear();
      lengthbackwards = 0;
      hlb.push_back(hlf.back());
      TVector3 lpos2(hits[hlb[0]].Position());
      TVector3 lastpos2 = lpos2;

      for (size_t inh=1;inh<hits.size();++inh)
	{
	  float dmin=0;
	  float jmin=-1;
	  for (size_t jh=0;jh<hits.size();++jh)
	    {
	      bool found = false;
	      for (size_t kh=0;kh<hlb.size();++kh)
		{
		  if (hlb[kh] == (int) jh)
		    {
		      found = true;
		      break;
		    }
		}
	      if (found) continue;   // skip if we've already assigned this hit on this track
	      TVector3 hpos(hits[jh].Position());
	      float d=(hpos-lpos2).Mag();
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
	  hlb.push_back(jmin);
	  lengthbackwards += (hpadd - lastpos2).Mag(); // add up length
	  lastpos2 = hpadd;
	}
      // replace our hit list with our newly sorted hit list.

      if (printlevel>2)
	{
	  for (size_t ihit=0; ihit<hits.size(); ++ihit)
	    {
	      printf("Sort compare: %5d %10.3f %10.3f %10.3f  %5d %10.3f %10.3f %10.3f\n",
		     (int) ihit,
		     hits[ihit].Position()[0],
		     hits[ihit].Position()[1],
		     hits[ihit].Position()[2],
		     hlb[ihit],
		     hits[hlb[ihit]].Position()[0],
		     hits[hlb[ihit]].Position()[1],
		     hits[hlb[ihit]].Position()[2]);
	    }
	}
}
