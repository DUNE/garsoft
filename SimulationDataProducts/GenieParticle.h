//
//  GenieParticle.h
//
//  Created by E. Brianne
//

#ifndef GAR_SIMULATIONDATAPRODUCTS_GenieParticle_h
#define GAR_SIMULATIONDATAPRODUCTS_GenieParticle_h

#include <string>

namespace gar {
  namespace sdp {

    class GenieParticle {
    public:
      GenieParticle();

#ifndef __GCCXML__

      GenieParticle(unsigned int interaction_index,
                    unsigned int index,
                    int pdg,
                    int status,
                    std::string name,
                    int firstmother,
                    int lastmother,
                    int firstdaughter,
                    int lastdaughter,
                    float px,
                    float py,
                    float pz,
                    float e,
                    float mass,
                    int rescatter)
        : fInteractionIndex(interaction_index)
        , fIndex(index)
        , fPdg(pdg)
        , fStatus(status)
        , fName(name)
        , fFirstMother(firstmother)
        , fLastMother(lastmother)
        , fFirstDaughter(firstdaughter)
        , fLastDaughter(lastdaughter)
        , fPx(px)
        , fPy(py)
        , fPz(pz)
        , fE(e)
        , fMass(mass)
        , fRescatterCode(rescatter)
      {}

      unsigned int InteractionIndex() const { return fInteractionIndex; }
      unsigned int Index() const { return fIndex; }
      int Pdg() const { return fPdg; }
      int Status() const { return fStatus; }
      std::string Name() const { return fName; }
      int FirstMother() const { return fFirstMother; }
      int LastMother() const { return fLastMother; }
      int FirstDaughter() const { return fFirstDaughter; }
      int LastDaughter() const { return fLastDaughter; }
      float Px() const { return fPx; }
      float Py() const { return fPy; }
      float Pz() const { return fPz; }
      float E() const { return fE; }
      float Mass() const { return fMass; }
      int RescatterCode() const { return fRescatterCode; }

#endif

    private:
      unsigned int fInteractionIndex;
      unsigned int fIndex;
      int fPdg;
      int fStatus;
      std::string fName;
      int fFirstMother;
      int fLastMother;
      int fFirstDaughter;
      int fLastDaughter;
      float fPx;
      float fPy;
      float fPz;
      float fE;
      float fMass;
      int fRescatterCode;
    };

  } // sdp
} // gar

#endif /* GAR_SIMULATIONDATAPRODUCTS_GenieParticle_h */
