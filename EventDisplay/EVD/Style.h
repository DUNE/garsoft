////////////////////////////////////////////////////////////////////////
/// \file Style.h
//
/// \version $Id: Style.h,v 1.1.1.1 2010/11/10 19:44:54 p-novaart Exp $
/// \author messier@indiana.edu
////////////////////////////////////////////////////////////////////////
#ifndef EVD_LINESTYLE_H
#define EVD_LINESTYLE_H
class TLine;

namespace gar {
  namespace evd {
    /// Parameters for drawing options. Allow a consistent style for
    /// drawing particle tracks
    class Style {
    public:
      static const char* LatexName(int pdgcode);
      static void FromPDG(TLine& line, int pdgcode);
      static int ColorFromPDG(int pdgcode);
      static int LineStyleFromPDG(int pdgcode);
      static int LineWidthFromPDG(int pdgcode);
    };
  }
}
#endif
