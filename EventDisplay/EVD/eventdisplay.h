/// \file    eventdisplay.h
/// \brief   Place to keep constants for event display
/// \author  brebel@fnal.gov
/// \version $Id: raw.h,v 1.2 2010/04/15 18:10:18 brebel Exp $
#ifndef EVENTDISPLAY_EVD_H
#define EVENTDISPLAY_EVD_H

///Event display classes
namespace gar {
  namespace evd {

    static const int kNCOLS = 14;
    static const int kColor[kNCOLS] = {kRed + 2,
                                       kGreen + 2,
                                       kBlue + 2,
                                       kYellow + 2,
                                       kMagenta - 9,
                                       kCyan - 6,
                                       8,
                                       29,
                                       30,
                                       38,
                                       40,
                                       41,
                                       42,
                                       46};
    static const int kColor2[kNCOLS] =
      {kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, 8, 29, 30, 38, 40, 41, 42, 46};

  }
}
#endif // EVENTDISPLAY_EVD_H
