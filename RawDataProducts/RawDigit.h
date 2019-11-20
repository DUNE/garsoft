/** ****************************************************************************
 * @file RawDigit.h
 * @brief Definition of basic raw digits
 * @author brebel@fnal.gov
 * @see  RawDigit.cxx raw.h
 * 
 * 
 * Changes:
 * 20141210 Gianluca Petrillo (petrillo@fnal.gov)
 *   data architecture revision changes:
 *   - fADC made private
 *   - removed constructor not assigning sample number
 *   - added flags interface
 * 
 * 2/19/2008 Mitch Soderberg
 * -modified RawDigit class slightly to be compatible with binary output of DAQ software,
 *  and to include #samples/channel explicity, instead of via sizeof() methods. 
 * 
 * ****************************************************************************/

#ifndef GAR_RAWDATA_RAWDIGIT_H
#define GAR_RAWDATA_RAWDIGIT_H

// C/C++ standard libraries
#include <stdint.h> // uint32_t
#include <cstdlib> // size_t
#include <vector>
#include "RawDataProducts/RawTypes.h"
#include "RtypesCore.h"  // ULong64_t

/// Raw data description and utilities
namespace gar {
  namespace raw {
    
    typedef uint32_t Channel_t;
    
    /**
     * @brief Collection of charge vs time digitized from a single readout channel
     *
     * This class hosts potentially compressed data.
     * It does not provide methods to uncompress it, not the same object can
     * become compressed/uncompressed or change compression type: to use a
     * compressed RawDigit, one has to create a new buffer, fill and use it:
     *
     *     raw::RawDigit::ADCvector_t ADCs(digits.Samples()); // fix the size!
     *     raw::Uncompress(digits.ADCs(), ADCs, digits.Compression());
     *
     * (remember that you have to provide raw::Uncompress() with a buffer large
     * enough to contain the uncompressed data).
     *
     * Removed; may need to add this back later.  The class provides some flags, defined in FlagIndices_t.
     * The construction of a RawDigit should be for example in the form:
     *
     *     raw::RawDigit::ADCvector_t ADCs;
     *     // ... fill the digits etc.
     *     raw::RawDigit saturatedDigit(
     *       channel, ADCs.size(), ADCs, raw::kNone,
     *       DefaultFlags | SaturationBit
     *       );
     *     raw::RawDigit unsaturatedDigit(
     *       channel, ADCs.size(), ADCs, raw::kNone,
     *       DefaultFlags & ~SaturationBit
     *       );
     *
     *
     */
    class RawDigit {
      
    public:
        /// Type representing a (compressed) vector of ADC counts
      typedef std::vector<short> ADCvector_t;
      
      /*
       // removed waiting for a real use for flags
       /// Type for the digit flags
       typedef std::bitset<16> Flags_t;
       */
      
        /// Default constructor: an empty raw digit with zeros put in for paraneters and an invalid channel
      RawDigit();
      
#ifndef __GCCXML__
    public:
      
      /*
       // removed waiting for a real use for flags
       typedef enum {
       fiSaturation,       ///< saturation flag: set if saturated at any time
       NGArSoftFlags = 8,  ///< GArSoft reserves flags up to here (this excluded)
       NFlagIndices = NGArSoftFlags ///< number of flags
       } FlagIndices_t; ///< type of index of flags
       */
      
      /**
       * @brief Constructor: sets all the fields
       * @param channel ID of the channel the digits were acquired from
       * @param samples number of ADC samples in the uncompressed collection
       * @param adclist list of ADC counts vs. time, compressed
       * @param compression compression algorithm used in adclist
       *
       * Data from the adclist is copied into the raw digits.
       * Pedestal is set to 0 by default.
       */
      RawDigit(Channel_t          channel,
               ULong64_t          samples,
               ADCvector_t const& adclist,
	       gar::raw::Compress_t compress,
	       ULong64_t          time);
      
      /**
       * @brief Constructor: sets all the fields
       * @param channel ID of the channel the digits were acquired from
       * @param samples number of ADC samples in the uncompressed collection
       * @param adclist list of ADC counts vs. time, compressed
       * @param compression compression algorithm used in adclist
       *
       * Data from the adclist is moved into the raw digits.
       * Pedestal is set to 0 by default.
       */
      RawDigit(Channel_t          channel,
               ULong64_t          samples,
               ADCvector_t&&      adclist,
	       gar::raw::Compress_t compress,
	       ULong64_t          time);
      
      /// Set pedestal and its RMS (the latter is 0 by default)
      void            SetPedestal(float ped, float sigma = 1.);
      
      ///@{
      ///@name Accessors
      
      /// Reference to the compressed ADC count vector
      const ADCvector_t& ADCs()     const;
      
      /// Number of elements in the compressed ADC sample vector
      size_t          NADC()        const;
      
      /// ADC vector element number i; no decompression is applied
      short           ADC(int i)    const;
      
      /// DAQ channel this raw data was read from
      Channel_t       Channel()     const;

      /// Compression algorithm selector
      gar::raw::Compress_t Compression() const;
      
      /// Number of samples in the uncompressed ADC data
      ULong64_t       Samples()     const;
      
      /// Pedestal level (ADC counts)
      /// @deprecated Might be removed soon
      float           Pedestal() const;
      
      /// TODO: RMS of the pedestal level?
      float           Sigma()    const;

      /// Timestmap
      ULong64_t       Time()    const;
      
      ///@}
      
      
      
#endif // !__GCCXML__
    private:
      std::vector<short> fADC;      ///< ADC readout per tick, before pedestal subtraction
      
      Channel_t          fChannel;  ///< channel number in the readout
      ULong64_t          fSamples;  ///< number of ticks of the clock
      
      float              fPedestal; ///< pedestal for this channel
      float              fSigma;    ///< sigma of the pedestal counts for this channel
      gar::raw::Compress_t         fCompression; ///< compression scheme used for the ADC vector
      ULong64_t          fTime;    ///< timestamp

    }; // class RawDigit
    
    
  } // namespace raw
} // gar

//------------------------------------------------------------------------------
//--- inline implementation
//---
#ifndef __GCCXML__

inline size_t           gar::raw::RawDigit::NADC()        const { return fADC.size();  }
inline short            gar::raw::RawDigit::ADC(int i)    const { return fADC.at(i);   }
inline gar::raw::RawDigit::ADCvector_t const& gar::raw::RawDigit::ADCs() const { return fADC; }
inline gar::raw::Channel_t gar::raw::RawDigit::Channel()  const { return fChannel;     }
inline ULong64_t        gar::raw::RawDigit::Samples()     const { return fSamples;     }
inline float            gar::raw::RawDigit::Pedestal()    const { return fPedestal;    }
inline float            gar::raw::RawDigit::Sigma()       const { return fSigma;       }
inline ULong64_t        gar::raw::RawDigit::Time()        const { return fTime;       }
inline gar::raw::Compress_t gar::raw::RawDigit::Compression() const  { return fCompression; }

#endif // !__GCCXML__

#endif // gar_RAWDATA_RAWDIGIT_H

////////////////////////////////////////////////////////////////////////
