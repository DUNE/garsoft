/// \file    raw.h
/// \brief   raw data utilities -- compression and decompression
/// \author  trj@fnal.gov
/// many thanks to Brian Rebel and Jonathan Insler

#ifndef GAR_RAWDATA_RAW_H
#define GAR_RAWDATA_RAW_H

#include "RawDataProducts/RawTypes.h"
#include <vector>

namespace gar {

  namespace raw {

    /**
   * @brief Uncompresses a raw data buffer
   * @param adc compressed buffer
   * @param uncompressed buffer to be filled with uncompressed data
   * @param compress type of compression in the adc buffer
   *
   * This function dispatches the uncompression to the correct uncompress
   * function according to compression type in compress.
   *
   * The uncompressed buffer *must* be already allocated with enough space
   * to store the full inflated adc data. Uncompressing raw::RawDigit can
   * be done as follows:
   *
   *     ADCvector_t uncompressed(digit.Samples(), 0);
   *     raw::Uncompress(digit.ADC(), uncompressed, digit.ADC());
   *
   *
   */

    // uncompressing Huffman-encoded data or other compression algs that do not need to backfill
    // with pedestal

    void Uncompress(const gar::raw::ADCvector_t& adc,
                    gar::raw::ADCvector_t& uncompressed,
                    gar::raw::Compress_t compress);

    // for filling in zero-suppressed data -- put pedestal in for the missing samples.

    void Uncompress(const gar::raw::ADCvector_t& adc,
                    gar::raw::ADCvector_t& uncompressed,
                    gar::raw::ADC_t pedestal,
                    gar::raw::Compress_t compress);

    int Compress(gar::raw::ADCvector_t& adc,
                 gar::raw::Compress_t compress,
                 gar::raw::ADC_t zerothreshold,
                 size_t ticksbefore,
                 size_t ticksafter);

    /**
   * @brief In-place compression of raw data buffer
   * @param adc buffer with uncompressed data
   * @param compress type of compression to be applied
   *
   * This function dispatches the compression to the function appropriate
   * for the specified compression type.
   * The resulting compressed data replaces the input buffer content, which is lost.
   * Compression is expected to reduce the size of the data, so that there is
   * in principle no need for reallocation of the input buffer, adc, to store
   * the result.
   */

    void Compress(gar::raw::ADCvector_t& adc, gar::raw::Compress_t compress);

    // duplicate of above?
    //int Compress(gar::raw::ADCvector_t   &adc,
    //            gar::raw::Compress_t    compress,
    //            gar::raw::ADC_t         zerothreshold,
    //		size_t                  ticksbefore,
    //		size_t                  ticksafter);

    void CompressHuffman(gar::raw::ADCvector_t& adc);

    int ZeroSuppression(gar::raw::ADCvector_t& adc,
                        gar::raw::ADC_t zerothreshold,
                        size_t ticksbefore,
                        size_t ticksafter);

    void ZeroUnsuppression(const gar::raw::ADCvector_t& adc, gar::raw::ADCvector_t& uncompressed);

    void UncompressHuffman(const gar::raw::ADCvector_t& adc, gar::raw::ADCvector_t& uncompressed);

  } // namespace raw

} // namespace gar

#endif // GAR_RAWDATA_RAW_H
