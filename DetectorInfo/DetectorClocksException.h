/**
 * \file DetectorClocksException.h
 *
 * \ingroup DetectorClocks
 *
 * \brief Class def header for exception classes in DetectorClocks data provider
 *
 * @author kazuhiro
 */

/** \addtogroup TimeService

    @{*/
#ifndef DETECTORCLOCKSEXCEPTION_H
#define DETECTORCLOCKSEXCEPTION_H

#include <exception>
#include <string>

namespace gar {
  namespace detinfo {
    /**
     \class DetectorClocksException
     Simple exception class for DetectorClocks
     */
    class DetectorClocksException : public std::exception {

    public:
      DetectorClocksException(std::string msg = "") : std::exception(), _msg(msg) {}

      virtual ~DetectorClocksException() throw(){};
      virtual const char* msg() const throw() { return _msg.c_str(); }

    private:
      std::string _msg;
    };

  }
} // gar

#endif
/** @} */ // end of doxygen group
