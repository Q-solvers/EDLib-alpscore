#ifndef EDLIB_ALPSCORE_STORAGE_H
#define EDLIB_ALPSCORE_STORAGE_H

#include <alps/params.hpp>

#include "edlib/alpscore/ParamsBridge.h"
#include "edlib/Storage.h"

namespace EDLib {
  namespace Storage {

    template <class Prec>
    class Storage : public ::edlib::Storage<Prec> {
    public:
      using ::edlib::Storage<Prec>::Storage;

#ifdef USE_MPI
      Storage(alps::params& p, MPI_Comm comm)
          : ::edlib::Storage<Prec>(EDLib::compat::from_alps(p), comm) {}
#endif
      explicit Storage(alps::params& p)
          : ::edlib::Storage<Prec>(EDLib::compat::from_alps(p)) {}
    };

  }
}

#endif
