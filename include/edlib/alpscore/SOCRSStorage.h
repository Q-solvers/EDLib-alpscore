#ifndef EDLIB_ALPSCORE_SOCRSSTORAGE_H
#define EDLIB_ALPSCORE_SOCRSSTORAGE_H

#include <alps/params.hpp>

#include "edlib/alpscore/Storage.h"
#include "edlib/alpscore/ParamsBridge.h"
#include "edlib/SOCRSStorage.h"

namespace EDLib {
  namespace Storage {

    template <class ModelType>
    class SOCRSStorage : public ::edlib::SOCRSStorage<ModelType> {
      using Base = ::edlib::SOCRSStorage<ModelType>;
    public:
      using typename Base::Model;
      using Base::Base;

#ifdef USE_MPI
      SOCRSStorage(alps::params& p, Model& m, MPI_Comm comm)
          : Base(EDLib::compat::from_alps(p), m, comm) {}
#endif
      SOCRSStorage(alps::params& p, Model& m)
          : Base(EDLib::compat::from_alps(p), m) {}
    };

  }
}

#endif
