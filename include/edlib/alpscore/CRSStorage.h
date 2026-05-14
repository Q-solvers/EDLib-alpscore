#ifndef EDLIB_ALPSCORE_CRSSTORAGE_H
#define EDLIB_ALPSCORE_CRSSTORAGE_H

#include <alps/params.hpp>

#include "edlib/alpscore/Storage.h"
#include "edlib/alpscore/ParamsBridge.h"
#include "edlib/CRSStorage.h"

namespace EDLib {
  namespace Storage {

    template <class ModelType>
    class CRSStorage : public ::edlib::CRSStorage<ModelType> {
      using Base = ::edlib::CRSStorage<ModelType>;
    public:
      using typename Base::Model;
      using Base::Base;

#ifdef USE_MPI
      CRSStorage(alps::params& p, Model& m, MPI_Comm comm)
          : Base(EDLib::compat::from_alps(p), m, comm) {}
#endif
      CRSStorage(alps::params& p, Model& m)
          : Base(EDLib::compat::from_alps(p), m) {}
    };

  }
}

#endif
