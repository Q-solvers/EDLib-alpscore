#ifndef EDLIB_ALPSCORE_SPINRESOLVEDSTORAGE_H
#define EDLIB_ALPSCORE_SPINRESOLVEDSTORAGE_H

#include <alps/params.hpp>

#include "edlib/alpscore/Storage.h"
#include "edlib/alpscore/ParamsBridge.h"
#include "edlib/SpinResolvedStorage.h"

namespace EDLib {
  namespace Storage {

    template <class ModelType>
    class SpinResolvedStorage : public ::edlib::SpinResolvedStorage<ModelType> {
      using Base = ::edlib::SpinResolvedStorage<ModelType>;
    public:
      using typename Base::Model;
      using typename Base::prec;
      using typename Base::Matrix;
      using Base::Base;

#ifdef USE_MPI
      SpinResolvedStorage(alps::params& p, Model& m, MPI_Comm comm)
          : Base(EDLib::compat::from_alps(p), m, comm) {}
#else
      SpinResolvedStorage(alps::params& p, Model& m)
          : Base(EDLib::compat::from_alps(p), m) {}
#endif
    };

  }
}

#endif
