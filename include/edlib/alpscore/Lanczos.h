#ifndef EDLIB_LANCZOS_H
#define EDLIB_LANCZOS_H

#include <alps/params.hpp>

#include "edlib/alpscore/MeshFactory.h"
#include "edlib/alpscore/ParamsBridge.h"
#include "edlib/core/Lanczos.h"

namespace EDLib {
  namespace gf {

    /**
     * Legacy alpscore-aware Lanczos shim. Templated on a legacy MeshFactory
     * (MatsubaraMeshFactory / RealFreqMeshFactory) plus its trailing Args
     * (statistics_type for matsubara). Forwards to the core
     * edlib::Lanczos<Hamiltonian, CoreMeshType>.
     */
    template <class Hamiltonian, class MeshFactory, typename... Args>
    class Lanczos : public ::edlib::Lanczos<Hamiltonian,
                                            typename MeshFactory::CoreMeshType> {
      using Base = ::edlib::Lanczos<Hamiltonian, typename MeshFactory::CoreMeshType>;
    public:
      using typename Base::precision;
      using Mesh = typename MeshFactory::CoreMeshType;

      Lanczos(alps::params& p, Hamiltonian& h, Args... args)
          : Base(EDLib::compat::from_alps(p), h,
                 MeshFactory::createCoreMesh(p, args...)) {}
    };

  }
}

#endif
