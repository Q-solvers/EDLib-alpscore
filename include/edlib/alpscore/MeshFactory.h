#ifndef EDLIB_MESHFACTORY_HPP
#define EDLIB_MESHFACTORY_HPP

#include <alps/gf/grid.hpp>
#include <alps/gf/mesh.hpp>
#include <alps/params.hpp>

#include "edlib/alpscore/MeshBridge.h"
#include "edlib/core/Mesh.h"

namespace EDLib {

  /// Legacy MeshFactory adapters. Each exposes both a MeshType (alps::gf::*
  /// mesh, for legacy public API) and a CoreMeshType (edlib::* mesh, used by
  /// the underlying core impl), plus matching createMesh / createCoreMesh.
  class MatsubaraMeshFactory {
  public:
    using MeshType     = alps::gf::matsubara_positive_mesh;
    using CoreMeshType = ::edlib::MatsubaraMesh;

    static MeshType createMesh(alps::params& p,
                               alps::gf::statistics::statistics_type s) {
      return MeshType(p["lanc.BETA"], p["lanc.NOMEGA"], s);
    }

    static CoreMeshType createCoreMesh(alps::params& p,
                                       alps::gf::statistics::statistics_type s) {
      return CoreMeshType(p["lanc.BETA"], p["lanc.NOMEGA"],
                          EDLib::compat::from_alps_statistics(s));
    }
  };

  class RealFreqMeshFactory {
  public:
    using MeshType     = alps::gf::real_frequency_mesh;
    using CoreMeshType = ::edlib::RealFreqMesh;

    static MeshType createMesh(alps::params& p) {
      alps::gf::grid::linear_real_frequency_grid g(
          p["lanc.EMIN"], p["lanc.EMAX"], p["lanc.NOMEGA"]);
      return MeshType(g);
    }

    static CoreMeshType createCoreMesh(alps::params& p) {
      return CoreMeshType(p["lanc.EMIN"], p["lanc.EMAX"], p["lanc.NOMEGA"]);
    }
  };

}

#endif
