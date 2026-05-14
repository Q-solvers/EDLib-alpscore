#ifndef EDLIB_COMPAT_MESHBRIDGE_H
#define EDLIB_COMPAT_MESHBRIDGE_H

#ifdef EDLIB_WITH_ALPSCORE

#include <alps/gf/grid.hpp>
#include <alps/gf/mesh.hpp>

#include "edlib/core/Mesh.h"

namespace EDLib {
  namespace compat {

    inline alps::gf::statistics::statistics_type
    to_alps_statistics(edlib::Statistics s) {
      return s == edlib::Statistics::Fermionic
          ? alps::gf::statistics::statistics_type::FERMIONIC
          : alps::gf::statistics::statistics_type::BOSONIC;
    }

    inline edlib::Statistics from_alps_statistics(alps::gf::statistics::statistics_type s) {
      return s == alps::gf::statistics::statistics_type::FERMIONIC
          ? edlib::Statistics::Fermionic
          : edlib::Statistics::Bosonic;
    }

    inline alps::gf::matsubara_positive_mesh
    to_alps_mesh(const edlib::MatsubaraMesh& m) {
      return alps::gf::matsubara_positive_mesh(
          m.beta(), m.extent(), to_alps_statistics(m.statistics()));
    }

    inline alps::gf::real_frequency_mesh
    to_alps_mesh(const edlib::RealFreqMesh& m) {
      alps::gf::grid::linear_real_frequency_grid g(m.emin(), m.emax(), m.extent());
      return alps::gf::real_frequency_mesh(g);
    }

  }
}

#endif // EDLIB_WITH_ALPSCORE
#endif
