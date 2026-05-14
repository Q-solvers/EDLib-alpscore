#ifndef EDLIB_ALPSCORE_COMMONUTILS_H
#define EDLIB_ALPSCORE_COMMONUTILS_H

// Legacy bridge for the alps::gf-mesh based freq_point overloads. New code
// should use edlib::freq_point with edlib::MatsubaraMesh / RealFreqMesh.

#include <complex>

#ifdef EDLIB_WITH_ALPSCORE
#include <alps/gf/mesh.hpp>
#endif

namespace EDLib {
  namespace common {

#ifdef EDLIB_WITH_ALPSCORE
    inline std::complex<double>
    freq_point(int index, const alps::gf::matsubara_positive_mesh& mesh, double /*beta*/) {
      return std::complex<double>(0.0, mesh.points()[index]);
    }

    inline std::complex<double>
    freq_point(int index, const alps::gf::real_frequency_mesh& mesh, double beta) {
      return std::complex<double>(mesh.points()[index], M_PI / beta);
    }
#endif

  }
}

#endif
