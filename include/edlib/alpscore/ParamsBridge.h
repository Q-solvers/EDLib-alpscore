#ifndef EDLIB_COMPAT_PARAMSBRIDGE_H
#define EDLIB_COMPAT_PARAMSBRIDGE_H

#ifdef EDLIB_WITH_ALPSCORE

#include <cstddef>
#include <alps/params.hpp>

#include "edlib/core/Parameters.h"

namespace EDLib {
  namespace compat {

    /**
     * Construct an edlib::Parameters from an alps::params instance defined by
     * EDLib::define_parameters(). Used by the legacy EDLib:: compat classes to
     * forward construction to the new edlib:: core.
     */
    inline edlib::Parameters from_alps(alps::params& p) {
      edlib::Parameters out;
      out.nsites             = p["NSITES"].as<int>();
      out.nspins             = p["NSPINS"].as<int>();

      // exists() = "user supplied this", as opposed to defined() = "declared
      // by define_parameters". arpack.SECTOR / arpack.NCV are declared without
      // defaults; reading them when unsupplied throws.
      if (p.exists("arpack.SECTOR")) {
        out.arpack_sector    = bool(p["arpack.SECTOR"]);
      }
      out.arpack_nev         = p["arpack.NEV"].as<int>();
      if (p.exists("arpack.NCV")) {
        out.arpack_ncv       = p["arpack.NCV"].as<int>();
      }

      out.storage_max_dim    = p["storage.MAX_DIM"].as<std::size_t>();
      out.storage_max_size   = p["storage.MAX_SIZE"].as<std::size_t>();
      out.eigenvalues_only   = (p["storage.EIGENVALUES_ONLY"].as<int>() != 0);
      out.spinstorage_orbital_number = p["spinstorage.ORBITAL_NUMBER"].as<int>();

      out.lanc_nomega        = p["lanc.NOMEGA"].as<int>();
      out.lanc_emin          = p["lanc.EMIN"].as<double>();
      out.lanc_emax          = p["lanc.EMAX"].as<double>();
      out.lanc_nlanc         = p["lanc.NLANC"].as<int>();
      out.lanc_beta          = p["lanc.BETA"].as<double>();
      out.lanc_boltzmann_cutoff = p["lanc.BOLTZMANN_CUTOFF"].as<double>();

      out.siam_norbitals     = p["siam.NORBITALS"].as<int>();
      return out;
    }

  }
}

#endif // EDLIB_WITH_ALPSCORE

#endif
