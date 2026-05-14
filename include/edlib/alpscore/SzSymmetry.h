#ifndef EDLIB_ALPSCORE_SZCOMBINATION_H
#define EDLIB_ALPSCORE_SZCOMBINATION_H

#include <array>
#include <string>
#include <vector>

#include <alps/hdf5.hpp>
#include <alps/params.hpp>

#include "edlib/alpscore/ParamsBridge.h"
#include "edlib/SzSymmetry.h"

namespace EDLib {
  namespace Symmetry {

    /**
     * Legacy alpscore-aware SzSymmetry shim. Forwards to edlib::SzSymmetry;
     * adds an alps::params ctor that pulls sector restrictions from the
     * INPUT_FILE HDF5 when arpack.SECTOR is set.
     */
    class SzSymmetry : public ::edlib::SzSymmetry {
    public:
      using ::edlib::SzSymmetry::SzSymmetry;
      using Sector = ::edlib::SzSymmetry::Sector;

      explicit SzSymmetry(alps::params& p)
          : ::edlib::SzSymmetry(EDLib::compat::from_alps(p), read_sectors(p)) {}

    private:
      static std::vector<std::array<int, 2>> read_sectors(alps::params& p) {
        if (!(p.exists("arpack.SECTOR") && bool(p["arpack.SECTOR"]))) return {};
        std::vector<std::vector<int>> sectors_list;
        std::string input = p["INPUT_FILE"];
        alps::hdf5::archive ar(input, "r");
        ar >> alps::make_pvp("sectors/values", sectors_list);
        ar.close();
        std::vector<std::array<int, 2>> out;
        out.reserve(sectors_list.size());
        for (const auto& v : sectors_list) out.push_back({v[0], v[1]});
        return out;
      }
    };

  }
}

#endif
