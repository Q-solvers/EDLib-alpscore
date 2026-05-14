#ifndef EDLIB_ALPSCORE_NSYMMETRY_H
#define EDLIB_ALPSCORE_NSYMMETRY_H

#include <string>
#include <vector>

#include <alps/hdf5.hpp>
#include <alps/params.hpp>

#include "edlib/alpscore/ParamsBridge.h"
#include "edlib/NSymmetry.h"

namespace EDLib {
  namespace Symmetry {

    /**
     * Legacy alpscore-aware NSymmetry shim. Forwards to edlib::NSymmetry;
     * adds an alps::params ctor that pulls sector restrictions from the
     * INPUT_FILE HDF5 when arpack.SECTOR is set.
     */
    class NSymmetry : public ::edlib::NSymmetry {
    public:
      using ::edlib::NSymmetry::NSymmetry;
      using Sector = ::edlib::NSymmetry::Sector;

      explicit NSymmetry(alps::params& p)
          : ::edlib::NSymmetry(EDLib::compat::from_alps(p), read_sectors(p)) {}

    private:
      static std::vector<int> read_sectors(alps::params& p) {
        if (!(p.exists("arpack.SECTOR") && bool(p["arpack.SECTOR"]))) return {};
        std::vector<std::vector<int>> sectors_list;
        std::string input = p["INPUT_FILE"];
        alps::hdf5::archive ar(input, "r");
        ar >> alps::make_pvp("sectors/values", sectors_list);
        ar.close();
        std::vector<int> out;
        out.reserve(sectors_list.size());
        for (const auto& v : sectors_list) out.push_back(v[0]);
        return out;
      }
    };

  }
}

#endif
