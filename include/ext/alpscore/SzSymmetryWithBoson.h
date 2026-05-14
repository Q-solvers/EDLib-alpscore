#ifndef EDLIB_NSYMMETRYWITHBOSON_H
#define EDLIB_NSYMMETRYWITHBOSON_H

#include <array>
#include <string>
#include <vector>

#include <alps/hdf5.hpp>
#include <alps/params.hpp>

#include <edlib/alpscore/HDF5Utils.h>
#include <edlib/alpscore/ParamsBridge.h>
#include <ext/core/SzSymmetryWithBoson.h>

namespace EDLib {
  namespace Symmetry {

    /**
     * Legacy alpscore-aware SzSymmetryWithBoson shim. Forwards to
     * edlib::ext::SzSymmetryWithBoson; adds an alps::params ctor that
     * reads sectors from HDF5 when arpack.SECTOR is set.
     */
    class SzSymmetryWithBoson : public ::edlib::ext::SzSymmetryWithBoson {
    public:
      using ::edlib::ext::SzSymmetryWithBoson::SzSymmetryWithBoson;
      using Sector = ::edlib::ext::SzSymmetryWithBoson::Sector;

      explicit SzSymmetryWithBoson(alps::params& p)
          : ::edlib::ext::SzSymmetryWithBoson(
                EDLib::compat::from_alps(p),
                p["NBBITS"].as<int>() * p["NBLEVEL"].as<int>(),
                read_sectors(p)) {}

    private:
      static std::vector<std::array<int, 2>> read_sectors(alps::params& p) {
        if (!(p.exists("arpack.SECTOR") && bool(p["arpack.SECTOR"]))) return {};
        std::vector<std::vector<int>> raw;
        std::string input = p["INPUT_FILE"];
        alps::hdf5::archive ar(input, "r");
        ar >> alps::make_pvp("sectors/values", raw);
        ar.close();
        std::vector<std::array<int, 2>> out;
        out.reserve(raw.size());
        for (const auto& v : raw) out.push_back({v[0], v[1]});
        return out;
      }
    };

  }

  // The legacy HDF5Utils specialisation for the bosonic-sector type lives
  // here for back-compat with old user code.
  namespace hdf5 {
    template <>
    inline void HDF5Utils<typename Symmetry::SzSymmetryWithBoson::Sector>::save(
        const typename Symmetry::SzSymmetryWithBoson::Sector& s,
        alps::hdf5::archive& ar, const std::string& path) {
      ar[path + "/nup"]   << s.nup();
      ar[path + "/ndown"] << s.ndown();
      ar[path + "/size"]  << s.size();
    }
  }
}

#endif
