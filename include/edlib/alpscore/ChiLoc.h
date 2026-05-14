#ifndef EDLIB_ALPSCORE_CHILOC_H
#define EDLIB_ALPSCORE_CHILOC_H

#include <array>
#include <complex>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include <alps/gf/gf.hpp>
#include <alps/hdf5.hpp>
#include <alps/params.hpp>

#include "edlib/alpscore/Lanczos.h"
#include "edlib/alpscore/MeshFactory.h"
#include "edlib/alpscore/ParamsBridge.h"
#include "edlib/ChiLoc.h"

namespace EDLib {
  namespace gf {

    // Re-export bosonic operators from the core under the legacy namespace.
    template <typename Prec> using BosonicOperator = ::edlib::BosonicOperator<Prec>;
    template <typename Prec> using SzOperator      = ::edlib::SzOperator<Prec>;
    template <typename Prec> using NOperator       = ::edlib::NOperator<Prec>;

    /**
     * Legacy alpscore-aware ChiLoc shim. Wraps edlib::ChiLoc; converts the
     * GF2 result into cached alps::gf::two_index_gf return types.
     */
    template <class Hamiltonian, class MeshFactory, typename... Args>
    class ChiLoc {
      using CoreMesh = typename MeshFactory::CoreMeshType;
      using AlpsMesh = typename MeshFactory::MeshType;
      using CoreImpl = ::edlib::ChiLoc<Hamiltonian, CoreMesh>;

    public:
      using precision = typename Hamiltonian::ModelType::precision;
      using Mesh      = AlpsMesh;
      using GF_TYPE   = alps::gf::two_index_gf<std::complex<double>, AlpsMesh,
                                                alps::gf::index_mesh>;

      ChiLoc(alps::params& p, Hamiltonian& h, Args... args)
          : _impl(EDLib::compat::from_alps(p), h,
                  MeshFactory::createCoreMesh(p, args...),
                  read_orbital_pairs(p)),
            _alps_mesh(MeshFactory::createMesh(p, args...)),
            _gf   (_alps_mesh, alps::gf::index_mesh(h.model().interacting_orbitals())),
            _gf_ij(_alps_mesh, alps::gf::index_mesh(h.model().interacting_orbitals()
                                                    * h.model().interacting_orbitals())),
            _type("Sz"),
            _ham(h) {}

      template <typename Op = SzOperator<precision>>
      void compute(const double* avg_ptr = nullptr) {
        _impl.template compute<Op>(avg_ptr);
        _type = Op().name();
        copy_gf2(_impl.G(),    _gf);
        copy_gf2(_impl.G_ij(), _gf_ij);
      }

      const GF_TYPE& G()    const { return _gf; }
      const GF_TYPE& G_ij() const { return _gf_ij; }

      void save(alps::hdf5::archive& ar, const std::string& path) const {
#ifdef USE_MPI
        int rank;
        MPI_Comm_rank(_ham.storage().comm(), &rank);
        if (rank != 0) return;
#endif
        if (_impl.G().shape(1) > 0) {
          _gf.save(ar, path + "/Chi" + _type + "_omega");
          std::ofstream f("Chi" + _type + "_omega");
          f << std::setprecision(14) << _gf;
          f.close();
        }
        if (_impl.G_ij().shape(1) > 0) {
          _gf_ij.save(ar, path + "/Chi" + _type + "_ij_omega");
          std::ofstream f("Chi_ij_" + _type + "_omega");
          f << std::setprecision(14) << _gf_ij;
          f.close();
        }
      }

    private:
      template <class CoreGF2>
      void copy_gf2(const CoreGF2& src, GF_TYPE& dst) {
        for (int iw = 0; iw < src.shape(0); ++iw) {
          for (int i = 0; i < src.shape(1); ++i) {
            dst(typename AlpsMesh::index_type(iw),
                alps::gf::index_mesh::index_type(i)) = src(iw, i);
          }
        }
      }

      static std::vector<std::array<int, 2>> read_orbital_pairs(alps::params& p) {
        std::vector<std::array<int, 2>> out;
        std::string input = p["INPUT_FILE"];
        alps::hdf5::archive ar(input, "r");
        if (ar.is_data("ChiLoc_orbitals/values")) {
          std::vector<std::vector<int>> raw;
          ar >> alps::make_pvp("ChiLoc_orbitals/values", raw);
          for (const auto& v : raw) out.push_back({v[0], v[1]});
        }
        ar.close();
        return out;
      }

      CoreImpl     _impl;
      AlpsMesh     _alps_mesh;
      GF_TYPE      _gf;
      GF_TYPE      _gf_ij;
      std::string  _type;
      Hamiltonian& _ham;
    };

  }
}

#endif
