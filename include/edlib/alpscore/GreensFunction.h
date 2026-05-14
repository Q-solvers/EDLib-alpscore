#ifndef EDLIB_ALPSCORE_GREENSFUNCTION_H
#define EDLIB_ALPSCORE_GREENSFUNCTION_H

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

#include "edlib/alpscore/ExecutionStatistic.h"
#include "edlib/alpscore/Lanczos.h"
#include "edlib/alpscore/MeshFactory.h"
#include "edlib/alpscore/GfBridge.h"
#include "edlib/alpscore/MeshBridge.h"
#include "edlib/alpscore/ParamsBridge.h"
#include "edlib/GreensFunction.h"

namespace EDLib {
  namespace gf {

    /**
     * Legacy alpscore-aware GreensFunction shim. Holds an
     * edlib::GreensFunction core impl. After compute(), bridges its edlib::GF3
     * results into cached alps::gf::three_index_gf objects so the legacy
     * accessors (G, G_ij, save, compute_selfenergy) keep returning alps types.
     */
    template <class Hamiltonian, class MeshFactory, typename... Args>
    class GreensFunction {
      using CoreMesh = typename MeshFactory::CoreMeshType;
      using AlpsMesh = typename MeshFactory::MeshType;
      using CoreImpl = ::edlib::GreensFunction<Hamiltonian, CoreMesh>;

    public:
      using precision = typename Hamiltonian::ModelType::precision;
      using Mesh      = AlpsMesh;
      using GF_TYPE   = alps::gf::three_index_gf<std::complex<double>, AlpsMesh,
                                                  alps::gf::index_mesh,
                                                  alps::gf::index_mesh>;

      GreensFunction(alps::params& p, Hamiltonian& h, Args... args)
          : _impl(EDLib::compat::from_alps(p), h,
                  MeshFactory::createCoreMesh(p, args...),
                  read_orbital_pairs(p)),
            _alps_mesh(MeshFactory::createMesh(p, args...)),
            _G_g    (_alps_mesh, alps::gf::index_mesh(h.model().interacting_orbitals()),
                     alps::gf::index_mesh(p["NSPINS"].as<int>())),
            _G_l    (_alps_mesh, alps::gf::index_mesh(h.model().interacting_orbitals()),
                     alps::gf::index_mesh(p["NSPINS"].as<int>())),
            _G      (_alps_mesh, alps::gf::index_mesh(h.model().interacting_orbitals()),
                     alps::gf::index_mesh(p["NSPINS"].as<int>())),
            _G_g_ij (_alps_mesh, alps::gf::index_mesh(h.model().interacting_orbitals() * h.model().interacting_orbitals()),
                     alps::gf::index_mesh(p["NSPINS"].as<int>())),
            _G_l_ij (_alps_mesh, alps::gf::index_mesh(h.model().interacting_orbitals() * h.model().interacting_orbitals()),
                     alps::gf::index_mesh(p["NSPINS"].as<int>())),
            _G_ij   (_alps_mesh, alps::gf::index_mesh(h.model().interacting_orbitals() * h.model().interacting_orbitals()),
                     alps::gf::index_mesh(p["NSPINS"].as<int>())),
            _ham(h) {}

      void compute() {
        _impl.compute();
        copy_gf(_impl.G_g(),    _G_g);
        copy_gf(_impl.G_l(),    _G_l);
        copy_gf(_impl.G(),      _G);
        copy_gf(_impl.G_g_ij(), _G_g_ij);
        copy_gf(_impl.G_l_ij(), _G_l_ij);
        copy_gf(_impl.G_ij(),   _G_ij);
      }

      const GF_TYPE& G_g()    const { return _G_g; }
      const GF_TYPE& G_l()    const { return _G_l; }
      const GF_TYPE& G()      const { return _G; }
      const GF_TYPE& G_g_ij() const { return _G_g_ij; }
      const GF_TYPE& G_l_ij() const { return _G_l_ij; }
      const GF_TYPE& G_ij()   const { return _G_ij; }

      void save(alps::hdf5::archive& ar, const std::string& path) const {
#ifdef USE_MPI
        int rank;
        MPI_Comm_rank(_ham.storage().comm(), &rank);
        if (rank != 0) return;
#endif
        if (_impl.G().shape(1) > 0) {
          _G.save(ar, path + "/G_omega");
          std::ofstream f("G_omega");
          f << std::setprecision(14) << _G;
          f.close();
        }
        if (_impl.G_ij().shape(1) > 0) {
          _G_ij.save(ar, path + "/G_ij_omega");
          std::ofstream f("G_ij_omega");
          f << std::setprecision(14) << _G_ij;
          f.close();
        }
      }

      /// Compute self-energy via Dyson's equation, write to archive, return
      /// the alps::gf-typed sigma.
      GF_TYPE compute_selfenergy(alps::hdf5::archive& ar, const std::string& path) {
        GF_TYPE bare (_G_ij.mesh1(), _G_ij.mesh2(), _G_ij.mesh3());
        GF_TYPE sigma(_G_ij.mesh1(), _G_ij.mesh2(), _G_ij.mesh3());
        _ham.model().bare_greens_function(bare, _impl.beta());
        bare.save(ar, path + "/G0_omega");
        {
          std::ofstream f("G0_omega");
          f << std::setprecision(14) << bare;
        }
        _ham.model().solve_dyson(bare, _G_ij, sigma);
        sigma.save(ar, path + "/Sigma_omega");
        {
          std::ofstream f("Sigma_omega");
          f << std::setprecision(14) << sigma;
        }
        return sigma;
      }

    private:
      template <class CoreGF3>
      void copy_gf(const CoreGF3& src, GF_TYPE& dst) {
        for (int iw = 0; iw < src.shape(0); ++iw) {
          typename AlpsMesh::index_type w(iw);
          for (int i = 0; i < src.shape(1); ++i) {
            for (int j = 0; j < src.shape(2); ++j) {
              dst(w, alps::gf::index_mesh::index_type(i),
                  alps::gf::index_mesh::index_type(j)) = src(iw, i, j);
            }
          }
        }
      }

      static std::vector<std::array<int, 2>> read_orbital_pairs(alps::params& p) {
        std::vector<std::array<int, 2>> out;
        std::string input = p["INPUT_FILE"];
        alps::hdf5::archive ar(input, "r");
        if (ar.is_data("GreensFunction_orbitals/values")) {
          std::vector<std::vector<int>> raw;
          ar >> alps::make_pvp("GreensFunction_orbitals/values", raw);
          for (const auto& v : raw) out.push_back({v[0], v[1]});
        }
        ar.close();
        return out;
      }

      CoreImpl     _impl;
      AlpsMesh     _alps_mesh;
      GF_TYPE      _G_g, _G_l, _G;
      GF_TYPE      _G_g_ij, _G_l_ij, _G_ij;
      Hamiltonian& _ham;
    };

  }
}

#endif
