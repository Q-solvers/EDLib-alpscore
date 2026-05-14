#ifndef EDLIB_HUBBARDMODEL_H
#define EDLIB_HUBBARDMODEL_H

#include <complex>
#include <string>
#include <vector>

#include <alps/gf/gf.hpp>
#include <alps/gf/mesh.hpp>
#include <alps/hdf5.hpp>
#include <alps/params.hpp>

#include <Eigen/Core>
#include <Eigen/LU>

#include "edlib/alpscore/CommonUtils.h"
#include "edlib/alpscore/ParamsBridge.h"
#include "edlib/core/HubbardModel.h"
#include "edlib/alpscore/FermionicModel.h"
#include "edlib/alpscore/SzSymmetry.h"

namespace EDLib {
  namespace Model {

    namespace Hubbard {
      template <typename Prec>
      using InnerState = ::edlib::hubbard::InnerState<Prec>;
    }

    /**
     * Legacy alpscore-aware HubbardModel shim. Reads hopping/U/μ/J/Hmag
     * from the INPUT_FILE HDF5 and forwards to the edlib::HubbardModel core.
     * Provides an alps::gf-taking bare_greens_function overload for legacy
     * GreensFunction callers.
     */
    template <typename Prec>
    class HubbardModel : public ::edlib::HubbardModel<Prec> {
      using Base = ::edlib::HubbardModel<Prec>;
    public:
      using typename Base::precision;
      using typename Base::SYMMETRY;
      using typename Base::St;
      using typename Base::Sector;

      using Base::Base;  // inherit (Parameters, ModelData) ctor

      explicit HubbardModel(alps::params& p)
          : Base(EDLib::compat::from_alps(p), read_model_data(p)) {}

      /// Legacy alps::gf solve_dyson kept here (and on SIAM) so the legacy
      /// GreensFunction shim's compute_selfenergy can still call
      /// _model.solve_dyson(...). The new code path uses edlib::solve_dyson.
      template <typename Mesh>
      void solve_dyson(
          const alps::gf::three_index_gf<std::complex<double>, Mesh,
              alps::gf::index_mesh, alps::gf::index_mesh>& bare_gf,
          const alps::gf::three_index_gf<std::complex<double>, Mesh,
              alps::gf::index_mesh, alps::gf::index_mesh>& G_ij,
          alps::gf::three_index_gf<std::complex<double>, Mesh,
              alps::gf::index_mesh, alps::gf::index_mesh>& sigma) const {
        const int Ns = this->orbitals();
        for (int iw = 0; iw < (int)bare_gf.mesh1().points().size(); ++iw) {
          typename Mesh::index_type w(iw);
          for (int is : bare_gf.mesh3().points()) {
            Eigen::MatrixXcd bare(Ns, Ns);
            Eigen::MatrixXcd bold(Ns, Ns);
            for (int im : bare_gf.mesh2().points()) {
              int I = im / Ns, J = im % Ns;
              bare(I, J) = bare_gf(w, alps::gf::index_mesh::index_type(im),
                                    alps::gf::index_mesh::index_type(is));
              bold(I, J) = G_ij   (w, alps::gf::index_mesh::index_type(im),
                                    alps::gf::index_mesh::index_type(is));
            }
            Eigen::MatrixXcd sigm = bare.inverse() - bold.inverse();
            for (int im : bare_gf.mesh2().points()) {
              int I = im / Ns, J = im % Ns;
              sigma(w, alps::gf::index_mesh::index_type(im),
                    alps::gf::index_mesh::index_type(is)) = sigm(I, J);
            }
          }
        }
      }

      /// Legacy overload that fills an alps::gf::three_index_gf.
      template <typename Mesh>
      void bare_greens_function(
          alps::gf::three_index_gf<std::complex<double>, Mesh,
              alps::gf::index_mesh, alps::gf::index_mesh>& bare_gf,
          double beta) const {
        const int Ns = this->orbitals();
        const int nw = static_cast<int>(bare_gf.mesh1().points().size());
        const auto& t   = this->hopping_matrix();
        const auto& mu  = this->chem_potential();
        const auto& eps = this->site_energy();
        for (int iw = 0; iw < nw; ++iw) {
          typename Mesh::index_type w(iw);
          for (int is : bare_gf.mesh3().points()) {
            Eigen::MatrixXcd G_inv = Eigen::MatrixXcd::Zero(Ns, Ns);
            for (int I = 0; I < Ns; ++I) {
              G_inv(I, I) = EDLib::common::freq_point(iw, bare_gf.mesh1(), beta)
                          + double(mu[I]) - double(eps[I][is]);
              for (int J = 0; J < Ns; ++J) {
                G_inv(I, J) += double(t[I][J]);
              }
            }
            G_inv = G_inv.inverse().eval();
            for (int im : bare_gf.mesh2().points()) {
              int I = im / Ns;
              int J = im % Ns;
              bare_gf(w, alps::gf::index_mesh::index_type(im),
                      alps::gf::index_mesh::index_type(is)) = G_inv(I, J);
            }
          }
        }
      }

    public:
      static typename Base::ModelData read_model_data(alps::params& p) {
        typename Base::ModelData b;
        const int N = p["NSITES"].as<int>();
        b.hopping.assign(N, std::vector<Prec>(N, Prec(0)));
        b.U.assign(N, Prec(0));
        b.mu.assign(N, Prec(0));
        b.magnetic_field.assign(N, Prec(0));
        b.exchange.assign(N, std::vector<Prec>(N, Prec(0)));

        std::string input = p["INPUT_FILE"];
        alps::hdf5::archive ar(input, "r");
        if (ar.is_data("magnetic_field/values")) {
          ar >> alps::make_pvp("magnetic_field/values", b.magnetic_field);
        }
        ar >> alps::make_pvp("hopping/values",            b.hopping);
        ar >> alps::make_pvp("interaction/values",        b.U);
        if (ar.is_data("exchange/values")) {
          ar >> alps::make_pvp("exchange/values",         b.exchange);
        }
        ar >> alps::make_pvp("chemical_potential/values", b.mu);
        ar.close();
        return b;
      }
    };

  }
}

#endif
