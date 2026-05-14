#ifndef EDLIB_ALPSCORE_SINGLEIMPURITYANDERSONMODEL_H
#define EDLIB_ALPSCORE_SINGLEIMPURITYANDERSONMODEL_H

#include <complex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <alps/gf/gf.hpp>
#include <alps/gf/mesh.hpp>
#include <alps/hdf5.hpp>
#include <alps/numeric/tensors.hpp>
#include <alps/params.hpp>

#include "edlib/alpscore/CommonUtils.h"
#include "edlib/alpscore/ParamsBridge.h"
#include "edlib/SingleImpurityAndersonModel.h"
#include "edlib/alpscore/FermionicModel.h"

namespace EDLib {
  namespace Model {

    namespace SingleImpurityAnderson {
      template <typename Prec>
      using InnerState              = ::edlib::siam::InnerState<Prec>;
      template <typename Prec>
      using InnerHybridizationState = ::edlib::siam::HybridisationInnerState<Prec>;
      template <typename Prec>
      using InnerInteractionState   = ::edlib::siam::InteractionInnerState<Prec>;
    }

    /**
     * Legacy alpscore-aware SIAM shim. Reads Vk/H0/Epsk/U from INPUT_FILE
     * and forwards to edlib::SingleImpurityAndersonModel. Exposes an
     * alps::gf-taking bare_greens_function for legacy GF callers.
     */
    template <typename Prec>
    class SingleImpurityAndersonModel : public ::edlib::SingleImpurityAndersonModel<Prec> {
      using Base = ::edlib::SingleImpurityAndersonModel<Prec>;
    public:
      using typename Base::precision;
      using typename Base::SYMMETRY;
      using typename Base::St;
      using typename Base::HSt;
      using typename Base::USt;
      using typename Base::Sector;
      using Base::Base;

      explicit SingleImpurityAndersonModel(alps::params& p)
          : Base(EDLib::compat::from_alps(p), read_model_data(p)) {}

      template <typename Mesh>
      void bare_greens_function(
          alps::gf::three_index_gf<std::complex<double>, Mesh,
              alps::gf::index_mesh, alps::gf::index_mesh>& bare_gf,
          double beta) const {
        const int nw = static_cast<int>(bare_gf.mesh1().points().size());
        const auto& Vk   = this->hybridisation();
        const auto& Epsk = this->bath_energies();
        const auto& H0   = this->H0();
        for (int iw = 0; iw < nw; ++iw) {
          typename Mesh::index_type w(iw);
          for (int im : bare_gf.mesh2().points()) {
            for (int is : bare_gf.mesh3().points()) {
              std::complex<double> delta(0, 0);
              for (int ik = 0; ik < static_cast<int>(Epsk.size()); ++ik) {
                delta += double(Vk[im][ik][is]) * double(Vk[im][ik][is])
                       / (EDLib::common::freq_point(iw, bare_gf.mesh1(), beta)
                          - double(Epsk[ik][is]));
              }
              bare_gf(w, alps::gf::index_mesh::index_type(im),
                      alps::gf::index_mesh::index_type(is))
                  = 1.0 / (EDLib::common::freq_point(iw, bare_gf.mesh1(), beta)
                           - double(H0[im][im][is]) - delta);
            }
          }
        }
      }

      static typename Base::ModelData read_model_data(alps::params& p) {
        typename Base::ModelData b;
        const int Ns = p["NSITES"].as<int>();
        const int ml = p["siam.NORBITALS"].as<int>();
        const int ms = p["NSPINS"].as<int>();
        if (ml > Ns) {
          throw std::invalid_argument(
              "SingleImpurityAndersonModel: siam.NORBITALS exceeds NSITES");
        }
        if (ms != 2) {
          throw std::invalid_argument("SingleImpurityAndersonModel: NSPINS must be 2");
        }

        b.Vk.assign(ml, std::vector<std::vector<Prec>>());
        b.H0.assign(ml, std::vector<std::vector<Prec>>(ml, std::vector<Prec>(ms, Prec(0))));

        std::string input = p["INPUT_FILE"];
        alps::hdf5::archive ar(input, "r");
        ar >> alps::make_pvp("ModelData/Epsk/values", b.Epsk);
        if ((int)b.Epsk.size() != Ns - ml) {
          throw std::invalid_argument(
              "SingleImpurityAndersonModel: bath levels + impurity orbitals != NSITES");
        }
        for (int im = 0; im < ml; ++im) {
          std::stringstream s;
          s << "ModelData/Vk_" << im << "/values";
          ar >> alps::make_pvp(s.str(), b.Vk[im]);
          s.str("");
          s << "H0_" << im << "/values";
          ar >> alps::make_pvp(s.str(), b.H0[im]);
        }
        Prec mu = Prec(0);
        ar >> alps::make_pvp("mu", mu);
        b.mu = mu;

        // Read U as an alps::numerics::tensor<Prec, 6>, copy into edlib::Gf<Prec,6>.
        alps::numerics::tensor<Prec, 6> U_alps;
        ar >> alps::make_pvp("interaction/values", U_alps);
        ar.close();
        if ((int)U_alps.shape()[2] != ml) {
          throw std::invalid_argument("SingleImpurityAndersonModel: U shape mismatch");
        }
        b.U = ::edlib::Gf<Prec, 6>({ms, ms, ml, ml, ml, ml});
        for (int s1 = 0; s1 < ms; ++s1)
          for (int s2 = 0; s2 < ms; ++s2)
            for (int i = 0; i < ml; ++i)
              for (int j = 0; j < ml; ++j)
                for (int k = 0; k < ml; ++k)
                  for (int l = 0; l < ml; ++l)
                    b.U(s1, s2, i, j, k, l) = U_alps(s1, s2, i, j, k, l);

        return b;
      }
    };

  }
}

#endif
