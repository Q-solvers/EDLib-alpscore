#ifndef EDLIB_EXT_ALPSCORE_HOLSTEINANDERSONMODEL_H
#define EDLIB_EXT_ALPSCORE_HOLSTEINANDERSONMODEL_H

#include <complex>
#include <string>
#include <vector>

#include <alps/gf/gf.hpp>
#include <alps/hdf5.hpp>
#include <alps/params.hpp>

#include <edlib/alpscore/CommonUtils.h>
#include <edlib/alpscore/SingleImpurityAndersonModel.h>
#include <edlib/alpscore/ParamsBridge.h>
#include <ext/alpscore/SzSymmetryWithBoson.h>
#include <ext/HolsteinAndersonModel.h>

namespace EDLib {
  namespace Ext {
    namespace Model {

      namespace HolsteinAnderson {
        template <typename Prec> using InnerState              = ::edlib::ext::holstein::InnerState<Prec>;
        template <typename Prec> using HybridisationInnerState = ::edlib::ext::holstein::HybridisationInnerState<Prec>;
        template <typename Prec> using BosonInnerState         = ::edlib::ext::holstein::BosonInnerState<Prec>;
      }

      /**
       * Legacy alpscore-aware HolsteinAndersonModel shim. Reads the full
       * fermionic+bosonic bath from INPUT_FILE and forwards to
       * edlib::ext::HolsteinAndersonModel.
       */
      template <typename Prec>
      class HolsteinAndersonModel : public ::edlib::ext::HolsteinAndersonModel<Prec> {
        using Base = ::edlib::ext::HolsteinAndersonModel<Prec>;
      public:
        using typename Base::precision;
        using typename Base::SYMMETRY;
        using typename Base::St;
        using typename Base::HSt;
        using typename Base::BSt;
        using typename Base::Sector;
        using Base::Base;

        explicit HolsteinAndersonModel(alps::params& p)
            : Base(EDLib::compat::from_alps(p), read_model_data(p)) {}

        template <typename Mesh>
        void bare_greens_function(
            alps::gf::three_index_gf<std::complex<double>, Mesh,
                alps::gf::index_mesh, alps::gf::index_mesh>& /*bare_gf*/,
            double /*beta*/) const {
          // Mirror of legacy: stub.
        }

        template <typename Mesh>
        void solve_dyson(
            const alps::gf::three_index_gf<std::complex<double>, Mesh,
                alps::gf::index_mesh, alps::gf::index_mesh>& /*bare_gf*/,
            const alps::gf::three_index_gf<std::complex<double>, Mesh,
                alps::gf::index_mesh, alps::gf::index_mesh>& /*G_ij*/,
            alps::gf::three_index_gf<std::complex<double>, Mesh,
                alps::gf::index_mesh, alps::gf::index_mesh>& /*sigma*/) const {
          // Not used by legacy ext path; provided for API parity.
        }

        static typename Base::ModelData read_model_data(alps::params& p) {
          typename Base::ModelData b;
          b.ep.nbbits  = p["NBBITS"].as<int>();
          b.ep.nblevel = p["NBLEVEL"].as<int>();
          b.ep.ml      = p["NORBITALS"].as<int>();

          std::string input = p["INPUT_FILE"];
          alps::hdf5::archive ar(input, "r");

          ar >> alps::make_pvp("ModelData/Vk/values",   b.Vk);
          ar >> alps::make_pvp("ModelData/Epsk/values", b.Epsk);
          ar >> alps::make_pvp("ModelData/w0/values",   b.w0);
          ar >> alps::make_pvp("ModelData/W/values",    b.W);
          ar >> alps::make_pvp("tk/values",        b.tk);
          ar >> alps::make_pvp("Eps0/values",      b.Eps);
          Prec mu = Prec(0);
          ar >> alps::make_pvp("mu", mu);
          b.mu = mu;
          ar >> alps::make_pvp("U", b.U);

          Prec avg = Prec(0);
          ar >> alps::make_pvp("AVG", avg);
          b.ep.avg = static_cast<double>(avg);

          ar.close();
          return b;
        }
      };

    }
  }
}

#endif
