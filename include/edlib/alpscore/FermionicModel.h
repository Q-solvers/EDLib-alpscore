#ifndef EDLIB_FERMIONICMODEL_H
#define EDLIB_FERMIONICMODEL_H

#include <complex>

#include <alps/gf/gf.hpp>
#include <alps/gf/mesh.hpp>
#include <alps/params.hpp>

#include <Eigen/Core>
#include <Eigen/LU>

#include "edlib/alpscore/ParamsBridge.h"
#include "edlib/core/CommonUtils.h"
#include "edlib/core/FermionicModel.h"

namespace EDLib {
  namespace Model {

    /**
     * Legacy alpscore-aware FermionicModel shim. Forwards to
     * edlib::FermionicModel; adds an alps::params ctor and an alps::gf
     * solve_dyson() for legacy GreensFunction callers (no Mesh bridge
     * needed since the matrix algebra is identical to edlib::solve_dyson).
     */
    class FermionicModel : public ::edlib::FermionicModel {
    public:
      using ::edlib::FermionicModel::FermionicModel;

      explicit FermionicModel(alps::params& p)
          : ::edlib::FermionicModel(EDLib::compat::from_alps(p)) {}

      template <typename Mesh>
      void solve_dyson(
          const alps::gf::three_index_gf<std::complex<double>, Mesh,
              alps::gf::index_mesh, alps::gf::index_mesh>& bare_gf,
          const alps::gf::three_index_gf<std::complex<double>, Mesh,
              alps::gf::index_mesh, alps::gf::index_mesh>& G_ij,
          alps::gf::three_index_gf<std::complex<double>, Mesh,
              alps::gf::index_mesh, alps::gf::index_mesh>& sigma) {
        for (int iw = 0; iw < (int)bare_gf.mesh1().points().size(); ++iw) {
          typename Mesh::index_type w(iw);
          for (int is : bare_gf.mesh3().points()) {
            Eigen::MatrixXcd bare(_Ns, _Ns);
            Eigen::MatrixXcd bold(_Ns, _Ns);
            Eigen::MatrixXcd sigm(_Ns, _Ns);
            for (int im : bare_gf.mesh2().points()) {
              int I = im / _Ns;
              int J = im % _Ns;
              bare(I, J) = bare_gf(w, alps::gf::index_mesh::index_type(im),
                                   alps::gf::index_mesh::index_type(is));
              bold(I, J) = G_ij(w, alps::gf::index_mesh::index_type(im),
                                alps::gf::index_mesh::index_type(is));
            }
            sigm = bare.inverse() - bold.inverse();
            for (int im : bare_gf.mesh2().points()) {
              int I = im / _Ns;
              int J = im % _Ns;
              sigma(w, alps::gf::index_mesh::index_type(im),
                    alps::gf::index_mesh::index_type(is)) = sigm(I, J);
            }
          }
        }
      }
    };

  }
}

#endif
