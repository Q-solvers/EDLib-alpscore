#ifndef EDLIB_ALPSCORE_COMPAT_GFBRIDGE_H
#define EDLIB_ALPSCORE_COMPAT_GFBRIDGE_H

#ifdef EDLIB_WITH_ALPSCORE

#include <complex>

#include <alps/gf/gf.hpp>
#include <alps/gf/mesh.hpp>

#include "edlib/alpscore/MeshBridge.h"
#include "edlib/Gf.h"
#include "edlib/Mesh.h"

namespace EDLib {
  namespace compat {

    /**
     * Convert an edlib::Gf<complex,3> on a frequency × index × index layout
     * into an alps::gf::three_index_gf. The supplied edlib mesh is bridged
     * to its alps equivalent for the first axis.
     */
    template <class EdlibMesh>
    auto to_alps_gf3(const edlib::Gf<std::complex<double>, 3>& g,
                     const EdlibMesh& mesh) {
      using AlpsMesh = decltype(to_alps_mesh(mesh));
      alps::gf::three_index_gf<std::complex<double>, AlpsMesh,
                               alps::gf::index_mesh, alps::gf::index_mesh>
          out(to_alps_mesh(mesh),
              alps::gf::index_mesh(g.shape(1)),
              alps::gf::index_mesh(g.shape(2)));
      for (int iw = 0; iw < g.shape(0); ++iw) {
        for (int i = 0; i < g.shape(1); ++i) {
          for (int j = 0; j < g.shape(2); ++j) {
            out(typename AlpsMesh::index_type(iw),
                alps::gf::index_mesh::index_type(i),
                alps::gf::index_mesh::index_type(j)) = g(iw, i, j);
          }
        }
      }
      return out;
    }

    /**
     * Convert an edlib::Gf<complex,2> (frequency × index) into an
     * alps::gf::two_index_gf<..., index_mesh>.
     */
    template <class EdlibMesh>
    auto to_alps_gf2(const edlib::Gf<std::complex<double>, 2>& g,
                     const EdlibMesh& mesh) {
      using AlpsMesh = decltype(to_alps_mesh(mesh));
      alps::gf::two_index_gf<std::complex<double>, AlpsMesh, alps::gf::index_mesh>
          out(to_alps_mesh(mesh), alps::gf::index_mesh(g.shape(1)));
      for (int iw = 0; iw < g.shape(0); ++iw) {
        for (int i = 0; i < g.shape(1); ++i) {
          out(typename AlpsMesh::index_type(iw),
              alps::gf::index_mesh::index_type(i)) = g(iw, i);
        }
      }
      return out;
    }

  }
}

#endif // EDLIB_WITH_ALPSCORE
#endif
