#ifndef EDLIB_ALPSCORE_CSRMATRIX_H
#define EDLIB_ALPSCORE_CSRMATRIX_H

#include "edlib/CRSMatrix.h"

namespace EDLib {
  namespace Storage {
    template <class Prec>
    using CRSMatrix = ::edlib::CRSMatrix<Prec>;
  }
}

#endif
