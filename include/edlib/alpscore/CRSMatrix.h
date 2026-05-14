#ifndef EDLIB_CSRMATRIX_H
#define EDLIB_CSRMATRIX_H

#include "edlib/core/CRSMatrix.h"

namespace EDLib {
  namespace Storage {
    template <class Prec>
    using CRSMatrix = ::edlib::CRSMatrix<Prec>;
  }
}

#endif
