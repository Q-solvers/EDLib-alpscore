#ifndef EDLIB_EIGENPAIR_H
#define EDLIB_EIGENPAIR_H

#include "edlib/core/EigenPair.h"

namespace EDLib {
  template <class Precision, class SectorType>
  using EigenPair = ::edlib::EigenPair<Precision, SectorType>;
}

#endif
