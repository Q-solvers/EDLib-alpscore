#ifndef EDLIB_EXECUTIONSTATISTIC_H
#define EDLIB_EXECUTIONSTATISTIC_H

#include "edlib/core/ExecutionStatistic.h"

namespace EDLib {
  namespace common {
    using ExecutionStatistic = ::edlib::ExecutionStatistic;
    inline ::edlib::ExecutionStatistic& statistics = ::edlib::statistics;
  }
}

#endif
