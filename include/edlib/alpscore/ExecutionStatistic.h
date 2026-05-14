#ifndef EDLIB_ALPSCORE_EXECUTIONSTATISTIC_H
#define EDLIB_ALPSCORE_EXECUTIONSTATISTIC_H

#include "edlib/ExecutionStatistic.h"

namespace EDLib {
  namespace common {
    using ExecutionStatistic = ::edlib::ExecutionStatistic;
    inline ::edlib::ExecutionStatistic& statistics = ::edlib::statistics;
  }
}

#endif
