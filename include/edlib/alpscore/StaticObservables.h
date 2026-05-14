#ifndef EDLIB_STATICOBSERVABLES_H
#define EDLIB_STATICOBSERVABLES_H

#include <string>

#include <alps/params.hpp>

#include "edlib/alpscore/EigenPair.h"
#include "edlib/alpscore/ParamsBridge.h"
#include "edlib/core/StaticObservables.h"

namespace EDLib {

  /**
   * Legacy alpscore-aware StaticObservables shim. Forwards to
   * edlib::StaticObservables. Static const std::string key names preserved
   * for callers that read them.
   */
  template <class Hamiltonian>
  class StaticObservables : public ::edlib::StaticObservables<Hamiltonian> {
    using Base = ::edlib::StaticObservables<Hamiltonian>;
  public:
    using typename Base::precision;
    using typename Base::sector;
    using Base::Base;

    explicit StaticObservables(alps::params& p)
        : Base(EDLib::compat::from_alps(p)) {}

    static const std::string _N_;
    static const std::string _N_UP_;
    static const std::string _N_DN_;
    static const std::string _M_;
    static const std::string _D_OCC_;
    static const std::string _N_EFF_;
    static const std::string _MI_MJ_;
    static const std::string _E_;
  };

  template <class H> const std::string StaticObservables<H>::_N_     = ::edlib::StaticObservables<H>::_N_;
  template <class H> const std::string StaticObservables<H>::_N_UP_  = ::edlib::StaticObservables<H>::_N_UP_;
  template <class H> const std::string StaticObservables<H>::_N_DN_  = ::edlib::StaticObservables<H>::_N_DN_;
  template <class H> const std::string StaticObservables<H>::_M_     = ::edlib::StaticObservables<H>::_M_;
  template <class H> const std::string StaticObservables<H>::_D_OCC_ = ::edlib::StaticObservables<H>::_D_OCC_;
  template <class H> const std::string StaticObservables<H>::_N_EFF_ = ::edlib::StaticObservables<H>::_N_EFF_;
  template <class H> const std::string StaticObservables<H>::_MI_MJ_ = ::edlib::StaticObservables<H>::_MI_MJ_;
  template <class H> const std::string StaticObservables<H>::_E_     = ::edlib::StaticObservables<H>::_E_;

}

#endif
