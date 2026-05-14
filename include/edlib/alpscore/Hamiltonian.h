#ifndef EDLIB_ALPSCORE_HAMILTONIAN_H
#define EDLIB_ALPSCORE_HAMILTONIAN_H

#include <alps/params.hpp>

#include "edlib/alpscore/CRSStorage.h"
#include "edlib/alpscore/HubbardModel.h"
#include "edlib/alpscore/SOCRSStorage.h"
#include "edlib/alpscore/SingleImpurityAndersonModel.h"
#include "edlib/alpscore/SpinResolvedStorage.h"
#include "edlib/alpscore/ParamsBridge.h"
#include "edlib/Hamiltonian.h"

namespace EDLib {

  /**
   * Legacy alpscore-aware Hamiltonian shim. Delegates to edlib::Hamiltonian;
   * adds an alps::params ctor that reads the ModelData via the legacy model's
   * static read_model_data helper.
   */
  template <class Storage>
  class Hamiltonian : public ::edlib::Hamiltonian<Storage> {
    using Base = ::edlib::Hamiltonian<Storage>;
  public:
    using typename Base::Model;
    using typename Base::ModelType;
    using typename Base::StorageType;
    using typename Base::prec;
    using Base::Base;

#ifdef USE_MPI
    Hamiltonian(alps::params& p, MPI_Comm comm)
        : Base(EDLib::compat::from_alps(p), Model::read_model_data(p), comm) {}
#endif
    explicit Hamiltonian(alps::params& p)
        : Base(EDLib::compat::from_alps(p), Model::read_model_data(p)) {}
  };

  using CSRHubbardHamiltonian         = Hamiltonian<Storage::CRSStorage         <Model::HubbardModel<double>>>;
  using SRSHubbardHamiltonian         = Hamiltonian<Storage::SpinResolvedStorage<Model::HubbardModel<double>>>;
  using SOCSRHubbardHamiltonian       = Hamiltonian<Storage::SOCRSStorage       <Model::HubbardModel<double>>>;

  using CSRHubbardHamiltonian_float   = Hamiltonian<Storage::CRSStorage         <Model::HubbardModel<float>>>;
  using SRSHubbardHamiltonian_float   = Hamiltonian<Storage::SpinResolvedStorage<Model::HubbardModel<float>>>;
  using SOCSRHubbardHamiltonian_float = Hamiltonian<Storage::SOCRSStorage       <Model::HubbardModel<float>>>;

  using CSRSIAMHamiltonian            = Hamiltonian<Storage::CRSStorage         <Model::SingleImpurityAndersonModel<double>>>;
  using CSRSIAMHamiltonian_float      = Hamiltonian<Storage::CRSStorage         <Model::SingleImpurityAndersonModel<float>>>;
  using SRSSIAMHamiltonian            = Hamiltonian<Storage::SpinResolvedStorage<Model::SingleImpurityAndersonModel<double>>>;
  using SRSSIAMHamiltonian_float      = Hamiltonian<Storage::SpinResolvedStorage<Model::SingleImpurityAndersonModel<float>>>;

}

#endif
