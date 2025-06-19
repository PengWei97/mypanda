#pragma once

#include "DerivativeMaterialInterface.h"
#include "GrainTrackerInterface.h"
#include "EBSDReaderMaterialProperty.h"

class DeformedGrainEBSDRho : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();

  DeformedGrainEBSDRho(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  const unsigned int _op_num;
  const std::vector<const VariableValue *> _vals;

  const Real _length_scale;
  const Real _time_scale;
  const Real _Elas_Mod;
  const Real _Burg_vec;
  const Real _stored_factor;
  const Real _execution_time;
  const Real _JtoeV;
  Real _beta;

  const GrainTrackerInterface & _grain_tracker;
  const EBSDReaderMaterialProperty & _GNDs_provider;

  MaterialProperty<Real> & _rho_eff;
  MaterialProperty<Real> & _feature_id;
  MaterialProperty<Real> & _grain_type_rho;
  std::vector<MaterialProperty<Real> *> _D_stored_energy;
};
