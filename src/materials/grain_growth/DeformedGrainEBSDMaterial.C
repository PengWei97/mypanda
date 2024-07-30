//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DeformedGrainEBSDMaterial.h"
#include "GrainTrackerInterface.h"

// Register the object with the application
registerMooseObject("mypandaApp", DeformedGrainEBSDMaterial);

// Define valid parameters for the class
InputParameters
DeformedGrainEBSDMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredCoupledVarWithAutoBuild("v", "var_name_base", "op_num", "Array of coupled variables");
  params.addParam<Real>("length_scale", 1.0e-9, "Length scale in meters (default is nm)");
  params.addParam<Real>("time_scale", 1.0e-6, "Time scale in seconds (default is microseconds)");
  params.addParam<Real>("Elas_Mod", 2.50e10, "Shear modulus in J/m^3");
  params.addParam<Real>("Burg_vec", 3.0e-10, "Length of Burgers vector in meters");
  params.addParam<Real>("stored_factor", 0.5, "Scaling factor in stored energy function");

  params.addParam<bool>("concurrent_recovery", false, "Consider concurrent recovery if true");
  params.addParam<Real>("rho_end_l2", 2.10e12, "Dislocation density after long-term concurrent recovery at level 2");
  params.addParam<Real>("rho_end_l3", 1.8e13, "Dislocation density after long-term concurrent recovery at level 3");
  params.addParam<Real>("rho_critical", 3.9e13, "Critical dislocation density of grains");
  params.addParam<Real>("a_rho_l2", 4.6e-4, "Evolution coefficient during medium time recovery at level 2");
  params.addParam<Real>("a_rho_l3", 6.0e-4, "Evolution coefficient during medium time recovery at level 3");

  params.addRequiredParam<UserObjectName>("grain_tracker", "GrainTracker UserObject to get values from");
  params.addRequiredParam<UserObjectName>("GNDs_provider", "GNDs provider for EBSD reader");
  return params;
}

DeformedGrainEBSDMaterial::DeformedGrainEBSDMaterial(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _op_num(coupledComponents("v")),
    _vals(coupledValues("v")),
    _length_scale(getParam<Real>("length_scale")),
    _time_scale(getParam<Real>("time_scale")),
    _Elas_Mod(getParam<Real>("Elas_Mod")),
    _Burg_vec(getParam<Real>("Burg_vec")),
    _stored_factor(getParam<Real>("stored_factor")),
    _JtoeV(6.24150974e18),

    _concurrent_recovery(getParam<bool>("concurrent_recovery")),
    _rho_end_l2(getParam<Real>("rho_end_l2")),
    _rho_end_l3(getParam<Real>("rho_end_l3")),
    _rho_critical(getParam<Real>("rho_critical")),
    _a_rho_l2(getParam<Real>("a_rho_l2")),
    _a_rho_l3(getParam<Real>("a_rho_l3")),

    _beta(declareProperty<Real>("beta")),
    _rho_eff(declareProperty<Real>("rho_eff")),
    _sumEtai2(declareProperty<Real>("sumEtai2")),
    _D_stored_energy(_op_num),
    
    _grain_tracker(getUserObject<GrainTrackerInterface>("grain_tracker")),
    _GNDs_provider(getUserObject<EBSDReader>("GNDs_provider"))
{
  if (_op_num == 0)
    paramError("op_num", "Model requires op_num > 0");

  // Loop over variables (ops)
  for (MooseIndex(_op_num) op_index = 0; op_index < _op_num; ++op_index)
    _D_stored_energy[op_index] = &declarePropertyDerivative<Real>(
      "stored_energy", coupledName("v", op_index));
}

void
DeformedGrainEBSDMaterial::computeQpProperties()
{
  Real SumEtai2 = 0.0;
  for (unsigned int i = 0; i < _op_num; ++i)
    SumEtai2 += (*_vals[i])[_qp] * (*_vals[i])[_qp];

  const auto & op_to_grains = _grain_tracker.getVarToFeatureVector(_current_elem->id());

  _rho_eff[_qp] = 0.0;
  for (MooseIndex(op_to_grains) op_index = 0; op_index < op_to_grains.size(); ++op_index)
  {
    auto grain_id = op_to_grains[op_index];
    if (op_to_grains[op_index] == FeatureFloodCount::invalid_id)
      continue;

    Real rho_i = getGNDsFromEBSD(grain_id);
    _rho_eff[_qp] += rho_i * (*_vals[op_index])[_qp] * (*_vals[op_index])[_qp];
  }

  _rho_eff[_qp] /= SumEtai2;
  _beta[_qp] = _stored_factor * _Elas_Mod * _Burg_vec * _Burg_vec * _JtoeV * _length_scale;

  // Calculate stored energy derivative
  for (MooseIndex(op_to_grains) op_index = 0; op_index < op_to_grains.size(); ++op_index)
  {
    Real & C_deriv = (*_D_stored_energy[op_index])[_qp];
    C_deriv = 0.0;

    auto grain_id = op_to_grains[op_index];
    if (grain_id == FeatureFloodCount::invalid_id)
      continue;

    C_deriv = (getGNDsFromEBSD(grain_id) - _rho_eff[_qp]) * _beta[_qp] / SumEtai2;
    C_deriv *= _length_scale * _length_scale; // Convert units
  }

  _sumEtai2[_qp] = SumEtai2;
  _rho_eff[_qp] *= _length_scale * _length_scale; // Convert units
}

Real
DeformedGrainEBSDMaterial::getGNDsFromEBSD(const unsigned int & grain_id)
{
  const auto & time_current = _fe_problem.time(); // Current simulation time in seconds
  Real rho_init = 2.0e15; // Initial dislocation density 1/m^2

  if (grain_id < _GNDs_provider.getGrainNum())
    rho_init = _GNDs_provider.getAvgData(grain_id)._custom[0]; // GNDs for each grain, 1/m^2
  
  Real rho_i = rho_init;
  if (_concurrent_recovery)
  {
    if (rho_init >= _rho_end_l3)
      rho_i = (rho_init - _rho_end_l3) * std::exp(-_a_rho_l3 * time_current) + _rho_end_l3; // Level 3 recovery
    else if (rho_init > _rho_end_l2 && rho_init < _rho_critical)
      rho_i = (rho_init - _rho_end_l2) * std::exp(-_a_rho_l2 * time_current) + _rho_end_l2; // Level 2 recovery
  }

  return rho_i;
}
