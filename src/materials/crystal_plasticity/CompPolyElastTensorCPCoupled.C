//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CompPolyElastTensorCPCoupled.h"
#include "RotationTensor.h"

registerMooseObject("mypandaApp", CompPolyElastTensorCPCoupled);

InputParameters
CompPolyElastTensorCPCoupled::validParams()
{
  InputParameters params = ComputeElasticityTensorBase::validParams();
  params.addClassDescription("Compute an elasticity tensor for crystal plasticity and phase field coupled model.");
  params.addRequiredParam<UserObjectName>(
      "grain_tracker", "Name of GrainTracker user object that provides RankFourTensors and euler angles");
  params.addRequiredCoupledVarWithAutoBuild(
      "v", "var_name_base", "op_num", "Array of coupled variables");
  params.addParam<Real>("length_scale", 1.0e-9, "Length scale of the problem, in meters");
  params.addParam<Real>("pressure_scale", 1.0e6, "Pressure scale of the problem, in pa");      
  return params;
}

CompPolyElastTensorCPCoupled::CompPolyElastTensorCPCoupled(const InputParameters & parameters)
  : ComputeElasticityTensorBase(parameters),
    _grain_tracker(getUserObject<GrainTrackerMateProp>("grain_tracker")),
    _op_num(coupledComponents("v")),
    _vals(coupledValues("v")),
    _elastic_tensor_gr(_op_num),
    _length_scale(getParam<Real>("length_scale")),
    _pressure_scale(getParam<Real>("pressure_scale")),
    _JtoeV(6.24150974e18),
    _Euler_angles_mat_prop_gr(_op_num),
    _crysrot_gr(_op_num),
    _D_elasticity_tensor_gr(_op_num)
{
  for (MooseIndex(_op_num) op_index = 0; op_index < _op_num; ++op_index)
  {
    _elastic_tensor_gr[op_index] = &declareProperty<RankFourTensor>(_base_name + 
        "elasticity_tensor_" + coupledName("v", op_index));
    _Euler_angles_mat_prop_gr[op_index] = &declareProperty<RealVectorValue>(
        "Euler_angles_" + coupledName("v", op_index));
    _crysrot_gr[op_index] = &declareProperty<RankTwoTensor>(_base_name + 
        "crysrot_" + coupledName("v", op_index));
    _D_elasticity_tensor_gr[op_index] = &declarePropertyDerivative<RankFourTensor>(_elasticity_tensor_name, coupledName("v", op_index));
  }
}

void
CompPolyElastTensorCPCoupled::initQpStatefulProperties()
{
  // Loop over variables (ops)
  for (MooseIndex(_op_num) op_index = 0; op_index < _op_num; ++op_index)
  {
    // Initialize Euler angles
    (*_Euler_angles_mat_prop_gr[op_index])[_qp] = RealVectorValue(0.0, 0.0, 0.0);

    // Initialize elasticity tensor
    RotationTensor R = RotationTensor((*_Euler_angles_mat_prop_gr[op_index])[_qp]);
    (*_crysrot_gr[op_index])[_qp] = R.transpose();

    (*_elastic_tensor_gr[op_index])[_qp].zero();
  }
}

void 
CompPolyElastTensorCPCoupled::computeQpElasticityTensor()
{
  // Get list of active order parameters from grain tracker
  const auto & op_to_grains = _grain_tracker.getVarToFeatureVector(_current_elem->id());

  Real sum_h = computeGrainContributions(op_to_grains);

  // Normalize the elasticity tensor
  const Real tol = 1.0e-10;
  sum_h = std::max(sum_h, tol);
  _elasticity_tensor[_qp] /= sum_h;

  // Update elasticity tensor derivatives
  updateElasticityTensorDerivatives(op_to_grains, sum_h);
}

Real 
CompPolyElastTensorCPCoupled::computeGrainContributions(const std::vector<unsigned int> & op_to_grains)
{
  Real sum_h = 0.0;
  for (MooseIndex(op_to_grains) op_index = 0; op_index < op_to_grains.size(); ++op_index)
  {
    auto grain_id = op_to_grains[op_index];
    if (grain_id == FeatureFloodCount::invalid_id)
      continue;

    // Interpolation factor for elasticity tensors
    Real h = (1.0 + std::sin(libMesh::pi * ((*_vals[op_index])[_qp] - 0.5))) / 2.0;

    // get Euler angles
    (*_Euler_angles_mat_prop_gr[op_index])[_qp] = _grain_tracker.updateEulerAngles(grain_id);

    // Initialize elasticity tensor
    RotationTensor R = RotationTensor((*_Euler_angles_mat_prop_gr[op_index])[_qp]);
    (*_crysrot_gr[op_index])[_qp] = R.transpose();

    (*_elastic_tensor_gr[op_index])[_qp] = _grain_tracker.getData(grain_id);

    sum_h += h;
  }
  return sum_h;
}

void CompPolyElastTensorCPCoupled::updateElasticityTensorDerivatives(const std::vector<unsigned int>& op_to_grains, Real sum_h)
{
  for (MooseIndex(_op_num) op_index = 0; op_index < _op_num; ++op_index)
    (*_D_elasticity_tensor_gr[op_index])[_qp].zero();

  for (MooseIndex(op_to_grains) op_index = 0; op_index < op_to_grains.size(); ++op_index)
  {
    auto grain_id = op_to_grains[op_index];
    if (grain_id == FeatureFloodCount::invalid_id)
      continue;

    Real dhdopi = libMesh::pi * std::cos(libMesh::pi * ((*_vals[op_index])[_qp] - 0.5)) / 2.0;
    RankFourTensor & C_deriv = (*_D_elasticity_tensor_gr[op_index])[_qp];

    C_deriv = (_grain_tracker.getData(grain_id) - _elasticity_tensor[_qp]) * dhdopi / sum_h;

    // Convert from XPa to eV/(xm)^3, where X is pressure scale and x is length scale
    C_deriv *= _JtoeV * (_length_scale * _length_scale * _length_scale) * _pressure_scale;
  }
}

