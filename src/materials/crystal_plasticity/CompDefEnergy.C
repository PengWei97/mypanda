#include "CompDefEnergy.h"

registerMooseObject("mypandaApp", CompDefEnergy);

InputParameters
CompDefEnergy::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("This material is used to compute the weigthed material properties for crystal plasticity model.");
  params.addParam<std::string>(
      "base_name",
      "Optional parameter that allows the user to define multiple crystal plasticity mechanisms");
  params.addRequiredParam<UserObjectName>(
      "grain_tracker", "Name of GrainTracker user object that provides RankFourTensors and euler angles");
  params.addRequiredCoupledVarWithAutoBuild(
      "v", "var_name_base", "op_num", "Array of coupled variables");
  params.addParam<Real>("length_scale", 1.0e-9, "Length scale of the problem, in meters");
  params.addParam<Real>("pressure_scale", 1.0e6, "Pressure scale of the problem, in pa");
  return params;
}

CompDefEnergy::CompDefEnergy(const InputParameters & parameters)
  : Material(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _grain_tracker(getUserObject<GrainTrackerMateProp>("grain_tracker")),
    _op_num(coupledComponents("v")),
    _vals(coupledValues("v")),
    _elastic_energy_gr(_op_num),
    _plastic_energy_gr(_op_num),
    _elastic_energy(declareProperty<Real>("elastic_energy")),
    _plastic_energy(declareProperty<Real>("plastic_energy")),
    _D_elastic_energy_gr(_op_num),
    _D_plastic_energy_gr(_op_num),
    _JtoeV(6.24150974e18),
    _length_scale(getParam<Real>("length_scale")),
    _pressure_scale(getParam<Real>("pressure_scale"))
{
  for (unsigned int i = 0; i < _op_num; i++)
  {
    _elastic_energy_gr[i] = &getMaterialProperty<Real>(_base_name + "elastic_energy_gr" + Moose::stringify(i));
    _plastic_energy_gr[i] = &getMaterialProperty<Real>(_base_name + "plastic_energy_gr" + Moose::stringify(i));

    MaterialPropertyName D_elastic_energy_name = derivativePropertyNameFirst("elastic_energy", coupledName("v", i));
    _D_elastic_energy_gr[i] = &declareProperty<Real>(D_elastic_energy_name);

    MaterialPropertyName D_plastic_energy_name = derivativePropertyNameFirst("plastic_energy", coupledName("v", i));
    _D_plastic_energy_gr[i] = &declareProperty<Real>(D_plastic_energy_name);
  }
}

void
CompDefEnergy::computeQpProperties()
{
  const auto & op_to_grains = _grain_tracker.getVarToFeatureVector(_current_elem->id());

  Real sum_h = 0.0;
  _elastic_energy[_qp] = 0.0;
  _plastic_energy[_qp] = 0.0;
  for (MooseIndex(op_to_grains) op_index = 0; op_index < op_to_grains.size(); ++op_index)
  {
    auto grain_id = op_to_grains[op_index];
    if (grain_id == FeatureFloodCount::invalid_id)
      continue;

    // Interpolation factor for material properties
    Real h = (1.0 + std::sin(libMesh::pi * ((*_vals[op_index])[_qp] - 0.5))) / 2.0;

    // Calculate the elastic energy
    _elastic_energy[_qp] += (*_elastic_energy_gr[op_index])[_qp] * h;
    _plastic_energy[_qp] += (*_plastic_energy_gr[op_index])[_qp] * h;

    sum_h += h;
  }

  const Real tol = 1.0e-10;
  sum_h = std::max(sum_h, tol);
  _elastic_energy[_qp] /= sum_h;
  _plastic_energy[_qp] /= sum_h;

  // initial elastic energy derivative with respect to the order parameter
  for (MooseIndex(_op_num) op_index = 0; op_index < _op_num; ++op_index)
  {
    (*_D_elastic_energy_gr[op_index])[_qp] = 0.0;
    (*_D_plastic_energy_gr[op_index])[_qp] = 0.0;
  }

  for (MooseIndex(op_to_grains) op_index = 0; op_index < op_to_grains.size(); ++op_index)
  {
    auto grain_id = op_to_grains[op_index];
    if (grain_id == FeatureFloodCount::invalid_id)
      continue;

    Real dhdopi = libMesh::pi * std::cos(libMesh::pi * ((*_vals[op_index])[_qp] - 0.5)) / 2.0;
    Real & fe_deriv = (*_D_elastic_energy_gr[op_index])[_qp];
    Real & fp_deriv = (*_D_plastic_energy_gr[op_index])[_qp];

    fe_deriv = ((*_elastic_energy_gr[op_index])[_qp] - _elastic_energy[_qp]) * dhdopi / sum_h;
    fp_deriv = ((*_plastic_energy_gr[op_index])[_qp] - _plastic_energy[_qp]) * dhdopi / sum_h;

    // Convert from XPa to eV/(xm)^3, where X is pressure scale and x is length scale;
    fe_deriv *= _JtoeV * (_length_scale * _length_scale * _length_scale) * _pressure_scale;    
  }
}