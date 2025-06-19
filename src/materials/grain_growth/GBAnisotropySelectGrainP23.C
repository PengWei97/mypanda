#include "GBAnisotropySelectGrainP23.h"

registerMooseObject("mypandaApp", GBAnisotropySelectGrainP23);

InputParameters
GBAnisotropySelectGrainP23::validParams()
{
  InputParameters params = GBAnisotropyMisoriAndTwin::validParams();
  params.addClassDescription("Apply grain boundary mobility scaling based on selected grain types.");
  params.addParam<bool>("is_select_grain_type", false, "Whether to apply scaling for selected grain types.");
  params.addRequiredParam<std::vector<Real>>("select_grain_type1", "Grain IDs for type 1 scaling.");
  params.addRequiredParam<std::vector<Real>>("select_grain_type2", "Grain IDs for type 2 scaling.");
  params.addRequiredParam<std::vector<Real>>("scaling_factors", "Two-element vector of scaling factors for type1 and type2 interactions.");
  return params;
}

GBAnisotropySelectGrainP23::GBAnisotropySelectGrainP23(const InputParameters & parameters)
  : GBAnisotropyMisoriAndTwin(parameters),
    _is_select_grain_type(getParam<bool>("is_select_grain_type")),
    _select_grain_type1(getParam<std::vector<Real>>("select_grain_type1")),
    _select_grain_type2(getParam<std::vector<Real>>("select_grain_type2")),
    _scaling_factors(getParam<std::vector<Real>>("scaling_factors"))
{
  mooseAssert(_scaling_factors.size() >= 2, "Scaling factors vector must have at least two elements.");
}

Real 
GBAnisotropySelectGrainP23::calculateGBEnergy(const MisorientationAngleData & misori_s)
{
  // Use base class implementation
  return GBAnisotropyMisoriAndTwin::calculateGBEnergy(misori_s);
}

Real 
GBAnisotropySelectGrainP23::calculateGBMobility(const MisorientationAngleData & misori_s)
{
  Real gbMob = GBAnisotropyMisoriAndTwin::calculateGBMobility(misori_s);

  if (!_is_select_grain_type)
    return gbMob;

  const bool i_in_type1 = isInList(_select_grain_type1, _current_grain_i);
  const bool j_in_type1 = isInList(_select_grain_type1, _current_grain_j);
  const bool i_in_type2 = isInList(_select_grain_type2, _current_grain_i);
  const bool j_in_type2 = isInList(_select_grain_type2, _current_grain_j);

  // 条件1: type1 只有一个
  if ((i_in_type1 && !j_in_type1) || (!i_in_type1 && j_in_type1))
    gbMob *= _scaling_factors[0];

  // 条件2: type2 只有一个
  else if ((i_in_type2 && !j_in_type2) || (!i_in_type2 && j_in_type2))
    gbMob *= _scaling_factors[1];

  // 条件3: 一个type1一个type2，保持gbMob不变
  return gbMob;
}

bool 
GBAnisotropySelectGrainP23::isInList(const std::vector<Real> & list, const Real val) const
{
  return std::find(list.begin(), list.end(), val) != list.end();
}
