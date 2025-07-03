#include "GBAnisotropy5SelectGrain.h"

registerMooseObject("mypandaApp", GBAnisotropy5SelectGrain);

InputParameters
GBAnisotropy5SelectGrain::validParams()
{
  InputParameters params = GBAnisotropyWtDeltaRho::validParams();
  params.addClassDescription("Grain boundary anisotropy with selective grain mobility adjustment");

  params.addParam<bool>("is_select_grains_for_high_mobility", false,
      "Flag to indicate if high mobility grains are selected");
  params.addParam<std::vector<unsigned int>>("select_grain_ids", {},
      "Grain IDs whose GB mobilities will be modified over time");
  params.addParam<std::vector<Real>>("amplification_factor", {},
      "Amplification factors corresponding to selected grain IDs");

  params.addParam<Real>("mobility_decay_rate", 2.0e-2,
      "Time decay factor applied to amplified mobility");
  return params;
}

GBAnisotropy5SelectGrain::GBAnisotropy5SelectGrain(const InputParameters & parameters)
  : GBAnisotropyWtDeltaRho(parameters),
    _is_select_grains_for_high_mobility(getParam<bool>("is_select_grains_for_high_mobility")),
    _mobility_decay_rate(getParam<Real>("mobility_decay_rate"))
{
  const auto & _selected_grain_ids = getParam<std::vector<unsigned int>>("select_grain_ids");
  const auto & _amplification_factors = getParam<std::vector<Real>>("amplification_factor");

  if (_selected_grain_ids.size() != _amplification_factors.size())
    mooseError("The size of 'select_grain_ids' and 'amplification_factor' must match.");

  for (std::size_t i = 0; i < _selected_grain_ids.size(); ++i)
    _selected_grain_mobility_factors[_selected_grain_ids[i]] = _amplification_factors[i];
}

Real 
GBAnisotropy5SelectGrain::calculateGBMobility(const MisorientationAngleData & misori_s)
{
  Real _gb_mobility = GBAnisotropyWtDeltaRho::calculateGBMobility(misori_s);

  if (!_is_select_grains_for_high_mobility)
    return _gb_mobility; // if not selecting grains, return the original mobility

  const Real _current_time = _fe_problem.time();
  auto _it_grain_i = _selected_grain_mobility_factors.find(_current_grain_i);
  auto _it_grain_j = _selected_grain_mobility_factors.find(_current_grain_j);

  const bool _is_grain_i_selected = (_it_grain_i != _selected_grain_mobility_factors.end());
  const bool _is_grain_j_selected = (_it_grain_j != _selected_grain_mobility_factors.end());

  if (_is_grain_i_selected != _is_grain_j_selected) // XOR：仅有一个晶粒在列表中
  {
    const Real _amplification_factor = _is_grain_i_selected
                                          ? _it_grain_i->second
                                          : _it_grain_j->second;

    const Real _adjustment_factor = (_amplification_factor - 1.0) *
                                      std::exp(-_mobility_decay_rate * _current_time) + 1.0;

    _gb_mobility *= _adjustment_factor;
  }

  return _gb_mobility;
}