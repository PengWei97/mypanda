#include "EBSDReaderMaterialProperty.h"

registerMooseObject("mypandaApp", EBSDReaderMaterialProperty);

InputParameters
EBSDReaderMaterialProperty::validParams()
{
  InputParameters params = EBSDReader::validParams();
  params.addClassDescription("EBSD reader for material property");
  params.addParam<bool>("is_concurrent_recovery", false, "Enable concurrent recovery mechanism.");
  params.addParam<Real>("rho_end", 2.10e12, "Dislocation density after long-term concurrent recovery");
  params.addParam<Real>("a_rho", 4.6e-4, "Evolution coefficient during medium time recovery");
  params.addParam<Real>("rho_default", 2.0e15, "Dislocation density at the beginning of simulation");
  return params;
}

EBSDReaderMaterialProperty::EBSDReaderMaterialProperty(const InputParameters & parameters)
  : EBSDReader(parameters),
    _is_concurrent_recovery(getParam<bool>("is_concurrent_recovery")),
    _rho_end(getParam<Real>("rho_end")),
    _a_rho(getParam<Real>("a_rho")),
    _rho_default(getParam<Real>("rho_default"))
{
}

Real
EBSDReaderMaterialProperty::getRhoInit(unsigned int grain_id) const
{
  // Return initial dislocation density for a given grain
  if (_custom_columns > 0 && grain_id < getGrainNum())
    return getAvgData(grain_id)._custom[0];
  else
    return _rho_default;

  return _rho_default;
}

Real
EBSDReaderMaterialProperty::getRhoWtTime(unsigned int grain_id) const
{
  // Step 1: Get initial dislocation density
  Real rho = getRhoInit(grain_id);

  // Step 2: Apply concurrent recovery if enabled
  if (_is_concurrent_recovery)
  {
    const Real time = _fe_problem.time();

    if (rho > _rho_end)
      rho = (rho - _rho_end) * std::exp(-_a_rho * time) + _rho_end;
    else
      rho = _rho_end;  // Already recovered
  }

  // Step 3: Clamp within physical limits and return
  return std::clamp(rho, _rho_end, _rho_default);
}