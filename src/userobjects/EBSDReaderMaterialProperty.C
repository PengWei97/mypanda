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
EBSDReaderMaterialProperty::getRhoWtTime(unsigned int grain_id) const
{
  // Initialize rho with the default value
  Real rho = _rho_default;

  // Check if custom columns are available and grain_id is valid
  if (_custom_columns > 0 && grain_id < getGrainNum())
    rho = getAvgData(grain_id)._custom[0];

  // Apply concurrent recovery mechanism if enabled
  if (_is_concurrent_recovery)
  {
    const Real time_current = _fe_problem.time();
    
    // Exponential decay if initial value is greater than the lower bound
    if (rho > _rho_end)
      rho = (rho - _rho_end) * std::exp(-_a_rho * time_current) + _rho_end;
    else
      rho = _rho_end;
  }

  // Clamp rho within the specified bounds
  return std::clamp(rho, _rho_end, _rho_default);
}

