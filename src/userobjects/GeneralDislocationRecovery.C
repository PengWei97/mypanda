#include "GeneralDislocationRecovery.h"

registerMooseObject("mypandaApp", GeneralDislocationRecovery);

InputParameters
GeneralDislocationRecovery::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params.addParam<bool>("is_concurrent_recovery", false, "Enable concurrent recovery model");
  params.addParam<Real>("a_rho1", 3.3e-3, "Evolution coefficient during medium time recovery");
  params.addParam<Real>("rho_end1", 2.10e12, "Dislocation density after long-term concurrent recovery");
  params.addParam<Real>("a_rho2", 3.3e-3, "Evolution coefficient during medium time recovery");
  params.addParam<Real>("rho_end2", 2.10e12, "Dislocation density after long-term concurrent recovery");
  params.addParam<Real>("rho_default", 2.0e15, "Dislocation density at the beginning of simulation");
  params.addRequiredParam<UserObjectName>("GNDs_provider", "GNDs provider for EBSD reader");
  return params;
}

GeneralDislocationRecovery::GeneralDislocationRecovery(const InputParameters & parameters)
  : GeneralUserObject(parameters),
    _is_concurrent_recovery(getParam<bool>("is_concurrent_recovery")),
    _rho_end1(getParam<Real>("rho_end1")),
    _rho_end2(getParam<Real>("rho_end2")),
    _a_rho1(getParam<Real>("a_rho1")),
    _a_rho2(getParam<Real>("a_rho2")),
    _rho_default(getParam<Real>("rho_default")),
    _GNDs_provider(getUserObject<EBSDReaderMaterialProperty>("GNDs_provider"))
{
}

Real
GeneralDislocationRecovery::getRhoWtTime(const unsigned int & grain_id, const Real & y_coord) const
{
  const Real time = _fe_problem.time(); // Get the current simulation time
  const Real rho_init = _GNDs_provider.getRhoInit(grain_id);
  Real rho_now = rho_init;

  if (!_is_concurrent_recovery)
    return rho_init;

  // High-y grain → Type 2
  if (y_coord > 100.0 && rho_init > 1.0e13)
  {
    if (rho_init > _rho_end2)
      rho_now = (rho_init - _rho_end2) * std::exp(-_a_rho2 * time) + _rho_end2;
    else
      rho_now = _rho_end2;
    return std::clamp(rho_now, _rho_end2, _rho_default);
  }

  // Normal grain → Type 1
  if (rho_init > _rho_end1)
    rho_now = (rho_init - _rho_end1) * std::exp(-_a_rho1 * time) + _rho_end1;
  else
    rho_now = _rho_end1;

  return std::clamp(rho_now, _rho_end1, _rho_default);
}

Real
GeneralDislocationRecovery::getRhoWtTime(const unsigned int & grain_id) const
{
  const Real time = _fe_problem.time(); // Get the current simulation time
  const Real rho_init = _GNDs_provider.getRhoInit(grain_id);
  Real rho_now = rho_init;

  if (!_is_concurrent_recovery)
    return rho_init;

  // Normal grain → Type 1
  if (rho_init > _rho_end1)
    rho_now = (rho_init - _rho_end1) * std::exp(-_a_rho1 * time) + _rho_end1;
  else
    rho_now = _rho_end1;

  return std::clamp(rho_now, _rho_end1, _rho_default);
}
