#include "DeformedGrainEBSDRho_ExtendR2.h"

// Register the object with the application
registerMooseObject("mypandaApp", DeformedGrainEBSDRho_ExtendR2);

// Define valid parameters for the class
InputParameters
DeformedGrainEBSDRho_ExtendR2::validParams()
{
  InputParameters params = DeformedGrainEBSDRho::validParams();
  params.addRequiredParam<std::vector<unsigned int>>("feature_ids", "Grain IDs for selected grains in EBSD");
  params.addRequiredParam<std::vector<Real>>("set_rhos", "Set rhos for the selected grains");
  params.addParam<Real>("threshold_DeltaRho", 20.0, "Threshold for delta_rho to apply the set_rhos");

  params.addParam<bool>("enable_concurrent_recovery", true, "Enable exponential concurrent recovery evolution");
  params.addParam<Real>("rho_default", 2.0e15, "Initial dislocation density");
  params.addParam<Real>("rho_end", 1.0e12, "Target dislocation density at long time");
  params.addParam<Real>("a_rho", 3.3e-4, "Evolution rate coefficient for recovery");
  
  return params;
}

// Constructor
DeformedGrainEBSDRho_ExtendR2::DeformedGrainEBSDRho_ExtendR2(const InputParameters & parameters)
  : DeformedGrainEBSDRho(parameters),
  _feature_ids(getParam<std::vector<unsigned int>>("feature_ids")),
  _set_rhos(getParam<std::vector<Real>>("set_rhos")),
  _threshold_DeltaRho(getParam<Real>("threshold_DeltaRho")),
  _enable_concurrent_recovery(getParam<bool>("enable_concurrent_recovery")),
  _rho_default(getParam<Real>("rho_default")),
  _rho_end(getParam<Real>("rho_end")),
  _a_rho(getParam<Real>("a_rho")),
  _delta_rho(getMaterialProperty<Real>("delta_rho"))
{
  // Sanity check: feature_ids and set_rhos must have the same length
  if (_feature_ids.size() != _set_rhos.size())
  mooseError("Length mismatch: 'feature_ids' and 'set_rhos' must have the same size.");

  // Construct a map for fast lookup
  for (std::size_t i = 0; i < _feature_ids.size(); ++i)
    _rho_map[_feature_ids[i]] = _set_rhos[i];
}

Real 
DeformedGrainEBSDRho_ExtendR2::getRhoWtTime(const unsigned int & grain_id) const
{
  auto feature_id = _GNDs_provider.getFeatureID(grain_id);

  // Check whether feature_id is in the list
  auto it = _rho_map.find(feature_id);
  if (it != _rho_map.end() && _delta_rho[_qp] > _threshold_DeltaRho)
    return computeRhoWithRecovery(it->second);

  // Get the time-evolved dislocation density for a specified grain
  return _GNDs_provider.getRhoWtTime(grain_id);
}

Real
DeformedGrainEBSDRho_ExtendR2::computeRhoWithRecovery(const Real & rho_init) const
{
  Real rho = rho_init;

  if (_enable_concurrent_recovery && rho > _rho_end)
  {
    const Real time_current = _fe_problem.time();
    rho = (_rho_default - _rho_end) * std::exp(-_a_rho * time_current) + _rho_end;
  }

  return std::clamp(rho, _rho_end, _rho_default);
}