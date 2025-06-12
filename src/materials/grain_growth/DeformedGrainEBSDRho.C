#include "DeformedGrainEBSDRho.h"
#include "GrainTrackerInterface.h"

// Register the object with the application
registerMooseObject("mypandaApp", DeformedGrainEBSDRho);

// Define valid parameters for the class
InputParameters
DeformedGrainEBSDRho::validParams()
{
  InputParameters params = DerivativeMaterialInterface<Material>::validParams();
  params.addRequiredCoupledVarWithAutoBuild("v", "var_name_base", "op_num", "Array of coupled variables");
  params.addParam<Real>("length_scale", 1.0e-9, "Length scale in meters (default is nm)");
  params.addParam<Real>("time_scale", 1.0e-6, "Time scale in seconds (default is microseconds)");
  params.addParam<Real>("Elas_Mod", 2.50e10, "Shear modulus in J/m^3");
  params.addParam<Real>("Burg_vec", 3.0e-10, "Length of Burgers vector in meters");
  params.addParam<Real>("stored_factor", 0.5, "Scaling factor in stored energy function");
  params.addParam<Real>("execution_time", 2.0, "Time to perform grain boundary anisotropy");

  params.addRequiredParam<UserObjectName>("grain_tracker", "GrainTracker UserObject to get values from");
  params.addRequiredParam<UserObjectName>("GNDs_provider", "GNDs provider for EBSD reader");
  return params;
}

DeformedGrainEBSDRho::DeformedGrainEBSDRho(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _op_num(coupledComponents("v")),
    _vals(coupledValues("v")),
    _length_scale(getParam<Real>("length_scale")),
    _time_scale(getParam<Real>("time_scale")),
    _Elas_Mod(getParam<Real>("Elas_Mod")),
    _Burg_vec(getParam<Real>("Burg_vec")),
    _stored_factor(getParam<Real>("stored_factor")),
    _execution_time(getParam<Real>("execution_time")),
    _JtoeV(6.24150974e18),
    _grain_tracker(getUserObject<GrainTrackerInterface>("grain_tracker")),
    _GNDs_provider(getUserObject<EBSDReaderMaterialProperty>("GNDs_provider")),

    _rho_eff(declareProperty<Real>("rho_eff")),
    _feature_id(declareProperty<Real>("feature_id")),
    _D_stored_energy(_op_num)
{
  if (_op_num == 0)
    paramError("op_num", "Model requires op_num > 0");

  // Loop over variables (ops)
  for (MooseIndex(_op_num) op_index = 0; op_index < _op_num; ++op_index)
    _D_stored_energy[op_index] = &declarePropertyDerivative<Real>(
      "stored_energy", coupledName("v", op_index));
}

void
DeformedGrainEBSDRho::computeQpProperties()
{
  // Step 1: Compute the sum of squared order parameters (SumEtai2)
  const Real SumEtai2 = std::accumulate(_vals.begin(), _vals.end(), 0.0,
                                        [this](Real sum, const VariableValue *val) {
                                          return sum + (*val)[_qp] * (*val)[_qp];
                                        });

  // Retrieve the mapping of order parameters to grain IDs
  const auto & op_to_grains = _grain_tracker.getVarToFeatureVector(_current_elem->id());
  _rho_eff[_qp] = 0.0;
  _feature_id[_qp] = 0.0;

  bool is_first_set_feature_id = true;
  bool is_first_set_delta_rho_2 = true;

  // Step 2: Compute effective dislocation density (rho_eff)
  for (const auto & grain_id : op_to_grains)
  {
    if (grain_id == FeatureFloodCount::invalid_id)
      continue;

    // Get order parameter value for current grain
    const unsigned int op_index = &grain_id - &op_to_grains[0]; // Calculate the index
    const Real op_value = (*_vals[op_index])[_qp];

    // Get the GNDs for the current grain only once
    const Real rho_i = getRhoWtTime(grain_id); // GNDs for each grain, 1/m^2
    _rho_eff[_qp] += rho_i * op_value * op_value; // rho_eff = sum(rho_i * eta_i^2)

    if (is_first_set_feature_id)
    {
      const unsigned int feature_id = _GNDs_provider.getFeatureID(grain_id);
      _feature_id[_qp] = static_cast<Real>(feature_id);
      is_first_set_feature_id = false;
    }  
  }

  // Normalize by the sum of squared order parameters
  _rho_eff[_qp] /= SumEtai2;

  // Step 3: Precompute beta factor
  const Real beta = _stored_factor * _Elas_Mod * _Burg_vec * _Burg_vec * 
                    _JtoeV * std::pow(_length_scale, 3);

  // Step 4: Check if the current simulation time is before the execution time
  const auto & time_current = _fe_problem.time(); // Current simulation time in seconds
  // If it's too early for the update, reset all derivatives to 0.0 and exit
  if (time_current < _execution_time)
  {
    // Use std::for_each to cleanly reset all derivatives
    std::for_each(_D_stored_energy.begin(), _D_stored_energy.end(),
                  [this](MaterialProperty<Real> *property) {
                    if (property) // Ensure pointer is valid
                      (*property)[_qp] = 0.0;
                  });
    return;
  }

  // Step 5: Compute stored energy derivatives
  for (MooseIndex(_op_num) op_index = 0; op_index < _op_num; ++op_index)
  {
    Real & C_deriv = (*_D_stored_energy[op_index])[_qp];
    C_deriv = 0.0;

    auto grain_id = op_to_grains[op_index];
    if (grain_id == FeatureFloodCount::invalid_id)
      continue;

    // Get the GNDs for the current grain (cached result)
    const Real rho_i = getRhoWtTime(grain_id); // GNDs for each grain, 1/m^2
    C_deriv = (rho_i - _rho_eff[_qp]) * beta / SumEtai2;
  }
}

Real 
DeformedGrainEBSDRho::getRhoWtTime(const unsigned int & grain_id) const
{
  const Real y_coord = _grain_tracker.getGrainCentroid(grain_id)(1);
  const Real rho = _GNDs_provider.getRhoInit(grain_id);

  // (Optional) Skip recovery for grains above a certain y-coordinate
  if (y_coord > 100.0 && rho > 1.0e13)
    return std::clamp(rho, 1.0e11, 2.0e15);

  // Get the time-evolved dislocation density for a specified grain
  return _GNDs_provider.getRhoWtTime(grain_id);
} 