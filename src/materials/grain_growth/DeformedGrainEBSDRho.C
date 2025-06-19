#include "DeformedGrainEBSDRho.h"

registerMooseObject("mypandaApp", DeformedGrainEBSDRho);

InputParameters
DeformedGrainEBSDRho::validParams()
{
  InputParameters params = DerivativeMaterialInterface<Material>::validParams();
  params.addRequiredCoupledVarWithAutoBuild("v", "var_name_base", "op_num", "Array of coupled variables");
  params.addParam<Real>("length_scale", 1.0e-9, "Length scale in meters");
  params.addParam<Real>("time_scale", 1.0e-6, "Time scale in seconds");
  params.addParam<Real>("Elas_Mod", 2.5e10, "Shear modulus in J/m^3");
  params.addParam<Real>("Burg_vec", 3.0e-10, "Burgers vector length in meters");
  params.addParam<Real>("stored_factor", 0.5, "Scaling factor for stored energy");
  params.addParam<Real>("execution_time", 2.0, "Simulation time to begin derivative evaluation");
  params.addRequiredParam<UserObjectName>("grain_tracker", "GrainTracker UserObject");
  params.addRequiredParam<UserObjectName>("GNDs_provider", "UserObject providing dislocation data");

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
    _grain_type_rho(declareProperty<Real>("grain_type_rho")),
    _D_stored_energy(_op_num)
{
  _beta = _stored_factor * _Elas_Mod * _Burg_vec * _Burg_vec * _JtoeV * std::pow(_length_scale, 3);

  if (_op_num == 0)
    paramError("op_num", "At least one order parameter (op) is required.");

  // Loop over variables (ops)
    for (MooseIndex(_op_num) op_index = 0; op_index < _op_num; ++op_index)
    _D_stored_energy[op_index] = &declarePropertyDerivative<Real>(
      "stored_energy", coupledName("v", op_index));
}

void
DeformedGrainEBSDRho::computeQpProperties()
{
  // Step 1: Compute sum(eta_i^2)
  const Real SumEtai2 = std::accumulate(_vals.begin(), _vals.end(), 0.0,
    [this](Real sum, const VariableValue * val) {
      return sum + (*val)[_qp] * (*val)[_qp];
    });

  if (SumEtai2 <= 1e-16)
  {
    _rho_eff[_qp] = 0.0;
    for (auto & deriv : _D_stored_energy)
      if (deriv)
        (*deriv)[_qp] = 0.0;
    return;
  }

  const auto & op_to_grains = _grain_tracker.getVarToFeatureVector(_current_elem->id());

  _rho_eff[_qp] = 0.0;
  _feature_id[_qp] = 0.0;
  _grain_type_rho[_qp] = 0.0;

  bool feature_set = false;

  for (unsigned int i = 0; i < op_to_grains.size(); ++i)
  {
    const auto grain_id = op_to_grains[i];
    if (grain_id == FeatureFloodCount::invalid_id)
      continue;

    const Real eta_i = (*_vals[i])[_qp];
    const Real y_coord = _grain_tracker.getGrainCentroid(grain_id)(1);
    const Real rho_i = _GNDs_provider.getRhoWtTime(grain_id, y_coord);

    _rho_eff[_qp] += rho_i * eta_i * eta_i;

    if (!feature_set)
    {
      _feature_id[_qp] = static_cast<Real>(_GNDs_provider.getFeatureID(grain_id));
      _grain_type_rho[_qp] = _GNDs_provider.isType2Grain(grain_id, y_coord) ? 2.0 : 1.0;
      feature_set = true;
    }
  }

  _rho_eff[_qp] /= SumEtai2;

  // Step 2: Check simulation time
  // If it's too early for the update, reset all derivatives to 0.0 and exit
  const auto & time_current = _fe_problem.time(); // Current simulation time in seconds
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

  // Step 3: Compute derivative of stored energy

  for (MooseIndex(_op_num) op_index = 0; op_index < _op_num; ++op_index)
  {
    Real & C_deriv = (*_D_stored_energy[op_index])[_qp];
    C_deriv = 0.0;

    const auto grain_id = op_to_grains[op_index];
    if (grain_id == FeatureFloodCount::invalid_id)
      continue;

    const Real y_coord = _grain_tracker.getGrainCentroid(grain_id)(1);
    const Real rho_i = _GNDs_provider.getRhoWtTime(grain_id, y_coord);

    C_deriv = (rho_i - _rho_eff[_qp]) * _beta / SumEtai2;
  }
}