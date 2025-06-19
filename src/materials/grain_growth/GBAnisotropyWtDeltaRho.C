#include "GBAnisotropyWtDeltaRho.h"

registerMooseObject("mypandaApp", GBAnisotropyWtDeltaRho);

InputParameters
GBAnisotropyWtDeltaRho::validParams()
{
  InputParameters params = GBAnisotropyMisoriAndTwin::validParams();
  params.addClassDescription("Grain boundary anisotropy considering Delta Rho and twin boundary effects");
  params.addRequiredParam<UserObjectName>("GNDs_provider", "EBSD reader for dislocation density");
  params.addParam<bool>("is_gb_mob_with_delta_rho", false, "Consider GB mobility with Delta Rho");
  params.addParam<Real>("trans_delta_rho", 90.0, "Critical delta rho for transformation considering GB mobility with delta rho"); 
  params.addParam<Real>("amplifier_factor", 1.0, "Amplification factor for GB mobility with delta rho");
  return params;
}

GBAnisotropyWtDeltaRho::GBAnisotropyWtDeltaRho(const InputParameters & parameters)
  : GBAnisotropyMisoriAndTwin(parameters),
    _GNDs_provider(getUserObject<EBSDReaderMaterialProperty>("GNDs_provider")),
    _is_gb_mob_with_delta_rho(getParam<bool>("is_gb_mob_with_delta_rho")),
    _trans_delta_rho(getParam<Real>("trans_delta_rho")),
    _amplifier_factor(getParam<Real>("amplifier_factor")),
    _delta_rho(declareProperty<Real>("delta_rho"))
{
}

void
GBAnisotropyWtDeltaRho::initOthersMaterialProperties()
{
  GBAnisotropyMisoriAndTwin::initOthersMaterialProperties();

  _delta_rho[_qp] = 0.0;
}

void
GBAnisotropyWtDeltaRho::calculateOthersMaterialProperties(const unsigned int & grain_i, const unsigned int & grain_j)
{
  // Efficiently calculate the absolute difference in dislocation density, scaled by length^2
  const Real y_coord_i = _grain_tracker.getGrainCentroid(grain_i)(1);
  const Real rho_i = _GNDs_provider.getRhoWtTime(grain_i, y_coord_i);
  const Real y_coord_j = _grain_tracker.getGrainCentroid(grain_j)(1);
  const Real rho_j = _GNDs_provider.getRhoWtTime(grain_j, y_coord_j);

  _delta_rho[_qp] = std::abs(rho_i - rho_j) * _length_scale * _length_scale;
}

Real 
GBAnisotropyWtDeltaRho::calculateGBMobility(const MisorientationAngleData & misori_s)
{
  // Step 1: Compute the base mobility using parent class method
  Real mob_ij = GBAnisotropyMisoriAndTwin::calculateGBMobility(misori_s);

  // Step 2: Apply Delta Rho influence if enabled
  if (_is_gb_mob_with_delta_rho)
    calculatedGBMobilityWtRho(_delta_rho[_qp], mob_ij);

  return mob_ij;
}

void
GBAnisotropyWtDeltaRho::calculatedGBMobilityWtRho(const Real & delta_rho, Real & mob_ij)
{
  // Step 1: Compute high mobility rate with amplification factor
  const Real mob_ij_high = mob_ij * _amplifier_factor;

  // Step 2: Efficient ternary-based computation for mobility update
  const Real exponent_term = std::exp(-_B * std::pow(delta_rho / _trans_delta_rho, _n));
  const Real mob_temp = (delta_rho <= _trans_delta_rho) ? mob_ij_high * (1 - exponent_term) : mob_ij_high;

  // Step 3: Update only if the calculated value is larger
  mob_ij = std::max(mob_ij, mob_temp);
}

