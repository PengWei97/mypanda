#include "GBAnisotropyWtDeltaRho3.h"

registerMooseObject("mypandaApp", GBAnisotropyWtDeltaRho3);

InputParameters
GBAnisotropyWtDeltaRho3::validParams()
{
  InputParameters params = GBAnisotropyMisoriAndTwin::validParams();
  params.addClassDescription("Grain boundary anisotropy considering Delta Rho and twin boundary effects");
  params.addRequiredParam<UserObjectName>("ebsd_reader", "EBSD data reader");
  params.addParam<bool>("is_gb_mob_with_delta_rho", false, "Consider GB mobility with Delta Rho");
  params.addRequiredParam<std::vector<unsigned int>>("spec_grains", "The vector of the specific grain IDs"); 
  params.addParam<Real>("amplifier_factor", 1.0, "Amplification factor for GB mobility with delta rho");
  params.addParam<Real>("threshold_for_sGrains", 50.0, "The threshold factor for delta_rho");
  return params;
}

GBAnisotropyWtDeltaRho3::GBAnisotropyWtDeltaRho3(const InputParameters & parameters)
  : GBAnisotropyMisoriAndTwin(parameters),
    _ebsd_reader(getUserObject<EBSDReaderMaterialProperty>("ebsd_reader")),
    _is_gb_mob_with_delta_rho(getParam<bool>("is_gb_mob_with_delta_rho")),
    _spec_grains(getParam<std::vector<unsigned int>>("spec_grains")),
    _amplifier_factor(getParam<Real>("amplifier_factor")),
    _threshold_for_sGrains(getParam<Real>("threshold_for_sGrains")),
    _delta_rho(declareProperty<Real>("delta_rho"))
{
}

void
GBAnisotropyWtDeltaRho3::initOthersMaterialProperties()
{
  GBAnisotropyMisoriAndTwin::initOthersMaterialProperties();

  _delta_rho[_qp] = 0.0;
}


void
GBAnisotropyWtDeltaRho3::calculateOthersMaterialProperties(const unsigned int & grain_i, const unsigned int & grain_j)
{
  // Efficiently calculate the absolute difference in dislocation density, scaled by length^2
  const Real rho_i = _ebsd_reader.getRhoWtTime(grain_i);
  const Real rho_j = _ebsd_reader.getRhoWtTime(grain_j);
  _delta_rho[_qp] = std::abs(rho_i - rho_j) * _length_scale * _length_scale;

  _grain_i = grain_i;
  _grain_j = grain_j;
}

Real 
GBAnisotropyWtDeltaRho3::calculateGBMobility(const MisorientationAngleData & misori_s)
{
  // Step 1: Compute the base mobility using parent class method
  Real mob_ij = GBAnisotropyMisoriAndTwin::calculateGBMobility(misori_s);

  // Step 2: Apply Delta Rho influence if enabled
  if (_is_gb_mob_with_delta_rho)
    calculatedGBMobilityWtRho(_delta_rho[_qp], mob_ij);

  return mob_ij;
}

void
GBAnisotropyWtDeltaRho3::calculatedGBMobilityWtRho(const Real & delta_rho, Real & mob_ij)
{
  if (delta_rho > _threshold_for_sGrains &&
    (std::find(_spec_grains.begin(), _spec_grains.end(), _grain_i) != _spec_grains.end() ||
     std::find(_spec_grains.begin(), _spec_grains.end(), _grain_j) != _spec_grains.end()))
{
    mob_ij *= _amplifier_factor;
}
}

