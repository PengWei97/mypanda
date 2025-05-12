#include "GBAnisotropyMisoriAndTwin.h"

registerMooseObject("mypandaApp", GBAnisotropyMisoriAndTwin);

InputParameters
GBAnisotropyMisoriAndTwin::validParams()
{
  InputParameters params = GBAnisotropyMisoriAng::validParams();
  params.addClassDescription("Grain boundary anisotropy based on crystal structure and misorientation with twin boundaries");
  params.addParam<Real>("TT1_sigma", 0.9, "Twin boundary energy for {10-12} tensile twin (type 1) based on MD, J/m^2");
  params.addParam<Real>("CT1_sigma", 0.9, "Twin boundary energy for {11-22} compresssion twin (type 1) based on MD, J/m^2");  
  params.addParam<Real>("TT1_mob", 2.5e-6, "Twin boundary mobility for {10-12} tensile twin (type 1) based on experiment, m^4/(J*s)");
  params.addParam<Real>("CT1_mob", 2.5e-6, "Twin boundary mobility for {11-22} compresssion twin (type 1) based on experiment, m^4/(J*s)");
  params.addParam<Real>("Sigma9_sigma", 0.9, "Twin boundary energy for Sigma 9 in FCC, J/m^2");
  params.addParam<Real>("Sigma3_sigma", 0.9, "Twin boundary energy for Sigma 3 based on MD, J/m^2");
  params.addParam<Real>("Sigma9_mob", 2.5e-6, "Twin boundary mobility for Sigma 9 in FCC based on experiment, m^4/(J*s)");
  params.addParam<Real>("Sigma3_mob", 2.5e-6, "Twin boundary mobility for Sigma 3 in FCC on experiment, m^4/(J*s)");

  return params;
}

GBAnisotropyMisoriAndTwin::GBAnisotropyMisoriAndTwin(const InputParameters & parameters)
  : GBAnisotropyMisoriAng(parameters),
    _TT1_sigma(getParam<Real>("TT1_sigma")),
    _CT1_sigma(getParam<Real>("CT1_sigma")),
    _TT1_mob(getParam<Real>("TT1_mob")),
    _CT1_mob(getParam<Real>("CT1_mob")),
    _Sigma9_sigma(getParam<Real>("Sigma9_sigma")),
    _Sigma3_sigma(getParam<Real>("Sigma3_sigma")),
    _Sigma9_mob(getParam<Real>("Sigma9_mob")),
    _Sigma3_mob(getParam<Real>("Sigma3_mob")),
    _twin_boundary_type(declareProperty<Real>("twin_boundary_type"))
{
}

void
GBAnisotropyMisoriAndTwin::initOthersMaterialProperties()
{
  _twin_boundary_type[_qp] = -2.0;
}

Real
GBAnisotropyMisoriAndTwin::calculateGBEnergy(const MisorientationAngleData & misori_s)
{
  // Default to general grain boundary (HAGB) energy
  Real gbSigma = GBAnisotropyMisoriAng::calculateGBEnergy(misori_s);

  // If not a twin boundary, return the calculated HAGB energy
  if (!misori_s._is_twin)
  {
    _twin_boundary_type[_qp] = 0.0;
    return gbSigma;
  }

  // Select the twin boundary energy based on the twin type
  const std::unordered_map<TwinType, std::pair<Real, Real>> twin_properties = {
    {TwinType::TT1_HCP, {_TT1_sigma, 1.0}},
    {TwinType::CT1_HCP, {_CT1_sigma, 2.0}},
    {TwinType::Sigma3_FCC, {_Sigma3_sigma, 3.0}},
    {TwinType::Sigma9_FCC, {_Sigma9_sigma, 9.0}}
  };

  auto it = twin_properties.find(misori_s._twin_type);
  if (it != twin_properties.end())
  {
    gbSigma = it->second.first; // update gbSigma to the twin boundary energy
    _twin_boundary_type[_qp] = it->second.second; // Twin boundary type (1.0 for TT1, 2.0 for CT1, 3.0 for Sigma3, 9.0 for Sigma9)
  }
  else
    _twin_boundary_type[_qp] = 0.0; // normal GB region, set to 0.0

  return gbSigma;
}

Real
GBAnisotropyMisoriAndTwin::calculateGBMobility(const MisorientationAngleData & misori_s)
{
  // Default to general grain boundary (HAGB) mobility
  Real gbMob = GBAnisotropyMisoriAng::calculateGBMobility(misori_s);

  // If not a twin boundary, return the calculated HAGB mobility
  if (!misori_s._is_twin)
    return gbMob;

  // Select the twin boundary mobility based on the twin type
  static const std::unordered_map<TwinType, Real> twin_mobility_map = {
    {TwinType::TT1_HCP, _TT1_mob},
    {TwinType::CT1_HCP, _CT1_mob},
    {TwinType::Sigma3_FCC, _Sigma3_mob},
    {TwinType::Sigma9_FCC, _Sigma9_mob}
  };

  return gbMob;
}