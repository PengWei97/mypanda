//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GBAnisotropyMisori.h"
#include "MooseMesh.h"
#include <cmath>

registerMooseObject("mypandaApp", GBAnisotropyMisori);

InputParameters
GBAnisotropyMisori::validParams()
{
  InputParameters params = GBAnisotropyMisoriBase::validParams();
  params.addClassDescription(
      "Computes necessary material properties for the anisotropic grain growth model");
  MooseEnum crystal_structures("FCC BCC HCP", "FCC");
  params.addParam<MooseEnum>("crystal_structure", crystal_structures, "The type of crystal structure");

  // Added parameters for twin boundary energy and mobility
  params.addParam<Real>("TT1_sigma",  0.9, "Twin boundary energy for {10-12} tensile twin (type 1) based on MD, J/m^2");
  params.addParam<Real>("CT1_sigma",  0.9, "Twin boundary energy for {11-22} compresssion twin (type 1) based on MD, J/m^2");  
  params.addParam<Real>("TT1_mob", 2.5e-6, "Twin boundary mobility for {10-12} tensile twin (type 1) based on experiment, m^4/(J*s)");
  params.addParam<Real>("CT1_mob", 2.5e-6, "Twin boundary mobility for {11-22} compresssion twin (type 1) based on experiment, m^4/(J*s)");
  params.addParam<Real>("Sigma9_sigma",  0.9, "Twin boundary energy for {10-12} tensile twin (type 1) based on MD, J/m^2");
  params.addParam<Real>("Sigma3_sigma",  0.9, "Twin boundary energy for {11-22} compresssion twin (type 1) based on MD, J/m^2");  
  params.addParam<Real>("Sigma9_mob", 2.5e-6, "Twin boundary mobility for {10-12} tensile twin (type 1) based on experiment, m^4/(J*s)");
  params.addParam<Real>("Sigma3_mob", 2.5e-6, "Twin boundary mobility for {11-22} compresssion twin (type 1) based on experiment, m^4/(J*s)");

  params.addRequiredParam<UserObjectName>(
      "grain_tracker", "Name of GrainTracker user object that provides Grain ID according to element ID");
  params.addRequiredParam<UserObjectName>("euler_angle_provider",
                                          "Name of Euler angle provider user object");
  params.addParam<bool>("gb_energy_anisotropy", false,
      "The GB energy anisotropy based on misorientation would be considered if true");
  params.addParam<bool>("gb_mobility_anisotropy", true,
      "The GB mobility anisotropy would be considered if true");

  return params;
}

GBAnisotropyMisori::GBAnisotropyMisori(const InputParameters & parameters)
  : GBAnisotropyMisoriBase(parameters),
    _crystal_structure(getParam<MooseEnum>("crystal_structure").getEnum<MisorientationAngleCalculator::CrystalType>()),

    _TT1_sigma(getParam<Real>("TT1_sigma")),
    _CT1_sigma(getParam<Real>("CT1_sigma")),
    _TT1_mob(getParam<Real>("TT1_mob")),
    _CT1_mob(getParam<Real>("CT1_mob")),
    _Sigma3_sigma(getParam<Real>("Sigma3_sigma")),
    _Sigma9_sigma(getParam<Real>("Sigma9_sigma")),
    _Sigma3_mob(getParam<Real>("Sigma3_mob")),
    _Sigma9_mob(getParam<Real>("Sigma9_mob")),

    _grain_tracker(getUserObject<GrainTracker>("grain_tracker")),
    _euler(getUserObject<EulerAngleProvider>("euler_angle_provider")),
    _gb_energy_anisotropy(getParam<bool>("gb_energy_anisotropy")),
    _gb_mobility_anisotropy(getParam<bool>("gb_mobility_anisotropy")),
    _misori_angle(declareProperty<Real>("misori_angle")),
    _twinning_type(declareProperty<Real>("twinning_type"))
{
}

void
GBAnisotropyMisori::computeGBProperties()
{
  auto & time_current = _fe_problem.time();

  _misori_angle[_qp] = 0.0;
  _twinning_type[_qp] = -1.0;

  // get the GB location based on the GrainTracker in the quadrature point
  const auto & op_to_grains = _grain_tracker.getVarToFeatureVector(_current_elem->id()); 
  std::vector<unsigned int> var_index, grain_id_index;

  for (MooseIndex(op_to_grains) op_index = 0; op_index < op_to_grains.size(); ++op_index)
  {
    auto grain_id = op_to_grains[op_index]; // grain id

    if (grain_id != FeatureFloodCount::invalid_id)
    {
      var_index.push_back(op_index);
      grain_id_index.push_back(grain_id);
    }
  }
   
  Real sigma_min = _GBsigma_HAGB, sigma_max = _GBsigma_HAGB;
  Real mob_min = _GBmob_HAGB, mob_max = _GBmob_HAGB;

  // When at grain boundaries or junction
  if (grain_id_index.size() > 1)
  {
    std::fill(_sigma.begin(), _sigma.end(), std::vector<Real>(_op_num, 0.0));
    std::fill(_mob.begin(), _mob.end(), std::vector<Real>(_op_num, 0.0));

    // Traverse to obtain the sigma_ij, mob_ij of the activation order parameters at the orthogonal point
    for (unsigned int i = 0; i < grain_id_index.size() - 1; ++i)
    {
      auto angles_i = _euler.getEulerAngles(grain_id_index[i]);
      for (unsigned int j = i+1; j < grain_id_index.size(); ++j)
      {
        auto angles_j = _euler.getEulerAngles(grain_id_index[j]);
        auto & _sigma_ij = _sigma[var_index[i]][var_index[j]];
        auto & _mob_ij = _mob[var_index[i]][var_index[j]];

        // calculate misorientation angle
        _misori_s = MisorientationAngleCalculator::calculateMisorientaion(angles_i, angles_j, _misori_s, _crystal_structure);

        _sigma_ij = _gb_energy_anisotropy ? calculatedGBEnergy(_misori_s) : _GBsigma_HAGB; // && (time_current > 10.0)
        _mob_ij = _gb_mobility_anisotropy ? calculatedGBMobility(_misori_s) : _GBmob_HAGB;

        _sigma[var_index[j]][var_index[i]] =  _sigma_ij;
        _mob[var_index[j]][var_index[i]] =  _mob_ij;

        // Set the twinning type
        if (i == 0 && j == 1)
        {
          _misori_angle[_qp] =  _misori_s._misor;
          _twinning_type[_qp] = determineTwinningType(_misori_s);
          sigma_min = sigma_max = _sigma_ij;
          mob_min = mob_max = _mob_ij;
        }

        sigma_min = std::min(sigma_min, _sigma_ij);
        sigma_max = std::max(sigma_max, _sigma_ij);
        mob_min = std::min(mob_min, _mob_ij);
        mob_max = std::max(mob_max, _mob_ij);
      }
    }

    fillSymmetricProperties(sigma_min, sigma_max, mob_min, mob_max);
  }
}



Real
GBAnisotropyMisori::calculatedGBEnergy(const MisorientationAngleData & misori_s)
{
  Real gbSigma = _GBsigma_HAGB;
  Real trans_misori_angle_HAGB = 15.0; // transition misorientation angle between low and high-angle grain boundary

  if ((_misori_s._misor <= 1.0) && !misori_s._is_twin)
    gbSigma = _GBsigma_HAGB * ((2.0 / trans_misori_angle_HAGB * (1 - std::log(2.0 / trans_misori_angle_HAGB))));    
  else if ((_misori_s._misor <= trans_misori_angle_HAGB) && !misori_s._is_twin)
    gbSigma = _GBsigma_HAGB * ((misori_s._misor / trans_misori_angle_HAGB * (1 - std::log(misori_s._misor / trans_misori_angle_HAGB))));
  else if (misori_s._is_twin)
  {
    switch (misori_s._twin_type)
    {
      case TwinType::TT1_HCP: gbSigma = _TT1_sigma; break;
      case TwinType::CT1_HCP: gbSigma = _CT1_sigma; break;
      case TwinType::Sigma3_FCC: gbSigma = _Sigma3_sigma; break;
      case TwinType::Sigma9_FCC: gbSigma = _Sigma9_sigma; break;
      default: break;
    }
  }

  return gbSigma; 
}
 
Real
GBAnisotropyMisori::calculatedGBMobility(const MisorientationAngleData & misori_s)
{
  // Initialize GB mobility
  Real gbMob = _GBmob_HAGB;

  // transition misorientation angle between low and high-angle grain boundary
  Real trans_misori_angle_HAGB = 15;

  // Equation constant
  Real B = 5;
  Real n = 4;
  
  if ((_misori_s._misor <= 1.0) && !misori_s._is_twin)
    gbMob = _GBmob_HAGB * ((1- std::exp(-B * std::pow( 2.0 / trans_misori_angle_HAGB, n)))); // Eq.8    
  else if ((_misori_s._misor <=  trans_misori_angle_HAGB) && !misori_s._is_twin )
    gbMob = _GBmob_HAGB * ((1- std::exp(-B * std::pow( misori_s._misor / trans_misori_angle_HAGB, n)))); // Eq.8
  else if (misori_s._is_twin)
  {
    switch (misori_s._twin_type)
    {
      case TwinType::TT1_HCP: gbMob = _TT1_mob; break;
      case TwinType::CT1_HCP: gbMob = _CT1_mob; break;
      case TwinType::Sigma3_FCC: gbMob = _Sigma3_mob; break;
      case TwinType::Sigma9_FCC: gbMob = _Sigma9_mob; break;
      default: break;
    }
  }

  return gbMob;
}

Real 
GBAnisotropyMisori::determineTwinningType(const MisorientationAngleData & misori_s)
{
  if (misori_s._is_twin)
  {
    switch (misori_s._twin_type)
    {
      case TwinType::TT1_HCP: return 1.0;
      case TwinType::CT1_HCP: return 2.0;
      case TwinType::Sigma3_FCC: return 3.0;
      case TwinType::Sigma9_FCC: return 4.0;
      default: return 0.0; // GB
    }
  }
  return 0.0; // GB
}

void 
GBAnisotropyMisori::fillSymmetricProperties(Real sigma_min, Real sigma_max, Real mob_min, Real mob_max)
{
  for (unsigned int i = 0; i < _op_num; ++i)
    for (unsigned int j = 0; j < _op_num; ++j)
    {
      if (_sigma[i][j] == 0.0)
        _sigma[i][j] = _sigma[j][i] = (sigma_max + sigma_min) / 2.0;

      if (_mob[i][j] == 0.0)
        _mob[i][j] = _mob[j][i] = (mob_max + mob_min) / 2.0;
    }
}