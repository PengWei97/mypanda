#include "GBAnisotropyWtGrainType.h"

registerMooseObject("mypandaApp", GBAnisotropyWtGrainType);

InputParameters
GBAnisotropyWtGrainType::validParams()
{
  InputParameters params = GBAnisotropyMisoriBase::validParams();
  params.addClassDescription(
    "Extends GBAnisotropyMisoriBase to account for grain-type-specific properties "
    "using data provided by an EBSD reader object.");
  params.addRequiredParam<UserObjectName>("ebsd_reader", 
    "Name of the EBSD reader user object that provides grain-specific data.");  
  params.addRequiredParam<UserObjectName>("grain_tracker", "GrainTracker user object");
  params.addRequiredParam<Real>("select_grain_type", "Selected grain type id.");
  params.addParam<bool>("gb_energy_anisotropy", false, "Consider GB energy anisotropy");
  params.addParam<bool>("gb_mobility_anisotropy", false, "Consider GB mobility anisotropy");
  params.addParam<Real>("execution_time", 2.0, "Time to perform grain boundary anisotropy");
  return params;
}

GBAnisotropyWtGrainType::GBAnisotropyWtGrainType(const InputParameters & parameters)
: GBAnisotropyMisoriBase(parameters),
  _ebsd_reader(getUserObject<EBSDReader>("ebsd_reader")),
  _grain_tracker(getUserObject<GrainTracker>("grain_tracker")),
  _select_grain_type(getParam<Real>("select_grain_type")),
  _gb_energy_anisotropy(getParam<bool>("gb_energy_anisotropy")),
  _gb_mobility_anisotropy(getParam<bool>("gb_mobility_anisotropy")),
  _execution_time(getParam<Real>("execution_time")),
  _grain_type(declareProperty<Real>("grain_type"))
{
}

void
GBAnisotropyWtGrainType::computeGBProperties()
{
  // Call the base class method to compute properties
  GBAnisotropyMisoriBase::computeGBProperties();

  const auto & op_to_grains = _grain_tracker.getVarToFeatureVector(_current_elem->id());
  std::vector<unsigned int> var_index, grain_ids;
  Real sum_h = 0.0;
  _grain_type[_qp] = 0.0;

  // Compute effective dislocation density (rho_eff)
  for (MooseIndex(op_to_grains) op_index = 0; op_index < op_to_grains.size(); ++op_index)
  {
    auto grain_id = op_to_grains[op_index];
    if (grain_id == FeatureFloodCount::invalid_id)
      continue;

    Real h = (1.0 + std::sin(libMesh::pi * ((*_vals[op_index])[_qp] - 0.5))) * 0.5;
    _grain_type[_qp] += std::round(_ebsd_reader.getAvgData(grain_id)._custom[0] * h);
    sum_h += h;

    var_index.push_back(op_index);
    grain_ids.push_back(grain_id);    
  }
  _grain_type[_qp] /= std::max(sum_h, 1e-10);

  // Set default values
  Real sigma_min = _GBsigma_HAGB, sigma_max = _GBsigma_HAGB;
  Real mob_min = _GBmob_HAGB, mob_max = _GBmob_HAGB;

  if (grain_ids.size() <= 1 || _fe_problem.time() <= _execution_time)
    return;

  Real init_sigma = _gb_energy_anisotropy ? 0.0 : _GBsigma_HAGB;
  Real init_mob = _gb_mobility_anisotropy ? 0.0 : _GBmob_HAGB;
  std::fill(_sigma.begin(), _sigma.end(), std::vector<Real>(_op_num, init_sigma));
  std::fill(_mob.begin(), _mob.end(), std::vector<Real>(_op_num, init_mob));

  if (_gb_energy_anisotropy || _gb_mobility_anisotropy)
  {
    Real sigma_min = _GBsigma_HAGB, sigma_max = _GBsigma_HAGB;
    Real mob_min = _GBmob_HAGB, mob_max = _GBmob_HAGB;
    computeSigmaAndMobility(var_index, grain_ids, sigma_min, sigma_max, mob_min, mob_max);
    fillSymmetricProperties(sigma_min, sigma_max, mob_min, mob_max);
  }
}

void
GBAnisotropyWtGrainType::computeSigmaAndMobility(const std::vector<unsigned int> & var_index,
                                                 const std::vector<unsigned int> & grain_ids,
                                                 Real & sigma_min, Real & sigma_max,
                                                 Real & mob_min, Real & mob_max)
{
  for (unsigned int i = 0; i < grain_ids.size() - 1; ++i)
    for (unsigned int j = i + 1; j < grain_ids.size(); ++j)
    {
      auto gi = grain_ids[i], gj = grain_ids[j];

      Real sigma = _gb_energy_anisotropy ? calculateGBEnergy(gi, gj) : _GBsigma_HAGB;
      Real mob = _gb_mobility_anisotropy ? calculateGBMobility(gi, gj) : _GBmob_HAGB;

      _sigma[var_index[i]][var_index[j]] = sigma;
      _sigma[var_index[j]][var_index[i]] = sigma;
      _mob[var_index[i]][var_index[j]] = mob;
      _mob[var_index[j]][var_index[i]] = mob;

      sigma_min = std::min(sigma_min, sigma);
      sigma_max = std::max(sigma_max, sigma);
      mob_min = std::min(mob_min, mob);
      mob_max = std::max(mob_max, mob);
    }
}

Real
GBAnisotropyWtGrainType::calculateGBEnergy(const Real & gi, const Real & gj)
{
  Real ti = std::round(_ebsd_reader.getAvgData(gi)._custom[0]);
  Real tj = std::round(_ebsd_reader.getAvgData(gj)._custom[0]);
  
  return (ti == _select_grain_type && tj == _select_grain_type) ? _GBsigma_HAGB * 0.1 : _GBsigma_HAGB;
}

Real
GBAnisotropyWtGrainType::calculateGBMobility(const Real & gi, const Real & gj)
{
  Real ti = std::round(_ebsd_reader.getAvgData(gi)._custom[0]);
  Real tj = std::round(_ebsd_reader.getAvgData(gj)._custom[0]);

  return (ti == _select_grain_type && tj == _select_grain_type) ? _GBmob_HAGB * 0.1 : _GBmob_HAGB;
}

void 
GBAnisotropyWtGrainType::fillSymmetricProperties(Real sigma_min, Real sigma_max, Real mob_min, Real mob_max)
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