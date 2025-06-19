#include "GBAnisotropyMisoriAng.h"

registerMooseObject("mypandaApp", GBAnisotropyMisoriAng);

InputParameters
GBAnisotropyMisoriAng::validParams()
{
  InputParameters params = GBAnisotropyMisoriBase::validParams();
  params.addClassDescription("Grain boundary anisotropy based on crystal structure and misorientation");
  MooseEnum crystal_structures("FCC BCC HCP", "FCC");
  params.addParam<MooseEnum>("crystal_structure", crystal_structures, "The type of crystal structure");
  params.addRequiredParam<UserObjectName>("grain_tracker", "GrainTracker user object");
  params.addRequiredParam<UserObjectName>("euler_angle_provider", "Euler angle provider user object");
  params.addParam<bool>("gb_energy_anisotropy", false, "Consider GB energy anisotropy");
  params.addParam<bool>("gb_mobility_anisotropy", false, "Consider GB mobility anisotropy");
  params.addParam<Real>("execution_time", 2.0, "Time to perform grain boundary anisotropy");
  return params;
}

GBAnisotropyMisoriAng::GBAnisotropyMisoriAng(const InputParameters & parameters)
: GBAnisotropyMisoriBase(parameters),
  _crystal_structure(getParam<MooseEnum>("crystal_structure").getEnum<MisorientationAngleCalculator::CrystalType>()),
  _execution_time(getParam<Real>("execution_time")),
  _grain_tracker(getUserObject<GrainTracker>("grain_tracker")),
  _euler(getUserObject<EulerAngleProvider>("euler_angle_provider")),
  _gb_energy_anisotropy(getParam<bool>("gb_energy_anisotropy")),
  _gb_mobility_anisotropy(getParam<bool>("gb_mobility_anisotropy")),
  _B(5),
  _n(4),
  _misori_angle(declareProperty<Real>("misori_angle"))
{
}

void
GBAnisotropyMisoriAng::computeGBProperties()
{
  auto & time_current = _fe_problem.time();

  _misori_angle[_qp] = 0.0;
  
  initOthersMaterialProperties();

  // get the GB location based on the GrainTracker in the quadrature point
  const auto & op_to_grains = _grain_tracker.getVarToFeatureVector(_current_elem->id()); 
  std::vector<unsigned int> var_index, grain_id_index;

  // Extract grain IDs and corresponding indices
  for (MooseIndex(op_to_grains) i = 0; i < op_to_grains.size(); ++i)
  {
    if (op_to_grains[i] != FeatureFloodCount::invalid_id)
    {
        var_index.push_back(i);
        grain_id_index.push_back(op_to_grains[i]);
    }
  }

  // Set default values
  Real sigma_min = _GBsigma_HAGB, sigma_max = _GBsigma_HAGB;
  Real mob_min = _GBmob_HAGB, mob_max = _GBmob_HAGB;

  // If at grain boundaries or junctions and the current time exceeds the execution time
  if (grain_id_index.size() > 1 && time_current > _execution_time)
  {
      // If there is anisotropy, initialize with 0.0; otherwise, use the HAGB values
      const Real init_sigma = _gb_energy_anisotropy ? 0.0 : _GBsigma_HAGB;
      const Real init_mob = _gb_mobility_anisotropy ? 0.0 : _GBmob_HAGB;

      // Create an initial vector filled with the appropriate values
      std::fill(_sigma.begin(), _sigma.end(), std::vector<Real>(_op_num, init_sigma));
      std::fill(_mob.begin(), _mob.end(), std::vector<Real>(_op_num, init_mob));

      // If anisotropy is enabled, perform additional computations
      if (_gb_energy_anisotropy || _gb_mobility_anisotropy)
      {
          // Compute the anisotropic grain boundary energy and mobility
          computeSigmaAndMobility(var_index, grain_id_index, sigma_min, sigma_max, mob_min, mob_max);

          fillSymmetricProperties(sigma_min, sigma_max, mob_min, mob_max);
      }
  }
}

void GBAnisotropyMisoriAng::computeSigmaAndMobility(const std::vector<unsigned int> & var_index,
                                                    const std::vector<unsigned int> & grain_id_index,
                                                    Real &sigma_min, Real & sigma_max,
                                                    Real &mob_min, Real & mob_max)
{
  for (unsigned int i = 0; i < grain_id_index.size() - 1; ++i)
  {
    _current_grain_i = grain_id_index[i];
    auto angles_i = _euler.getEulerAngles(_current_grain_i);
    for (unsigned int j = i + 1; j < grain_id_index.size(); ++j)
    {
        _current_grain_j = grain_id_index[j];
        auto angles_j = _euler.getEulerAngles(_current_grain_j);
        _misori_s = MisorientationAngleCalculator::calculateMisorientaion(angles_i, angles_j, _misori_s, _crystal_structure);

        _misori_angle[_qp] = _misori_s._misor;

        calculateOthersMaterialProperties(grain_id_index[i], grain_id_index[j]);

        Real sigma = _gb_energy_anisotropy ? calculateGBEnergy(_misori_s) : _GBsigma_HAGB;
        Real mobility = _gb_mobility_anisotropy ? calculateGBMobility(_misori_s) : _GBmob_HAGB;

        _sigma[var_index[i]][var_index[j]] = sigma;
        _sigma[var_index[j]][var_index[i]] = sigma;
        _mob[var_index[i]][var_index[j]] = mobility;
        _mob[var_index[j]][var_index[i]] = mobility;

        sigma_min = std::min(sigma_min, sigma);
        sigma_max = std::max(sigma_max, sigma);
        mob_min = std::min(mob_min, mobility);
        mob_max = std::max(mob_max, mobility);
      }
    }
}

Real
GBAnisotropyMisoriAng::calculateGBEnergy(const MisorientationAngleData & misori_s)
{
  const Real trans_misori_angle_HAGB = 15.0;
  const Real misori_angle = misori_s._misor;

  if (misori_angle <= 1.0)
      return _GBsigma_HAGB * (2.0 / trans_misori_angle_HAGB * (1 - std::log(2.0 / trans_misori_angle_HAGB)));
  else if (misori_angle <= trans_misori_angle_HAGB)
      return _GBsigma_HAGB * (misori_angle / trans_misori_angle_HAGB * (1 - std::log(misori_angle / trans_misori_angle_HAGB)));
      
  return _GBsigma_HAGB;
}
 
Real
GBAnisotropyMisoriAng::calculateGBMobility(const MisorientationAngleData & misori_s)
{
  const Real misori_angle = misori_s._misor;
  const Real trans_misori_angle_HAGB = 15.0;

  if (misori_angle <= 1.0)
    return _GBmob_HAGB * ((1-std::exp(-_B*std::pow(2.0/trans_misori_angle_HAGB,_n))));
  else if (misori_angle <= trans_misori_angle_HAGB)
    return _GBmob_HAGB * (1 - std::exp(-_B * std::pow(misori_angle / trans_misori_angle_HAGB, _n)));

  return _GBmob_HAGB;
}

void 
GBAnisotropyMisoriAng::fillSymmetricProperties(Real sigma_min, Real sigma_max, Real mob_min, Real mob_max)
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