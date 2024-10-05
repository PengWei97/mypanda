#include "GBAnisotropyWtStoredEnergy.h"

registerMooseObject("mypandaApp", GBAnisotropyWtStoredEnergy);

InputParameters
GBAnisotropyWtStoredEnergy::validParams()
{
  InputParameters params = GBAnisotropyMisori::validParams();
  params.addClassDescription("TODO");
  params.addParam<bool>("stored_energy_mobility", false, "The GB mobility related to stored energy would be considered if true");
  params.addRequiredParam<UserObjectName>("ebsd_reader", "User object name for the EBSD reader");  
  return params;
}

GBAnisotropyWtStoredEnergy::GBAnisotropyWtStoredEnergy(const InputParameters & parameters)
  : GBAnisotropyMisori(parameters),
  _stored_energy_mobility(getParam<bool>("stored_energy_mobility")),
  _ebsd_reader(getUserObject<EBSDReader>("ebsd_reader"))
{
}

void 
GBAnisotropyWtStoredEnergy::computeSigmaAndMobility(const std::vector<unsigned int> & var_index,
                                                    const std::vector<unsigned int> & grain_id_index)
{ 
  for (unsigned int i = 0; i < grain_id_index.size() - 1; ++i)
  {
    auto & grain_i = grain_id_index[i];
    auto angles_i = _euler.getEulerAngles(grain_i);
    auto rho_i = _ebsd_reader.getAvgData(grain_i)._custom[0]; // GNDs for each grain, 1/m^2

    for (unsigned int j = i + 1; j < grain_id_index.size(); ++j)
    {
      auto & grain_j = grain_id_index[j];
      auto angles_j = _euler.getEulerAngles(grain_j);
      auto rho_j = _ebsd_reader.getAvgData(grain_j)._custom[0]; // GNDs for each grain, 1/m^2

      auto & _sigma_ij = _sigma[var_index[i]][var_index[j]];
      auto & _mob_ij = _mob[var_index[i]][var_index[j]];

      // Calculate misorientation angle
      _misori_s = MisorientationAngleCalculator::calculateMisorientaion(angles_i, angles_j, _misori_s, _crystal_structure);

      // TODO
      Real delta_rho = std::abs(rho_i - rho_j) * _length_scale * _length_scale;
      _delta_rho[_qp] = delta_rho;

      // Compute sigma_ij
      _sigma_ij = _gb_energy_anisotropy ? calculatedGBEnergy(_misori_s) : _GBsigma_HAGB;

      // Compute mob_ij
      _mob_ij = _gb_mobility_anisotropy ? calculatedGBMobility(_misori_s) : _GBmob_HAGB;

      // Compute mob_ij with rho
      if (_stored_energy_mobility)
       calculatedGBMobilityWtRho(delta_rho, _mob_ij);

      _sigma[var_index[j]][var_index[i]] = _sigma_ij;
      _mob[var_index[j]][var_index[i]] = _mob_ij;
    }
  }
}

void
GBAnisotropyWtStoredEnergy::calculatedGBMobilityWtRho(const Real & delta_rho, Real & mob_ij)
{
  Real mob_ij_high = mob_ij * 10.0;
  Real mob_temp = mob_ij;

  Real trans_misori_rho = 90; // 1/(\miu m)^2

  // Equation constant
  Real B = 5;
  Real n = 4;

  if (delta_rho <= trans_misori_rho)
    mob_temp = mob_ij_high * ((1- std::exp(-B * std::pow( delta_rho / trans_misori_rho, n)))); // Eq.8
  else
    mob_temp = mob_ij_high;

  if (mob_temp > mob_ij)
    mob_ij = mob_temp;
}
