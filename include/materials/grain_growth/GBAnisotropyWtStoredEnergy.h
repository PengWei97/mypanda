#pragma once

#include "GBAnisotropyMisori.h"
#include "EBSDReader.h"

/**
 * TODO
 */
class GBAnisotropyWtStoredEnergy : public GBAnisotropyMisori
{
public:
  static InputParameters validParams();

  GBAnisotropyWtStoredEnergy(const InputParameters & parameters);

protected:  
  // Calculates the sigma_ij and mob_ij values
  virtual void computeSigmaAndMobility(const std::vector<unsigned int> & var_index,
                                       const std::vector<unsigned int> & grain_id_index) override;

  virtual void calculatedGBMobilityWtRho(const Real & delta_rho, Real & mob_ij);

  const bool _stored_energy_mobility;
  const EBSDReader & _ebsd_reader;
  MaterialProperty<Real> & _delta_rho;
};