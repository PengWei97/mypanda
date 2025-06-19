#pragma once

#include "GBAnisotropyMisoriAndTwin.h"
#include "EBSDReaderMaterialProperty.h"

/**
 * GBAnisotropyWtDeltaRho:
 * A class to calculate grain boundary (GB) mobility considering 
 * the difference in dislocation density (Delta Rho) and twin boundary effects.
 */
class GBAnisotropyWtDeltaRho : public GBAnisotropyMisoriAndTwin
{
public:
  static InputParameters validParams();

  GBAnisotropyWtDeltaRho(const InputParameters & parameters);

protected:
  // initial _delta_rho
  virtual void initOthersMaterialProperties() override;
  
  // override to compute the other material properties based on grain id
  virtual void calculateOthersMaterialProperties(const unsigned int & grain_i, const unsigned int & grain_j) override;

  // Override to compute GB mobility with consideration of Delta Rho
  virtual Real calculateGBMobility(const MisorientationAngleData & misori_s) override;

  // Compute GB mobility adjustment based on dislocation density difference
  virtual void calculatedGBMobilityWtRho(const Real & delta_rho, Real & mob_ij);

  // Reference to the EBSD data reader
  const EBSDReaderMaterialProperty & _GNDs_provider;

  // Flag indicating if GB mobility calculation considers Delta Rho
  const bool _is_gb_mob_with_delta_rho;

  // Critical misorientation angle for transformation considering dislocation density
  const Real _trans_delta_rho;

  // Amplification factor for enhancing GB mobility
  const Real _amplifier_factor;

  // Material property to store the dislocation density difference (Delta Rho)
  MaterialProperty<Real> & _delta_rho;
};